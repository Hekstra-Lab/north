import numpy as np 
import pandas as pd
import gspread
import gspread_dataframe as gd
import argparse
import re
import tempfile
from multiprocessing import cpu_count
from io import StringIO
from os import listdir,devnull
from os.path import isdir,exists
from subprocess import call,STDOUT
from time import sleep
from oauth2client.service_account import ServiceAccountCredentials


description = """North is a command line server for indexing test diffraction images. Results are uploaded with some basic statistics to a Google Spreadsheet"""

parser = argparse.ArgumentParser(description=description)
parser.add_argument("data_directory", help="Directory containing test shots.", type=str)
parser.add_argument("database_filename", help="Location of the HDF5 database file saved by North. This can be a new file or North can be initialized from a pre-existing database.", type=str)
parser.add_argument("reference_geo_file", help="Reference DIALS experimental geometry file", type=str)
parser.add_argument("user_credentials", help="JSON file with Google API access credentials", type=str)
parser.add_argument("results_filename", help="Google sheets results filename. This needs to be created in advance by the user")
parser.add_argument("-t", "--temp", help="Temporary directory for running dials. The specified directory will not be cleaned on exit.",  type=str, default=None)
parser.add_argument("-c", "--sigma-cutoff", help="I over Sigma(I) cutoff for the resolution estimate in Ångstroms", type=float, default=1.5)
parser.add_argument("-p", "--nproc", help="Number of processors to use for dials programs.", default=cpu_count()-1,  type=int)
parser.add_argument("-w", "--wait-time", help="Wait time in seconds between Google API requests.",  type=int)
parser.add_argument("--space-group-number", help="Space group number for indexing. Defaults to 1 (P1).", default=1,  type=int)
parser.add_argument("-e", "--email", help="Email(s) of users you want to access this spreadsheet.", type=str, nargs='+')

parser = parser.parse_args()

if parser.temp is not None:
    tmpdir_name = parser.temp
else:
    tmpdir = tempfile.TemporaryDirectory()
    tmpdir_name = tmpdir.name

datadir       = parser.data_directory
reference_geo = parser.reference_geo_file
sigma_cutoff  = parser.sigma_cutoff
user_credential_file = parser.user_credentials
results_filename = parser.results_filename
nproc = parser.nproc
waittime = parser.wait_time
space_group = parser.space_group_number


if datadir[-1] != "/":
    datadir += "/"

scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']
credentials = ServiceAccountCredentials.from_json_keyfile_name(user_credential_file,scope)
gc = gspread.authorize(credentials)
sh = gc.create(results_filename)

if parser.email is not None:
    for addr in parser.email:
        sh.share(addr, perm_type='user', role='writer')

dbFN = parser.database_filename

#Global variables are evil
if exists(dbFN):
    database = pd.read_hdf(dbFN)
else:
    database = None

def check_for_new_data():
    directories = [datadir + i for i in listdir(datadir) if isdir(datadir + i)]
    directories = [i for i in directories if re.match(r".*base[0-9]*$", i, re.DOTALL) is not None]
    if database is not None:
       old = set(database['directory'])
    else:
       old = set()
    new = set(directories) - old
    return new

def to_arr(string):
    return np.array([float(i) for i in re.sub(r"[^0-9\.\,]", "", string).split(',')])

def parse_dials_output(directory):
    if directory[-1] != '/':
        directory += '/'
    inFN = directory + 'indexed.expt'
    print(f"Parsing {inFN} ...")

    entry = None
    if exists(inFN):
        text = open(inFN).read()
        A =  to_arr(re.search(r"(?<=real_space_a\":).*?\]", text, flags=re.DOTALL).group()  )
        B =  to_arr(re.search(r"(?<=real_space_b\":).*?\]", text, flags=re.DOTALL).group()  )
        C =  to_arr(re.search(r"(?<=real_space_c\":).*?\]", text, flags=re.DOTALL).group()  )
        G =  to_arr(re.search(r"(?<=rotation_axis\":).*?\]", text, flags=re.DOTALL).group() )
         
        A,B,C,G = np.array(A),np.array(B),np.array(C),np.array(G)
         
        
        entry = pd.DataFrame()
        entry["A1"]= [A[0]]
        entry["A2"]= [A[1]]
        entry["A3"]= [A[2]]
        entry["B1"]= [B[0]]
        entry["B2"]= [B[1]]
        entry["B3"]= [B[2]]
        entry["C1"]= [C[0]]
        entry["C2"]= [C[1]]
        entry["C3"]= [C[2]]
        
        A_offset = np.rad2deg(np.arccos(np.dot(A / np.linalg.norm(A, 2), G)))
        B_offset = np.rad2deg(np.arccos(np.dot(B / np.linalg.norm(B, 2), G)))
        C_offset = np.rad2deg(np.arccos(np.dot(C / np.linalg.norm(C, 2), G)))

        A_offset = np.min(np.abs([180-A_offset, A_offset]))
        B_offset = np.min(np.abs([180-B_offset, B_offset]))
        C_offset = np.min(np.abs([180-C_offset, C_offset]))
        entry["Vertical Axis"] = ['A', 'B', 'C'][np.argmin([A_offset, B_offset, C_offset])]
        entry["Alignment Error"] = [np.min([A_offset, B_offset, C_offset])]

        entry["A to gonio axis"]= [A_offset]
        entry["B to gonio axis"]= [B_offset]
        entry["C to gonio axis"]= [C_offset]

        inFN = directory + 'dials.index.log'
        text = open(inFN).read()
        cell = re.search(r"(?<=Unit cell:).*$", text, re.MULTILINE).group()
        cell = re.sub(r'[\(\)]', '', cell)
        cell = np.array(cell.split(','), dtype=np.float)
        a,b,c,alpha,beta,gamma = cell
        entry["a"]= a
        entry["b"]= b
        entry["c"]= c
        entry["alpha"]= alpha
        entry["beta"]= beta
        entry["gamma"]= gamma

        inFN = directory + 'dials.index.log'
        text = open(inFN).read()
        table = re.findall(r'(?<=% indexed).*?(?=\n\n)', text, re.DOTALL)[-1] 
        table =  re.sub(r'[-|]', '', table).strip()
        num_indexed = sum([int(i.split()[1]) for i in table.split('\n')])
        num_unindexed = sum([int(i.split()[2]) for i in table.split('\n')])
        entry['Percent Indexed'] = 100. * num_indexed / (num_indexed + num_unindexed)

    inFN   = directory + 'dials.integrate.log'
    reflFN = directory + 'integrated.refl'
    if exists(inFN) and exists(reflFN):
        print(f"Parsing {inFN} ...")
        text = open(inFN).read()
        table = re.search(r'(?<=Summary vs resolution).*?\n\n', text, re.DOTALL).group()
        table = re.sub(r'[|-]', '', table).strip()
        df = pd.read_csv(StringIO(table), skiprows=3, delim_whitespace=True, usecols=[1, 10], names=['d_min', 'SNR'])
        #estimate = df[df['SNR'] > sigma_cutoff]['d_min'].min()

    inFN   = directory + 'dials.resolutionizer.log'
    if exists(inFN) and exists(reflFN):
        text = open(inFN).read()
        estimate = float(re.search('(?<=Resolution I/sig:)[0-9\.\s]*$', text, re.MULTILINE).group())
        entry["Resolution Estimate"] = estimate

    return entry


def process_data(directory):
    process_dir= tmpdir_name + directory
    print(f"Processing directory {directory} in {process_dir}...")
    call(['mkdir', '-p', process_dir])

    FNULL = open(devnull, 'w')
    #import data
    command = f"dials.import reference_geometry={reference_geo} image_range=1,5 {directory}/*.cbf"
    call(command.split(), cwd=process_dir, stdout=FNULL, stderr=STDOUT)

    #find spots
    command = f"dials.find_spots imported.expt nproc={nproc}"
    call(command.split(), cwd=process_dir, stdout=FNULL, stderr=STDOUT)

    #index
    command= f"dials.index imported.expt strong.refl space_group={space_group}"
    call(command.split(), cwd=process_dir, stdout=FNULL, stderr=STDOUT)

    #refine
    command="dials.refine scan_varying=False indexed.expt indexed.refl"
    call(command.split(), cwd=process_dir, stdout=FNULL, stderr=STDOUT)

    #integrate
    command=f"dials.integrate refined.expt refined.refl nproc={nproc}"
    call(command.split(), cwd=process_dir, stdout=FNULL, stderr=STDOUT)

    #resolution estimate
    command=f"dials.resolutionizer integrated.expt integrated.refl"
    call(command.split(), cwd=process_dir, stdout=FNULL, stderr=STDOUT)

    FNULL.close()

    db_entry = parse_dials_output(process_dir)
    if db_entry is None:
        db_entry = pd.DataFrame()
    db_entry['directory'] = directory
    db_entry['base'] = int(re.search(r'(?<=base)[0-9]*', directory).group())
    return db_entry

def push_databse_to_sheets():
    #This is bad control flow
    #TODO: instantiate the database properly and sensibly
    gc = gspread.authorize(credentials)
    if database is not None and gc is not None and 'Vertical Axis' in database:
        workbook = gc.open(results_filename)
        sh1 = workbook.get_worksheet(0)
        sh2 = workbook.get_worksheet(1)
        if sh2 is None:
            sh2 = workbook.add_worksheet(title="Sheet2", rows=500, cols=20)

        gd.set_with_dataframe(sh1, database[["base", "Vertical Axis", "Alignment Error", "Percent Indexed", "Resolution Estimate", "a","b","c","alpha","beta","gamma"]])
        gd.set_with_dataframe(sh2, database)

while True:
    new_directories = check_for_new_data()
    for directory in new_directories:
        print(f"Working on directory: {directory}")
        db_entry = process_data(directory)
        if db_entry is None:
            print(f"Failed processing {directory}.")
        else:
            database = pd.concat((database, db_entry))
            database.sort_values(by='base', inplace=True)           
            database.to_hdf(dbFN, 'database')
        push_databse_to_sheets()
    sleep(waittime)
