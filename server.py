import numpy as np 
import pandas as pd
import re
from io import StringIO
from os import listdir,devnull
from os.path import isdir,exists
from subprocess import call,STDOUT
from time import sleep

datadir = "/home/userbmc/hekstra_201911/"
workdir= "/home/userbmc/processing_201911/"
reference_geo="/home/userbmc/processing/reference_geometry.expt"
sigma_cutoff = 1.5

if not exists(workdir):
    call(['mkdir', '-p', workdir])

space_group=1

waittime = 0
nproc=10


#Global variables are evil
database = None


if datadir[-1] != "/":
    datadir += "/"
if workdir[-1] != "/":
    workdir += "/"

dbFN = workdir + "database.h5"

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
        entry["A to gonio axis"]= [np.rad2deg(np.arccos(np.dot(A / np.linalg.norm(A, 2), G)))]
        entry["B to gonio axis"]= [np.rad2deg(np.arccos(np.dot(B / np.linalg.norm(B, 2), G)))]
        entry["C to gonio axis"]= [np.rad2deg(np.arccos(np.dot(C / np.linalg.norm(C, 2), G)))]

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



    
    inFN   = directory + 'dials.integrate.log'
    reflFN = directory + 'integrated.refl'
    if exists(inFN) and exists(reflFN):
        print(f"Parsing {inFN} ...")
        text = open(inFN).read()
        table = re.search(r'(?<=Summary vs resolution).*?\n\n', text, re.DOTALL).group()
        table = re.sub(r'[|-]', '', table).strip()
        df = pd.read_csv(StringIO(table), skiprows=3, delim_whitespace=True, usecols=[1, 10], names=['d_min', 'SNR'])
        estimate = df[df['SNR'] > sigma_cutoff]['d_min'].min()
        entry["Resolution Estimate"] = estimate

    return entry


def process_data(directory):
    process_dir="/tmp/hekstra/directory" + directory
    print(f"Processing directory {directory} in {process_dir}...")
    call(['mkdir', '-p', process_dir])

    FNULL = open(devnull, 'w')
    #import data
    command = f"dials.import invert_rotation_axis=True reference_geometry={reference_geo} image_range=1,5 {directory}/*.cbf"
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

    FNULL.close()

    db_entry = parse_dials_output(process_dir)
    if db_entry is None:
        db_entry = pd.DataFrame()
    db_entry['directory'] = directory
    db_entry['base'] = int(re.search(r'(?<=base)[0-9]*', '/home/userbmc/hekstra_201911/base223').group())
    return db_entry

while True:
    new_directories = check_for_new_data()
    for directory in new_directories:
        print(f"Working on directory: {directory}")
        db_entry = process_data(directory)
        if db_entry is None:
             print(f"Failed processing {directory}.")
        else:
             database = pd.concat((database, db_entry))
             database.to_hdf(dbFN, 'database')
        print(database)
    sleep(waittime)