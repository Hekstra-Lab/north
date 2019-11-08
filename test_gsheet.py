import pandas as pd
import gspread as gs
import gspread_dataframe as gd
from oauth2client.service_account import ServiceAccountCredentials

scope = ['https://docs.google.com', 'https://googleapis.com/auth/drive']

credentials = ServiceAccountCredentials('/home/userbmc/aps_bmc_screening.json', scope)

gc = gs.authorize(credentials)

sheet = gc.open("Screening Results").sheet1
