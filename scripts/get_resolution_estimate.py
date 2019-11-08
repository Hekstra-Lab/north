inFN = "testfiles/dials.integrate.log"

import numpy as np
import pandas as pd
from io import StringIO
import re

sigma_cutoff = 1.5

text = open(inFN).read()


table = re.search(r'(?<=Summary vs resolution).*?\n\n', text, re.DOTALL).group()
table = re.sub(r'[|-]', '', table).strip()

df = pd.read_csv(StringIO(table), skiprows=3, delim_whitespace=True, usecols=[1, 10], names=['d_min', 'SNR'])

print(df[df['SNR'] > sigma_cutoff]['d_min'].min())


