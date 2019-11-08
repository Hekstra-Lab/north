inFN = "testfiles/dials.index.log"
import numpy as np
import re


text = open(inFN).read()


cell = re.search(r"(?<=Unit cell:).*$", text, re.MULTILINE).group()
print(cell)
cell = re.sub(r'[\(\)]', '', cell)
cell = np.array(cell.split(','), dtype=np.float)
print(cell)

a,b,c,alpha,beta,gamma = cell
