inFN = "testfiles/indexed.expt"

import numpy as np
import re


def to_arr(string):
  print(string)
  return np.array([float(i) for i in re.sub(r"[^0-9\.\,]", "", string).split(',')])

text = open(inFN).read()
A =  to_arr(re.search(r"(?<=real_space_a\":).*?\]", text, flags=re.DOTALL).group()  )
B =  to_arr(re.search(r"(?<=real_space_b\":).*?\]", text, flags=re.DOTALL).group()  )
C =  to_arr(re.search(r"(?<=real_space_c\":).*?\]", text, flags=re.DOTALL).group()  )
G =  to_arr(re.search(r"(?<=rotation_axis\":).*?\]", text, flags=re.DOTALL).group() )

A,B,C,G = np.array(A),np.array(B),np.array(C),np.array(G)

print(f"A: {np.rad2deg(np.arccos(np.dot(A / np.linalg.norm(A, 2), G)))} degrees")
print(f"B: {np.rad2deg(np.arccos(np.dot(B / np.linalg.norm(B, 2), G)))} degrees")
print(f"C: {np.rad2deg(np.arccos(np.dot(C / np.linalg.norm(C, 2), G)))} degrees")
