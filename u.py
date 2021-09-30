import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

u = pd.read_csv("u_projected.txt", header=None)
df = pd.DataFrame(index=np.arange((u.size)/31), columns=np.arange(1))

n = 15

for i in range(0, 129):
	df.iloc[i, 0] = u.iloc[i*31 + n, 0]
	
df.to_csv("u_final.txt", header=None, index=False)
