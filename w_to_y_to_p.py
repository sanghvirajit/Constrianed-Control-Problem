import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

x = 0.5

w = pd.read_csv("w_final.txt", header=None);
df = pd.DataFrame(index=np.arange((w.size)/62), columns=np.arange(1))

n = 15

for i in range(0, 129):
	df.iloc[i, 0] = w.iloc[i*62 + n, 0]
	
df.to_csv("y_final_experimental_{}.txt".format(x), header=None, index=False)

n = 46

for i in range(0, 129):
	df.iloc[i, 0] = w.iloc[i*62 + n, 0]

df.to_csv("p_final_experimental_{}.txt".format(x), header=None, index=False)
