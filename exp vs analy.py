import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

mu = 1;
x = 0.03125;

T = np.linspace(0, 1, 81)

#PLOT Y
#y1 = pd.read_csv("y_final_experimental_{}.txt".format(x), header=None)
y2 = []

for j in range(len(T)):
	t = T[j]
	y = np.around(math.sin(2 * t * math.pi), decimals=5) * np.around(math.sin(x * math.pi), decimals=5)
	y2.append(y)

fig = plt.gcf()
plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

#plt.plot(T, y1, "k+")
#plt.plot(T, y1, label="y exp")

#plt.plot(T, y2, "kx")
plt.plot(T, y2, label="y analytical")

plt.xlabel("Time step")
plt.ylabel("y")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical Y_{}.png".format(x))
plt.show()

#PLOT P
#p1 = pd.read_csv("p_final_experimental_{}.txt".format(x), header=None)
P = []

for j in range(len(T)):
	t = T[j]	
	p = np.around(math.sin(2 * t * math.pi), decimals=5) * np.around(math.sin(x * math.pi), decimals=5)
	P.append(p)
	
U = []

for j in range(len(T)):
	t = T[j]	
	u = (-1.0 / mu) * np.around(math.sin(2 * t * math.pi), decimals=5) * np.around(math.sin(x * math.pi), decimals=5)
	U.append(u)
		
fig = plt.gcf()

plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

#plt.plot(t, p1, "k+")
#plt.plot(T, p1, label="p exp")

#plt.plot(T, P, "kx")
plt.plot(T, P, label="p analytical")
plt.plot(T, U, label="u analytical")

plt.xlabel("Time step")
plt.ylabel("p")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical P_{}.png".format(x))
plt.show()


