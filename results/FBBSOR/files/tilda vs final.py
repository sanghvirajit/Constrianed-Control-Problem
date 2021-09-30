import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import math

mu = 1;
x = 0.5;

T = np.linspace(0, 1, 129)

#PLOT Y
z1 = []
y1 = []
u1 = []

for j in range(len(T)):
	t = T[j]
	y = (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5)
	y1.append(y)
		
for j in range(len(T)):
	t = T[j]
	z = (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5) * (1 + mu * math.pi**4) - 2 * mu * np.around(math.sin(x * math.pi), decimals=5)
	z1.append(z)

y2 = pd.read_csv("y_tilda_experimental_{}.txt".format(x), header=None)
y3 = pd.read_csv("y_tilda_vector_experimental_{}.txt".format(x), header=None)
y4 = pd.read_csv("y_final_experimental_{}.txt".format(x), header=None)

fig = plt.gcf()
plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

#plt.plot(t, y3, "k+")
plt.plot(T, y2, label="y tilda")
plt.plot(T, y3, label="y tilda vector")
plt.plot(T, y4, label="y final")

#plt.plot(t, y2, "kx")
plt.plot(T, y1, label="y analytical")
#plt.plot(T, z1, label="z")

plt.xlabel("Time step")
plt.ylabel("y")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical Y_{}.png".format(x))
plt.show()

#PLOT P
P = []

for j in range(len(T)):
	t = T[j]	
	p = -1 * mu * ( (-1.0) * 2 * (1-t) * np.around(math.sin(x * math.pi), decimals=5) + math.pi**2 * (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5))
	P.append(p)

p2 = pd.read_csv("p_tilda_experimental_{}.txt".format(x), header=None)
p3 = pd.read_csv("p_tilda_vector_experimental_{}.txt".format(x), header=None)
p4 = pd.read_csv("p_final_experimental_{}.txt".format(x), header=None)

fig = plt.gcf()

plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

plt.plot(T, p2, label="p tilda")
plt.plot(T, p3, label="p tilda vector")
plt.plot(T, p4, label="p final")

plt.plot(T, P, label="p analytical")

plt.xlabel("Time step")
plt.ylabel("p")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical P_{}.png".format(x))
plt.show()

U = []

for j in range(len(T)):
	t = T[j]	
	u = ( (-1.0) * 2 * (1-t) * np.around(math.sin(x * math.pi), decimals=5) + math.pi**2 * (1-t)**2 * np.around(math.sin(x * math.pi), decimals=5))
	U.append(u)
	
u_final = pd.read_csv("u_final.txt", header=None)

fig = plt.gcf()

plt.plot([], [], ' ', label="parameter v={}".format(mu))
plt.plot([], [], ' ', label="parameter x={}".format(x))

plt.plot(T, u_final, label="u final")
plt.plot(T, U, label="u analytical")

plt.xlabel("Time step")
plt.ylabel("u")

plt.legend()
plt.grid()

fig.savefig("exp vs analytical u_{}.png".format(x))
plt.show()

