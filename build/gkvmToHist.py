import numpy as np
import matplotlib.pyplot as plt

data = np.genfromtxt('./build/gkvm.dat', delimiter=' ')

eng = np.array(data[:, 0])
nue = np.array(data[:, 1])
p2 = np.array(data[:, 2])
p3 = np.array(data[:, 3])
nueBar = np.array(data[:, 4])
p5 = np.array(data[:, 5])
p6 = np.array(data[:, 6])

nux = p2+p3+p5+p6

with open('./build/gkvmHist.txt', 'w') as file:
    file.write('[')
    for j in range(len(eng)):
        file.write(" " + str(round(eng[j]*1000,1)) + ",")
    file.write("]\n[")
    for k in range(len(nue)):
        file.write(" " + str(nue[k]) + ",")
    file.write(']')

plt.figure()
plt.plot(eng, nue, 'r')
plt.plot(eng, nueBar, 'k')
plt.plot(eng, nux, 'g')
plt.show()