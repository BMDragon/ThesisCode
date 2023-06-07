import numpy as np
import matplotlib.pyplot as plt
import math

fermiCorrect = .1
gtCorrect = .1

minima = np.load('./build/Skymaps/stdevRot' + str(int(10*fermiCorrect)) + str(int(10*gtCorrect)) + '.npy')
tdc = np.array([0.0, 0.6, -0.2])
trueDir = np.array([math.atan2(tdc[0], tdc[1]), math.acos(tdc[2]/np.linalg.norm(tdc))])

def angleDist(th, ph):
    theta = trueDir[0]
    phi = trueDir[1]
    ph = -ph + math.pi/2
    angle = math.cos(ph)*math.cos(phi) + math.sin(ph)*math.sin(phi)*math.cos(th-theta)
    return math.acos(angle)*180/math.pi

angles = np.array([angleDist(*x) for x in minima])
print(np.mean(angles), np.std(angles))

pointSize = 16
transparancy = 0.9
fig = plt.figure(figsize=(8,6))
ax = fig.add_subplot(projection = "aitoff")
plt.scatter(minima[:,0], minima[:,1], c='b', label='minima', s=pointSize, alpha=transparancy)
ax.grid()
plt.scatter(trueDir[0], -trueDir[1]+math.pi/2, c='r',label='truth', marker='x', s=pointSize)
plt.title('Pointing resolution\nFermi correctness ratio: ' + str(fermiCorrect) + '\n     GT correctness ratio: ' + str(gtCorrect))
plt.legend(loc=1)
#plt.show()
plt.savefig('./build/Skymaps/standardDevRot' + str(int(10*fermiCorrect)) + str(int(10*gtCorrect)))