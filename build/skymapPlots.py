import numpy as np
import matplotlib.pyplot as plt
import math

data = np.load('./build/skymapDataRot.npy', allow_pickle=True)
numBins = data[0]
maxNuE = data[1]
phiF = data[2]
thetaF = data[3]
phiGT = data[4]
thetaGT = data[5]
latRes = data[6]
longRes = data[7]
likelihoods = data[8]

z = np.where(likelihoods == np.min(likelihoods))
minY = z[0][0]*math.pi/latRes - math.pi/2
minX = 2*z[1][0]*math.pi/longRes - math.pi

print(z[0], z[1])
print(minX, minY)

tdc = np.array([0.0, 0.6, -0.2])
trueDir = np.array([math.atan2(tdc[0], tdc[1]), math.acos(tdc[2]/np.linalg.norm(tdc))])
def angleDist(th, ph):
    sigTh = math.pi/180
    sigPh = math.pi/180
    theta = trueDir[0]
    phi = trueDir[1]
    ph = -ph + math.pi/2
    angle = math.cos(ph)*math.cos(phi) + math.sin(ph)*math.sin(phi)*math.cos(th-theta)
    dfdph = -math.sin(ph)*math.cos(phi) + math.cos(ph)*math.sin(phi)*math.cos(th-theta)
    dfdth = -math.sin(ph)*math.sin(phi)*math.sin(th-theta)
    uncCos = math.sqrt(dfdph**2 * sigPh**2 + dfdth**2 * sigTh**2)
    unc = (np.abs(math.acos(angle)-math.acos(angle-uncCos)) + np.abs(math.acos(angle)-math.acos(angle+uncCos)))/2
    return math.acos(angle)*180/math.pi, unc*180/math.pi

print(angleDist(minX, minY))
print(trueDir[0], -trueDir[1]+math.pi/2)

pointSize = 8
transparancy = 0.4
fig1 = plt.figure(figsize=(8,6))
ax1 = fig1.add_subplot(projection = "aitoff")
plt.scatter(thetaF, -phiF+math.pi/2, c='b', label='Fermi', s=pointSize, alpha=transparancy)
plt.scatter(thetaGT, -phiGT+math.pi/2, c='r', label='Gamow-Teller', s=pointSize, alpha=transparancy)
ax1.grid()
plt.scatter(trueDir[0], -trueDir[1]+math.pi/2, c='k', marker='x', s=5*pointSize, label='truth')
plt.legend(loc=4)
plt.title('Distribution of charged-current interactions\n ')
plt.savefig('./build/Skymaps/eventDistRot')

y=np.linspace(-math.pi/2, math.pi/2, latRes)
x=np.linspace(-math.pi, math.pi, longRes)
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(projection = "aitoff")
plt.pcolormesh(x, y, likelihoods, cmap='viridis_r', shading='gouraud')
plt.colorbar()
ax.grid()
plt.scatter(minX, minY, c='k', s=4*pointSize, label='minimum')
plt.scatter(trueDir[0], -trueDir[1]+math.pi/2, c='r', marker='x', s=4*pointSize, label='truth')
plt.legend(loc=1)
plt.title('Supernova negative log-likelihood skymap\n')
plt.savefig('./build/Skymaps/likelihoodRot')