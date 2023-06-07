import numpy as np
from scipy.stats import chi2
import matplotlib.pyplot as plt
import math

fermiCorrect = 0.9
gtCorrect = 0.9

data0 = np.load('./build/Wilks/wilksRot' + str(int(100*fermiCorrect)) + str(int(100*gtCorrect)) + '.npy')
data = chi2.cdf(data0, 2)
latRes = 90
longRes = 180

minim = np.where(data == 0)
minX = 2*minim[1][0]*math.pi/longRes - math.pi
minY = minim[0][0]*math.pi/latRes - math.pi/2

tdc = np.array([0.0, 0.6, -0.2])
trueDir = np.array([math.atan2(tdc[0], tdc[1]), math.acos(tdc[2]/np.linalg.norm(tdc))])

print(np.max(data), np.min(data))

x=np.linspace(-math.pi, math.pi, longRes)
y=np.linspace(-math.pi/2, math.pi/2, latRes)
fig = plt.figure(figsize=(12,9))
ax = fig.add_subplot(projection = "aitoff")
plt.contour(x, y, data, levels=[0.95], extend='both', colors=['red'])
plt.pcolormesh(x, y, data0, cmap='viridis_r', shading='nearest')
plt.colorbar()
ax.grid()
plt.scatter(trueDir[0], -trueDir[1]+math.pi/2, c='k', marker='x')
plt.scatter(minX, minY, c='k')
plt.title('Fermi correctness ratio: ' + str(fermiCorrect) + '\nGT correctness ratio: ' + str(gtCorrect))
plt.show()
