import numpy as np
import matplotlib.pyplot as plt
import math

longRes = 360
latRes = 180

cosines = np.zeros((latRes, longRes))

for lat in range(latRes):
    for long in range(longRes):
        ph = lat*math.pi/latRes
        th = 2*long*math.pi/longRes
        angle = math.cos(ph)*math.cos(0) + math.sin(ph)*math.sin(0)*math.cos(th)
        cosines[lat][long] = angle

y=np.linspace(-math.pi/2, math.pi/2, latRes)
x=np.linspace(-math.pi, math.pi, longRes)
plt.figure(1)
plt.subplot(projection='aitoff')
plt.pcolormesh(x, y, cosines, cmap='viridis', shading='gouraud')
plt.colorbar()
plt.show()
