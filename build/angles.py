import numpy as np
import matplotlib.pyplot as plt
import math

cutoff = 30
numBins = (cutoff-5)*2+1
ePhi = np.array([])
eAz = np.array([])
gPhi = np.array([])
gAz = np.array([])
nuE = np.array([])
mapGamma = np.zeros((numBins+1, numBins))
mapE = np.zeros((numBins+1, numBins))

# 2 MeV to 60 MeV, 0.5 binning, 1000 events per bin
# 2d plot of energy vs cos(angle)
#file = open("build/histE.ascii")
file = open(r"build/flatE.ascii")
#file = open("build/triangleE.ascii")
line = file.readline()
line = file.readline()

def getAngle(v1, v2):
    dot = 0
    mag1 = 0; mag2 = 0
    for i in range(len(v1)):
        dot += v1[i]*v2[i]
        mag1 += v1[i]**2; mag2 += v2[i]**2
    return dot/(mag1**0.5 * mag2**0.5)

counter = 1
while line:
    spl = line.split()
    ni = int(spl[0]); nf = int(spl[1])
    line = file.readline()
    counter += 1
    spl = line.split()
    assert int(spl[0]) == 12, 'Initial particle not neutrino'
    v1 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
    eng = float(spl[1])
    if eng > cutoff: break
    nuE = np.append(nuE, eng)
    print(eng, counter)
    line = file.readline()
    counter += 1

    for j in range(nf):
        line = file.readline()
        counter += 1
        spl = line.split()
        if int(spl[0]) == 11:
            v2 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
            angle = getAngle(v1, v2)
            ePhi = np.append(ePhi, angle)
            cross = np.cross(v1, v2)
            phase = math.pi/2
            if cross[1] < 0: phase += math.pi
            eAz = np.append(eAz, np.arctan(cross[1]/cross[0])+phase)
            mapE[round(angle*numBins)][round((eng-5)*2)] += 1
        elif int(spl[0]) == 22:
            v2 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
            angle = getAngle(v1, v2)
            gPhi = np.append(gPhi, angle)
            cross = np.cross(v1, v2)
            phase = math.pi/2
            if cross[1] < 0: phase += math.pi
            gAz = np.append(gAz, np.arctan(cross[1]/cross[0])+phase)
            mapGamma[round(angle*numBins)][round((eng-5)*2)] += 1
    line = file.readline()
    counter += 1

plt.figure(1)
plt.hist(np.arccos(ePhi), bins=numBins)
plt.title('Angular distribution of electrons produced in MARLEY simulation')
plt.xlabel('angle')
plt.ylabel('count')

plt.figure(2)
plt.hist(np.arccos(gPhi), bins=numBins)
plt.title('Angular distribution of photons produced in MARLEY simulation')
plt.xlabel('angle')
plt.ylabel('count')

plt.figure(3)
plt.hist(ePhi, bins=numBins)
plt.title('Angular distribution of electrons produced in MARLEY simulation')
plt.xlabel('cos(angle)')
plt.ylabel('count')

plt.figure(4)
plt.hist(gPhi, bins=numBins)
plt.title('Angular distribution of photons produced in MARLEY simulation')
plt.xlabel('cos(angle)')
plt.ylabel('count')

plt.figure(5)
plt.hist(eAz, bins=numBins)
plt.title('Azimuthal distribution of electrons produced in MARLEY simulation')
plt.xlabel('angle')
plt.ylabel('count')

plt.figure(6)
plt.hist(gAz, bins=numBins)
plt.title('Azimuthal distribution of photons produced in MARLEY simulation')
plt.xlabel('angle')
plt.ylabel('count')

plt.figure(7)
plt.hist(nuE, bins=numBins)
plt.title(r'Energy distribution of incoming $\nu_e$')
plt.xlabel('energy (MeV)')
plt.ylabel('count')

y=np.linspace(0., 1., numBins+1)
x=np.arange(5., cutoff+0.1, 0.5)
plt.figure(8)
plt.pcolormesh(x, y, mapE, cmap='viridis', shading='nearest')
plt.title('Correspondence between neutrino energy and \n electron angular distribution')
plt.xlabel('energy (MeV)')
plt.ylabel('cos(angle)')
plt.colorbar()

plt.figure(9)
plt.pcolormesh(x, y, mapGamma, cmap='viridis', shading='nearest')
plt.title('Correspondence between neutrino energy and \n photon angular distribution')
plt.xlabel('energy (MeV)')
plt.ylabel('cos(angle)')
plt.colorbar()
plt.show()

file.close()