import numpy as np
import matplotlib.pyplot as plt

fermiE = np.array([])
gamowTellerE = np.array([])
numBins = 100
maxNuE = 101
ePhiF = np.array([])
ePhiGT = np.array([])
eEngF = np.array([])
eEngGT = np.array([])
nuEF = np.array([])
nuEGT = np.array([])
nuETot = np.array([])
mapEF = np.zeros((numBins+1, maxNuE))
mapEGT = np.zeros((numBins+1, maxNuE))
location = './'

fileTag = 'BigE'

file = open(location + "gkvmBig.ascii")
line = file.readline()
line = file.readline()

def getAngle(v1, v2): # Returns cos(angle)
    dot = 0
    mag1 = 0; mag2 = 0
    for i in range(len(v1)):
        dot += v1[i]*v2[i]
        mag1 += v1[i]**2; mag2 += v2[i]**2
    return dot/(mag1**0.5 * mag2**0.5)

counter = 1
while line:
    fermi = False
    spl = line.split()
    ni = int(spl[0]); nf = int(spl[1]); ex = float(spl[2])
    assert ni == 2, 'More than 2 initial particles in collision'

    if round(ex, 4) == 4.3837:
        fermi=True
        fermiE = np.append(fermiE, ex)
    else: gamowTellerE = np.append(gamowTellerE, ex)

    line = file.readline()
    counter += 1
    spl = line.split()
    assert int(spl[0]) == 12, 'Initial particle not neutrino'
    v1 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
    eng = float(spl[1]) # neutrino energy
    if fermi: nuEF = np.append(nuEF, eng)
    else: nuEGT = np.append(nuEGT, eng)
    nuETot = np.append(nuETot, eng)
    if counter%1e5 < 10: print(eng, counter)
    line = file.readline()
    counter += 1

    for j in range(nf):
        line = file.readline()
        counter += 1
        spl = line.split()
        if int(spl[0]) == 11:
            v2 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
            angle = getAngle(v1, v2)
            energy = float(spl[1]) # electron energy
            if fermi:
                ePhiF = np.append(ePhiF, angle)
                eEngF = np.append(eEngF, energy)
                mapEF[round(angle*numBins/2+numBins/2)][round(energy)] += 1
            else:
                ePhiGT = np.append(ePhiGT, angle)
                eEngGT = np.append(eEngGT, energy)
                mapEGT[round(angle*numBins/2+numBins/2)][round(energy)] += 1

    line = file.readline()
    counter += 1

plt.figure(1)
plt.hist([gamowTellerE, fermiE], bins=300, color=['r','b'], label=['Gamow-Teller','Fermi'])
plt.title('Fermi vs Gamow-Teller Transitions')
plt.xlabel('Transition Energy (MeV)')
plt.ylabel('count')
plt.legend()

plt.figure(2)
valueEPhiF, binsEPhiF, patches = plt.hist(ePhiF, bins=numBins)
plt.title('Angular distribution of electrons produced in MARLEY simulation\n(Only Fermi transitions)')
plt.xlabel('cos(angle)')
plt.ylabel('count')
plt.savefig(location + 'Figures/FermiAngles' + fileTag)

plt.figure(3)
valueEPhiGT, binsEPhiGT, patches = plt.hist(ePhiGT, bins=numBins)
plt.title('Angular distribution of electrons produced in MARLEY simulation\n(Only Gamow-Teller transitions)')
plt.xlabel('cos(angle)')
plt.ylabel('count')
plt.savefig(location + 'Figures/GTAngles' + fileTag)

plt.figure(4)
plt.hist(eEngF, bins=numBins)
plt.title('Total energy of electrons produced in MARLEY simulation\n(Only Fermi transitions)')
plt.xlabel('energy (MeV)')
plt.ylabel('count')
plt.savefig(location + 'Figures/FermiElecEnergyDist' + fileTag)

plt.figure(5)
plt.hist(eEngGT, bins=numBins)
plt.title('Total energy of electrons produced in MARLEY simulation\n(Only Gamow-Teller transitions)')
plt.xlabel('energy (MeV)')
plt.ylabel('count')
plt.savefig(location + 'Figures/GTElecEnergyDist' + fileTag)

plt.figure(6)
plt.hist(nuEF, bins=numBins)
plt.title(r'Energy distribution of incoming $\nu_e$'+'\n(Only Fermi transitions)')
plt.xlabel('energy (MeV)')
plt.ylabel('count')

plt.figure(7)
plt.hist(nuEGT, bins=numBins)
plt.title(r'Energy distribution of incoming $\nu_e$'+'\n(Only Gamow-Teller transitions)')
plt.xlabel('energy (MeV)')
plt.ylabel('count')

plt.figure(8)
valueNuEEng, binsNuEEng, patches = plt.hist(nuETot, bins=numBins)
plt.title(r'Total Energy distribution of incoming $\nu_e$')
plt.xlabel('energy (MeV)')
plt.ylabel('count')

y=np.linspace(-1., 1., numBins+1)
x=np.arange(0., maxNuE, 1.)
plt.figure(9)
fermiQuadmesh = plt.pcolormesh(x, y, mapEF, cmap='viridis', shading='nearest')
plt.title('Correspondence between neutrino energy and electron \n angular distribution (Only Fermi transitions)')
plt.xlabel('energy (MeV)')
plt.ylabel('cos(angle)')
plt.colorbar()
plt.savefig(location + 'Figures/FermiHeatmap' + fileTag)

plt.figure(10)
gtQuadmesh = plt.pcolormesh(x, y, mapEGT, cmap='viridis', shading='nearest')
plt.title('Correspondence between neutrino energy and electron \n angular distribution (Only Gamow-Teller transitions)')
plt.xlabel('energy (MeV)')
plt.ylabel('cos(angle)')
plt.colorbar()
plt.savefig(location + 'Figures/GTHeatmap' + fileTag)

saveData = np.array([numBins, maxNuE, fermiQuadmesh, gtQuadmesh], dtype=object)
np.save(location + 'histData', saveData)