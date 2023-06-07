import numpy as np
import matplotlib.pyplot as plt

particles = {}
gammas = {}
gammaE = np.array([])
elecE = np.array([])
numToSym = {'11':r'e$^-$', '22':r'$\gamma$', '1000190400':r'$^{40}$K', '1000180400':r'$^{40}$Ar', 
            '2112':'n', '2212':'p', '12':r'$\nu_e$', '1000190390':r'$^{39}$K', 
            '1000180390':r'$^{39}$Ar', '1000170360':r'$^{36}$Cl', '1000020040':r'$\alpha$',
            '1000180380':r'$^{38}$Ar', '1000010020':r'D', '1000170350':r'$^{35}$Cl'}

file = open("build/histE.ascii")
line = file.readline()
line = file.readline()

while line:
    spl = line.split()
    ni = int(spl[0]); nf = int(spl[1])

    for j in range(ni):
        line = file.readline()
    
    numGamma = 0
    for k in range(nf):
        line = file.readline()
        spl = line.split()
        sym = numToSym[spl[0]]
        print(sym)
        if sym in particles.keys():
            particles[sym] = particles[sym]+1
        else: particles[sym] = 1
        if spl[0] == '22':
            numGamma += 1
            gammaE = np.append(gammaE, float(spl[1]))
        if spl[0] == '11':
            elecE = np.append(elecE, float(spl[1]))
    
    if numGamma in gammas.keys():
        gammas[numGamma] = gammas[numGamma]+1
    else: gammas[numGamma] = 1

    line = file.readline()

plt.figure(1)
plt.bar(particles.keys(), particles.values())
plt.title('Particles produced in Marley simulation')
plt.xlabel('particle type')
plt.ylabel('counts')

plt.figure(2)
plt.hist(gammaE, bins=20)
plt.title('Energy distribution of photons')
plt.xlabel('energy (MeV)')

plt.figure(3)
plt.bar(gammas.keys(), gammas.values())
plt.title('Number of photons produced per event')

plt.figure(4)
plt.hist(elecE, bins=20)
plt.title('Energy distribution of electrons')
plt.xlabel('energy (MeV)')
plt.show()

file.close()