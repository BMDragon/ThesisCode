import numpy as np
import matplotlib.pyplot as plt
import math

filepath = './build/'
data = np.load(filepath + 'histData.npy', allow_pickle=True)
numBins = data[0]
maxNuE = data[1]
fermiQuadMesh = data[2].get_array().reshape(numBins+1, maxNuE)
gtQuadMesh = data[3].get_array().reshape(numBins+1, maxNuE)

# For each electron:
#   find the angle and energy
#   for each point in sky:
#       calculate negative log likelihood electron was from there
#       add this value to previous values for same point
#   

fermiNormConst = sum(np.array([sum(fermiQuadMesh[:,x]) for x in range(len(fermiQuadMesh[0]))]))
gtNormConst = sum(np.array([sum(gtQuadMesh[:,x]) for x in range(len(gtQuadMesh[0]))]))

thetaF = np.array([])
phiF = np.array([])
thetaGT = np.array([])
phiGT = np.array([])

longRes = 360
latRes = 180

likelihoods = np.zeros((latRes, longRes))

def addNLL(theta, phi, energy, isFermi):
    if isFermi: qm = fermiQuadMesh; const = fermiNormConst
    else: qm = gtQuadMesh; const = gtNormConst
    for lat in range(latRes):
        for long in range(longRes):
            ph = lat*math.pi/latRes
            th = 2*long*math.pi/longRes
            angle = math.cos(ph)*math.cos(phi) + math.sin(ph)*math.sin(phi)*math.cos(th-theta)
            lh = qm[round(angle*numBins/2+numBins/2)][round(energy)]/const
            if lh <= 0: lh = 1e-10
            likelihoods[lat][long] += -1*math.log(lh)

file = open(filepath + "gkvm.ascii")
line = file.readline()
line = file.readline()

counter = 1
while line:
    fermi = False
    spl = line.split()
    ni = int(spl[0]); nf = int(spl[1]); ex = float(spl[2])
    assert ni == 2, 'More than 2 initial particles in collision'

    if round(ex, 4) == 4.3837:
        fermi=True

    line = file.readline()
    counter += 1
    spl = line.split()
    assert int(spl[0]) == 12, 'Initial particle not neutrino'
    v1 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
    line = file.readline()
    counter += 1

    for j in range(nf):
        line = file.readline()
        counter += 1
        spl = line.split()
        # For electrons produced
        if int(spl[0]) == 11:            
            eng = float(spl[1]) # electron energy
            v2 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
            th = math.atan2(v2[0], v2[1])
            ph = math.acos(v2[2]/np.linalg.norm(v2))
            if fermi:
                thetaF = np.append(thetaF, th)
                phiF = np.append(phiF, ph)
            else:
                thetaGT = np.append(thetaGT, th)
                phiGT = np.append(phiGT, ph)
            addNLL(th, ph, eng, fermi)
    line = file.readline()
    counter += 1

saveData = np.array([numBins, maxNuE, phiF, thetaF, phiGT, thetaGT, 
                     latRes, longRes, likelihoods], dtype=object)
np.save(filepath + 'skymapData', saveData)