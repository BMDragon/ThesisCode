import numpy as np
import matplotlib.pyplot as plt
import math
import sys

rand = np.random.default_rng(3141)

fermiCorrect = 1
gtCorrect = 1
params = sys.argv
if len(params) >= 3:
    fCorrect = float(params[1])
    gCorrect = float(params[2])
    assert 0. <= fCorrect <= 1. and 0. <= gCorrect <= 1., \
        'Argument out of range. Please input a number in the closed interval [0., 1.]'
    fermiCorrect = fCorrect
    gtCorrect = gCorrect
elif len(params) == 2:
    correct = float(params[1])
    assert 0. <= correct <= 1., \
        'Argument out of range. Please input a number in the closed interval [0., 1.]'
    fermiCorrect = correct
    gtCorrect = correct

filepath = './'
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

longRes = 180
latRes = 90
eventsPerSN = 3000
numSN = 100

minima = np.zeros((numSN, 2))
angles = np.array([])
likelihoods = np.zeros((latRes, longRes))

def addNLL(theta, phi, energy, isFermi, randNum):
    if (isFermi and randNum < fermiCorrect) or (not isFermi and randNum >= gtCorrect): 
        qm = fermiQuadMesh; const = fermiNormConst
    else: qm = gtQuadMesh; const = gtNormConst
    for lat in range(latRes):
        for long in range(longRes):
            ph = lat*math.pi/latRes
            th = 2*long*math.pi/longRes
            angle = math.cos(ph)*math.cos(phi) + math.sin(ph)*math.sin(phi)*math.cos(th-theta)
            lh = qm[round(angle*numBins/2+numBins/2)][round(energy)]/const
            if lh <= 0: lh = 1e-10
            likelihoods[lat][long] += -1*math.log(lh)

tdc = np.array([0.0, 0.6, -0.2])
trueDir = np.array([math.atan2(tdc[0], tdc[1]), math.acos(tdc[2]/np.linalg.norm(tdc))])
def pointToAngle(long, lat):
    th = long
    ph = -lat+math.pi/2
    theta = trueDir[0]
    phi = -trueDir[1]+math.pi/2
    angle = math.cos(ph)*math.cos(phi) + math.sin(ph)*math.sin(phi)*math.cos(th-theta)
    return math.acos(angle)

file = open(filepath + "stdevDataRot.ascii")
line = file.readline()
line = file.readline()

eventsCounter = 0
snCounter = 0
while line and snCounter < numSN:
    fermi = False
    spl = line.split()
    ni = int(spl[0]); nf = int(spl[1]); ex = float(spl[2])
    assert ni == 2, 'More than 2 initial particles in collision'
    eventsCounter += 1

    if round(ex, 4) == 4.3837:
        fermi=True

    line = file.readline()
    spl = line.split()
    assert int(spl[0]) == 12, 'Initial particle not neutrino'
    v1 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
    line = file.readline()

    for j in range(nf):
        line = file.readline()
        spl = line.split()
        # For electrons produced
        if int(spl[0]) == 11:            
            eng = float(spl[1]) # electron energy
            v2 = np.array([float(spl[2]), float(spl[3]), float(spl[4])])
            th = math.atan2(v2[0], v2[1])
            ph = math.acos(v2[2]/np.linalg.norm(v2))
            rn = rand.random()
            addNLL(th, ph, eng, fermi, rn)

    if eventsCounter == eventsPerSN:
        z = np.where(likelihoods == np.min(likelihoods))
        minY = z[0][0]*math.pi/latRes - math.pi/2
        minX = 2*z[1][0]*math.pi/longRes - math.pi
        minima[snCounter][0] = minX
        minima[snCounter][1] = minY
        print(minX, minY)
        angles = np.append(angles, pointToAngle(minX, minY))
        snCounter += 1
        eventsCounter = 0
        likelihoods = np.zeros((latRes, longRes))
    
    line = file.readline()

np.save('./Skymaps/stdevRot' + str(int(10*fermiCorrect)) + str(int(10*gtCorrect)), minima)
