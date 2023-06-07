import numpy as np
import matplotlib.pyplot as plt

l_data = np.load("./build/l_data.npy", allow_pickle=True)
data = [dict(point) for point in l_data]

ar40 = np.zeros(200)

for i in range(len(data)):
    ar40 += data[i]['nue_Ar40']

with open('./build/histEnergies.txt', 'w') as file:
    eng = data[0]['Energy']
    file.write('[')
    for j in range(len(eng)):
        file.write(" " + str(eng[j]*1000) + ",")
    file.write("]\n[")
    for k in range(len(ar40)):
        file.write(" " + str(ar40[k]) + ",")
    file.write(']')
