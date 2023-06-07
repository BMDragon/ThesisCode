import subprocess

repeat = 100

for i in range(repeat):
    subprocess.run('./marley ./stdevData.js', shell=True)

