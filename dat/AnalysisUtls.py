import numpy as np
from statistics import mean

def ReadEnergies(energyFile):
    lines = []
    with open(energyFile, 'r') as eF:
        lines = eF.readlines()
        lines = [int(x.strip()) for x in lines]
    return np.array(lines)

def AverageEnergies(energyArray, binSize):
    binNum=int(len(energyArray)/binSize)
    energyArray=np.array(energyArray, dtype=float)
    meanEnergies=[]
    for i in range(0, binNum):
        low=i*binSize
        hi=(i+1) * binSize
        meanEnergies.append(mean(energyArray[low:hi]))
    return meanEnergies

