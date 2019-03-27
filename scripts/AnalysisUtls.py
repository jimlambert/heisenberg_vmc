import numpy as np
from statistics import mean

def ReshapeData(array):
    reshapedData = []
    for i in range(0, len(array[0])):
        ccol = []
        for j in range(0, len(array)):
            ccol.append(array[j][i])
        reshapedData.append(np.array(ccol))
    return np.array(reshapedData)


def ReadEnergies(energyFile):
    lines = []
    with open(energyFile, 'r') as eF:
        lines = eF.readlines()
        lines = [int(x.strip()) for x in lines]
    return np.array(lines)


def ReadMeasurements(measFile):
    head = []
    data = []
    with open(measFile, 'r') as mF:
        head = mF.readline().strip().split()
        head = [str(x) for x in head]
        next(mF)
        for line in mF:
            dline = line.strip().split()
            dline[0] = float(dline[0])
            data.append(dline)
    return head, data


def ReadParamsDat(varparFile):
    head = []
    data = []
    with open(varparFile, 'r') as vF:
        head = vF.readline().strip().split()
        head = [str(x) for x in head]
        next(vF)
        for line in vF:
            cline = line.strip().split()
            cline = [float(x) for x in cline]
            data.append(np.array(cline))
    data=ReshapeData(data)
    return head, data


def AverageEnergies(energyArray, binSize):
    binNum=int(len(energyArray)/binSize)
    energyArray=np.array(energyArray, dtype=float)
    meanEnergies=[]
    for i in range(0, binNum):
        low=i*binSize
        hi=(i+1) * binSize
        meanEnergies.append(mean(energyArray[low:hi]))
    return meanEnergies

