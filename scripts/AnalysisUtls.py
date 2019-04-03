import numpy as np
import glob
import uncertainties
from statistics import mean, stdev

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


def ReadMeasurements(measFile, varFileType):
    head = []
    datArr = []
    with open(measFile, 'r') as mF:
        # if this data was extracted during a variational step, its step
        # identifier is the first number in the file        
        if varFileType:
            vstep = mF.readline().strip().split()
            vstep = int(vstep[0])
            datArr.append(vstep)
            # read file names
            head = mF.readline().strip().split()
            head = [str(x) for x in head]
            next(mF)
            values = []
            for line in mF:
                rawLine=line.strip().split()
                datLine=[]
                for ele in rawLine:
                    ele=ele.strip("(").strip(")").split(",")
                    if len(ele)==1:
                        datLine.append(float(ele[0]))
                    else:
                        datLine.append(float(ele[0])+float(ele[1])*1.0j) 
                values.append(datLine)
            datArr.append(values)
        else:
            # read file names
            head = mF.readline().strip().split()
            head = [str(x) for x in head]
            next(mF)
            for line in mF:
                rawLine=line.strip().split()
                datLine=[]
                for ele in rawLine:
                    ele=ele.strip("(").strip(")").split(",")
                    if len(ele)==1:
                        datLine.append(float(ele[0]))
                    else:
                        datLine.append(float(ele[0])+float(ele[1])*1.0j) 
                datArr.append(datLine)
    return head, datArr


def ProcessDataFiles(namePattern, varFileType):
    dataFiles = glob.glob(namePattern)
    varSteps = []
    header = []
    obsAve = []
    obsErr = []
    if varFileType:
        for dataFile in dataFiles:
            cHeader, cData = ReadMeasurements(dataFile, varFileType)
            if not header:
                header=cHeader
            cData[1] = ReshapeData(cData[1])
            varSteps.append(cData[0])
            vals = []
            errs = []
            for data in cData[1]:
                vals.append(mean(data.real))
                errs.append(stdev(data.real))
            obsAve.append(vals)
            obsErr.append(errs)
    else:
        for dataFile in dataFiles:
            cHeader, cData = ReadMeasurements(dataFile, varFileType)
            if not header:
                header=cHeader
            cData = ReshapeData(cData)
            vals = []
            errs = []
            for data in cData:
                vals.append(mean(data.real))
                errs.append(stdev(data.real))
            obsAve.append(vals)
            obsErr.append(errs)
    obsAve = ReshapeData([x for _,x in sorted(zip(varSteps, obsAve))])
    obsErr = ReshapeData([x for _,x in sorted(zip(varSteps, obsErr))])
    return header, obsAve, obsErr 

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

