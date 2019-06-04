import numpy as np
import glob
import uncertainties
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import AutoMinorLocator
from statistics import mean, stdev, variance

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
                errs.append(np.sqrt(variance(data.real)/len(data)))
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
                errs.append(np.sqrt(variance(data.real)/len(data)))
            obsAve.append(vals)
            obsErr.append(errs)
    obsAve = ReshapeData([x for _,x in sorted(zip(varSteps, obsAve))])
    obsErr = ReshapeData([x for _,x in sorted(zip(varSteps, obsErr))])
    return header, obsAve, obsErr 


def ReadParamsDat(varparFile):
    """
    Read in an undefined number of variational parameters from a data file
    along with the names of each variational parameter.
    """
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


def PlotVarParams(head, data, defaultWidth, fs, ls, ms, x, y,
                  figTitle):
    """
    Plot the variational paramters in file.
    """
    nPar = len(data)
    minorLocator=AutoMinorLocator()
    nSteps = np.arange(0, len(data[0]))
    width = defaultWidth    # width of plots is set by the user
    height = int(np.ceil(nPar / width))  # height is set by default
    # create subplot
    fig, ax = plt.subplots(nrows=height, ncols=width)
    fig.set_size_inches(x, y)
    fig.suptitle(figTitle, fontsize=fs, y=1.08)
    gs1=gridspec.GridSpec(height, width)
    for i in range(0, height):
        for j in range(0, width):
            currPar = (i*(width)) + j
            # check that we haven't exceeded number of plots
            if currPar >= nPar:
                ax[i][j].axis('off')
            else:
                ax[i][j].plot(nSteps, data[currPar], '-o', markersize=ms)
                ax[i][j].set_xlabel("optimization step", fontsize=ls)
                ax[i][j].set_ylabel(head[currPar], fontsize=ls)
                ax[i][j].xaxis.set_minor_locator(minorLocator)
                ax[i][j].tick_params(axis='both', which='both', width=1)
                ax[i][j].tick_params(axis='both', which='major', labelsize=ls)
                ax[i][j].grid(which='both', axis='both')
                ax[i][j].grid(which='major', axis='both', linestyle='-',
                              linewidth=0.3)
                ax[i][j].grid(which='minor', axis='both', linestyle='--',
                              linewidth=0.1)
    fig.tight_layout()
    return fig, ax


def AverageEnergies(energyArray, binSize):
    binNum=int(len(energyArray)/binSize)
    energyArray=np.array(energyArray, dtype=float)
    meanEnergies=[]
    for i in range(0, binNum):
        low=i*binSize
        hi=(i+1) * binSize
        meanEnergies.append(mean(energyArray[low:hi]))
    return meanEnergies

