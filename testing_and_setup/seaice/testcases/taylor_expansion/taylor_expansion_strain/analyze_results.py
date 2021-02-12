from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import math
import sys

#gridTypes = ["hex","quad"]
gridTypes = ["hex"]
#operatorMethods = ["wachspress","pwl","weak"]
operatorMethods = ["wachspress","pwl"]
#velocityTypes = ["linear","quadratic","cubic","quartic"]
velocityTypes = ["quadratic"]
velocityTypeOrders = [1,2,3,4]

lineColors = {"wachspress":"black","pwl":"grey","weak":"red"}

labels = {"wachspress":"Wachs.","pwl":"PWL","weak":"Weak"}
lineStyle = {"hex":"solid","quad":"dashed"}
markerStyle = {"hex":"v","quad":"^"}

#-------------------------------------------------------------------------------

def analyze_results_var(axes, operatorMethod, gridType):

    # grid
    filenamein = "grid_var_%s.nc" %(gridType)

    filein = Dataset(filenamein,"r")

    iVertexCenter = filein.iVertexCenter

    vertexDegree = len(filein.dimensions["vertexDegree"])
    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1
    verticesOnCell = filein.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    filein.close()

    strain11varAvgDiff = []
    strain22varAvgDiff = []
    strain12varAvgDiff = []

    for velocityType in velocityTypes:

        print("    velocityType: ", velocityType)

        # ic
        filenamein = "ic_var_%s_%s.nc" %(gridType, velocityType)

        filein = Dataset(filenamein,"r")

        strain11VertexAnalytical = filein.variables["strain11VertexAnalytical"][:]
        strain22VertexAnalytical = filein.variables["strain22VertexAnalytical"][:]
        strain12VertexAnalytical = filein.variables["strain12VertexAnalytical"][:]

        filein.close()

        # sim
        filenamein = "./output_%s_%s_%s/output.2000.nc" %(operatorMethod, gridType, velocityType)

        filein = Dataset(filenamein,"r")

        strain11var = filein.variables["strain11var"][0,:,:]
        strain22var = filein.variables["strain22var"][0,:,:]
        strain12var = filein.variables["strain12var"][0,:,:]

        filein.close()

        strain11varAvg = 0.0
        strain22varAvg = 0.0
        strain12varAvg = 0.0

        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertexCenter,iCellOnVertex]
            for iVertexOnCell in range(nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                if (iVertex == iVertexCenter):
                    strain11varAvg += strain11var[iCell,iVertexOnCell]
                    strain22varAvg += strain22var[iCell,iVertexOnCell]
                    strain12varAvg += strain12var[iCell,iVertexOnCell]
                    print("         {0:f} {1:f} {2:f}".format(strain11var[iCell,iVertexOnCell] / float(vertexDegree), strain22var[iCell,iVertexOnCell] / float(vertexDegree), strain12var[iCell,iVertexOnCell] / float(vertexDegree)))

        strain11varAvg /= float(vertexDegree)
        strain22varAvg /= float(vertexDegree)
        strain12varAvg /= float(vertexDegree)

        strain11varAvgDiff.append(math.fabs(strain11varAvg - strain11VertexAnalytical[iVertexCenter]))
        strain22varAvgDiff.append(math.fabs(strain22varAvg - strain22VertexAnalytical[iVertexCenter]))
        strain12varAvgDiff.append(math.fabs(strain12varAvg - strain12VertexAnalytical[iVertexCenter]))

        print("      {0:e} {1:e} {2:e}".format(strain11varAvg, strain11VertexAnalytical[iVertexCenter], strain11varAvg-strain11VertexAnalytical[iVertexCenter]))
        print("      {0:e} {1:e} {2:e}".format(strain22varAvg, strain22VertexAnalytical[iVertexCenter], strain22varAvg-strain22VertexAnalytical[iVertexCenter]))
        print("      {0:e} {1:e} {2:e}".format(strain12varAvg, strain12VertexAnalytical[iVertexCenter], strain12varAvg-strain12VertexAnalytical[iVertexCenter]))

    #axes[0].semilogy(velocityTypeOrders, strain11varAvgDiff, label="%s %s" %(labels[operatorMethod],gridType), linestyle=lineStyle[gridType], color=lineColors[operatorMethod], marker=markerStyle[gridType])
    #axes[1].semilogy(velocityTypeOrders, strain22varAvgDiff, label="%s %s" %(labels[operatorMethod],gridType), linestyle=lineStyle[gridType], color=lineColors[operatorMethod], marker=markerStyle[gridType])
    #axes[2].semilogy(velocityTypeOrders, strain12varAvgDiff, label="%s %s" %(labels[operatorMethod],gridType), linestyle=lineStyle[gridType], color=lineColors[operatorMethod], marker=markerStyle[gridType])

#-------------------------------------------------------------------------------

def analyze_results_weak(axes, gridType):

    strain11weakDiff = []
    strain22weakDiff = []
    strain12weakDiff = []

    for velocityType in velocityTypes:

        print("    velocityType: ", velocityType)

        # ic
        filenamein = "ic_weak_%s_%s.nc" %(gridType, velocityType)

        filein = Dataset(filenamein,"r")

        strain11CellAnalytical = filein.variables["strain11CellAnalytical"][:]
        strain22CellAnalytical = filein.variables["strain22CellAnalytical"][:]
        strain12CellAnalytical = filein.variables["strain12CellAnalytical"][:]

        filein.close()

        # sim
        filenamein = "./output_weak_%s_%s/output.2000.nc" %(gridType, velocityType)

        filein = Dataset(filenamein,"r")

        strain11weak = filein.variables["strain11weak"][0,:][0]
        strain22weak = filein.variables["strain22weak"][0,:][0]
        strain12weak = filein.variables["strain12weak"][0,:][0]

        filein.close()

        strain11weakDiff.append(math.fabs(strain11weak - strain11CellAnalytical[0]))
        strain22weakDiff.append(math.fabs(strain22weak - strain22CellAnalytical[0]))
        strain12weakDiff.append(math.fabs(strain12weak - strain12CellAnalytical[0]))

        print("      {0:e} {1:e} {2:e}".format(strain11weak, strain11CellAnalytical[0], strain11weak-strain11CellAnalytical[0]))
        print("      {0:e} {1:e} {2:e}".format(strain22weak, strain22CellAnalytical[0], strain22weak-strain22CellAnalytical[0]))
        print("      {0:e} {1:e} {2:e}".format(strain12weak, strain12CellAnalytical[0], strain12weak-strain12CellAnalytical[0]))

    #axes[0].semilogy(velocityTypeOrders, strain11weakDiff, label="%s %s" %("Weak",gridType), linestyle=lineStyle[gridType], color=lineColors["weak"], marker=markerStyle[gridType])
    #axes[1].semilogy(velocityTypeOrders, strain22weakDiff, label="%s %s" %("Weak",gridType), linestyle=lineStyle[gridType], color=lineColors["weak"], marker=markerStyle[gridType])
    #axes[2].semilogy(velocityTypeOrders, strain12weakDiff, label="%s %s" %("Weak",gridType), linestyle=lineStyle[gridType], color=lineColors["weak"], marker=markerStyle[gridType])

#-------------------------------------------------------------------------------

def analyze_results():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(3, 1, figsize=(5,7))

    for gridType in gridTypes:

        print("Grid type: ", gridType)

        for operatorMethod in operatorMethods:

            print("  Operator Method: ", operatorMethod)

            if (operatorMethod in ["wachspress","pwl"]):
                analyze_results_var(axes, operatorMethod, gridType)
            else:
                analyze_results_weak(axes, gridType)

    #axes[0].legend()

    #axes[0].set_title(r'$\epsilon_{11}$')
    #axes[1].set_title(r'$\epsilon_{22}$')
    #axes[2].set_title(r'$\epsilon_{12}$')

    #plt.show()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    analyze_results()
