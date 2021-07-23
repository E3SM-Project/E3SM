from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math

fig, axes = plt.subplots()


subcycleNumbers = [120,240,480,960,1920,3840,7680]
operatorMethods = ["wachspress","pwl","weak"]

Lx = 1280000.0


iMethod = 0
for operatorMethod in operatorMethods:
    unegs = []
    for subcycleNumber in subcycleNumbers:

        # data in
        filenameIn = "./output_hex_%s_%i/output.2000.nc" %(operatorMethod,subcycleNumber)
        filein = Dataset(filenameIn, "r")

        nCells = len(filein.dimensions["nCells"])
        nVertices = len(filein.dimensions["nVertices"])
        vertexDegree = len(filein.dimensions["vertexDegree"])
        nTimes = len(filein.dimensions["Time"])

        cellsOnVertex = filein.variables["cellsOnVertex"][:]
        cellsOnVertex -= 1

        xVertex = filein.variables["xVertex"][:]
        yVertex = filein.variables["yVertex"][:]

        xCell = filein.variables["xCell"][:]
        yCell = filein.variables["yCell"][:]

        uVelocity = filein.variables["uVelocity"][-1,:]
        vVelocity = filein.variables["vVelocity"][-1,:]

        uVelocities = filein.variables["uVelocity"][:,:]

        filein.close()
        
        xmin = np.amin(xVertex)
        xmax = np.amax(xVertex)
        ymin = np.amin(yVertex)
        ymax = np.amax(yVertex)

        uneg = 0.0
        for iVertex in range(0,nVertices):
            if (math.fabs(yVertex[iVertex] - 508068.236886871) < 1e-8 and
                xVertex[iVertex] >= 0.0 and xVertex[iVertex] <= Lx):
                u = uVelocities[-1,iVertex]
                if (u < 0.0):
                    uneg += u
        unegs.append(math.fabs(uneg))

    axes.loglog(subcycleNumbers, unegs, label=operatorMethod)


axes.legend()
axes.set_xlabel("subcycle number")
axes.set_ylabel("total negative velocity")
plt.savefig("1D_velocity_negu.png",dpi=300)
