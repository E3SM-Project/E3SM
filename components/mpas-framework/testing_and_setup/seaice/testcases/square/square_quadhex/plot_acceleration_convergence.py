from netCDF4 import Dataset
import math
import numpy as np
import matplotlib.pyplot as plt

gridSize = 640
#gridTypes = ["quad","hex"]
gridTypes = ["quad"]
operatorMethods = ["wachspress","pwl","weak"]
#subcycleNumbers = [120,240,480,960,1920,3840,7680,15360,30720]
subcycleNumbers = [120,240,480,960,1920,3840,7680]

linestyles = {"quad":"solid",
              "hex":"dashed"}
colors = {"wachspress":"red",
          "pwl":"blue",
          "weak":"green"}

for gridType in gridTypes:

    linestyle = linestyles[gridType]

    for operatorMethod in operatorMethods:

        color = colors[operatorMethod]

        for subcycleNumber in subcycleNumbers:

            filenameIn = "./output_%s_%s_%i_%i/output.2000.nc" %(gridType, operatorMethod, gridSize, subcycleNumber)

            filein = Dataset(filenameIn,"r")

            nTimes = len(filein.dimensions["Time"])

            uAcceleration = filein.variables["uAcceleration"][:]
            uVelocity = filein.variables["uVelocity"][:]

            filein.close()

            rmss = []

            for iTime in range(0,nTimes):

                rms = math.sqrt(np.sum(np.power(uAcceleration[iTime,:],2)) / float(np.shape(uAcceleration[iTime,:])[0]))

                rmss.append(rms)

            plt.semilogy(rmss, label="%s %s %i" %(gridType, operatorMethod, subcycleNumber), linestyle=linestyle, color=color)


plt.xlabel("Time")
plt.ylabel("RMS acceleration")
#plt.legend()
plt.savefig("acceleration_convergence.png",dpi=300)
