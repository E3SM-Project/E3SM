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

for gridType in gridTypes:

    for operatorMethod in operatorMethods:

        rmss = []

        for subcycleNumber in subcycleNumbers:

            filenameIn = "./output_%s_%s_%i_%i/output.2000.nc" %(gridType, operatorMethod, gridSize, subcycleNumber)

            filein = Dataset(filenameIn,"r")

            uAcceleration = filein.variables["uAcceleration"][-1,:]
            uVelocity = filein.variables["uVelocity"][-1,:]

            filein.close()

            rms = math.sqrt(np.sum(np.power(uAcceleration,2)) / float(np.shape(uAcceleration)[0]))

            rmss.append(rms)

        plt.loglog(subcycleNumbers, rmss, label="%s %s" %(gridType,operatorMethod))

plt.xlabel("subcycle number per time step")
plt.ylabel("RMS acceleration")
plt.legend()
plt.savefig("acceleration_scaling.png",dpi=300)
