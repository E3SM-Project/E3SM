from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def aggregate_tracers():

    advectionMethods = ["IR","upwind"]

    icTypes = ["cosine_bell","slotted_cylinder"]

    gridSizes = [2562, 10242, 40962, 163842]

    for advectionMethod in advectionMethods:

        print("Advection method: ", advectionMethod)

        for icType in icTypes:

            print("  IC type: ", icType)

            for gridSize in gridSizes:

                print("    Gridsize: ", gridSize)

                filename = "./output_%s_%s_%i/output.2000.nc" %(advectionMethod,icType,gridSize)
                filein = Dataset(filename,"a")

                nTimes = len(filein.dimensions["Time"])
                nCells = len(filein.dimensions["nCells"])
                nCategories = len(filein.dimensions["nCategories"])

                iceAreaCategory = filein.variables["iceAreaCategory"][:]

                iceAreaCell = np.sum(iceAreaCategory,axis=(2,3))

                try:
                    var = filein.createVariable("iceAreaCell","d",dimensions=["nCells"])
                except:
                    var = filein.variables["iceAreaCell"]
                var[:] = iceAreaCell[:]

                filein.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    aggregate_tracers()
