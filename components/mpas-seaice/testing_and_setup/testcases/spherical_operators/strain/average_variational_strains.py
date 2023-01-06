from netCDF4 import Dataset
import numpy as np

#-------------------------------------------------------------------------------

def average_variational_strains():

    gridSizes = [2562, 10242, 40962, 163842]

    operatorMethods = ["wachspress","pwl","weakwachs"]

    for operatorMethod in operatorMethods:

        print("Operator Method: ", operatorMethod)

        for gridSize in gridSizes:

            print("  Gridsize: ", gridSize)
    
            filenameModify = "./output_%s_%i/output.2000.nc" %(operatorMethod, gridSize)
            fileModify = Dataset(filenameModify,"a")

            nVertices = len(fileModify.dimensions["nVertices"])
            nCells = len(fileModify.dimensions["nCells"])
            vertexDegree = len(fileModify.dimensions["vertexDegree"])
            nTimes = len(fileModify.dimensions["Time"])

            cellsOnVertex = fileModify.variables["cellsOnVertex"][:]
            cellVerticesAtVertex = fileModify.variables["cellVerticesAtVertex"][:]

            cellsOnVertex[:] -= 1
            cellVerticesAtVertex[:] -= 1

            strain11var = fileModify.variables["strain11var"][:]
            strain22var = fileModify.variables["strain22var"][:]
            strain12var = fileModify.variables["strain12var"][:]

            strain11varAvg = np.zeros((nTimes,nVertices))
            strain22varAvg = np.zeros((nTimes,nVertices))
            strain12varAvg = np.zeros((nTimes,nVertices))

            for iTime in range(0, nTimes):

                for iVertex in range(0, nVertices):

                    strain11avg = 0.0
                    strain22avg = 0.0
                    strain12avg = 0.0
                    nCellsSum = 0

                    for iVertexDegree in range(0,vertexDegree):

                        iCell = cellsOnVertex[iVertex,iVertexDegree]

                        if (iCell <= nCells-1):

                            iVertexOnCell = cellVerticesAtVertex[iVertex,iVertexDegree]

                            strain11avg = strain11avg + strain11var[iTime,iCell,iVertexOnCell]
                            strain22avg = strain22avg + strain22var[iTime,iCell,iVertexOnCell]
                            strain12avg = strain12avg + strain12var[iTime,iCell,iVertexOnCell]
                            nCellsSum = nCellsSum + 1

                    strain11varAvg[iTime,iVertex] = strain11avg / float(nCellsSum)
                    strain22varAvg[iTime,iVertex] = strain22avg / float(nCellsSum)
                    strain12varAvg[iTime,iVertex] = strain12avg / float(nCellsSum)

                    #print(iTime,iVertex,strain11varAvg[iTime,iVertex],strain11avg,float(nCellsSum))

            try:
                var = fileModify.createVariable("strain11varAvg","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain11varAvg"]
            var[:] = strain11varAvg[:]

            try:
                var = fileModify.createVariable("strain22varAvg","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain22varAvg"]
            var[:] = strain22varAvg[:]

            try:
                var = fileModify.createVariable("strain12varAvg","d",dimensions=["Time","nVertices"])
            except:
                var = fileModify.variables["strain12varAvg"]
            var[:] = strain12varAvg[:]

            fileModify.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    average_variational_strains()
