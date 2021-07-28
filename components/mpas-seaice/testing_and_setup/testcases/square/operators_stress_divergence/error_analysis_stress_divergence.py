from netCDF4 import Dataset
import numpy as np

#-----------------------------------------------------------------------

def error_analysis_stress_divergence():

    # taylor series terms:
    # 0: f(x0,y0)
    # 1: df/dx(x0,y0) \Delta x
    # 2: df/dy(x0,y0) \Delta y
    # 3: 1/2 d2f/dx2(x0,y0) \Delta x^2
    # 4: 1/2 d2f/dy2(x0,y0) \Delta y^2
    # 5: d2f/dxdy(x0,y0) \Delta x \Delta y
    nTerms = 6
    dxPow = [0,1,0,2,1,0]
    dyPow = [0,0,1,0,1,2]
    coeff = [1.0,1.0,1.0,0.5,1.0,0.5]
    terms = ["f","df/dx","df/dy","d2f/dx2","d2f/dxdy","d2f/dy2"]


    fileout = open("error_report_stress_divergence.txt","w")


    gridTypes = ["hex","quad"]
    grids = {"hex" :"0082x0094",
             "quad":"0080x0080"}
    for gridType in gridTypes:

        print("GridType: ", gridType)
        fileout.write("GridType: %s" %(gridType))


        # grid file
        filegrid = Dataset("grid_%s_%s.nc" %(gridType,grids[gridType]),"r")

        nVertices = len(filegrid.dimensions["nVertices"])
        vertexDegree = len(filegrid.dimensions["vertexDegree"])

        nEdgesOnCell = filegrid.variables["nEdgesOnCell"][:]
        cellsOnVertex = filegrid.variables["cellsOnVertex"][:]
        cellsOnVertex[:] = cellsOnVertex[:] - 1
        cellsOnEdge = filegrid.variables["cellsOnEdge"][:]
        cellsOnEdge[:] = cellsOnEdge[:] - 1
        verticesOnCell = filegrid.variables["verticesOnCell"][:]
        verticesOnCell[:] = verticesOnCell[:] - 1
        verticesOnEdge = filegrid.variables["verticesOnEdge"][:]
        verticesOnEdge[:] = verticesOnEdge[:] - 1
        edgesOnVertex = filegrid.variables["edgesOnVertex"][:]
        edgesOnVertex[:] = edgesOnVertex[:] - 1
        edgesOnCell = filegrid.variables["edgesOnCell"][:]
        edgesOnCell[:] = edgesOnCell[:] - 1

        xVertex = filegrid.variables["xVertex"][:]
        yVertex = filegrid.variables["yVertex"][:]
        xCell = filegrid.variables["xCell"][:]
        yCell = filegrid.variables["yCell"][:]
        dcEdges = filegrid.variables["dcEdge"][:]
        dvEdge = filegrid.variables["dvEdge"][:]
        areaCell = filegrid.variables["areaCell"][:]
        areaTriangle = filegrid.variables["areaTriangle"][:]

        filegrid.close()


        # variational
        print("  Variational:")
        fileout.write("  Variational:\n")

        iVertexTest = 1123
        print("  iVertexTest: ",iVertexTest)
        fileout.write("  iVertexTest: %i\n" %(iVertexTest))
        for iCellOnVertex in range(0, vertexDegree):
            print("  iCell: ", iCellOnVertex, cellsOnVertex[iVertexTest,iCellOnVertex])
            fileout.write("  iCell: %i %i\n" %(iCellOnVertex, cellsOnVertex[iVertexTest,iCellOnVertex]))

        iEdge = edgesOnVertex[iVertexTest,0]
        dcEdge = dcEdges[iEdge]
        print("  dcEdge: ", dcEdge)
        fileout.write("  dcEdge: %f\n" %(dcEdge))


        basises = ["wachspress","pwl"]
        for basis in basises:

            print("   ",basis)
            fileout.write("   %s\n" %(basis))

            # variational
            filein = Dataset("./output_%s_%s_%s/output.2000.nc" %(gridType,basis,grids[gridType]),"r")

            basisIntegralsU = filein.variables["basisIntegralsU"][:]
            basisIntegralsV = filein.variables["basisIntegralsV"][:]
            cellVerticesAtVertex = filein.variables["cellVerticesAtVertex"][:]
            cellVerticesAtVertex[:] = cellVerticesAtVertex[:] - 1

            filein.close()

            # gradient taylor terms
            dfdx = np.zeros(nTerms) # ideal: (0,1,0,0,0,0)
            dfdy = np.zeros(nTerms) # ideal: (0,0,1,0,0,0)

            for iTerm in range(0,nTerms):

                # loop over surrounding cells
                for iSurroundingCell in range(0,vertexDegree):

                    # get the cell number of this cell
                    iCell = cellsOnVertex[iVertexTest,iSurroundingCell]

                    # get the vertexOnCell number of the iVertex velocity point from cell iCell
                    iVelocityVertex = cellVerticesAtVertex[iVertexTest,iSurroundingCell]

                    # loop over the vertices of the surrounding cell
                    for iStressVertex in range(0,nEdgesOnCell[iCell]):

                        iVertex = verticesOnCell[iCell,iStressVertex]

                        dx = xVertex[iVertex] - xVertex[iVertexTest]
                        dy = yVertex[iVertex] - yVertex[iVertexTest]

                        fTerm = coeff[iTerm] * pow(dx,dxPow[iTerm]) * pow(dy,dyPow[iTerm])

                        dfdx[iTerm] -= (fTerm * basisIntegralsU[iCell,iVelocityVertex,iStressVertex]) / areaTriangle[iVertexTest]
                        dfdy[iTerm] -= (fTerm * basisIntegralsV[iCell,iVelocityVertex,iStressVertex]) / areaTriangle[iVertexTest]

            dfdxAvg = np.mean(dfdx,axis=0)
            dfdyAvg = np.mean(dfdy,axis=0)
            print("          # Term       df/dx                  df/dy")
            print("          ---------------------------------------------------------")
            fileout.write("          # Term       df/dx                  df/dy\n")
            fileout.write("          ---------------------------------------------------------\n")
            for iTerm in range(0,nTerms):
                print("      %5i %-8s: % 20.15e % 20.15e" %(iTerm, terms[iTerm], dfdx[iTerm], dfdy[iTerm]))
                fileout.write("      %5i %-8s: % 20.15e % 20.15e\n" %(iTerm, terms[iTerm], dfdx[iTerm], dfdy[iTerm]))

        # weak
        print("  Weak:")
        fileout.write("  Weak:\n")
        iVertexTest = 1123

        # weak
        filein = Dataset("./output_%s_weak_%s/output.2000.nc" %(gridType,grids[gridType]),"r")

        normalVectorTriangle = filein.variables["normalVectorTriangle"][:]

        filein.close()


        dfdx = np.zeros(nTerms) # ideal: (0,1,0,0,0,0)
        dfdy = np.zeros(nTerms) # ideal: (0,0,1,0,0,0)

        for iTerm in range(0,nTerms):

            for iVertexDegree in range(0,vertexDegree):

                # interpolated edge velocity
                iEdge = edgesOnVertex[iVertexTest,iVertexDegree]

                stressEdge = 0.0

                for iCellOnEdge in range(0,2):

                    iCell = cellsOnEdge[iEdge,iCellOnEdge]

                    dx = xCell[iCell] - xVertex[iVertexTest]
                    dy = yCell[iCell] - yVertex[iVertexTest]
                    fTerm = coeff[iTerm] * pow(dx,dxPow[iTerm]) * pow(dy,dyPow[iTerm])

                    stressEdge += fTerm

                stressEdge *= 0.5

                dfdx[iTerm] += (stressEdge * normalVectorTriangle[iVertexTest,iVertexDegree,0] * dcEdges[iEdge]) / areaTriangle[iVertexTest]
                dfdy[iTerm] += (stressEdge * normalVectorTriangle[iVertexTest,iVertexDegree,1] * dcEdges[iEdge]) / areaTriangle[iVertexTest]



        print("          # Term       df/dx                  df/dy")
        print("          ---------------------------------------------------------")
        fileout.write("          # Term       df/dx                  df/dy\n")
        fileout.write("          ---------------------------------------------------------\n")
        for iTerm in range(0,nTerms):
            print("      %5i %-8s: % 20.15e % 20.15e" %(iTerm, terms[iTerm], dfdx[iTerm], dfdy[iTerm]))
            fileout.write("      %5i %-8s: % 20.15e % 20.15e\n" %(iTerm, terms[iTerm], dfdx[iTerm], dfdy[iTerm]))


    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    error_analysis_stress_divergence()
