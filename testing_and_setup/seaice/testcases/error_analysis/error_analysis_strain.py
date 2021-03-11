from netCDF4 import Dataset
import numpy as np
from math import pow, sqrt

#-------------------------------------------------------------------------------

def error_analysis_strain():

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


    fileout = open("error_report_strain.txt","w")


    gridTypes = ["hex"]
    grids = {"hex" :"0082x0094"}
    for gridType in gridTypes:

        print("GridType: ", gridType)
        fileout.write("GridType: %s" %(gridType))


        # grid file
        #filegrid = Dataset("grid_%s_%s.nc" %(gridType,grids[gridType]),"r")
        filegrid = Dataset("grid.nc","r")

        nVertices = len(filegrid.dimensions["nVertices"])
        vertexDegree = len(filegrid.dimensions["vertexDegree"])

        nEdgesOnCell = filegrid.variables["nEdgesOnCell"][:]
        cellsOnVertex = filegrid.variables["cellsOnVertex"][:]
        cellsOnVertex[:] = cellsOnVertex[:] - 1
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
        kiteAreasOnVertex = filegrid.variables["kiteAreasOnVertex"][:]

        filegrid.close()



        # variational
        print("  Variational:")
        fileout.write("  Variational:\n")

        iVertexTest = 8
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

            print("  ",basis)
            fileout.write("  %s\n" %(basis))

            # variational
            #filein = Dataset("./output_%s_%s_%s/output.2000.nc" %(gridType,basis,grids[gridType]),"r")
            filein = Dataset("./output_%s_%s/output.2000.nc" %(gridType,basis),"r")

            basisGradientU = filein.variables["basisGradientU"][:]
            basisGradientV = filein.variables["basisGradientV"][:]

            filein.close()

            # gradient taylor terms
            dfdx = np.zeros((vertexDegree,nTerms)) # ideal: (0,1,0,0,0,0)
            dfdy = np.zeros((vertexDegree,nTerms)) # ideal: (0,0,1,0,0,0)

            for iCellOnVertex in range(0, vertexDegree):
                iCell = cellsOnVertex[iVertexTest, iCellOnVertex]

                # find correct vertex
                for iGradientVertex in range(0, nEdgesOnCell[iCell]):
                    iVertex1 = verticesOnCell[iCell, iGradientVertex]
                    if (iVertex1 == iVertexTest):

                        # loop over basis functions
                        for iBasisVertex in range(0, nEdgesOnCell[iCell]):

                            iVertex = verticesOnCell[iCell, iBasisVertex]
                            dx = xVertex[iVertex] - xVertex[iVertexTest]
                            dy = yVertex[iVertex] - yVertex[iVertexTest]

                            for iTerm in range(0,nTerms):
                                fTerm = coeff[iTerm] * pow(dx,dxPow[iTerm]) * pow(dy,dyPow[iTerm])
                                dfdx[iCellOnVertex,iTerm] += fTerm * basisGradientU[iCell,iGradientVertex,iBasisVertex]
                                dfdy[iCellOnVertex,iTerm] += fTerm * basisGradientV[iCell,iGradientVertex,iBasisVertex]


            for iCellOnVertex in range(0, vertexDegree):
                print("    iCellOnVertex: ",iCellOnVertex)
                print("          # Term       df/dx                  df/dy")
                print("          ---------------------------------------------------------")
                fileout.write("    iCellOnVertex: %i\n" %(iCellOnVertex))
                fileout.write("          # Term       df/dx                  df/dy\n")
                fileout.write("          ---------------------------------------------------------\n")
                for iTerm in range(0,nTerms):
                    print("      %5i %-8s: % 20.15e % 20.15e" %(iTerm, terms[iTerm], dfdx[iCellOnVertex,iTerm], dfdy[iCellOnVertex,iTerm]))
                    fileout.write("      %5i %-8s: % 20.15e % 20.15e\n" %(iTerm, terms[iTerm], dfdx[iCellOnVertex,iTerm], dfdy[iCellOnVertex,iTerm]))

            print("    Avg:")
            fileout.write("    Avg:\n")

            dfdxAvg = np.zeros(nTerms)
            dfdyAvg = np.zeros(nTerms)

            for iTerm in range(0,nTerms):
                denominator = 0.0
                for iCellOnVertex in range(0, vertexDegree):
                    iCell = cellsOnVertex[iVertexTest, iCellOnVertex]
                    #dfdxAvg[iTerm] += dfdx[iCellOnVertex,iTerm]
                    #dfdyAvg[iTerm] += dfdy[iCellOnVertex,iTerm]
                    #denominator += 1.0
                    #dfdxAvg[iTerm] += dfdx[iCellOnVertex,iTerm] * areaCell[iCell]
                    #dfdyAvg[iTerm] += dfdy[iCellOnVertex,iTerm] * areaCell[iCell]
                    #denominator += areaCell[iCell]
                    dfdxAvg[iTerm] += dfdx[iCellOnVertex,iTerm] * kiteAreasOnVertex[iVertexTest,iCellOnVertex]
                    dfdyAvg[iTerm] += dfdy[iCellOnVertex,iTerm] * kiteAreasOnVertex[iVertexTest,iCellOnVertex]
                    denominator += kiteAreasOnVertex[iVertexTest,iCellOnVertex]
                dfdxAvg[iTerm] /= denominator
                dfdyAvg[iTerm] /= denominator

            print("          # Term       df/dx                  df/dy")
            print("          ---------------------------------------------------------")
            fileout.write("          # Term       df/dx                  df/dy\n")
            fileout.write("          ---------------------------------------------------------\n")
            for iTerm in range(0,nTerms):
                print("      %5i %-8s: % 20.15e % 20.15e" %(iTerm, terms[iTerm], dfdxAvg[iTerm], dfdyAvg[iTerm]))
                fileout.write("      %5i %-8s: % 20.15e % 20.15e\n" %(iTerm, terms[iTerm], dfdxAvg[iTerm], dfdyAvg[iTerm]))


        # weak
        print("  Weak:")
        fileout.write("  Weak:\n")
        iCellTest = 0

        # weak
        #filein = Dataset("./output_%s_weak_%s/output.2000.nc" %(gridType,grids[gridType]),"r")
        filein = Dataset("./output_%s_weak/output.2000.nc" %(gridType),"r")

        normalVectorPolygon = filein.variables["normalVectorPolygon"][:]

        filein.close()

        # gradient taylor terms
        dfdx = np.zeros(nTerms) # ideal: (0,1,0,0,0,0)
        dfdy = np.zeros(nTerms) # ideal: (0,0,1,0,0,0)

        for iTerm in range(0,nTerms):

            for iEdgeOnCell in range(0,nEdgesOnCell[iCellTest]):

                # interpolated edge velocity
                iEdge = edgesOnCell[iCellTest,iEdgeOnCell]

                velocityEdge = 0.0

                for iVertexOnEdge in range(0,2):

                    iVertex = verticesOnEdge[iEdge,iVertexOnEdge]

                    dx = xVertex[iVertex] - xCell[iCellTest]
                    dy = yVertex[iVertex] - yCell[iCellTest]
                    fTerm = coeff[iTerm] * pow(dx,dxPow[iTerm]) * pow(dy,dyPow[iTerm])

                    velocityEdge += fTerm

                velocityEdge *= 0.5

                # summation over edges
                dfdx[iTerm] += (velocityEdge * normalVectorPolygon[iCellTest,iEdgeOnCell,0] * dvEdge[iEdge]) / areaCell[iCellTest]
                dfdy[iTerm] += (velocityEdge * normalVectorPolygon[iCellTest,iEdgeOnCell,1] * dvEdge[iEdge]) / areaCell[iCellTest]

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

    error_analysis_strain()
