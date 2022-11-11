from netCDF4 import Dataset
import numpy as np
import math
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def create_boundary_info_grid(nx, plot=False):

    gridname = "grid_%4.4i.nc" %(nx)

    print(gridname)

    gridFile = Dataset(gridname,"r")

    Lx = gridFile.Lx
    Ly = gridFile.Ly

    #nx = gridFile.nx
    ny = gridFile.ny

    dcx = gridFile.dcx
    dcy = gridFile.dcy

    a = dcx / math.sqrt(3.0) # side of hexagon

    nVertices = len(gridFile.dimensions["nVertices"])

    xVertex = gridFile.variables["xVertex"][:]
    yVertex = gridFile.variables["yVertex"][:]

    nCells = len(gridFile.dimensions["nCells"])

    xCell = gridFile.variables["xCell"][:]
    yCell = gridFile.variables["yCell"][:]

    verticesOnCell = gridFile.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    nEdgesOnCell = gridFile.variables["nEdgesOnCell"][:]

    gridFile.close()




    pairType = np.zeros(nVertices,dtype="i")
    pairGlobalIndex = np.zeros(nVertices,dtype="i")
    solveVelocity = np.zeros(nVertices,dtype="i")
    solveStress = np.zeros(nCells,dtype="i")

    # pairType
    # 0: Interior
    # 1: On boundary / Exterior (both set to zero)
    # 2: Reflected
    # 3: Corner
    VELOCITY_BOUNDARY_NONE = 0
    VELOCITY_BOUNDARY_PERIODIC = 1
    VELOCITY_BOUNDARY_REVERSE = 2
    VELOCITY_BOUNDARY_ZERO = 3
    
    e = 1.0

    # set vertices on or out of boundary to type 1
    for iVertex in range(0,nVertices):
        if (xVertex[iVertex] <= 0.0 + e or
            xVertex[iVertex] >= Lx  - e or
            yVertex[iVertex] <= 0.0 + e or
            yVertex[iVertex] >= Ly  - e): pairType[iVertex] = VELOCITY_BOUNDARY_ZERO

    # now find cells that are on the boundaries but not corners
    # bottom
    #for iCell in range(0,nCells):
    #    if (xCell[iCell] > 0.0 - e and
    #        xCell[iCell] < Lx  + e and
    #        yCell[iCell] > 0.0 - e and
    #        yCell[iCell] < 0.0 + e):

    #        for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
    #            iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

    #            for iVertexOnCell2 in range(0,nEdgesOnCell[iCell]):
    #                iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

    #                if (math.fabs(xVertex[iVertex1] - xVertex[iVertex2]) < e and
    #                    yVertex[iVertex1] < yVertex[iVertex2] - e):

    #                    #pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
    #                    pairType[iVertex1] = VELOCITY_BOUNDARY_PERIODIC
    #                    pairGlobalIndex[iVertex1] = iVertex2

    # bottom
    for iCell in range(0,nCells):
        if (xCell[iCell] > 0.0 - e and
            xCell[iCell] < Lx  + e and
            yCell[iCell] > 0.0 - e and
            yCell[iCell] < 0.0 + e):

            for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

                if (yVertex[iVertex1] > yCell[iCell]):

                    for iVertex2 in range(0,nVertices):

                        if (math.fabs(xVertex[iVertex1] - xVertex[iVertex2]) < e and
                            math.fabs(yVertex[iVertex1] - yVertex[iVertex2] + 2.0*dcy) < e):

                            #pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
                            pairType[iVertex1] = VELOCITY_BOUNDARY_PERIODIC
                            pairGlobalIndex[iVertex1] = iVertex2

    # top
    #for iCell in range(0,nCells):
    #    if (xCell[iCell] > 0.0 - e and
    #        xCell[iCell] < Lx  + e and
    #        yCell[iCell] > Ly  - e and
    #        yCell[iCell] < Ly  + e):

    #        for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
    #            iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

    #            for iVertexOnCell2 in range(0,nEdgesOnCell[iCell]):
    #                iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

    #                if (math.fabs(xVertex[iVertex1] - xVertex[iVertex2]) < e and
    #                    yVertex[iVertex1] > yVertex[iVertex2] + e):

    #                    #pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
    #                    pairType[iVertex1] = VELOCITY_BOUNDARY_PERIODIC
    #                    pairGlobalIndex[iVertex1] = iVertex2

    # top
    for iCell in range(0,nCells):
        if (xCell[iCell] > 0.0 - e and
            xCell[iCell] < Lx  + e and
            yCell[iCell] > Ly  - e and
            yCell[iCell] < Ly  + e):

            for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

                if (yVertex[iVertex1] < yCell[iCell]):

                    for iVertex2 in range(0,nVertices):

                        if (math.fabs(xVertex[iVertex1] - xVertex[iVertex2]) < e and
                            math.fabs(yVertex[iVertex1] - yVertex[iVertex2] - 2.0*dcy) < e):

                            #pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
                            pairType[iVertex1] = VELOCITY_BOUNDARY_PERIODIC
                            pairGlobalIndex[iVertex1] = iVertex2

    # left
    for iCell in range(0,nCells):
        if (yCell[iCell] > 0.0 - e and
            yCell[iCell] < Ly  + e and
            xCell[iCell] > 0.0 - e and
            xCell[iCell] < 0.0 + e):

            for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

                for iVertexOnCell2 in range(0,nEdgesOnCell[iCell]):
                    iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

                    if (math.fabs(yVertex[iVertex1] - yVertex[iVertex2]) < e and
                        xVertex[iVertex1] < xVertex[iVertex2] - e):

                        pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
                        pairGlobalIndex[iVertex1] = iVertex2

    # right
    for iCell in range(0,nCells):
        if (yCell[iCell] > 0.0 - e and
            yCell[iCell] < Ly  + e and
            xCell[iCell] > Lx  - e and
            xCell[iCell] < Lx  + e):

            for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

                for iVertexOnCell2 in range(0,nEdgesOnCell[iCell]):
                    iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

                    if (math.fabs(yVertex[iVertex1] - yVertex[iVertex2]) < e and
                        xVertex[iVertex1] > xVertex[iVertex2] + e):

                        pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
                        pairGlobalIndex[iVertex1] = iVertex2

    # bottom left
    for iCell in range(0,nCells):
        if (xCell[iCell] > 0.0 - e and
            xCell[iCell] < 0.0 + e and
            yCell[iCell] > 0.0 - e and
            yCell[iCell] < 0.0 + e):

            for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

                if (xVertex[iVertex1] >= xCell[iCell] - dcx * 0.5 - e and
                    xVertex[iVertex1] <= xCell[iCell] - dcx * 0.5 + e and
                    yVertex[iVertex1] >= yCell[iCell] - a   * 0.5 - e and
                    yVertex[iVertex1] <= yCell[iCell] - a   * 0.5 + e):

                    for iVertexOnCell2 in range(0,nEdgesOnCell[iCell]):
                        iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

                        if (xVertex[iVertex2] >= xVertex[iVertex1] + dcx - e and
                            xVertex[iVertex2] <= xVertex[iVertex1] + dcx + e and
                            yVertex[iVertex2] >= yVertex[iVertex1] + a   - e and
                            yVertex[iVertex2] <= yVertex[iVertex1] + a   + e):

                            #pairType[iVertex1] = VELOCITY_BOUNDARY_PERIODIC
                            pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
                            pairGlobalIndex[iVertex1] = iVertex2

    # bottom right
    for iCell in range(0,nCells):
        if (xCell[iCell] > Lx  - e and
            xCell[iCell] < Lx  + e and
            yCell[iCell] > 0.0 - e and
            yCell[iCell] < 0.0 + e):

            for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

                if (xVertex[iVertex1] >= xCell[iCell] + dcx * 0.5 - e and
                    xVertex[iVertex1] <= xCell[iCell] + dcx * 0.5 + e and
                    yVertex[iVertex1] >= yCell[iCell] - a   * 0.5 - e and
                    yVertex[iVertex1] <= yCell[iCell] - a   * 0.5 + e):

                    for iVertexOnCell2 in range(0,nEdgesOnCell[iCell]):
                        iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

                        if (xVertex[iVertex2] >= xVertex[iVertex1] - dcx - e and
                            xVertex[iVertex2] <= xVertex[iVertex1] - dcx + e and
                            yVertex[iVertex2] >= yVertex[iVertex1] + a   - e and
                            yVertex[iVertex2] <= yVertex[iVertex1] + a   + e):

                            #pairType[iVertex1] = VELOCITY_BOUNDARY_PERIODIC
                            pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
                            pairGlobalIndex[iVertex1] = iVertex2

    # top left
    for iCell in range(0,nCells):
        if (xCell[iCell] > 0.0 - e and
            xCell[iCell] < 0.0 + e and
            yCell[iCell] > Ly  - e and
            yCell[iCell] < Ly  + e):

            for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

                if (xVertex[iVertex1] >= xCell[iCell] - dcx * 0.5 - e and
                    xVertex[iVertex1] <= xCell[iCell] - dcx * 0.5 + e and
                    yVertex[iVertex1] >= yCell[iCell] + a   * 0.5 - e and
                    yVertex[iVertex1] <= yCell[iCell] + a   * 0.5 + e):

                    for iVertexOnCell2 in range(0,nEdgesOnCell[iCell]):
                        iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

                        if (xVertex[iVertex2] >= xVertex[iVertex1] + dcx - e and
                            xVertex[iVertex2] <= xVertex[iVertex1] + dcx + e and
                            yVertex[iVertex2] >= yVertex[iVertex1] - a   - e and
                            yVertex[iVertex2] <= yVertex[iVertex1] - a   + e):

                            #pairType[iVertex1] = VELOCITY_BOUNDARY_PERIODIC
                            pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
                            pairGlobalIndex[iVertex1] = iVertex2

    # top right
    for iCell in range(0,nCells):
        if (xCell[iCell] > Lx  - e and
            xCell[iCell] < Lx  + e and
            yCell[iCell] > Ly  - e and
            yCell[iCell] < Ly  + e):

            for iVertexOnCell1 in range(0,nEdgesOnCell[iCell]):
                iVertex1 = verticesOnCell[iCell,iVertexOnCell1]

                if (xVertex[iVertex1] >= xCell[iCell] + dcx * 0.5 - e and
                    xVertex[iVertex1] <= xCell[iCell] + dcx * 0.5 + e and
                    yVertex[iVertex1] >= yCell[iCell] + a   * 0.5 - e and
                    yVertex[iVertex1] <= yCell[iCell] + a   * 0.5 + e):

                    for iVertexOnCell2 in range(0,nEdgesOnCell[iCell]):
                        iVertex2 = verticesOnCell[iCell,iVertexOnCell2]

                        if (xVertex[iVertex2] >= xVertex[iVertex1] - dcx - e and
                            xVertex[iVertex2] <= xVertex[iVertex1] - dcx + e and
                            yVertex[iVertex2] >= yVertex[iVertex1] - a   - e and
                            yVertex[iVertex2] <= yVertex[iVertex1] - a   + e):

                            #pairType[iVertex1] = VELOCITY_BOUNDARY_PERIODIC
                            pairType[iVertex1] = VELOCITY_BOUNDARY_REVERSE
                            pairGlobalIndex[iVertex1] = iVertex2

    for iVertex in range(0,nVertices):
        if (pairType[iVertex] == VELOCITY_BOUNDARY_NONE):
            solveVelocity[iVertex] = 1
        else:
            solveVelocity[iVertex] = 0

    for iCell in range(0,nCells):
        if (xCell[iCell] > 0.0 - e and
            xCell[iCell] < Lx  + e and
            yCell[iCell] > 0.0 - e and
            yCell[iCell] < Ly  + e):
            solveStress[iCell] = 1
        else:
            solveStress[iCell] = 0


    # fortran indices
    pairGlobalIndex[:] = pairGlobalIndex[:] + 1


    boundaryname = "boundary_%4.4i.nc" %(nx)

    boundaryFile = Dataset(boundaryname,"w",format="NETCDF3_CLASSIC")

    boundaryFile.createDimension("nVertices",nVertices)
    boundaryFile.createDimension("nCells",nCells)

    var = boundaryFile.createVariable("vertexBoundarySource","i",dimensions=["nVertices"])
    var[:] = pairGlobalIndex[:]

    var = boundaryFile.createVariable("vertexBoundaryType","i",dimensions=["nVertices"])
    var[:] = pairType[:]

    var = boundaryFile.createVariable("solveVelocitySpecialBoundaries","i",dimensions=["nVertices"])
    var[:] = solveVelocity[:]

    var = boundaryFile.createVariable("solveStressSpecialBoundaries","i",dimensions=["nCells"])
    var[:] = solveStress[:]

    boundaryFile.close()

    if (plot):

        plt.scatter(xVertex, yVertex, c=pairType, s=3, cmap='jet')
        plt.colorbar()

        rect = plt.Rectangle((0.0, 0.0), Lx, Ly, fill=False, edgecolor="black", linewidth=0.5)
        plt.gca().add_patch(rect)

        for iVertex1 in range(0,nVertices):
            if (pairType[iVertex1] == 2 or pairType[iVertex1] == 3):
                iVertex2 = pairGlobalIndex[iVertex1]
                plt.arrow(xVertex[iVertex1],yVertex[iVertex1],xVertex[iVertex2]-xVertex[iVertex1],yVertex[iVertex2]-yVertex[iVertex1], length_includes_head=True, head_width=400, head_length=500)
        
        plt.show()

#-------------------------------------------------------------------------------

def create_boundary_info():

    #nxs = [51,101,201,401,801,1601]
    nxs = [81]

    for nx in nxs:
        create_boundary_info_grid(nx)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_boundary_info()
