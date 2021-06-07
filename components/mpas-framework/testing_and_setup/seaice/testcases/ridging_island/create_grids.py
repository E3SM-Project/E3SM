from netCDF4 import Dataset
import math
import os
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def in_box(x, y, x1, x2, y1, y2):

    if (x > x1 and x <= x2 and y > y1 and y <= y2):
        return True
    else:
        return False

#-------------------------------------------------------------------------------

def create_grid(gridfileOut, icfileOut, specialBoundariedFileOut, island=True):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    dc = 10000.0
    nx = 100
    ny = 100
    vertexDegree = 4

    fileGrid = Dataset("grid_in.nc","w",format="NETCDF3_CLASSIC")

    fileGrid.on_a_sphere = "NO"

    nCells = nx * ny
    nVertices = (nx+1) * (ny+1)

    fileGrid.createDimension("nCells", nCells)
    fileGrid.createDimension("nVertices", nVertices)
    fileGrid.createDimension("vertexDegree", vertexDegree)

    xCell = np.zeros(nCells)
    yCell = np.zeros(nCells)
    zCell = np.zeros(nCells)

    for j in range(0,ny):
        for i in range(0,nx):
            iCell = i + j * nx
            xCell[iCell] = (float(i) + 0.5) * dc
            yCell[iCell] = (float(j) + 0.5) * dc

    xVertex = np.zeros(nVertices)
    yVertex = np.zeros(nVertices)
    zVertex = np.zeros(nVertices)
    cellsOnVertex = np.zeros((nVertices,vertexDegree),dtype="i")

    for j in range(0,ny+1):
        for i in range(0,nx+1):
            iVertex = i + j * (nx+1)
            xVertex[iVertex] = float(i) * dc
            yVertex[iVertex] = float(j) * dc

            ic = i - 1 ; jc = j - 1
            iCell = ic + jc * nx
            if (ic >= 0 and jc >= 0):
                cellsOnVertex[iVertex,0] = iCell
            else:
                cellsOnVertex[iVertex,0] = -1

            ic = i ; jc = j - 1
            iCell = ic + jc * nx
            if (ic < nx and jc >= 0):
                cellsOnVertex[iVertex,1] = iCell
            else:
                cellsOnVertex[iVertex,1] = -1

            ic = i ; jc = j
            iCell = ic + jc * nx
            if (ic < nx and jc < ny):
                cellsOnVertex[iVertex,2] = iCell
            else:
                cellsOnVertex[iVertex,2] = -1

            ic = i - 1 ; jc = j
            iCell = ic + jc * nx
            if (ic >= 0 and jc < ny):
                cellsOnVertex[iVertex,3] = iCell
            else:
                cellsOnVertex[iVertex,3] = -1

            #print(iVertex+1,cellsOnVertex[iVertex,:]+1)

    var = fileGrid.createVariable("xCell","d",dimensions=["nCells"])
    var[:] = xCell[:]
    var = fileGrid.createVariable("yCell","d",dimensions=["nCells"])
    var[:] = yCell[:]
    var = fileGrid.createVariable("zCell","d",dimensions=["nCells"])
    var[:] = zCell[:]

    var = fileGrid.createVariable("xVertex","d",dimensions=["nVertices"])
    var[:] = xVertex[:]
    var = fileGrid.createVariable("yVertex","d",dimensions=["nVertices"])
    var[:] = yVertex[:]
    var = fileGrid.createVariable("zVertex","d",dimensions=["nVertices"])
    var[:] = zVertex[:]

    cellsOnVertex[:] = cellsOnVertex[:] + 1
    var = fileGrid.createVariable("cellsOnVertex","i",dimensions=["nVertices","vertexDegree"])
    var[:] = cellsOnVertex[:]

    fileGrid.close()

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_in.nc grid_conv.nc" %(mpas_tools_dir))

    filein = Dataset("grid_conv.nc","a")

    nCells = len(filein.dimensions["nCells"])

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    cullCell = np.zeros(nCells, dtype="i")

    if (island):
        for iCell in range(0,nCells):
            x = xCell[iCell]
            y = yCell[iCell]
            if (    in_box(x,y,400000.0,600000.0,400000.0,600000.0) and
                not in_box(x,y,400000.0,550000.0,400000.0,550000.0)):
                cullCell[iCell] = 1

    var = filein.createVariable("cullCell","i",dimensions=["nCells"])
    var[:] = cullCell[:]

    filein.close()

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasCellCuller.x grid_conv.nc grid_culled_island.nc" %(mpas_tools_dir))

    cmd = "mv grid_culled_island.nc %s" %(gridfileOut)
    os.system(cmd)

    # grid
    filein = Dataset(gridfileOut,"a")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])

    print("%s: nCells: %i" %(gridfileOut,nCells))

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    verticesOnCell = filein.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    filein.close()

    # plot
    plot = True
    if (plot):

        filediag = open("grid_diag.txt","w")

        for iCell in range(0,nCells):
            for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
                iVertex = verticesOnCell[iCell,iVertexOnCell]
                filediag.write("%f %f\n" %(xVertex[iVertex],yVertex[iVertex]))
            iVertex = verticesOnCell[iCell,0]
            filediag.write("%f %f\n" %(xVertex[iVertex],yVertex[iVertex]))
            filediag.write("\n")

        filediag.close()

    # initial conditions
    fileic = Dataset(icfileOut,"w",format="NETCDF3_CLASSIC")

    fileic.createDimension("nCells",nCells)

    var = fileic.createVariable("uAirVelocity","d",dimensions=["nCells"])
    var[:] = 10.0

    var = fileic.createVariable("vAirVelocity","d",dimensions=["nCells"])
    var[:] = 10.0

    fileic.close()

    # special boundaries
    TRACER_BOUNDARY_ZERO = 1

    tracerBoundaryType = np.zeros(nCells,dtype="i")

    for iCell in range(0,nCells):
        if (xCell[iCell] == 995000.0 or yCell[iCell] == 995000.0):
            tracerBoundaryType[iCell] = TRACER_BOUNDARY_ZERO

    fileboundary = Dataset(specialBoundariedFileOut,"w",format="NETCDF3_CLASSIC")

    fileic.createDimension("nCells",nCells)

    var = fileic.createVariable("tracerBoundaryType","i",dimensions=["nCells"])
    var[:] = tracerBoundaryType[:]

    fileboundary.close()

#-------------------------------------------------------------------------------

def create_grids():

    create_grid("grid_island.nc", "ic_island.nc", "special_boundaries_island.nc", True)
    create_grid("grid_no_island.nc", "ic_no_island.nc", "special_boundaries_no_island.nc", False)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_grids()
