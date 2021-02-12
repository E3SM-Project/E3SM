from netCDF4 import Dataset
import math
import os
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def create_grid(gridType):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    fileout = open("namelist.input","w")

    line = "&periodic_grid\n"
    fileout.write(line)

    line = "   nx = 10\n"
    fileout.write(line)

    line = "   ny = 10\n"
    fileout.write(line)

    line = "   dc = 1\n"
    fileout.write(line)

    line = "   nVertLevels = 1,\n"
    fileout.write(line)

    line = "   nTracers = 1,\n"
    fileout.write(line)

    line = "   nproc = 1,\n"
    fileout.write(line)

    line = "/\n"
    fileout.write(line)

    fileout.close()

    if (gridType == "hex"):
        os.system("%s/mesh_tools/periodic_hex/periodic_grid" %(mpas_tools_dir))
    elif (gridType == "quad"):
        os.system("%s/mesh_tools/periodic_quad/periodic_grid" %(mpas_tools_dir))
    os.system("python %s/mesh_tools/periodic_hex/mark_periodic_boundaries_for_culling.py -f grid.nc" %(mpas_tools_dir))
    os.system("%s/mesh_tools/mesh_conversion_tools/MpasCellCuller.x grid.nc grid_culled.nc" %(mpas_tools_dir))

    filein = Dataset("grid_culled.nc","a")

    nCells = len(filein.dimensions["nCells"])
    vertexDegree = len(filein.dimensions["vertexDegree"])

    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    iVertex = 46

    xC = xVertex[iVertex]
    yC = yVertex[iVertex]

    cullCell = np.ones(nCells,dtype="i")
    for iCellOnVertex in range(0,vertexDegree):
        iCell = cellsOnVertex[iVertex,iCellOnVertex]
        print(iCell)
        cullCell[iCell] = 0

    var = filein.createVariable("cullCell","i",dimensions=["nCells"])
    var[:] = cullCell[:]

    filein.close()

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasCellCuller.x grid_culled.nc grid_culled2.nc" %(mpas_tools_dir))

    filein = Dataset("grid_culled2.nc","a")

    xCell = filein.variables["xCell"]
    yCell = filein.variables["yCell"]
    xCell[:] -= xC
    yCell[:] -= yC

    xVertex = filein.variables["xVertex"]
    yVertex = filein.variables["yVertex"]
    xVertex[:] -= xC
    yVertex[:] -= yC

    xEdge = filein.variables["xEdge"]
    yEdge = filein.variables["yEdge"]
    xEdge[:] -= xC
    yEdge[:] -= yC

    nVertices = len(filein.dimensions["nVertices"])
    vertexDegree = len(filein.dimensions["vertexDegree"])
    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    for iVertex in range(0,nVertices):
        isCenterVertex = True
        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            if (iCell == -1):
                isCenterVertex = False
        if (isCenterVertex):
            iVertexCenter = iVertex

    filein.iVertexCenter = iVertexCenter

    filein.close()

    gridname = "grid_var_%s.nc" %(gridType)

    cmd = "mv grid_culled2.nc %s" %(gridname)
    os.system(cmd)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_grid("hex")
    create_grid("quad")
