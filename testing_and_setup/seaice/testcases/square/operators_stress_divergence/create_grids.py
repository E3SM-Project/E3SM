from netCDF4 import Dataset
import math
import os
import matplotlib.pyplot as plt
import numpy as np

#-------------------------------------------------------------------------------

def create_grid_hex(nx, ny, dc):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    Lx = float(nx) * dc
    Ly = float(ny) * dc

    gridname = "grid_hex_%4.4ix%4.4i.nc" %(nx,ny)

    fileout = open("namelist.input","w")

    line = "&periodic_grid\n"
    fileout.write(line)

    line = "   nx = %i,\n" %(nx+2)
    fileout.write(line)

    line = "   ny = %i,\n" %(ny+2)
    fileout.write(line)

    line = "   dc = %f,\n" %(dc)
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

    os.system("%s/mesh_tools/periodic_hex/periodic_grid" %(mpas_tools_dir))
    os.system("python %s/mesh_tools/periodic_hex/mark_periodic_boundaries_for_culling.py -f grid.nc" %(mpas_tools_dir))
    os.system("%s/mesh_tools/mesh_conversion_tools/MpasCellCuller.x grid.nc grid_culled.nc" %(mpas_tools_dir))

    filein = Dataset("grid_culled.nc","a")

    filein.Lx = Lx
    filein.Ly = Ly

    filein.nx = nx
    filein.ny = ny

    filein.dc = dc

    filein.close()

    cmd = "mv grid_culled.nc %s" %(gridname)
    os.system(cmd)

#-------------------------------------------------------------------------------

def create_grid_quad(nx, ny, dc):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    Lx = float(nx) * dc
    Ly = float(ny) * dc

    gridname = "grid_quad_%4.4ix%4.4i.nc" %(nx,ny)

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

    os.system("%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_in.nc %s" %(mpas_tools_dir, gridname))

    filein = Dataset(gridname,"a")

    nCells = len(filein.dimensions["nCells"])
    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    verticesOnCell = filein.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1
    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    filein.Lx = Lx
    filein.Ly = Ly

    filein.nx = nx
    filein.ny = ny

    filein.dc = dc

    filein.close()

    # plot
    filenameplot = "grid_quad_%4.4ix%4.4i.txt" %(nx,ny)
    fileplot = open(filenameplot,"w")

    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            fileplot.write("%f %f\n" %(xVertex[iVertex],yVertex[iVertex]))
        iVertex = verticesOnCell[iCell,0]
        fileplot.write("%f %f\n" %(xVertex[iVertex],yVertex[iVertex]))
        fileplot.write("\n")

    fileplot.close()

#-------------------------------------------------------------------------------

def create_grids():

    nGrid = 4

    # hex
    dc = 0.0125
    nx = 82
    ny = 94

    nxs = []
    nys = []
    dcs = []
    for i in range(0,nGrid):
        nxs.append(nx)
        nys.append(ny)
        dcs.append(dc)
        nx = nx*2
        ny = ny*2
        dc = dc/2

    for nx, ny, dc in zip(nxs, nys, dcs):
        create_grid_hex(nx, ny, dc)

    # quad
    dc = 0.0125
    nx = 80
    ny = 80

    nxs = []
    nys = []
    dcs = []
    for i in range(0,nGrid):
        nxs.append(nx)
        nys.append(ny)
        dcs.append(dc)
        nx = nx*2
        ny = ny*2
        dc = dc/2

    for nx, ny, dc in zip(nxs, nys, dcs):
        create_grid_quad(nx, ny, dc)

#-------------------------------------------------------------

if __name__ == "__main__":

    create_grids()
