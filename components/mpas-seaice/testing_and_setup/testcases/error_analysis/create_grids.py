from shutil import copyfile
from shapely.geometry import Polygon
import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import os

#-------------------------------------------------------------------------------

def create_grids():

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']


    # copy file
    src = "../square/operators_strain/grid_hex_0082x0094.nc"
    dst = "grid_in.nc"
    copyfile(src, dst)


    # cull cells
    iVertexTest = 1123

    filein = Dataset("grid_in.nc","a")

    nCells = len(filein.dimensions["nCells"])
    vertexDegree = len(filein.dimensions["vertexDegree"])

    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    cullCell = np.ones(nCells)

    for iCellOnVertex in range(0,vertexDegree):
        iCell = cellsOnVertex[iVertexTest,iCellOnVertex]
        cullCell[iCell] = 0

    var = filein.createVariable("cullCell","i",dimensions=["nCells"])
    var[:] = cullCell[:]

    filein.close()

    cmd = "%s/mesh_tools/mesh_conversion_tools/MpasCellCuller.x grid_in.nc grid_culled.nc" %(mpas_tools_dir)
    os.system(cmd)


    # create input file
    filein = Dataset("grid_culled.nc","a")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])
    vertexDegree = len(filein.dimensions["vertexDegree"])

    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1
    verticesOnCell = filein.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    dcEdge = filein.variables["dcEdge"][:]
    dxy = dcEdge[0] * 0.1

    for iVertex in range(0,nVertices):
        isCentre = True
        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            if (iCell == -1):
                isCentre = False
        if (isCentre):
            iVertexTest = iVertex

    xVertex = filein.variables["xVertex"]
    yVertex = filein.variables["yVertex"]

    xC = xVertex[iVertexTest]
    yC = yVertex[iVertexTest]

    for iVertex in range(0,nVertices):
        if (iVertex != iVertexTest):
            dx = np.random.uniform(low=-dxy, high=dxy)
            dy = np.random.uniform(low=-dxy, high=dxy)
            xVertex[iVertex] += dx
            yVertex[iVertex] += dy

    xCell = filein.variables["xCell"]
    yCell = filein.variables["yCell"]

    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            vertices.append((xVertex[iVertex],yVertex[iVertex]))
        cell = Polygon(vertices)
        xCell[iCell] = cell.centroid.x
        yCell[iCell] = cell.centroid.y

    xCell[:] -= xC
    yCell[:] -= yC

    xVertex[:] -= xC
    yVertex[:] -= yC

    filein.close()

    # make converter input file
    filein = Dataset("grid_culled.nc","r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])
    vertexDegree = len(filein.dimensions["vertexDegree"])

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]
    zVertex = filein.variables["zVertex"][:]

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]
    zCell = filein.variables["zCell"][:]

    cellsOnVertex = filein.variables["cellsOnVertex"][:]

    filein.close()

    fileout = Dataset("grid_converter.nc","w",format="NETCDF3_CLASSIC")

    fileout.on_a_sphere = "NO"

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("nVertices",nVertices)
    fileout.createDimension("vertexDegree",vertexDegree)

    var = fileout.createVariable("xCell","d",dimensions=["nCells"])
    var[:] = xCell[:]
    var = fileout.createVariable("yCell","d",dimensions=["nCells"])
    var[:] = yCell[:]
    var = fileout.createVariable("zCell","d",dimensions=["nCells"])
    var[:] = zCell[:]

    var = fileout.createVariable("xVertex","d",dimensions=["nVertices"])
    var[:] = xVertex[:]
    var = fileout.createVariable("yVertex","d",dimensions=["nVertices"])
    var[:] = yVertex[:]
    var = fileout.createVariable("zVertex","d",dimensions=["nVertices"])
    var[:] = zVertex[:]

    var = fileout.createVariable("cellsOnVertex","i",dimensions=["nVertices","vertexDegree"])
    var[:] = cellsOnVertex[:]

    fileout.close()


    # create mesh
    cmd = "%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_converter.nc grid.nc" %(mpas_tools_dir)
    os.system(cmd)


    # plot mesh
    filein = Dataset("grid.nc","r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])
    vertexDegree = len(filein.dimensions["vertexDegree"])

    nEdgesOnCell = filein.variables["nEdgesOnCell"][:]
    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1
    verticesOnCell = filein.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]
    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]
    xEdge = filein.variables["xEdge"][:]
    yEdge = filein.variables["yEdge"][:]

    for iCellOnVertex in range(0,vertexDegree):
        iCell = cellsOnVertex[iVertexTest,iCellOnVertex]
        x = []
        y = []
        for iVertexOnCell in range(nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            x.append(xVertex[iVertex])
            y.append(yVertex[iVertex])
        iVertex = verticesOnCell[iCell,0]
        x.append(xVertex[iVertex])
        y.append(yVertex[iVertex])

        plt.plot(x,y)

    plt.scatter(xCell,yCell)
    plt.scatter(xEdge,yEdge)

    plt.savefig("grid.png")

    filein.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_grids()
