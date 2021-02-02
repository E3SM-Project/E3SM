from netCDF4 import Dataset
import numpy as np
import math
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def calc_testcase_ic():

    filenameGrid = "square_mesh_80x80_culled.nc"
    filenameIC = "ic_quad.nc"
    specialBoundaries = "True"
    dy = 16000.0

    filein = Dataset(filenameGrid, "r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    filein.close()

    xMin = np.amin(xVertex)
    yMin = np.amin(yVertex)

    xMax = np.amax(xVertex)
    yMax = np.amax(yVertex)

    Lx = np.amax(xVertex) - np.amin(xVertex)
    Ly = np.amax(yVertex) - np.amin(yVertex)


    fileout = Dataset(filenameIC, "w", format="NETCDF3_CLASSIC")

    fileout.createDimension("nCells", nCells)
    fileout.createDimension("nVertices", nVertices)
    fileout.createDimension("nCategories", 1)
    fileout.createDimension("ONE", 1)

    uAirVelocityVar      = fileout.createVariable("uAirVelocity",      "d", dimensions=["nCells"])
    vAirVelocityVar      = fileout.createVariable("vAirVelocity",      "d", dimensions=["nCells"])
    uOceanVelocityVar    = fileout.createVariable("uOceanVelocity",    "d", dimensions=["nCells"])
    vOceanVelocityVar    = fileout.createVariable("vOceanVelocity",    "d", dimensions=["nCells"])
    iceAreaCategoryVar   = fileout.createVariable("iceAreaCategory",   "d", dimensions=["nCells", "nCategories", "ONE"])
    iceVolumeCategoryVar = fileout.createVariable("iceVolumeCategory", "d", dimensions=["nCells", "nCategories", "ONE"])
    fVertexVar           = fileout.createVariable("fVertex",           "d", dimensions=["nVertices"])

    uAirVelocityVar[:] = 1.0
    vAirVelocityVar[:] = 0.0
    uOceanVelocityVar[:] = 0.0
    vOceanVelocityVar[:] = 0.0
    fVertexVar[:] = 0.0

    for iCell in range(0,nCells):
        x = xCell[iCell] - xMin
        iceAreaCategoryVar[iCell] = x / Lx
        iceVolumeCategoryVar[iCell] = 2.0 * iceAreaCategoryVar[iCell]

    fileout.close()

    if (specialBoundaries):

        fileout = Dataset("special_boundaries.nc", "w", format="NETCDF3_CLASSIC")

        fileout.createDimension("nVertices", nVertices)

        vertexBoundaryTypeVar = fileout.createVariable("vertexBoundaryType", "i", dimensions=["nVertices"])
        vertexBoundarySourceVar = fileout.createVariable("vertexBoundarySource", "i", dimensions=["nVertices"])

        vertexBoundaryType = np.zeros(nVertices,dtype="i")
        vertexBoundarySource = np.zeros(nVertices,dtype="i")

        for iVertex in range(0,nVertices):
            if (math.fabs(yVertex[iVertex] - yMax) < 1.0):
                for iVertex2 in range(0,nVertices):
                    if (math.fabs(yVertex[iVertex2] - (yMin + dy)) < 1.0 and
                        math.fabs(xVertex[iVertex2] - xVertex[iVertex]) < 1.0):
                        vertexBoundaryType[iVertex] = 1
                        vertexBoundarySource[iVertex] = iVertex2+1
            if (math.fabs(yVertex[iVertex] - yMin) < 1.0):
                for iVertex2 in range(0,nVertices):
                    if (math.fabs(yVertex[iVertex2] - (yMax - dy)) < 1.0 and
                        math.fabs(xVertex[iVertex2] - xVertex[iVertex]) < 1.0):
                        vertexBoundaryType[iVertex] = 1
                        vertexBoundarySource[iVertex] = iVertex2+1

        vertexBoundaryTypeVar[:] = vertexBoundaryType[:]
        vertexBoundarySourceVar[:] = vertexBoundarySource[:]

        fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    calc_testcase_ic()
