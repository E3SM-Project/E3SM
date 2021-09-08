from netCDF4 import Dataset
import numpy as np
from math import sin, cos, pi, radians

#-------------------------------------------------------------

def velocities_strains_stress_divergences(x, y):

    A = 2.56#np.random.uniform(2,4)
    B = 2.56#np.random.uniform(2,4)
    C = 2.56#np.random.uniform(2,4)
    D = 2.56#np.random.uniform(2,4)

    Lx = 1.0
    Ly = 1.0

    u = sin((2.0 * pi * x * A) / Lx) * sin((2.0 * pi * y * B) / Ly)
    v = sin((2.0 * pi * x * C) / Lx) * sin((2.0 * pi * y * D) / Ly)

    dudx = ((2.0 * pi * A) / Lx) * cos((2.0 * pi * x * A) / Lx) * sin((2.0 * pi * y * B) / Ly)
    dudy = ((2.0 * pi * B) / Ly) * sin((2.0 * pi * x * A) / Lx) * cos((2.0 * pi * y * B) / Ly)

    dvdx = ((2.0 * pi * C) / Lx) * cos((2.0 * pi * x * C) / Lx) * sin((2.0 * pi * y * D) / Ly)
    dvdy = ((2.0 * pi * D) / Ly) * sin((2.0 * pi * x * C) / Lx) * cos((2.0 * pi * y * D) / Ly)

    d2udx2  = -((2.0 * pi * A) / Lx)*((2.0 * pi * A) / Lx) * sin((2.0 * pi * x * A) / Lx) * sin((2.0 * pi * y * B) / Ly)
    d2udy2  = -((2.0 * pi * B) / Ly)*((2.0 * pi * B) / Ly) * sin((2.0 * pi * x * A) / Lx) * sin((2.0 * pi * y * B) / Ly)
    d2udxdy =  ((2.0 * pi * A) / Lx)*((2.0 * pi * B) / Ly) * cos((2.0 * pi * x * A) / Lx) * cos((2.0 * pi * y * B) / Ly)

    d2vdx2  = -((2.0 * pi * C) / Lx)*((2.0 * pi * C) / Lx) * sin((2.0 * pi * x * C) / Lx) * sin((2.0 * pi * y * D) / Ly)
    d2vdy2  = -((2.0 * pi * D) / Ly)*((2.0 * pi * D) / Ly) * sin((2.0 * pi * x * C) / Lx) * sin((2.0 * pi * y * D) / Ly)
    d2vdxdy =  ((2.0 * pi * C) / Lx)*((2.0 * pi * D) / Ly) * cos((2.0 * pi * x * C) / Lx) * cos((2.0 * pi * y * D) / Ly)

    e11 = dudx
    e22 = dvdy
    e12 = 0.5 * (dudy + dvdx)

    de11dx = d2udx2

    de12dy = 0.5 * (d2udy2  + d2vdxdy)
    de12dx = 0.5 * (d2udxdy + d2vdx2)

    de22dy = d2vdy2

    divu = de11dx + de12dy
    divv = de12dx + de22dy

    return u, v, e11, e22, e12, divu, divv

#-------------------------------------------------------------

def create_ic(gridfile, icfile):

    # load grid file
    grid = Dataset(gridfile, "r")

    nCells = len(grid.dimensions["nCells"])
    nVertices = len(grid.dimensions["nVertices"])
    maxEdges = len(grid.dimensions["maxEdges"])

    nEdgesOnCell = grid.variables["nEdgesOnCell"][:]
    verticesOnCell = grid.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    xCell = grid.variables["xCell"][:]
    yCell = grid.variables["yCell"][:]

    xVertex = grid.variables["xVertex"][:]
    yVertex = grid.variables["yVertex"][:]

    grid.close()

    xMin = np.amin(xVertex)
    xMax = np.amax(xVertex)
    yMin = np.amin(yVertex)
    yMax = np.amax(yVertex)

    # calculate output variables
    uVelocity = np.empty(nVertices)
    vVelocity = np.empty(nVertices)

    stressDivergenceU = np.empty(nVertices)
    stressDivergenceV = np.empty(nVertices)

    strain11Vertex = np.zeros(nVertices)
    strain22Vertex = np.zeros(nVertices)
    strain12Vertex = np.zeros(nVertices)

    solveVelocity = np.ones(nVertices,dtype="i")
    solveVelocityPrevious = np.ones(nVertices,dtype="i")
    solveStress = np.ones(nCells,dtype="i")

    for iVertex in range(0,nVertices):

        x = (xVertex[iVertex] - xMin)
        y = (yVertex[iVertex] - yMin)

        u, v, e11, e22, e12, divu, divv = velocities_strains_stress_divergences(x, y)

        uVelocity[iVertex] = u
        vVelocity[iVertex] = v

        strain11Vertex[iVertex] = e11
        strain22Vertex[iVertex] = e22
        strain12Vertex[iVertex] = e12

        stressDivergenceU[iVertex] = divu
        stressDivergenceV[iVertex] = divv

    strain11Cell = np.zeros(nCells)
    strain22Cell = np.zeros(nCells)
    strain12Cell = np.zeros(nCells)

    for iCell in range(0, nCells):

        x = (xCell[iCell] - xMin)
        y = (yCell[iCell] - yMin)

        u, v, e11, e22, e12, divu, divv = velocities_strains_stress_divergences(x, y)

        strain11Cell[iCell] = e11
        strain22Cell[iCell] = e22
        strain12Cell[iCell] = e12

    stress11var = np.zeros((nCells, maxEdges))
    stress22var = np.zeros((nCells, maxEdges))
    stress12var = np.zeros((nCells, maxEdges))

    for iCell in range(0, nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            stress11var[iCell,iVertexOnCell] = strain11Vertex[iVertex]
            stress22var[iCell,iVertexOnCell] = strain22Vertex[iVertex]
            stress12var[iCell,iVertexOnCell] = strain12Vertex[iVertex]


    # create output file
    fileOut = Dataset(icfile, "w", format="NETCDF3_CLASSIC")

    fileOut.createDimension("nVertices", nVertices)
    fileOut.createDimension("nCells", nCells)
    fileOut.createDimension("maxEdges", maxEdges)

    var = fileOut.createVariable("uVelocity","d",dimensions=["nVertices"])
    var[:] = uVelocity[:]

    var = fileOut.createVariable("vVelocity","d",dimensions=["nVertices"])
    var[:] = vVelocity[:]

    var = fileOut.createVariable("solveVelocity","i",dimensions=["nVertices"])
    var[:] = solveVelocity[:]

    var = fileOut.createVariable("solveVelocityPrevious","i",dimensions=["nVertices"])
    var[:] = solveVelocityPrevious[:]

    var = fileOut.createVariable("solveStress","i",dimensions=["nCells"])
    var[:] = solveStress[:]

    var = fileOut.createVariable("strain11VertexAnalytical","d",dimensions=["nVertices"])
    var[:] = strain11Vertex[:]

    var = fileOut.createVariable("strain22VertexAnalytical","d",dimensions=["nVertices"])
    var[:] = strain22Vertex[:]

    var = fileOut.createVariable("strain12VertexAnalytical","d",dimensions=["nVertices"])
    var[:] = strain12Vertex[:]

    var = fileOut.createVariable("strain11CellAnalytical","d",dimensions=["nCells"])
    var[:] = strain11Cell[:]

    var = fileOut.createVariable("strain22CellAnalytical","d",dimensions=["nCells"])
    var[:] = strain22Cell[:]

    var = fileOut.createVariable("strain12CellAnalytical","d",dimensions=["nCells"])
    var[:] = strain12Cell[:]

    var = fileOut.createVariable("stress11var","d",dimensions=["nCells","maxEdges"])
    var[:] = stress11var[:]

    var = fileOut.createVariable("stress22var","d",dimensions=["nCells","maxEdges"])
    var[:] = stress22var[:]

    var = fileOut.createVariable("stress12var","d",dimensions=["nCells","maxEdges"])
    var[:] = stress12var[:]

    var = fileOut.createVariable("stress11weak","d",dimensions=["nCells"])
    var[:] = strain11Cell[:]

    var = fileOut.createVariable("stress22weak","d",dimensions=["nCells"])
    var[:] = strain22Cell[:]

    var = fileOut.createVariable("stress12weak","d",dimensions=["nCells"])
    var[:] = strain12Cell[:]

    var = fileOut.createVariable("stressDivergenceUAnalytical","d",dimensions=["nVertices"])
    var[:] = stressDivergenceU[:]

    var = fileOut.createVariable("stressDivergenceVAnalytical","d",dimensions=["nVertices"])
    var[:] = stressDivergenceV[:]

    fileOut.close()

#-------------------------------------------------------------

def create_ics():

    gridTypes = ["hex","quad"]

    grids = {"hex": ["0082x0094",
                     "0164x0188",
                     "0328x0376",
                     "0656x0752"],
             "quad":["0080x0080",
                     "0160x0160",
                     "0320x0320",
                     "0640x0640"]}

    for gridType in gridTypes:
        for grid in grids[gridType]:

            gridfile = "grid_%s_%s.nc" %(gridType,grid)
            icfile = "ic_%s_%s.nc" %(gridType,grid)

            create_ic(gridfile, icfile)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ics()
