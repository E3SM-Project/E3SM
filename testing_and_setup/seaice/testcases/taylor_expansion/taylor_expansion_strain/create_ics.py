from netCDF4 import Dataset
import numpy as np
from math import sin, cos, pi, pow

#-------------------------------------------------------------

def velocities_linear(x,y):

    u = x
    v = y

    dudx = 1.0
    dudy = 0.0

    dvdx = 0.0
    dvdy = 1.0

    d2udx2  = 0.0
    d2udy2  = 0.0
    d2udxdy = 0.0

    d2vdx2  = 0.0
    d2vdy2  = 0.0
    d2vdxdy = 0.0

    return u, v, dudx, dudy, dvdx, dvdy, d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy

#-------------------------------------------------------------

def velocities_quadratic(x,y):

    u = pow(x, 2)
    v = pow(y, 2)

    dudx = 2.0 * x
    dudy = 0.0

    dvdx = 0.0
    dvdy = 2.0 * y

    d2udx2  = 2.0
    d2udy2  = 0.0
    d2udxdy = 0.0

    d2vdx2  = 0.0
    d2vdy2  = 2.0
    d2vdxdy = 0.0

    return u, v, dudx, dudy, dvdx, dvdy, d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy

#-------------------------------------------------------------

def velocities_cubic(x,y):

    u = pow(x, 3)
    v = pow(y, 3)

    dudx = 3.0 * pow(x, 2)
    dudy = 0.0

    dvdx = 0.0
    dvdy = 3.0 * pow(y, 2)

    d2udx2  = 6.0 * x
    d2udy2  = 0.0
    d2udxdy = 0.0

    d2vdx2  = 0.0
    d2vdy2  = 6.0 * y
    d2vdxdy = 0.0

    return u, v, dudx, dudy, dvdx, dvdy, d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy

#-------------------------------------------------------------

def velocities_quartic(x,y):

    u = pow(x, 4)
    v = pow(y, 4)

    dudx = 4.0 * pow(x, 3)
    dudy = 0.0

    dvdx = 0.0
    dvdy = 4.0 * pow(y, 3)

    d2udx2  = 12.0 * pow(x, 2)
    d2udy2  = 0.0
    d2udxdy = 0.0

    d2vdx2  = 0.0
    d2vdy2  = 12.0 * pow(y, 2)
    d2vdxdy = 0.0

    return u, v, dudx, dudy, dvdx, dvdy, d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy

#-------------------------------------------------------------

def velocities_strains_stress_divergences(x, y, velType):

    if (velType == "linear"):
        u, v, dudx, dudy, dvdx, dvdy, d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy = \
            velocities_linear(x,y)
    elif (velType == "quadratic"):
        u, v, dudx, dudy, dvdx, dvdy, d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy = \
            velocities_quadratic(x,y)
    elif (velType == "cubic"):
        u, v, dudx, dudy, dvdx, dvdy, d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy = \
            velocities_cubic(x,y)
    elif (velType == "quartic"):
        u, v, dudx, dudy, dvdx, dvdy, d2udx2, d2udy2, d2udxdy, d2vdx2, d2vdy2, d2vdxdy = \
            velocities_quartic(x,y)

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

def create_ic(gridfile, icfile, velType):

    # load grid file
    grid = Dataset(gridfile, "r")

    nCells = len(grid.dimensions["nCells"])
    nVertices = len(grid.dimensions["nVertices"])

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

        u, v, e11, e22, e12, divu, divv = velocities_strains_stress_divergences(x, y, velType)
        
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

        u, v, e11, e22, e12, divu, divv = velocities_strains_stress_divergences(x, y, velType)

        strain11Cell[iCell] = e11
        strain22Cell[iCell] = e22
        strain12Cell[iCell] = e12


    # create output file
    fileOut = Dataset(icfile, "w", format="NETCDF3_CLASSIC")

    fileOut.velType = velType

    fileOut.createDimension("nVertices", nVertices)
    fileOut.createDimension("nCells", nCells)

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

    var = fileOut.createVariable("stressDivergenceUAnalytical","d",dimensions=["nVertices"])
    var[:] = stressDivergenceU[:]

    var = fileOut.createVariable("stressDivergenceVAnalytical","d",dimensions=["nVertices"])
    var[:] = stressDivergenceV[:]

    fileOut.close()

#-------------------------------------------------------------

def create_ics():

    operatorTypes = ["var","weak"]
    gridTypes = ["hex","quad"]
    velocityTypes = ["linear","quadratic","cubic","quartic"]
    
    for operatorType in operatorTypes:
        for gridType in gridTypes:
            for velocityType in velocityTypes:

                gridfile = "grid_%s_%s.nc" %(operatorType,gridType)
                icfile = "ic_%s_%s_%s.nc" %(operatorType,gridType,velocityType)
                
                create_ic(gridfile, icfile, velocityType)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ics()
