from netCDF4 import Dataset
import numpy as np
import math

#-------------------------------------------------------------

def wind_velocity(x, y, Lx, Ly):

    t = 0.0#21600.0
    tau = 4.0 * 24.0 * 3600.0

    u = 5.0 + (math.sin((2.0 * math.pi * t) / tau) - 3.0) * math.sin((2.0 * math.pi * x) / Lx) * math.sin((math.pi * y) / Ly)
    v = 5.0 + (math.sin((2.0 * math.pi * t) / tau) - 3.0) * math.sin((2.0 * math.pi * y) / Ly) * math.sin((math.pi * x) / Lx)
    #u = 5.0
    #v = 0.0

    return u, v

#-------------------------------------------------------------

def ocean_currents(x, y, Lx, Ly):

    u =  0.1 * ((2.0 * y - Ly) / Ly)
    v = -0.1 * ((2.0 * x - Lx) / Lx)
    #u = 0.0
    #v = 0.0

    return u, v

#-------------------------------------------------------------

def ice_concentration(x, y, Lx, Ly):

    conc = max(min(x / Lx, 1.0),0.0)

    return conc

#-------------------------------------------------------------

def create_ic(gridfile, icfile):

    # load grid file
    grid = Dataset(gridfile, "r")

    xCell = grid.variables["xCell"][:]
    yCell = grid.variables["yCell"][:]

    xVertex = grid.variables["xVertex"][:]
    yVertex = grid.variables["yVertex"][:]

    nCells = len(grid.dimensions["nCells"])
    nVertices = len(grid.dimensions["nVertices"])

    # calculate output variables
    uAirVelocity = np.empty(nCells)
    vAirVelocity = np.empty(nCells)

    uOceanVelocity = np.empty(nCells)
    vOceanVelocity = np.empty(nCells)

    iceConcentration = np.empty((nCells,1,1))
    iceVolume        = np.empty((nCells,1,1))

    fVertex = np.empty((nVertices,1,1))

    xmin = np.amin(xVertex)
    xmax = np.amax(xVertex)
    Lx = xmax - xmin

    print(xmin, xmax, Lx)

    ymin = np.amin(yVertex)
    ymax = np.amax(yVertex)
    Ly = ymax - ymin

    print(ymin, ymax, Ly)

    for iCell in range(0,nCells):

        x = xCell[iCell] - xmin
        y = yCell[iCell] - ymin

        uAirVelocity[iCell],   vAirVelocity[iCell]   = wind_velocity (x, y, Lx, Ly)
        uOceanVelocity[iCell], vOceanVelocity[iCell] = ocean_currents(x, y, Lx, Ly)

        iceConcentration[iCell,0,0] = ice_concentration(x, y, Lx, Ly)

        iceVolume[iCell,0,0] = 2.0 * iceConcentration[iCell,0,0]

    for iVertex in range(0,nVertices):
        fVertex[iVertex] = 1.46e-4

    # create output file
    ic = Dataset(icfile, "w", format="NETCDF3_64BIT")
    # format: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET' or 'NETCDF3_64BIT_DATA

    ic.createDimension("nCells", nCells)
    ic.createDimension("nVertices", nVertices)
    ic.createDimension("nCategories", 1)
    ic.createDimension("ONE", 1)
    ic.createDimension("Time", None)

    uAirVelocityVar = ic.createVariable("uAirVelocity", 'd', ('nCells'))
    vAirVelocityVar = ic.createVariable("vAirVelocity", 'd', ('nCells'))
    uAirVelocityVar[:] = uAirVelocity[:]
    vAirVelocityVar[:] = vAirVelocity[:]

    uOceanVelocityVar = ic.createVariable("uOceanVelocity", 'd', ('nCells'))
    vOceanVelocityVar = ic.createVariable("vOceanVelocity", 'd', ('nCells'))
    uOceanVelocityVar[:] = uOceanVelocity[:]
    vOceanVelocityVar[:] = vOceanVelocity[:]

    iceConcentrationVar = ic.createVariable("iceAreaCategory", 'd', ('nCells', 'nCategories', 'ONE'))
    iceConcentrationVar[:,:,:] = iceConcentration[:,:,:]

    iceVolumeVar = ic.createVariable("iceVolumeCategory", 'd', ('nCells', 'nCategories', 'ONE'))
    iceVolumeVar[:,:,:] = iceVolume[:,:,:]

    fVertexVar = ic.createVariable("fVertex", 'd', ('nVertices'))
    fVertexVar[:] = fVertex[:]

    ic.close()
    grid.close()

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
