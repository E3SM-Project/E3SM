from netCDF4 import Dataset
import numpy as np
import math

#-------------------------------------------------------------

def wind_velocity(x, y, Lx, Ly):

    t = 21600.0
    tau = 4.0 * 24.0 * 3600.0

    u = 5.0 + (math.sin((2.0 * math.pi * t) / tau) - 3.0) * math.sin((2.0 * math.pi * x) / Lx) * math.sin((math.pi * y) / Ly)
    v = 5.0 + (math.sin((2.0 * math.pi * t) / tau) - 3.0) * math.sin((2.0 * math.pi * y) / Ly) * math.sin((math.pi * x) / Lx)
    #u = 1.0
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

def create_ic_grid(n, plot=False):

    # load grid file
    gridfile = "grid_%4.4i.nc" %(n)
    icfile   = "ic_%4.4i.nc" %(n)

    grid = Dataset(gridfile, "r")

    Lx = grid.Lx
    Ly = grid.Ly

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


    xcentre = 0.5 * Lx
    ycentre = 0.5 * Ly
    e = 1.0

    for iCell in range(0,nCells):

        xCorrection = Lx / 2.0 - xcentre
        yCorrection = Ly / 2.0 - ycentre

        uAirVelocity[iCell],   vAirVelocity[iCell]   = wind_velocity( xCell[iCell], yCell[iCell], Lx, Ly)
        uOceanVelocity[iCell], vOceanVelocity[iCell] = ocean_currents(xCell[iCell], yCell[iCell], Lx, Ly)

        if (xCell[iCell] > 0.0 - e and
            xCell[iCell] < Lx  + e and
            yCell[iCell] > 0.0 - e and
            yCell[iCell] < Ly  + e):
            iceConcentration[iCell,0,0] = ice_concentration(xCell[iCell], yCell[iCell], Lx, Ly)
            iceVolume[iCell,0,0] = 2.0 * iceConcentration[iCell,0,0]

    for iVertex in range(0,nVertices):
        fVertex[iVertex] = 1.46e-4

    # create output file

    ic = Dataset(icfile, "w", format="NETCDF3_64BIT")
    # format: 'NETCDF4', 'NETCDF4_CLASSIC', 'NETCDF3_CLASSIC', 'NETCDF3_64BIT_OFFSET' or 'NETCDF3_64BIT_DATA

    if (plot):
        for dname, the_dim in grid.dimensions.iteritems():
            ic.createDimension(dname, len(the_dim) if not the_dim.isunlimited() else None)

        for v_name, varin in grid.variables.iteritems():
            outVar = ic.createVariable(v_name, varin.datatype, varin.dimensions)
            outVar[:] = varin[:]

    else:
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

def create_ic():

    #ns = [51,101,201,401,801,1601]
    ns = [51]

    for n in ns:

        create_ic_grid(n)

#-------------------------------------------------------------

if __name__ == "__main__":

    create_ic()
