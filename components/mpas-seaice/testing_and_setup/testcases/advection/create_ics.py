from netCDF4 import Dataset
import math, os, sys

#--------------------------------------------------------------------

def create_ic_file(res, icType):

    # grid in
    gridFilename = "x1.%s.grid.nc" %(res)
    gridFile = Dataset(gridFilename, "r")

    nCells = len(gridFile.dimensions["nCells"])
    nVertices = len(gridFile.dimensions["nVertices"])

    xCell = gridFile.variables["xCell"][:]
    yCell = gridFile.variables["yCell"][:]
    zCell = gridFile.variables["zCell"][:]

    latVertex = gridFile.variables["latVertex"][:]

    gridFile.close()

    # ice out
    icFilename = "ic_%s_%s.nc" %(icType, res)
    icFile = Dataset(icFilename, "w", format="NETCDF3_CLASSIC")

    icFile.createDimension("nCells", nCells)
    icFile.createDimension("nVertices", nVertices)
    icFile.createDimension("nCategories", size=1)
    icFile.createDimension("ONE", size=1)

    uVelocity = icFile.createVariable("uVelocity", 'd', dimensions=("nVertices"))
    vVelocity = icFile.createVariable("vVelocity", 'd', dimensions=("nVertices"))

    days = 120.0
    seconds = days * 24.0 * 3600.0
    radius = 6371229.0
    uVelocityEquator = (2.0 * math.pi * radius) / (seconds)
    print("uVelocityEquator: ",uVelocityEquator)
    for iVertex in range(0,nVertices):
        uVelocity[iVertex] = uVelocityEquator * math.cos(latVertex[iVertex])
        vVelocity[iVertex] = 0.0

    iceAreaCell   = icFile.createVariable("iceAreaCell",   'd', dimensions=("nCells"))
    iceVolumeCell = icFile.createVariable("iceVolumeCell", 'd', dimensions=("nCells"))

    iceAreaCategory   = icFile.createVariable("iceAreaCategory",   'd', dimensions=("nCells","nCategories","ONE"))
    iceVolumeCategory = icFile.createVariable("iceVolumeCategory", 'd', dimensions=("nCells","nCategories","ONE"))

    iceAreaCell[:]   = 0.0
    iceVolumeCell[:] = 0.0

    if (icType == "slotted_cylinder"):

        circleRadius = 0.5

        for iCell in range(0,nCells):

            r = math.sqrt(math.pow(zCell[iCell],2) + math.pow(xCell[iCell],2))

            if (r < circleRadius and yCell[iCell] > 0.0):

                iceAreaCell[iCell]   = 1.0
                iceVolumeCell[iCell] = 1.0

        for iCell in range(0,nCells):

            if (math.fabs(xCell[iCell]) < 1.0/12.0 and zCell[iCell] > -2.0/6.0):

                iceAreaCell[iCell]   = 0.0
                iceVolumeCell[iCell] = 0.0

    elif (icType == "cosine_bell_volume"):

        circleRadius = 1.0/3.0

        for iCell in range(0,nCells):

            r = math.sqrt(math.pow(zCell[iCell],2) + math.pow(xCell[iCell],2))

            if (r < circleRadius and yCell[iCell] > 0.0):

                iceAreaCell[iCell]   = 1.0

                iceVolumeCell[iCell] = 0.5 * (1.0 + math.cos((math.pi * r) / circleRadius))

    elif (icType == "cosine_bell"):

        circleRadius = 1.0/3.0

        for iCell in range(0,nCells):

            r = math.sqrt(math.pow(zCell[iCell],2) + math.pow(xCell[iCell],2))

            if (r < circleRadius and yCell[iCell] > 0.0):

                iceAreaCell[iCell] = 0.5 * (1.0 + math.cos((math.pi * r) / circleRadius))

                iceVolumeCell[iCell]   = 1.0

    elif (icType == "cylinder"):

        circleRadius = 0.5

        for iCell in range(0,nCells):

            r = math.sqrt(math.pow(zCell[iCell],2) + math.pow(xCell[iCell],2))

            if (r < circleRadius and yCell[iCell] > 0.0):

                iceAreaCell[iCell]   = 1.0
                iceVolumeCell[iCell] = 1.0

    iceAreaCategory[:,0]   = iceAreaCell[:]
    iceVolumeCategory[:,0] = iceVolumeCell[:]

    icFile.close()

#--------------------------------------------------------------------

def create_ics():

    reses = ["2562","10242","40962","163842"]

    icTypes = ["cosine_bell","slotted_cylinder"]

    for icType in icTypes:

        print("icType: ", icType)

        for res in reses:

            print("  Res: ", res)

            create_ic_file(res, icType)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_ics()
