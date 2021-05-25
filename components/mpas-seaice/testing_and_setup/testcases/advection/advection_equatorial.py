from netCDF4 import Dataset
import math
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl

#---------------------------------------------------------

def get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, array):

    # find initial point
    radius = 6371229.0
    xStart = radius
    yStart = 0.0
    zStart = 0.0

    minDistance = 1e30
    iCellStart = -1
    for iCell in range(0,nCells):

        distance = math.sqrt(math.pow((xCell[iCell] - xStart),2) + \
                             math.pow((yCell[iCell] - yStart),2) + \
                             math.pow((zCell[iCell] - zStart),2))

        if (distance < minDistance):
            minDistance = distance
            iCellStart = iCell


    iCell = iCellStart

    lats = []
    lons = []
    y = []
    deltaLon = 0.0

    while True:

        lats.append(latCell[iCell])
        lons.append(lonCell[iCell])
        y.append(array[iCell])

        iCellNext = -1
        nextSmallestLatitude = 1e30

        for iCellOnCell in range(0,nEdgesOnCell[iCell]):

            iCellNeighbour = cellsOnCell[iCell,iCellOnCell] - 1

            lonNeighbour = lonCell[iCellNeighbour]

            if (lonNeighbour < lonCell[iCell] - math.pi):
                lonNeighbour = lonNeighbour + 2.0 * math.pi
            if (lonNeighbour > lonCell[iCell] + math.pi):
                lonNeighbour = lonNeighbour - 2.0 * math.pi

            if (lonNeighbour > lonCell[iCell]):

                if (math.fabs(latCell[iCellNeighbour]) < nextSmallestLatitude):

                    nextSmallestLatitude = math.fabs(latCell[iCellNeighbour])
                    iCellNext = iCellNeighbour
                    lonNeighbourAll = lonNeighbour

        deltaLon = deltaLon + (lonNeighbourAll - lonCell[iCell])

        if (deltaLon > 2.0 * math.pi):
            break

        iCell = iCellNext

    points = zip(lons,y)
    sorted_points = sorted(points)
    lons = [point[0] for point in sorted_points]
    y    = [point[1] for point in sorted_points]

    return lons, y

#---------------------------------------------------------

def advection_equatorial():

    #res = "2562"
    #res = "10242"
    res = "40962"
    #res = "163842"

    experiment1 = "cosine_bell"
    experiment2 = "slotted_cylinder"

    iTime = -1

    # cosine bell IR
    filename = "./output_IR_%s_%s/output.2000.nc" %(experiment1,res)

    print(filename)
    fileMPAS = Dataset(filename,"r")

    nCells = len(fileMPAS.dimensions["nCells"])

    latCell = fileMPAS.variables["latCell"][:]
    lonCell = fileMPAS.variables["lonCell"][:]

    xCell = fileMPAS.variables["xCell"][:]
    yCell = fileMPAS.variables["yCell"][:]
    zCell = fileMPAS.variables["zCell"][:]

    cellsOnCell = fileMPAS.variables["cellsOnCell"][:]
    nEdgesOnCell = fileMPAS.variables["nEdgesOnCell"][:]

    iceArea_CB_Initial   = fileMPAS.variables["iceAreaCategory"][0,:,0,0]
    iceVolume_CB_Initial = fileMPAS.variables["iceVolumeCategory"][0,:,0,0]
    iceThickness_CB_Initial = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_CB_Initial[iCell] > 0.0):
            iceThickness_CB_Initial[iCell] = iceVolume_CB_Initial[iCell] / iceArea_CB_Initial[iCell]

    iceArea_CB_IR   = fileMPAS.variables["iceAreaCategory"][iTime,:,0,0]
    iceVolume_CB_IR = fileMPAS.variables["iceVolumeCategory"][iTime,:,0,0]
    iceThickness_CB_IR = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_CB_IR[iCell] > 0.0):
            iceThickness_CB_IR[iCell] = iceVolume_CB_IR[iCell] / iceArea_CB_IR[iCell]

    fileMPAS.close()

    # cosine bell upwind
    filename = "./output_upwind_%s_%s/output.2000.nc" %(experiment1,res)

    print(filename)
    fileMPAS = Dataset(filename,"r")

    iceArea_CB_upwind   = fileMPAS.variables["iceAreaCategory"][iTime,:,0,0]
    iceVolume_CB_upwind = fileMPAS.variables["iceVolumeCategory"][iTime,:,0,0]
    iceThickness_CB_upwind = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_CB_upwind[iCell] > 0.0):
            iceThickness_CB_upwind[iCell] = iceVolume_CB_upwind[iCell] / iceArea_CB_upwind[iCell]

    fileMPAS.close()

    # slotted cylinder IR
    filename = "./output_IR_%s_%s/output.2000.nc" %(experiment2,res)

    print(filename)
    fileMPAS = Dataset(filename,"r")

    iceArea_SC_Initial   = fileMPAS.variables["iceAreaCategory"][0,:,0,0]
    iceVolume_SC_Initial = fileMPAS.variables["iceVolumeCategory"][0,:,0,0]
    iceThickness_SC_Initial = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_SC_Initial[iCell] > 0.0):
            iceThickness_SC_Initial[iCell] = iceVolume_SC_Initial[iCell] / iceArea_SC_Initial[iCell]

    iceArea_SC_IR   = fileMPAS.variables["iceAreaCategory"][iTime,:,0,0]
    iceVolume_SC_IR = fileMPAS.variables["iceVolumeCategory"][iTime,:,0,0]
    iceThickness_SC_IR = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_SC_IR[iCell] > 0.0):
            iceThickness_SC_IR[iCell] = iceVolume_SC_IR[iCell] / iceArea_SC_IR[iCell]

    print(np.amin(iceArea_SC_Initial), np.amax(iceArea_SC_Initial), np.amin(iceVolume_SC_Initial), np.amax(iceVolume_SC_Initial))
    print(np.amin(iceArea_SC_IR), np.amax(iceArea_SC_IR), np.amin(iceVolume_SC_IR), np.amax(iceVolume_SC_IR))

    fileMPAS.close()

    # slotted cylinder upwind
    filename = "./output_upwind_%s_%s/output.2000.nc" %(experiment2,res)

    print(filename)
    fileMPAS = Dataset(filename,"r")

    iceArea_SC_upwind   = fileMPAS.variables["iceAreaCategory"][iTime,:,0,0]
    iceVolume_SC_upwind = fileMPAS.variables["iceVolumeCategory"][iTime,:,0,0]
    iceThickness_SC_upwind = np.zeros(nCells)
    for iCell in range(0,nCells):
        if (iceArea_SC_upwind[iCell] > 0.0):
            iceThickness_SC_upwind[iCell] = iceVolume_SC_upwind[iCell] / iceArea_SC_upwind[iCell]

    fileMPAS.close()

    lons, yiceArea_CB_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_CB_Initial)
    lons, yiceArea_CB_IR      = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_CB_IR)
    lons, yiceArea_CB_upwind  = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_CB_upwind)

    lons, yiceArea_SC_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_SC_Initial)
    lons, yiceArea_SC_IR      = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_SC_IR)
    lons, yiceArea_SC_upwind  = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceArea_SC_upwind)

    lons, yiceThickness_CB_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_CB_Initial)
    lons, yiceThickness_CB_IR      = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_CB_IR)
    lons, yiceThickness_CB_upwind  = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_CB_upwind)

    lons, yiceThickness_SC_Initial = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_SC_Initial)
    lons, yiceThickness_SC_IR      = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_SC_IR)
    lons, yiceThickness_SC_upwind  = get_equatorial_array(nCells, latCell, lonCell, xCell, yCell, zCell, cellsOnCell, nEdgesOnCell, iceThickness_SC_upwind)

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(1, 2)
    fig.set_size_inches(7.2, 3)
    #fig, axes = plt.subplots(2, 2)
    #fig.set_size_inches(7.2, 5)

    # area
    axes[0].plot(lons, yiceArea_CB_Initial, linewidth=0.5, color="black", ls="-")
    axes[0].plot(lons, yiceArea_CB_IR,      linewidth=0.5, color="black", ls="--")
    axes[0].plot(lons, yiceArea_CB_upwind,  linewidth=0.5, color="black", ls=":")

    axes[0].legend(["Initial","IR","Upwind"],frameon=False,fontsize=8)
    axes[0].set_xlabel("Longitude (radians)")
    axes[0].set_ylabel("Equatorial ice concentration")
    axes[0].set_xlim([0.5, math.pi-0.5])
    axes[0].text(0.07, 0.9, "(a)", verticalalignment='bottom', horizontalalignment='right',transform=axes[0].transAxes, fontsize=8)


    axes[1].plot(lons, yiceArea_SC_Initial, linewidth=0.5, color="black", ls="-")
    axes[1].plot(lons, yiceArea_SC_IR,      linewidth=0.5, color="black", ls="--")
    axes[1].plot(lons, yiceArea_SC_upwind,  linewidth=0.5, color="black", ls=":")

    #axes[1].legend(["Initial","IR","Upwind"],frameon=False,fontsize=8)
    axes[1].set_xlabel("Longitude (radians)")
    #axes[1].set_ylabel("Equatorial ice concentration")
    axes[1].set_xlim([0.5, math.pi-0.5])
    axes[1].text(0.07, 0.9, "(b)", verticalalignment='bottom', horizontalalignment='right',transform=axes[1].transAxes, fontsize=8)




    # volume
    #axes[1,0].plot(lons, yiceThickness_CB_Initial, linewidth=0.5, color="black", ls="-")
    #axes[1,0].plot(lons, yiceThickness_CB_IR,      linewidth=0.5, color="black", ls="--")
    #axes[1,0].plot(lons, yiceThickness_CB_upwind,  linewidth=0.5, color="black", ls=":")

    #axes[1,0].legend(["Initial","IR","Upwind"],frameon=False,fontsize=8)
    #axes[1,0].set_xlabel("Longitude (radians)")
    #axes[1,0].set_ylabel("Equatorial ice thickness (m)")
    #axes[1,0].set_xlim([0.5, math.pi-0.5])
    #axes[1,0].text(0.07, 0.9, "(c)", verticalalignment='bottom', horizontalalignment='right',transform=axes[1,0].transAxes, fontsize=8)


    #axes[1,1].plot(lons, yiceThickness_SC_Initial, linewidth=0.5, color="black", ls="-")
    #axes[1,1].plot(lons, yiceThickness_SC_IR,      linewidth=0.5, color="black", ls="--")
    #axes[1,1].plot(lons, yiceThickness_SC_upwind,  linewidth=0.5, color="black", ls=":")

    #axes[1,0].legend(["Initial","IR","Upwind"],frameon=False,fontsize=8)
    #axes[1,1].set_xlabel("Longitude (radians)")
    #axes[1,1].set_ylabel("Equatorial ice volume (m)")
    #axes[1,1].set_xlim([0.5, math.pi-0.5])
    #axes[1,1].text(0.07, 0.9, "(d)", verticalalignment='bottom', horizontalalignment='right',transform=axes[1,1].transAxes, fontsize=8)


    plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
    #plt.savefig("advection_equatorial.eps")
    plt.savefig("advection_equatorial.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    advection_equatorial()
