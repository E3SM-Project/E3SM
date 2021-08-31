from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime, timedelta

#-------------------------------------------------------------------------------

def plot_thicknesses():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    timesPlot = ["2000-01-02_00:00:00",
                 "2000-01-06_00:00:00",
                 "2000-01-31_00:00:00"]

    labels = {"2000-01-02_00:00:00":"Day 1",
              "2000-01-06_00:00:00":"Day 5",
              "2000-01-31_00:00:00":"Day 30"}

    lineColours = {"2000-01-02_00:00:00":"black",
                   "2000-01-06_00:00:00":"black",
                   "2000-01-31_00:00:00":"black"}

    lineStyles = {"2000-01-02_00:00:00":"dashdot",
                  "2000-01-06_00:00:00":"dashed",
                  "2000-01-31_00:00:00":"solid"}

    filein = Dataset("./output/output.2000.nc","r")

    nTimes = len(filein.dimensions["Time"])
    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])

    times = filein.variables["xtime"][:]

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    iceVolumeCell = filein.variables["iceVolumeCell"][:]
    icePressure = filein.variables["icePressure"][:]
    uVelocity = filein.variables["uVelocity"][:]

    filein.close()

    # thickness
    fig, axis = plt.subplots()

    yc = 245000.0

    for timePlot in timesPlot:
        for iTime in range(0,nTimes):

            timeStr = b''.join(times[iTime]).decode('UTF-8').replace('\x00','')
            if (timeStr == timePlot):

                filenameTxt = "thickness_%s.txt" %(timePlot)
                fileout = open(filenameTxt,"w")

                x = []
                h = []
                for iCell in range(0,nCells):
                    if (yCell[iCell] == yc):
                        x.append(xCell[iCell]/1000.0)
                        h.append(iceVolumeCell[iTime,iCell])
                        fileout.write("%f %f\n" %(xCell[iCell]/1000.0,iceVolumeCell[iTime,iCell]))

                fileout.close()

                axis.plot(x,h,label=labels[timePlot], color=lineColours[timePlot], ls=lineStyles[timePlot], markersize=5.0)

                break

    axis.legend(frameon=False, handlelength=4)
    axis.set_xlabel("Position (km)")
    axis.set_ylabel("Ice thickness (m)")
    axis.set_xlim((0.0,1000.0))
    axis.set_ylim((0.0,5.0))

    plt.tight_layout()
    plt.savefig("thicknesses.png", dpi=400)
    plt.cla()
    plt.close(fig)

    # pressure
    fig, axis = plt.subplots()

    yc = 245000.0

    for timePlot in timesPlot:
        for iTime in range(0,nTimes):

            timeStr = b''.join(times[iTime]).decode('UTF-8').replace('\x00','')
            if (timeStr == timePlot):

                x = []
                p = []
                for iCell in range(0,nCells):
                    if (yCell[iCell] == yc):
                        x.append(xCell[iCell]/1000.0)
                        p.append(icePressure[iTime,iCell])

                break

        axis.plot(x,p,label=labels[timePlot], color=lineColours[timePlot], ls=lineStyles[timePlot], markersize=5.0)

    axis.legend(frameon=False, handlelength=4)
    axis.set_xlabel("Position (km)")
    axis.set_ylabel("Pressure")
    axis.set_xlim((0.0,1000.0))
    axis.set_ylim((0.0,None))

    plt.tight_layout()
    plt.savefig("pressure.png", dpi=400)
    plt.cla()
    plt.close(fig)

    # velocity
    fig, axis = plt.subplots()

    yv = 240000.0

    for timePlot in timesPlot:
        for iTime in range(0,nTimes):

            timeStr = b''.join(times[iTime]).decode('UTF-8').replace('\x00','')
            if (timeStr == timePlot):

                x = []
                v = []
                for iVertex in range(0,nVertices):
                    if (yVertex[iVertex] == yv):
                        x.append(xVertex[iVertex]/1000.0)
                        v.append(uVelocity[iTime,iVertex])

                break

        axis.plot(x,v,label=labels[timePlot], color=lineColours[timePlot], ls=lineStyles[timePlot], markersize=5.0)

    axis.legend(frameon=False, handlelength=4)
    axis.set_xlabel("Position (km)")
    axis.set_ylabel("Velocity (m/s)")
    axis.set_xlim((0.0,1000.0))
    axis.set_ylim((0.0,None))

    plt.tight_layout()
    plt.savefig("velocity.png", dpi=400)
    plt.cla()
    plt.close(fig)

    # mean ice position vs time
    initialIcePosition = 500000.0

    filein = Dataset("./output/output.2000.nc","r")

    nTimes = len(filein.dimensions["Time"])
    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])

    times = filein.variables["xtime"][:]

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    iceVolumeCell = filein.variables["iceVolumeCell"][:]
    uVelocity = filein.variables["uVelocity"][:]

    filein.close()


    Time0 = b''.join(times[0,0:19]).decode("utf-8").rstrip('\x00')
    Time0 = datetime.strptime(Time0, '%Y-%m-%d_%H:%M:%S')


    timesDays = []
    meanIcePositions = []
    velocities = []

    for iTime in range(0,nTimes):

        Time = b''.join(times[iTime,0:19]).decode("utf-8").rstrip('\x00')
        Time = datetime.strptime(Time, '%Y-%m-%d_%H:%M:%S')

        TimeDiff = Time - Time0
        timesDays.append(TimeDiff / timedelta(days=1))

        meanIcePosition = 0.0
        denominator = 0.0
        for iCell in range(0,nCells):
            if (yCell[iCell] == yc):
                meanIcePosition += xCell[iCell] * iceVolumeCell[iTime,iCell]
                denominator += iceVolumeCell[iTime,iCell]

        meanIcePosition /= denominator
        meanIcePosition -= initialIcePosition
        meanIcePosition /= 1000.0

        meanIcePositions.append(meanIcePosition)


    fileout = open("meanIcePositions.txt","w")
    for timeDay, meanIcePosition in zip(timesDays,meanIcePositions):
        fileout.write("%f %f\n" %(timeDay, meanIcePosition))
    fileout.close()

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_thicknesses()
