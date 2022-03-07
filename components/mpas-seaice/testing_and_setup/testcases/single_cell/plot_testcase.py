from netCDF4 import Dataset
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def plot_testcase():

    cm = 1/2.54  # centimeters in inches
    plt.rc('font', family="Times New Roman")
    plt.rc('mathtext',fontset="stix")
    SMALL_SIZE = 8
    MEDIUM_SIZE = 8
    BIGGER_SIZE = 8
    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    fig, axis = plt.subplots(figsize=(8*cm,7*cm))

    filein = Dataset("./output/output.2000.nc","r")

    iceAreaCell = filein.variables["iceAreaCell"][:]
    iceVolumeCell = filein.variables["iceVolumeCell"][:]
    snowVolumeCell = filein.variables["snowVolumeCell"][:]
    surfaceTemperatureCell = filein.variables["surfaceTemperatureCell"][:]

    filein.close()
    
    axis.plot(surfaceTemperatureCell,color="green")
    axis.set_ylabel("Temperature (C)")
    axis.set_xlabel("Time step")
    axis.set_xlim(0,4700)
    axis.set_ylim(None,0)
    axis.set_title("MPAS_Seaice single cell")

    axis2 = axis.twinx()

    axis2.plot(iceVolumeCell,color="red")
    axis2.plot(snowVolumeCell,color="blue")
    axis2.set_ylabel("Thickness (m)")
    axis2.set_ylim(0,None)

    plt.tight_layout()
    plt.savefig("single_cell.eps")
    plt.savefig("single_cell.png",dpi=300)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_testcase()
