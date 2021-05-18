from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

#-------------------------------------------------------------------------------

def strain_stress_divergence_hist():

    # grid
    fileGrid = Dataset("x1.40962.grid.nc","r")

    nVertices = len(fileGrid.dimensions["nVertices"])

    latVertex = fileGrid.variables["latVertex"][:]

    fileGrid.close()


    # ic
    fileIC = Dataset("ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    print("Stress divergence: ",
          np.amin(stressDivergenceUAnalytical), np.amax(stressDivergenceUAnalytical),
          np.amin(stressDivergenceVAnalytical), np.amax(stressDivergenceVAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./output_wachspress_alt_40962/output.2000.nc","r")

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("Wachs: ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./output_pwl_alt_40962/output.2000.nc","r")

    stressDivergenceUPWL = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWL = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLDiff = (stressDivergenceUPWL - stressDivergenceUAnalytical)
    stressDivergenceVPWLDiff = (stressDivergenceVPWL - stressDivergenceVAnalytical)

    print("PWL:   ",
          np.amin(stressDivergenceUPWLDiff), np.amax(stressDivergenceUPWLDiff),
          np.amin(stressDivergenceVPWLDiff), np.amax(stressDivergenceVPWLDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./output_weak_40962/output.2000.nc","r")

    stressDivergenceUWeak = fileWeak.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeak = fileWeak.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakDiff = (stressDivergenceUWeak - stressDivergenceUAnalytical)
    stressDivergenceVWeakDiff = (stressDivergenceVWeak - stressDivergenceVAnalytical)

    print("Weak:  ",
          np.amin(stressDivergenceUWeakDiff), np.amax(stressDivergenceUWeakDiff),
          np.amin(stressDivergenceVWeakDiff), np.amax(stressDivergenceVWeakDiff))

    fileWeak.close()

    # WeakWachs
    fileWeakWachs = Dataset("./output_weakwachs_40962/output.2000.nc","r")

    stressDivergenceUWeakWachs = fileWeakWachs.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeakWachs = fileWeakWachs.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakWachsDiff = (stressDivergenceUWeakWachs - stressDivergenceUAnalytical)
    stressDivergenceVWeakWachsDiff = (stressDivergenceVWeakWachs - stressDivergenceVAnalytical)

    print("WeakWachs:  ",
          np.amin(stressDivergenceUWeakWachsDiff), np.amax(stressDivergenceUWeakWachsDiff),
          np.amin(stressDivergenceVWeakWachsDiff), np.amax(stressDivergenceVWeakWachsDiff))

    fileWeakWachs.close()


    # histograms
    stressDivergenceUWachDiffHist = []
    stressDivergenceUPWLDiffHist  = []
    stressDivergenceUWeakDiffHist = []

    maxValue = 0.0
    for iVertex in range(0,nVertices):

        if (latVertex[iVertex] > math.radians(20.0)):

            stressDivergenceUWachDiffHist.append(math.fabs(stressDivergenceUWachDiff[iVertex]))
            stressDivergenceUPWLDiffHist.append(math.fabs(stressDivergenceUPWLDiff[iVertex]))
            stressDivergenceUWeakDiffHist.append(math.fabs(stressDivergenceUWeakDiff[iVertex]))

            maxValue = max(math.fabs(stressDivergenceUWachDiff[iVertex]),maxValue)
            maxValue = max(math.fabs(stressDivergenceUPWLDiff[iVertex]),maxValue)


    mpl.rc('text', usetex=True)
    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rcParams['axes.linewidth'] = 0.5

    plt.figure(figsize=(3.74016, 3))

    plt.hist(stressDivergenceUWachDiffHist, 50, range=[0.0,1.0], histtype='step', lw=1, color='blue',  label='Wachspress')
    plt.hist(stressDivergenceUPWLDiffHist,  50, range=[0.0,1.0], histtype='step', lw=1, color='red',   label='PWL')
    plt.hist(stressDivergenceUWeakDiffHist, 50, range=[0.0,1.0], histtype='step', lw=1, color='green', label='Weak')

    plt.yscale('log', nonpositive='clip')

    plt.xlabel("Error")
    plt.ylabel("Frequency")

    plt.legend(["Wachspress","PWL","Weak"], frameon=False, fontsize=8)

    plt.xlim([0,1.0])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    plt.savefig("strain_stress_divergence_hist.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    strain_stress_divergence_hist()
