from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

#-------------------------------------------------------------------------------

def stress_divergence_hist():

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
    fileWach = Dataset("./output_wachspress_40962/output.2000.nc","r")

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("Wachs:     ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./output_pwl_40962/output.2000.nc","r")

    stressDivergenceUPWL = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWL = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLDiff = (stressDivergenceUPWL - stressDivergenceUAnalytical)
    stressDivergenceVPWLDiff = (stressDivergenceVPWL - stressDivergenceVAnalytical)

    print("PWL:       ",
          np.amin(stressDivergenceUPWLDiff), np.amax(stressDivergenceUPWLDiff),
          np.amin(stressDivergenceVPWLDiff), np.amax(stressDivergenceVPWLDiff))

    filePWL.close()

    # Wachspress alt
    fileWachAlt = Dataset("./output_wachspress_alt_40962/output.2000.nc","r")

    stressDivergenceUWachAlt = fileWachAlt.variables["stressDivergenceU"][0,:]
    stressDivergenceVWachAlt = fileWachAlt.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachAltDiff = (stressDivergenceUWachAlt - stressDivergenceUAnalytical)
    stressDivergenceVWachAltDiff = (stressDivergenceVWachAlt - stressDivergenceVAnalytical)

    print("Wachs Alt: ",
          np.amin(stressDivergenceUWachAltDiff), np.amax(stressDivergenceUWachAltDiff),
          np.amin(stressDivergenceVWachAltDiff), np.amax(stressDivergenceVWachAltDiff))

    fileWachAlt.close()

    # PWL alt
    filePWLAlt = Dataset("./output_pwl_alt_40962/output.2000.nc","r")

    stressDivergenceUPWLAlt = filePWLAlt.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWLAlt = filePWLAlt.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLAltDiff = (stressDivergenceUPWLAlt - stressDivergenceUAnalytical)
    stressDivergenceVPWLAltDiff = (stressDivergenceVPWLAlt - stressDivergenceVAnalytical)

    print("PWL Alt:   ",
          np.amin(stressDivergenceUPWLAltDiff), np.amax(stressDivergenceUPWLAltDiff),
          np.amin(stressDivergenceVPWLAltDiff), np.amax(stressDivergenceVPWLAltDiff))

    filePWLAlt.close()

    # Weak
    fileWeak = Dataset("./output_weak_40962/output.2000.nc","r")

    stressDivergenceUWeak = fileWeak.variables["stressDivergenceU"][0,:]
    stressDivergenceVWeak = fileWeak.variables["stressDivergenceV"][0,:]

    stressDivergenceUWeakDiff = (stressDivergenceUWeak - stressDivergenceUAnalytical)
    stressDivergenceVWeakDiff = (stressDivergenceVWeak - stressDivergenceVAnalytical)

    print("Weak:      ",
          np.amin(stressDivergenceUWeakDiff), np.amax(stressDivergenceUWeakDiff),
          np.amin(stressDivergenceVWeakDiff), np.amax(stressDivergenceVWeakDiff))

    fileWeak.close()


    # histograms
    stressDivergenceUWachDiffHist    = []
    stressDivergenceUPWLDiffHist     = []
    stressDivergenceUWachAltDiffHist = []
    stressDivergenceUPWLAltDiffHist  = []
    stressDivergenceUWeakDiffHist    = []

    for iVertex in range(0,nVertices):

        if (latVertex[iVertex] > math.radians(20.0)):

            stressDivergenceUWachDiffHist.append(math.fabs(stressDivergenceUWachDiff[iVertex]))
            stressDivergenceUPWLDiffHist.append(math.fabs(stressDivergenceUPWLDiff[iVertex]))
            stressDivergenceUWachAltDiffHist.append(math.fabs(stressDivergenceUWachAltDiff[iVertex]))
            stressDivergenceUPWLAltDiffHist.append(math.fabs(stressDivergenceUPWLAltDiff[iVertex]))
            stressDivergenceUWeakDiffHist.append(math.fabs(stressDivergenceUWeakDiff[iVertex]))


    mpl.rc('text', usetex=True)
    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rcParams['axes.linewidth'] = 0.5

    plt.figure(figsize=(3.74016, 3))

    plt.hist(stressDivergenceUWachDiffHist,    50, range=[0.0,1.0], histtype='step', lw=1, color='red', label='Wachspress')
    plt.hist(stressDivergenceUPWLDiffHist,     50, range=[0.0,1.0], histtype='step', lw=1, color='blue',  label='PWL')
    plt.hist(stressDivergenceUWachAltDiffHist, 50, range=[0.0,1.0], histtype='step', lw=1, color='darkturquoise', label='Wachs. Alt')
    plt.hist(stressDivergenceUPWLAltDiffHist,  50, range=[0.0,1.0], histtype='step', lw=1, color='darkorange',  label='PWL Alt')
    plt.hist(stressDivergenceUWeakDiffHist,    50, range=[0.0,1.0], histtype='step', lw=1, color='green',   label='Weak')

    plt.yscale('log', nonpositive='clip')

    plt.xlabel("Error")
    plt.ylabel("Frequency")

    plt.legend(["Wachs.","PWL","Wachs. Alt","PWL Alt","Weak"], frameon=False, fontsize=8)

    plt.xlim([0,1.0])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)

    plt.savefig("stress_divergence_hist.png",dpi=400)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    stress_divergence_hist()
