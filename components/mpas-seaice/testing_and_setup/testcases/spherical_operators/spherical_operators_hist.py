from netCDF4 import Dataset
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib as mpl

#-------------------------------------------------------------------------------

def strain_hist(axis):

    # grid
    fileGrid = Dataset("./strain/x1.40962.grid.nc","r")

    nCells = len(fileGrid.dimensions["nCells"])
    nVertices = len(fileGrid.dimensions["nVertices"])
    vertexDegree = len(fileGrid.dimensions["vertexDegree"])
    maxEdges = len(fileGrid.dimensions["maxEdges"])

    vertexDegreeArr = np.zeros(nVertices,dtype="i")
    vertexDegreeArr[:] = vertexDegree

    nEdgesOnCell = fileGrid.variables["nEdgesOnCell"][:]

    cellsOnVertex = fileGrid.variables["cellsOnVertex"][:]
    cellsOnVertex[:] = cellsOnVertex[:] - 1

    verticesOnCell = fileGrid.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    latVertex = fileGrid.variables["latVertex"][:]
    latCell = fileGrid.variables["latCell"][:]

    xCell = fileGrid.variables["xCell"][:]
    yCell = fileGrid.variables["yCell"][:]
    zCell = fileGrid.variables["zCell"][:]

    xVertex = fileGrid.variables["xVertex"][:]
    yVertex = fileGrid.variables["yVertex"][:]
    zVertex = fileGrid.variables["zVertex"][:]

    fileGrid.close()

    # ic
    fileIC = Dataset("./strain/ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    strain11VertexAnalytical = fileIC.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileIC.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileIC.variables["strain12VertexAnalytical"][:]

    strain11CellAnalytical = fileIC.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileIC.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileIC.variables["strain12CellAnalytical"][:]

    strain11VarAnalytical = np.zeros((nCells, maxEdges))
    strain22VarAnalytical = np.zeros((nCells, maxEdges))
    strain12VarAnalytical = np.zeros((nCells, maxEdges))
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            strain11VarAnalytical[iCell,iVertexOnCell] = strain11VertexAnalytical[iVertex]
            strain22VarAnalytical[iCell,iVertexOnCell] = strain22VertexAnalytical[iVertex]
            strain12VarAnalytical[iCell,iVertexOnCell] = strain12VertexAnalytical[iVertex]

    print("Strain: ",
          np.amin(strain11VertexAnalytical), np.amax(strain11VertexAnalytical),
          np.amin(strain22VertexAnalytical), np.amax(strain22VertexAnalytical),
          np.amin(strain12VertexAnalytical), np.amax(strain12VertexAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./strain/output_wachspress_40962/output.2000.nc","r")

    strain11varWachspress = fileWach.variables["strain11var"][0,:]
    strain22varWachspress = fileWach.variables["strain22var"][0,:]
    strain12varWachspress = fileWach.variables["strain12var"][0,:]

    strain11varWachspressDiff = strain11varWachspress - strain11VarAnalytical
    strain22varWachspressDiff = strain22varWachspress - strain22VarAnalytical
    strain12varWachspressDiff = strain12varWachspress - strain12VarAnalytical

    print("Wachs: ",
          np.amin(strain11varWachspressDiff), np.amax(strain11varWachspressDiff),
          np.amin(strain22varWachspressDiff), np.amax(strain22varWachspressDiff),
          np.amin(strain12varWachspressDiff), np.amax(strain12varWachspressDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./strain/output_pwl_40962/output.2000.nc","r")

    strain11varPWL = filePWL.variables["strain11var"][0,:]
    strain22varPWL = filePWL.variables["strain22var"][0,:]
    strain12varPWL = filePWL.variables["strain12var"][0,:]

    strain11varPWLDiff = strain11varPWL - strain11VarAnalytical
    strain22varPWLDiff = strain22varPWL - strain22VarAnalytical
    strain12varPWLDiff = strain12varPWL - strain12VarAnalytical

    print("PWL:   ",
          np.amin(strain11varPWLDiff), np.amax(strain11varPWLDiff),
          np.amin(strain22varPWLDiff), np.amax(strain22varPWLDiff),
          np.amin(strain12varPWLDiff), np.amax(strain12varPWLDiff))

    filePWL.close()

    # Weak
    fileWeak = Dataset("./strain/output_weak_40962/output.2000.nc","r")

    strain11weakWeak = fileWeak.variables["strain11weak"][0,:]
    strain22weakWeak = fileWeak.variables["strain22weak"][0,:]
    strain12weakWeak = fileWeak.variables["strain12weak"][0,:]

    strain11weakWeakDiff = strain11weakWeak - strain11CellAnalytical
    strain22weakWeakDiff = strain22weakWeak - strain22CellAnalytical
    strain12weakWeakDiff = strain12weakWeak - strain12CellAnalytical

    print("Weak:  ",
          np.amin(strain11weakWeakDiff), np.amax(strain11weakWeakDiff),
          np.amin(strain22weakWeakDiff), np.amax(strain22weakWeakDiff),
          np.amin(strain12weakWeakDiff), np.amax(strain12weakWeakDiff))

    fileWeak.close()

    # Weak Wachs
    fileWeakWachs = Dataset("./strain/output_weakwachs_40962/output.2000.nc","r")

    strain11weakWachs = fileWeakWachs.variables["strain11varAvg"][0,:]
    strain22weakWachs = fileWeakWachs.variables["strain22varAvg"][0,:]
    strain12weakWachs = fileWeakWachs.variables["strain12varAvg"][0,:]

    strain11weakWachsDiff = strain11weakWachs - strain11VertexAnalytical
    strain22weakWachsDiff = strain22weakWachs - strain22VertexAnalytical
    strain12weakWachsDiff = strain12weakWachs - strain12VertexAnalytical

    print("WeakWachs:  ",
          np.amin(strain11weakWachsDiff), np.amax(strain11weakWachsDiff),
          np.amin(strain22weakWachsDiff), np.amax(strain22weakWachsDiff),
          np.amin(strain12weakWachsDiff), np.amax(strain12weakWachsDiff))

    fileWeakWachs.close()


    # histograms
    strain11UWachDiffHist = []
    strain11UPWLDiffHist  = []
    strain11UWeakDiffHist = []

    # variational
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            iVertex = verticesOnCell[iCell,iVertexOnCell]

            if (latVertex[iVertex] > math.radians(20.0)):

                strain11UWachDiffHist.append(math.fabs(strain11varWachspressDiff[iCell,iVertexOnCell]))
                strain11UPWLDiffHist.append(math.fabs(strain11varPWLDiff[iCell,iVertexOnCell]))

    # weak
    for iCell in range(0,nCells):

        if (latCell[iCell] > math.radians(20.0)):

            strain11UWeakDiffHist.append(math.fabs(strain11weakWeakDiff[iCell]))

    # plot
    axis.hist(strain11UWachDiffHist, 50, range=[0.0,0.125], histtype='step', lw=1, color='red',   label='Wachs.')
    axis.hist(strain11UPWLDiffHist,  50, range=[0.0,0.125], histtype='step', lw=1, color='blue',  label='PWL')
    axis.hist(strain11UWeakDiffHist, 50, range=[0.0,0.125], histtype='step', lw=1, color='green', label='Weak')

    axis.set_yscale('log', nonpositive='clip')
    #axis.set_xscale('log', nonpositive='clip')

    axis.set_title(r'(a) $\dot{\epsilon}_{11}$')
    axis.set_xlabel("Error")
    axis.set_ylabel("Frequency")

    axis.legend(["Wachs.","PWL","Weak"], frameon=False, fontsize=8)

    axis.set_xlim([0,0.125])

#-------------------------------------------------------------------------------

def stress_divergence_hist(axis):

    # grid
    fileGrid = Dataset("./stress_divergence/x1.40962.grid.nc","r")

    nVertices = len(fileGrid.dimensions["nVertices"])

    latVertex = fileGrid.variables["latVertex"][:]

    fileGrid.close()


    # ic
    fileIC = Dataset("./stress_divergence/ic_40962.nc","r")

    uVelocity = fileIC.variables["uVelocity"][:]
    vVelocity = fileIC.variables["vVelocity"][:]

    stressDivergenceUAnalytical = fileIC.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileIC.variables["stressDivergenceVAnalytical"][:]

    print("Stress divergence: ",
          np.amin(stressDivergenceUAnalytical), np.amax(stressDivergenceUAnalytical),
          np.amin(stressDivergenceVAnalytical), np.amax(stressDivergenceVAnalytical))

    fileIC.close()

    # Wachspress
    fileWach = Dataset("./stress_divergence/output_wachspress_40962/output.2000.nc","r")

    stressDivergenceUWach = fileWach.variables["stressDivergenceU"][0,:]
    stressDivergenceVWach = fileWach.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachDiff = (stressDivergenceUWach - stressDivergenceUAnalytical)
    stressDivergenceVWachDiff = (stressDivergenceVWach - stressDivergenceVAnalytical)

    print("Wachs:     ",
          np.amin(stressDivergenceUWachDiff), np.amax(stressDivergenceUWachDiff),
          np.amin(stressDivergenceVWachDiff), np.amax(stressDivergenceVWachDiff))

    fileWach.close()

    # PWL
    filePWL = Dataset("./stress_divergence/output_pwl_40962/output.2000.nc","r")

    stressDivergenceUPWL = filePWL.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWL = filePWL.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLDiff = (stressDivergenceUPWL - stressDivergenceUAnalytical)
    stressDivergenceVPWLDiff = (stressDivergenceVPWL - stressDivergenceVAnalytical)

    print("PWL:       ",
          np.amin(stressDivergenceUPWLDiff), np.amax(stressDivergenceUPWLDiff),
          np.amin(stressDivergenceVPWLDiff), np.amax(stressDivergenceVPWLDiff))

    filePWL.close()

    # Wachspress alt
    fileWachAlt = Dataset("./stress_divergence/output_wachspress_alt_40962/output.2000.nc","r")

    stressDivergenceUWachAlt = fileWachAlt.variables["stressDivergenceU"][0,:]
    stressDivergenceVWachAlt = fileWachAlt.variables["stressDivergenceV"][0,:]

    stressDivergenceUWachAltDiff = (stressDivergenceUWachAlt - stressDivergenceUAnalytical)
    stressDivergenceVWachAltDiff = (stressDivergenceVWachAlt - stressDivergenceVAnalytical)

    print("Wachs Alt: ",
          np.amin(stressDivergenceUWachAltDiff), np.amax(stressDivergenceUWachAltDiff),
          np.amin(stressDivergenceVWachAltDiff), np.amax(stressDivergenceVWachAltDiff))

    fileWachAlt.close()

    # PWL alt
    filePWLAlt = Dataset("./stress_divergence/output_pwl_alt_40962/output.2000.nc","r")

    stressDivergenceUPWLAlt = filePWLAlt.variables["stressDivergenceU"][0,:]
    stressDivergenceVPWLAlt = filePWLAlt.variables["stressDivergenceV"][0,:]

    stressDivergenceUPWLAltDiff = (stressDivergenceUPWLAlt - stressDivergenceUAnalytical)
    stressDivergenceVPWLAltDiff = (stressDivergenceVPWLAlt - stressDivergenceVAnalytical)

    print("PWL Alt:   ",
          np.amin(stressDivergenceUPWLAltDiff), np.amax(stressDivergenceUPWLAltDiff),
          np.amin(stressDivergenceVPWLAltDiff), np.amax(stressDivergenceVPWLAltDiff))

    filePWLAlt.close()

    # Weak
    fileWeak = Dataset("./stress_divergence/output_weak_40962/output.2000.nc","r")

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


    # plot
    axis.hist(stressDivergenceUWachDiffHist,    50, range=[0.0,1.0], histtype='step', lw=1, color='red',           label='Wachs.')
    axis.hist(stressDivergenceUPWLDiffHist,     50, range=[0.0,1.0], histtype='step', lw=1, color='blue',          label='PWL')
    axis.hist(stressDivergenceUWachAltDiffHist, 50, range=[0.0,1.0], histtype='step', lw=1, color='darkturquoise', label='Wachs. Alt')
    axis.hist(stressDivergenceUPWLAltDiffHist,  50, range=[0.0,1.0], histtype='step', lw=1, color='darkorange',    label='PWL Alt')
    axis.hist(stressDivergenceUWeakDiffHist,    50, range=[0.0,1.0], histtype='step', lw=1, color='green',         label='Weak')

    axis.set_yscale('log', nonpositive='clip')

    axis.set_title(r'(b) $(\nabla \cdot \sigma)_u$')
    axis.set_xlabel("Error")
    axis.set_ylabel("Frequency")

    axis.legend(["Wachs.","PWL","Wachs. Alt","PWL Alt","Weak"], frameon=False, fontsize=8)

    axis.set_xlim([0,1.0])

#-------------------------------------------------------------------------------

def spherical_operators_hist():

    mpl.rc('font', family='Times New Roman', size=8)
    mpl.rc('text', usetex=True)
    mpl.rcParams['axes.linewidth'] = 0.5

    fig, axes = plt.subplots(1, 2, figsize=(7.3,3))

    strain_hist(axes[0])

    stress_divergence_hist(axes[1])

    plt.tight_layout(pad=0.5, w_pad=0.5, h_pad=0.5)
    plt.savefig("spherical_operators_hist.png",dpi=400)
    plt.savefig("spherical_operators_hist.eps")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    spherical_operators_hist()
