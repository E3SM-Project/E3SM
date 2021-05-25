from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import math, sys
import numpy as np

#---------------------------------------------------------------

def cm2inch(value):
    return value/2.54

#------------------------------------------------------------

def generate_ellipse(nEllipse):

    xEllipse = []
    yEllipse = []

    theta = math.pi / 4.0

    a = math.sqrt(2.0) * 0.5
    b = a / 2.0

    for i in range(0,nEllipse+1):

        t = (i / float(nEllipse)) * 2.0 * math.pi

        x = a * math.cos(t)
        y = b * math.sin(t)

        xp = x * math.cos(theta) - y * math.sin(theta)
        yp = x * math.sin(theta) + y * math.cos(theta)
    
        xp = xp - 0.5
        yp = yp - 0.5

        xEllipse.append(xp)
        yEllipse.append(yp)

    return xEllipse, yEllipse

#------------------------------------------------------------

def plot_yield_curve(axes, sig1, sig2, title, subfigureLabel, addAxisLabels=True):

    axes.scatter(sig1, sig2, s=0.01, color="black", zorder=5, rasterized=False, marker=',')
    #axes.scatter(sig1, sig2, s=0.1, color="black", zorder=5, rasterized=True, marker=',')
    
    axes.plot((0.0,0.0), (-1.2,0.2), '--', color="grey", zorder=2, linewidth=0.5)
    axes.plot((-1.2,0.2), (0.0,0.0), '--', color="grey", zorder=3, linewidth=0.5)
    
    xEllipse, yEllipse = generate_ellipse(100)

    axes.tick_params('both', length=3)
    
    axes.plot(xEllipse,yEllipse, ':', color="grey", zorder=4, linewidth=0.5)

    axes.set_xlim([-1.2,0.2])
    axes.set_ylim([-1.2,0.2])

    if (not addAxisLabels):
        axes.get_xaxis().set_ticklabels([])
        axes.get_yaxis().set_ticklabels([])

    axes.set_aspect('equal')

    axes.set_title(title, fontsize=8)

    if (addAxisLabels):
        axes.set_xlabel(r"$\sigma_1$", fontsize=8)
        axes.set_ylabel(r"$\sigma_2$", fontsize=8)

    axes.text(0.12, 0.9, subfigureLabel, verticalalignment='bottom', horizontalalignment='right',transform=axes.transAxes, fontsize=8)
        
#------------------------------------------------------------

def get_cice_stresses(filename):

    fileCICE = Dataset(filename, "r")

    sig1in = fileCICE.variables["sig1"][0,:,:]
    sig2in = fileCICE.variables["sig2"][0,:,:]

    fileCICE.close()

    sig1 = []
    sig2 = []
    
    for j in range(0,84):
        for i in range(0,84):
            if (sig1in[j,i] < 1e20):
                sig1.append(sig1in[j,i])
                sig2.append(sig2in[j,i])

    return sig1, sig2

#------------------------------------------------------------

def get_mpas_variational_stresses(filename, iTime):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])

    nEdgesOnCell = fileMPAS.variables["nEdgesOnCell"][:]
    
    principalStress1 = fileMPAS.variables["principalStress1Var"][iTime,0:nCells-1,:]
    principalStress2 = fileMPAS.variables["principalStress2Var"][iTime,0:nCells-1,:]

    sig1 = []
    sig2 = []

    for iCell in range(0,nCells-1):
        for iVertexOnCell in range(0,nEdgesOnCell[iCell]):
            if (principalStress1[iCell,iVertexOnCell] < 1e20):
                sig1.append(principalStress1[iCell,iVertexOnCell])
                sig2.append(principalStress2[iCell,iVertexOnCell])
        
    fileMPAS.close()

    return sig1, sig2

#------------------------------------------------------------

def get_mpas_weak_stresses(filename, iTime):

    fileMPAS = Dataset(filename, "r")

    nCells = len(fileMPAS.dimensions["nCells"])
    
    principalStress1 = fileMPAS.variables["principalStress1Weak"][iTime,0:nCells-1]
    principalStress2 = fileMPAS.variables["principalStress2Weak"][iTime,0:nCells-1]

    sig1 = []
    sig2 = []

    for iCell in range(0,nCells-1):
        if (principalStress1[iCell] < 1e20):
            sig1.append(principalStress1[iCell])
            sig2.append(principalStress2[iCell])
        
    fileMPAS.close()

    return sig1, sig2

#------------------------------------------------------------

iTime = 3

mpl.rc('font',family='Times New Roman', size=8)
#mpl.rc('text', usetex=True)
mpl.rcParams['axes.linewidth'] = 0.5

sig1CICE,         sig2CICE         = get_cice_stresses("/home/akt/Work/CICE/rundirs/rundir_square/rundir/history/iceh_inst.1997-01-01-14400.nc")

sig1MPASQuadWach, sig2MPASQuadWach = get_mpas_variational_stresses("./output_quad_wachspress/output.2000.nc", iTime)
sig1MPASQuadPWL,  sig2MPASQuadPWL  = get_mpas_variational_stresses("./output_quad_pwl/output.2000.nc", iTime)
sig1MPASQuadWeak, sig2MPASQuadWeak = get_mpas_weak_stresses("./output_quad_weak/output.2000.nc", iTime)

sig1MPASHexWach,  sig2MPASHexWach = get_mpas_variational_stresses("./output_hex_wachspress/output.2000.nc", iTime)
sig1MPASHexPWL,   sig2MPASHexPWL  = get_mpas_variational_stresses("./output_hex_pwl/output.2000.nc", iTime)
sig1MPASHexWeak,  sig2MPASHexWeak = get_mpas_weak_stresses("./output_hex_weak/output.2000.nc", iTime)

fig, axes = plt.subplots(2, 4)
fig.set_size_inches(7.48031, 3.5)

plot_yield_curve(axes[0,0], sig1CICE,         sig2CICE,         "CICE",                 "(a)")

plot_yield_curve(axes[0,1], sig1MPASQuadWach, sig2MPASQuadWach, "Quad MPAS Wachspress", "(b)", False)
plot_yield_curve(axes[0,2], sig1MPASQuadPWL,  sig2MPASQuadPWL,  "Quad MPAS PWL",        "(c)", False)
plot_yield_curve(axes[0,3], sig1MPASQuadWeak, sig2MPASQuadWeak, "Quad MPAS Weak",       "(d)", False)

plot_yield_curve(axes[1,1], sig1MPASHexWach,  sig2MPASHexWach,  "Hex MPAS Wachspress",  "(e)", False)
plot_yield_curve(axes[1,2], sig1MPASHexPWL,   sig2MPASHexPWL,   "Hex MPAS PWL",         "(f)", False)
plot_yield_curve(axes[1,3], sig1MPASHexWeak,  sig2MPASHexWeak,  "Hex MPAS Weak",        "(g)", False)

axes[1,0].get_xaxis().set_visible(False)
axes[1,0].get_yaxis().set_visible(False)
axes[1,0].axis('off')


plt.tight_layout(pad=0.2, w_pad=0.2, h_pad=0.2)
plt.savefig("principal_stresses.eps")
#plt.savefig("principal_stresses.pdf", dpi=400)
#plt.savefig("principal_stresses.png",bbox_inches="tight",dpi=2000)






