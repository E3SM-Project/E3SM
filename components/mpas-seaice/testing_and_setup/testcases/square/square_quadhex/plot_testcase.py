from netCDF4 import Dataset
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
import numpy as np
from mpl_toolkits.axes_grid1 import make_axes_locatable

#-------------------------------------------------------------------------------

def plot_testcase_diff(testDir, baseDir, strainType, filenameOut):

    if (strainType not in ["varvar","weakvar"]):
        raise RuntimeError("Unknown strainType: %s" %(strainType))

    iTime = -1

    # base dir
    filein = Dataset("./%s/output.2000.nc" %(baseDir),"r")

    uVelocityBase = filein.variables["uVelocity"][iTime,:]
    vVelocityBase = filein.variables["vVelocity"][iTime,:]

    stressDivergenceUBase = filein.variables["stressDivergenceU"][iTime,:]
    stressDivergenceVBase = filein.variables["stressDivergenceV"][iTime,:]

    if (strainType == "varvar"):

        principalStress1Base = filein.variables["principalStress1Var"][iTime,:,2]
        principalStress2Base = filein.variables["principalStress2Var"][iTime,:,2]

        strain11varBase = filein.variables["strain11var"][iTime,:]
        strain22varBase = filein.variables["strain22var"][iTime,:]
        strain12varBase = filein.variables["strain12var"][iTime,:]

    filein.close()

    # test dir
    filein = Dataset("./%s/output.2000.nc" %(testDir),"r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])
    vertexDegree = len(filein.dimensions["vertexDegree"])
    maxEdges = len(filein.dimensions["maxEdges"])

    cellsOnVertex = filein.variables["cellsOnVertex"][:]
    cellsOnVertex -= 1

    verticesOnCell = filein.variables["verticesOnCell"][:]
    verticesOnCell -= 1

    edgesOnCell = filein.variables["edgesOnCell"][:]
    edgesOnCell -= 1

    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]

    xEdge = filein.variables["xEdge"][:]
    yEdge = filein.variables["yEdge"][:]

    uVelocityTest = filein.variables["uVelocity"][iTime,:]
    vVelocityTest = filein.variables["vVelocity"][iTime,:]

    stressDivergenceUTest = filein.variables["stressDivergenceU"][iTime,:]
    stressDivergenceVTest = filein.variables["stressDivergenceV"][iTime,:]

    if (strainType == "varvar"):

        principalStress1 = filein.variables["principalStress1Var"][iTime,:,2]
        principalStress2 = filein.variables["principalStress2Var"][iTime,:,2]

        strain11varTest = filein.variables["strain11var"][iTime,:]
        strain22varTest = filein.variables["strain22var"][iTime,:]
        strain12varTest = filein.variables["strain12var"][iTime,:]

    elif (strainType == "weakvar"):

        strain11weakTest = filein.variables["strain11weak"][iTime,:]
        strain22weakTest = filein.variables["strain22weak"][iTime,:]
        strain12weakTest = filein.variables["strain12weak"][iTime,:]

        principalStress1 = filein.variables["principalStress1Weak"][iTime,:]
        principalStress2 = filein.variables["principalStress2Weak"][iTime,:]

    filein.close()

    # differences
    uVelocity = uVelocityTest[:]# - uVelocityBase[:]
    vVelocity = vVelocityTest[:]# - vVelocityBase[:]

    stressDivergenceU = stressDivergenceUTest[:]# - stressDivergenceUBase[:]
    stressDivergenceV = stressDivergenceVTest[:]# - stressDivergenceVBase[:]

    if (strainType == "varvar"):
        strain11var = strain11varTest[:] - strain11varBase[:]
        strain22var = strain22varTest[:] - strain22varBase[:]
        strain12var = strain12varTest[:] - strain12varBase[:]
    elif (strainType == "weakvar"):
        strain11weak = strain11weakTest[:]
        strain22weak = strain22weakTest[:]
        strain12weak = strain12weakTest[:]


    xmin = np.amin(xVertex) - 10000.0
    xmax = np.amax(xVertex) + 10000.0
    ymin = np.amin(yVertex) - 10000.0
    ymax = np.amax(yVertex) + 10000.0

    # get vertex patches list
    patchesVertex = []
    deletedVertices = []
    for iVertex in range(0,nVertices):
        vertices = []
        useVertex = True
        for iCellOnVertex in range(0,vertexDegree):
            iCell = cellsOnVertex[iVertex,iCellOnVertex]
            if (iCell != -1):
                vertices.append((xCell[iCell],yCell[iCell]))
            else:
                useVertex = False
        if (useVertex):
            patchesVertex.append(Polygon(vertices,True))
        else:
            deletedVertices.append(iVertex)
    deletedVertices = np.array(deletedVertices)

    # get strain patches list
    patchesStrain = []
    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,maxEdges):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            iEdgeOnCell1 = (iVertexOnCell - 1)
            if (iEdgeOnCell1 < 0): iEdgeOnCell1 += maxEdges
            iEdgeOnCell2 = iVertexOnCell
            #iEdgeOnCell1 = iVertexOnCell
            #iEdgeOnCell2 = iVertexOnCell + 1
            #if (iEdgeOnCell2 > maxEdges-1): iEdgeOnCell2 -= maxEdges
            iEdge1 = edgesOnCell[iCell,iEdgeOnCell1]
            iEdge2 = edgesOnCell[iCell,iEdgeOnCell2]
            vertices = [(xCell[iCell],yCell[iCell]),
                        (xEdge[iEdge1],yEdge[iEdge1]),
                        (xVertex[iVertex],yVertex[iVertex]),
                        (xEdge[iEdge2],yEdge[iEdge2])]
            patchesStrain.append(Polygon(vertices,True))

    # patches cell
    patchesCell = []
    for iCell in range(0,nCells):
        vertices = []
        for iVertexOnCell in range(0,maxEdges):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            vertices.append((xVertex[iVertex],yVertex[iVertex]))
        patchesCell.append(Polygon(vertices,True))

    uVelocity = np.delete(uVelocity, deletedVertices)
    vVelocity = np.delete(vVelocity, deletedVertices)
    stressDivergenceU = np.delete(stressDivergenceU, deletedVertices)
    stressDivergenceV = np.delete(stressDivergenceV, deletedVertices)
    principalStress1 = np.where(principalStress1 < 1e20, principalStress1, 0.0)
    principalStress2 = np.where(principalStress2 < 1e20, principalStress2, 0.0)

    # patch collections
    pcUVelocity = PatchCollection(patchesVertex, cmap=plt.get_cmap("jet"))
    pcUVelocity.set_array(uVelocity)

    pcVVelocity = PatchCollection(patchesVertex, cmap=plt.get_cmap("jet"))
    pcVVelocity.set_array(vVelocity)

    pcStressDivergenceU = PatchCollection(patchesVertex, cmap=plt.get_cmap("jet"))
    pcStressDivergenceU.set_array(stressDivergenceU)

    pcStressDivergenceV = PatchCollection(patchesVertex, cmap=plt.get_cmap("jet"))
    pcStressDivergenceV.set_array(stressDivergenceV)

    strainMin = None#-5e-7
    strainMax = None# 5e-7

    if (strainType == "varvar"):

        pcStrain11 = PatchCollection(patchesStrain, cmap=plt.get_cmap("jet"))
        pcStrain11.set_array(strain11var.flatten())
        pcStrain11.set_clim([strainMin, strainMax])

        pcStrain22 = PatchCollection(patchesStrain, cmap=plt.get_cmap("jet"))
        pcStrain22.set_array(strain22var.flatten())
        pcStrain22.set_clim([strainMin, strainMax])

        pcStrain12 = PatchCollection(patchesStrain, cmap=plt.get_cmap("jet"))
        pcStrain12.set_array(strain12var.flatten())
        pcStrain12.set_clim([strainMin, strainMax])

    elif (strainType == "weakvar"):

        pcStrain11 = PatchCollection(patchesCell, cmap=plt.get_cmap("jet"))
        pcStrain11.set_array(strain11weak.flatten())
        pcStrain11.set_clim([strainMin, strainMax])

        pcStrain22 = PatchCollection(patchesCell, cmap=plt.get_cmap("jet"))
        pcStrain22.set_array(strain22weak.flatten())
        pcStrain22.set_clim([strainMin, strainMax])

        pcStrain12 = PatchCollection(patchesCell, cmap=plt.get_cmap("jet"))
        pcStrain12.set_array(strain12weak.flatten())
        pcStrain12.set_clim([strainMin, strainMax])


    # plot
    fig, axes = plt.subplots(2,4,figsize=(10,5))

    axes[0,0].add_collection(pcUVelocity)
    axes[0,0].set_xlim((xmin,xmax))
    axes[0,0].set_ylim((ymin,ymax))
    axes[0,0].set_aspect('equal')
    axes[0,0].set_xticks([])
    axes[0,0].set_yticks([])
    axes[0,0].set_title("uVelocity")
    divider = make_axes_locatable(axes[0,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcUVelocity, cax=cax)

    axes[0,1].add_collection(pcVVelocity)
    axes[0,1].set_xlim((xmin,xmax))
    axes[0,1].set_ylim((ymin,ymax))
    axes[0,1].set_aspect('equal')
    axes[0,1].set_xticks([])
    axes[0,1].set_yticks([])
    axes[0,1].set_title("vVelocity")
    divider = make_axes_locatable(axes[0,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcVVelocity, cax=cax)

    axes[0,2].add_collection(pcStressDivergenceU)
    axes[0,2].set_xlim((xmin,xmax))
    axes[0,2].set_ylim((ymin,ymax))
    axes[0,2].set_aspect('equal')
    axes[0,2].set_xticks([])
    axes[0,2].set_yticks([])
    axes[0,2].set_title("stressDivU")
    divider = make_axes_locatable(axes[0,2])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcStressDivergenceU, cax=cax)

    axes[0,3].add_collection(pcStressDivergenceV)
    axes[0,3].set_xlim((xmin,xmax))
    axes[0,3].set_ylim((ymin,ymax))
    axes[0,3].set_aspect('equal')
    axes[0,3].set_xticks([])
    axes[0,3].set_yticks([])
    axes[0,3].set_title("stressDivV")
    divider = make_axes_locatable(axes[0,3])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcStressDivergenceV, cax=cax)

    axes[1,0].add_collection(pcStrain11)
    axes[1,0].set_xlim((xmin,xmax))
    axes[1,0].set_ylim((ymin,ymax))
    axes[1,0].set_aspect('equal')
    axes[1,0].set_xticks([])
    axes[1,0].set_yticks([])
    axes[1,0].set_title("strain11")
    divider = make_axes_locatable(axes[1,0])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcStrain11, cax=cax)

    axes[1,1].add_collection(pcStrain22)
    axes[1,1].set_xlim((xmin,xmax))
    axes[1,1].set_ylim((ymin,ymax))
    axes[1,1].set_aspect('equal')
    axes[1,1].set_xticks([])
    axes[1,1].set_yticks([])
    axes[1,1].set_title("strain22")
    divider = make_axes_locatable(axes[1,1])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcStrain22, cax=cax)

    axes[1,2].add_collection(pcStrain12)
    axes[1,2].set_xlim((xmin,xmax))
    axes[1,2].set_ylim((ymin,ymax))
    axes[1,2].set_aspect('equal')
    axes[1,2].set_xticks([])
    axes[1,2].set_yticks([])
    axes[1,2].set_title("strain12")
    divider = make_axes_locatable(axes[1,2])
    cax = divider.append_axes('right', size='5%', pad=0.05)
    fig.colorbar(pcStrain12, cax=cax)

    axes[1,3].scatter(principalStress1, principalStress2, s=1, c="k")
    axes[1,3].set_aspect('equal')
    axes[1,3].set_xticks([])
    axes[1,3].set_yticks([])



    #axes[1,1].axis('off')


    plt.tight_layout()
    plt.savefig(filenameOut,dpi=1200)

#-------------------------------------------------------------------------------

def plot_testcase():

    plot_testcase_diff("output_hex_pwl_0082x0094_120",
                       "output_hex_wachspress_0082x0094_120",
                       "varvar",
                       "square_pwl_wachs_diff.png")

    plot_testcase_diff("output_hex_pwlavg_0082x0094_120",
                       "output_hex_wachsavg_0082x0094_120",
                       "varvar",
                       "square_pwlavg_wachsavg_diff.png")

    plot_testcase_diff("output_hex_weak_0082x0094_120",
                       "output_hex_wachspress_0082x0094_120",
                       "weakvar",
                       "square_weak_wachs_diff.png")

    plot_testcase_diff("output_hex_weakwachs_0082x0094_120",
                       "output_hex_wachspress_0082x0094_120",
                       "varvar",
                       "square_weakwachs_wachs_diff.png")

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    plot_testcase()
