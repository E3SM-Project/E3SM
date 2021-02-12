from netCDF4 import Dataset
import numpy as np

velocityType = "xquadratic"
gridType = "hex"

operatorMethods = ["wachspress","pwl"]

for operatorMethod in operatorMethods:

    print(operatorMethod)

    # grid
    filegridname = "grid_var_%s.nc" %(gridType)

    filegrid = Dataset(filegridname,"r")

    nCells = len(filegrid.dimensions["nCells"])
    maxEdges = len(filegrid.dimensions["maxEdges"])
    verticesOnCell = filegrid.variables["verticesOnCell"][:]
    verticesOnCell[:] = verticesOnCell[:] - 1

    iVertexCenter = filegrid.iVertexCenter

    filegrid.close()

    # ic
    fileicname = "ic_var_%s_%s.nc" %(gridType, velocityType)

    fileic = Dataset(fileicname,"r")

    strain11VertexAnalytical = fileic.variables["strain11VertexAnalytical"][:]
    strain22VertexAnalytical = fileic.variables["strain22VertexAnalytical"][:]
    strain12VertexAnalytical = fileic.variables["strain12VertexAnalytical"][:]
    strain11CellAnalytical = fileic.variables["strain11CellAnalytical"][:]
    strain22CellAnalytical = fileic.variables["strain22CellAnalytical"][:]
    strain12CellAnalytical = fileic.variables["strain12CellAnalytical"][:]
    stressDivergenceUAnalytical = fileic.variables["stressDivergenceUAnalytical"][:]
    stressDivergenceVAnalytical = fileic.variables["stressDivergenceVAnalytical"][:]

    fileic.close()

    # sim
    fileinname = "./output_%s_%s_%s/output.2000.nc" %(operatorMethod, gridType, velocityType)

    filein = Dataset(fileinname,"r")

    strain11var = filein.variables["strain11var"][0,:]
    strain22var = filein.variables["strain22var"][0,:]
    strain12var = filein.variables["strain12var"][0,:]
    stressDivergenceU = filein.variables["stressDivergenceU"][0,:]
    stressDivergenceV = filein.variables["stressDivergenceV"][0,:]

    filein.close()

    # diffs
    strain11varDiff = np.zeros((nCells, maxEdges))
    strain22varDiff = np.zeros((nCells, maxEdges))
    strain12varDiff = np.zeros((nCells, maxEdges))

    for iCell in range(0,nCells):
        for iVertexOnCell in range(0,maxEdges):
            iVertex = verticesOnCell[iCell,iVertexOnCell]
            strain11varDiff[iCell,iVertexOnCell] = \
                strain11var[iCell,iVertexOnCell] - strain11VertexAnalytical[iVertex]
            strain22varDiff[iCell,iVertexOnCell] = \
                strain22var[iCell,iVertexOnCell] - strain22VertexAnalytical[iVertex]
            strain12varDiff[iCell,iVertexOnCell] = \
                strain12var[iCell,iVertexOnCell] - strain12VertexAnalytical[iVertex]

            #print("  ", iCell, iVertexOnCell, strain22varDiff[iCell,iVertexOnCell], \
            #      strain22var[iCell,iVertexOnCell], strain22VertexAnalytical[iVertex])

    #print("  strain11: ", np.mean(strain11varDiff), np.std(strain11varDiff))
    #print("  strain22: ", np.mean(strain22varDiff), np.std(strain22varDiff))
    #print("  strain12: ", np.mean(strain12varDiff), np.std(strain12varDiff))

    print(stressDivergenceU[iVertexCenter], stressDivergenceUAnalytical[iVertexCenter], \
          stressDivergenceU[iVertexCenter] - stressDivergenceUAnalytical[iVertexCenter])
    print(stressDivergenceV[iVertexCenter], stressDivergenceVAnalytical[iVertexCenter], \
          stressDivergenceV[iVertexCenter] - stressDivergenceVAnalytical[iVertexCenter])
