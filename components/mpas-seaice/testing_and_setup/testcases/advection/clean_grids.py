from netCDF4 import Dataset
import os


nGrids = [2562,10242,40962,163842]

mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

for nGrid in nGrids:

    filein = Dataset("x1.%i.grid.nc" %(nGrid), "r")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])
    vertexDegree = len(filein.dimensions["vertexDegree"])

    xCell = filein.variables["xCell"][:]
    yCell = filein.variables["yCell"][:]
    zCell = filein.variables["zCell"][:]
    xVertex = filein.variables["xVertex"][:]
    yVertex = filein.variables["yVertex"][:]
    zVertex = filein.variables["zVertex"][:]
    cellsOnVertex = filein.variables["cellsOnVertex"][:]

    filein.close()

    fileout = Dataset("grid_in.%i.nc" %(nGrid),"w",format="NETCDF3_CLASSIC")

    fileout.ON_A_SPHERE = "YES"
    fileout.sphere_radius = 1.0

    fileout.createDimension("nCells",nCells)
    fileout.createDimension("nVertices",nVertices)
    fileout.createDimension("vertexDegree",vertexDegree)

    var = fileout.createVariable("xCell","d",dimensions=["nCells"])
    var[:] = xCell[:]
    var = fileout.createVariable("yCell","d",dimensions=["nCells"])
    var[:] = yCell[:]
    var = fileout.createVariable("zCell","d",dimensions=["nCells"])
    var[:] = zCell[:]

    var = fileout.createVariable("xVertex","d",dimensions=["nVertices"])
    var[:] = xVertex[:]
    var = fileout.createVariable("yVertex","d",dimensions=["nVertices"])
    var[:] = yVertex[:]
    var = fileout.createVariable("zVertex","d",dimensions=["nVertices"])
    var[:] = zVertex[:]

    var = fileout.createVariable("cellsOnVertex","i",dimensions=["nVertices","vertexDegree"])
    var[:] = cellsOnVertex[:]

    fileout.close()

    cmd = "%s/mesh_tools/mesh_conversion_tools/MpasMeshConverter.x grid_in.%i.nc grid.%i.nc" %(mpas_tools_dir, nGrid,nGrid)
    os.system(cmd)
    
