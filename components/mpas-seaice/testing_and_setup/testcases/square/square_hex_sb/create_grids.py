from netCDF4 import Dataset
import math
import os
import matplotlib.pyplot as plt

#-------------------------------------------------------------------------------

def create_grid(nx, ny, Lx, Ly, plot=False):

    mpas_tools_dir = os.environ['MPAS_TOOLS_DIR']

    dcx = Lx / float(nx-1)
    dcy = 0.5 * math.sqrt(3.0) * dcx

    print(dcx, dcy, nx, ny)

    gridname = "grid_%4.4i.nc" %(nx)

    fileout = open("namelist.input","w")

    line = "&periodic_grid\n"
    fileout.write(line)

    nxUse = nx + 4
    if ((nxUse % 2) != 0): nxUse += 1
    line = "   nx = %i,\n" %(nxUse)
    fileout.write(line)

    nyUse = ny + 4
    if ((nyUse % 2) != 0): nyUse += 1
    line = "   ny = %i,\n" %(nyUse)
    fileout.write(line)

    line = "   dc = %f,\n" %(dcx)
    fileout.write(line)

    line = "   nVertLevels = 1,\n"
    fileout.write(line)

    line = "   nTracers = 1,\n"
    fileout.write(line)

    line = "   nproc = 2, 4, 8,\n"
    fileout.write(line)

    line = "/\n"
    fileout.write(line)

    fileout.close()

    os.system("%s/mesh_tools/periodic_hex/periodic_grid" %(mpas_tools_dir))
    os.system("python %s/mesh_tools/periodic_hex/mark_periodic_boundaries_for_culling.py -f grid.nc" %(mpas_tools_dir))
    os.system("%s/mesh_tools/mesh_conversion_tools/MpasCellCuller.x grid.nc grid_culled.nc" %(mpas_tools_dir))

    dcyInt = 6*int(round(dcy / 6.0))

    filein = Dataset("grid_culled.nc","a")

    nCells = len(filein.dimensions["nCells"])
    nVertices = len(filein.dimensions["nVertices"])
    nEdges = len(filein.dimensions["nEdges"])

    filein.Lx = Lx
    filein.Ly = Ly
    #filein.Ly = float((ny-1) * dcyInt)

    filein.nx = nx
    filein.ny = ny

    filein.dcx = dcx
    filein.dcy = dcy
    #filein.dcy = float(dcyInt)

    xCell = filein.variables["xCell"]
    xCell[:] -= 2.5 * dcx

    yCell = filein.variables["yCell"]
    yCell[:] -= 3.0 * dcy
    #for iCell in range(0,nCells):
    #    yCell[iCell] = float((int(round((6.0 * yCell[iCell]) / dcy)) * dcyInt) / 6)

    xVertex = filein.variables["xVertex"]
    xVertex[:] -= 2.5 * dcx

    yVertex = filein.variables["yVertex"]
    yVertex[:] -= 3.0 * dcy
    #for iVertex in range(0,nVertices):
    #    yVertex[iVertex] = float((int(round((6.0 * yVertex[iVertex]) / dcy)) * dcyInt) / 6)

    xEdge = filein.variables["xEdge"]
    xEdge[:] -= 2.5 * dcx

    yEdge = filein.variables["yEdge"]
    yEdge[:] -= 3.0 * dcy
    #for iEdge in range(0,nEdges):
    #    yEdge[iEdge] = float((int(round((6.0 * xEdge[iEdge]) / dcy)) * dcyInt) / 6)

    filein.close()

    # plot
    if (plot):
        filein = Dataset("grid_culled.nc","r")

        xCell = filein.variables["xCell"][:]
        yCell = filein.variables["yCell"][:]

        xVertex = filein.variables["xVertex"][:]
        yVertex = filein.variables["yVertex"][:]

        xEdge = filein.variables["xEdge"][:]
        yEdge = filein.variables["yEdge"][:]

        plt.scatter(xCell[:],yCell[:],s=3.0)
        plt.scatter(xVertex[:],yVertex[:],s=3.0)
        plt.scatter(xEdge[:],yEdge[:],s=3.0)

        rect = plt.Rectangle((0.0, 0.0), Lx, Ly, fill=False, edgecolor="red", linewidth=2.5)
        plt.gca().add_patch(rect)

        plt.show()

        filein.close()

    cmd = "mv grid_culled.nc %s" %(gridname)
    os.system(cmd)

#-------------------------------------------------------------------------------

def create_grids():

    Lx = 80000.0
    nx = 51
    #Lx = 1280000.0
    #nx = 81

    #nGrid = 8
    nGrid = 1

    dcx = Lx / float(nx-1)
    dcy = 0.5 * math.sqrt(3.0) * dcx

    print(dcx, dcy)

    ny = Lx / dcy + 1
    ny = int(math.ceil(ny))
    if ((ny % 2) == 0): ny -= 1
    Ly = float(ny-1) * dcy

    print(Lx, Ly)


    nxs = []
    nys = []
    for i in range(0,nGrid):
        nxs.append(nx)
        nys.append(ny)
        nx = 2*nx - 1
        ny = 2*ny - 1


    for nx, ny in zip(nxs, nys):
        create_grid(nx, ny, Lx, Ly)

#-------------------------------------------------------------------------------

if __name__ == "__main__":

    create_grids()
