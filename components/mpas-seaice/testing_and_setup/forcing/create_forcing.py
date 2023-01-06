from netCDF4 import Dataset
import numpy as np
import scipy.sparse as sparse

#-------------------------------------------------------------------------------

def create_scrip_grid_file(filenameScrip, nGridSize, nGridCorners, gridRank, gridDims, gridCenterLat, gridCenterLon, gridImask, gridCornerLat, gridCornerLon, title):

    fileScrip = Dataset(filenameScrip,"w",format="NETCDF3_CLASSIC")

    # dimensions
    fileScrip.createDimension("grid_size", nGridSize)
    fileScrip.createDimension("grid_corners", nGridCorners)
    fileScrip.createDimension("grid_rank", gridRank)

    # variables
    var = fileScrip.createVariable("grid_dims","i",dimensions=["grid_rank"])
    var[:] = gridDims[:]

    var = fileScrip.createVariable("grid_center_lat","d",dimensions=["grid_size"])
    var.units = "radians"
    var[:] = gridCenterLat[:]

    var = fileScrip.createVariable("grid_center_lon","d",dimensions=["grid_size"])
    var.units = "radians"
    var[:] = gridCenterLon[:]

    var = fileScrip.createVariable("grid_imask","i",dimensions=["grid_size"])
    var[:] = gridImask[:]

    var = fileScrip.createVariable("grid_corner_lat","d",dimensions=["grid_size","grid_corners"])
    var.units = "radians"
    var[:] = gridCornerLat[:]

    var = fileScrip.createVariable("grid_corner_lon","d",dimensions=["grid_size","grid_corners"])
    var.units = "radians"
    var[:] = gridCornerLon[:]

    # attributes
    fileScrip.title = title

    fileScrip.close()

#-------------------------------------------------------------------------------

def get_mpas_grid_info(filenameMPASGrid):

    fileGrid = Dataset(filenameMPASGrid,"r")

    nCells    = len(fileGrid.dimensions["nCells"])
    maxEdges  = len(fileGrid.dimensions["maxEdges"])
    nVertices = len(fileGrid.dimensions["nVertices"])

    latCell   = fileGrid.variables["latCell"][:]
    lonCell   = fileGrid.variables["lonCell"][:]
    latVertex = fileGrid.variables["latVertex"][:]
    lonVertex = fileGrid.variables["lonVertex"][:]

    nEdgesOnCell   = fileGrid.variables["nEdgesOnCell"][:]
    verticesOnCell = fileGrid.variables["verticesOnCell"][:]

    fileGrid.close()

    return nCells, maxEdges, nVertices, latCell, lonCell, latVertex, lonVertex, nEdgesOnCell, verticesOnCell

#-------------------------------------------------------------------------------

def create_scrip_file_MPAS(filenameMPASGrid, filenameScrip):

    nCells, maxEdges, nVertices, centerLat, centerLon, latVertex, lonVertex, nEdgesOnCell, verticesOnCell = get_mpas_grid_info(filenameMPASGrid)

    cornerLat = np.zeros((nCells,maxEdges))
    cornerLon = np.zeros((nCells,maxEdges))

    # create the corner arrays
    for iCell in range(0,nCells):

       # fill vertices we have
       for iVertexOnCell in range(0,nEdgesOnCell[iCell]):

          iVertex = verticesOnCell[iCell,iVertexOnCell] - 1

          cornerLat[iCell,iVertexOnCell] = latVertex[iVertex]
          cornerLon[iCell,iVertexOnCell] = lonVertex[iVertex]

       # fill up maxEdges corners with repeated vertex
       for iVertexOnCell in range(nEdgesOnCell[iCell], maxEdges):

          cornerLat[iCell,iVertexOnCell] = cornerLat[iCell,nEdgesOnCell[iCell]-1]
          cornerLon[iCell,iVertexOnCell] = cornerLon[iCell,nEdgesOnCell[iCell]-1]

    # create scrip file
    gridDims = np.array([nCells])
    gridImask = np.ones(nCells,dtype="i")

    create_scrip_grid_file(filenameScrip, nCells, maxEdges, 1, gridDims, centerLat, centerLon, gridImask, cornerLat, cornerLon, "MPAS")

#-------------------------------------------------------------------------------

def write_scrip_in_file(srcTitle):

    scripFile = open("scrip_in","w")

    scripFile.write("&remap_inputs\n")
    scripFile.write("    num_maps = 1\n")
    scripFile.write("    grid1_file = 'remap_grid_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    grid2_file = 'remap_grid_MPAS_tmp.nc'\n")
    scripFile.write("    interp_file1 = 'remap_%s_to_MPAS_tmp.nc'\n" %(srcTitle))
    scripFile.write("    interp_file2 = 'remap_MPAS_to_%s_tmp.nc'\n" %(srcTitle))
    scripFile.write("    map1_name = '%s to MPAS bilinear mapping'\n" %(srcTitle))
    scripFile.write("    map2_name = 'MPAS to %s bilinear mapping'\n" %(srcTitle))
    scripFile.write("    map_method = 'bilinear'\n")
    scripFile.write("    normalize_opt = 'frac'\n")
    scripFile.write("    output_opt = 'scrip'\n")
    scripFile.write("    restrict_type = 'latitude'\n")
    scripFile.write("    num_srch_bins = 90 \n")
    scripFile.write("    luse_grid1_area = .false.\n")
    scripFile.write("    luse_grid2_area = .false.\n")
    scripFile.write("/\n")

    scripFile.close()

#-------------------------------------------------------------------------------

def create_output_times(inputTimesPerYear, year):

    daysInMonth = [31,28,31,30,31,30,31,31,30,31,30,31]

    xtimes = []

    if (inputTimesPerYear == 1460):

        minute = 0
        second = 0

        for iMonth in range(0,12):
            for iDay in range(0,daysInMonth[iMonth]):
                for iSixHours in range(0,4):

                    month = iMonth + 1
                    day = iDay + 1
                    hour = (iSixHours + 1) * 6

                    timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
                    xtimes.append(timeStr)

    elif (inputTimesPerYear == 12):

        day = 15
        hour = 0
        minute = 0
        second = 0

        for iMonth in range(0,12):

            month = iMonth + 1

            timeStr = "%4.4i-%2.2i-%2.2i_%2.2i:%2.2i:%2.2i" %(year,month,day,hour,minute,second)
            xtimes.append(timeStr)

    return xtimes

#-------------------------------------------------------------------------------

def get_remapping_data(filenameRemapping):

    fileRemap = Dataset(filenameRemapping,"r")

    numLinks    = len(fileRemap.dimensions["num_links"])
    numWgts     = len(fileRemap.dimensions["num_wgts"])
    srcGridSize = len(fileRemap.dimensions["src_grid_size"])
    dstGridSize = len(fileRemap.dimensions["dst_grid_size"])

    srcAddress  = fileRemap.variables["src_address"][:]
    dstAddress  = fileRemap.variables["dst_address"][:]
    remapMatrix = fileRemap.variables["remap_matrix"][:]

    # fortran indices
    srcAddress[:] = srcAddress[:] - 1
    dstAddress[:] = dstAddress[:] - 1

    fileRemap.close()

    # covert to python sparse arrays
    remapMatrixSparse = sparse.coo_matrix((remapMatrix[:,0], (dstAddress, srcAddress)), shape=(dstGridSize,srcGridSize))
    remapMatrixSparse = remapMatrixSparse.tocsr()

    return remapMatrixSparse, dstGridSize

#-------------------------------------------------------------------------------
