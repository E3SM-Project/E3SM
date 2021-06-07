from netCDF4 import Dataset
import numpy as np

filenameOut = "basis_fields.nc"

filenameWacshpress = "output_hex_wachspress_0082x0094_120/output.2000.nc"
filenamePWL        = "output_hex_pwl_0082x0094_120/output.2000.nc"

fileWachspress = Dataset(filenameWacshpress,"r")

nCells = len(fileWachspress.dimensions["nCells"])
maxEdges = len(fileWachspress.dimensions["maxEdges"])

interiorCell = fileWachspress.variables["interiorCell"][0,:]

basisIntegralsUWachspress = fileWachspress.variables["basisIntegralsU"][:]
basisIntegralsVWachspress = fileWachspress.variables["basisIntegralsV"][:]
basisIntegralsMetricWachspress = fileWachspress.variables["basisIntegralsMetric"][:]
basisGradientUWachspress = fileWachspress.variables["basisGradientU"][:]
basisGradientVWachspress = fileWachspress.variables["basisGradientV"][:]

fileWachspress.close()

filePWL = Dataset(filenamePWL,"r")

basisIntegralsUPWL = filePWL.variables["basisIntegralsU"][:]
basisIntegralsVPWL = filePWL.variables["basisIntegralsV"][:]
basisIntegralsMetricPWL = filePWL.variables["basisIntegralsMetric"][:]
basisGradientUPWL = filePWL.variables["basisGradientU"][:]
basisGradientVPWL = filePWL.variables["basisGradientV"][:]

filePWL.close()


basisIntegralsU = np.zeros((nCells, maxEdges, maxEdges))
basisIntegralsV = np.zeros((nCells, maxEdges, maxEdges))
basisIntegralsMetric = np.zeros((nCells, maxEdges, maxEdges))
basisGradientU = np.zeros((nCells, maxEdges, maxEdges))
basisGradientV = np.zeros((nCells, maxEdges, maxEdges))

for iCell in range(0,nCells):

    if (interiorCell[iCell] == 1):
        # interior

        #basisIntegralsU[iCell,:,:] = basisIntegralsUWachspress[iCell,:,:]
        #basisIntegralsV[iCell,:,:] = basisIntegralsVWachspress[iCell,:,:]
        #basisIntegralsMetric[iCell,:,:] = basisIntegralsMetricWachspress[iCell,:,:]
        #basisGradientU[iCell,:,:] = basisGradientUWachspress[iCell,:,:]
        #basisGradientV[iCell,:,:] = basisGradientVWachspress[iCell,:,:]

        basisIntegralsU[iCell,:,:] = basisIntegralsUPWL[iCell,:,:]
        basisIntegralsV[iCell,:,:] = basisIntegralsVPWL[iCell,:,:]
        basisIntegralsMetric[iCell,:,:] = basisIntegralsMetricPWL[iCell,:,:]
        basisGradientU[iCell,:,:] = basisGradientUPWL[iCell,:,:]
        basisGradientV[iCell,:,:] = basisGradientVPWL[iCell,:,:]

    else:
        # boundary

        #basisIntegralsU[iCell,:,:] = basisIntegralsUWachspress[iCell,:,:]
        #basisIntegralsV[iCell,:,:] = basisIntegralsVWachspress[iCell,:,:]
        #basisIntegralsMetric[iCell,:,:] = basisIntegralsMetricWachspress[iCell,:,:]
        #basisGradientU[iCell,:,:] = basisGradientUWachspress[iCell,:,:]
        #basisGradientV[iCell,:,:] = basisGradientVWachspress[iCell,:,:]

        basisIntegralsU[iCell,:,:] = basisIntegralsUPWL[iCell,:,:]
        basisIntegralsV[iCell,:,:] = basisIntegralsVPWL[iCell,:,:]
        basisIntegralsMetric[iCell,:,:] = basisIntegralsMetricPWL[iCell,:,:]
        basisGradientU[iCell,:,:] = basisGradientUPWL[iCell,:,:]
        basisGradientV[iCell,:,:] = basisGradientVPWL[iCell,:,:]


fileOut = Dataset(filenameOut,"w",format="NETCDF3_CLASSIC")

fileOut.createDimension("nCells",nCells)
fileOut.createDimension("maxEdges",maxEdges)

var = fileOut.createVariable("basisIntegralsU","d",dimensions=["nCells", "maxEdges", "maxEdges"])
var[:] = basisIntegralsU[:]

var = fileOut.createVariable("basisIntegralsV","d",dimensions=["nCells", "maxEdges", "maxEdges"])
var[:] = basisIntegralsV[:]

var = fileOut.createVariable("basisIntegralsMetric","d",dimensions=["nCells", "maxEdges", "maxEdges"])
var[:] = basisIntegralsMetric[:]

var = fileOut.createVariable("basisGradientU","d",dimensions=["nCells", "maxEdges", "maxEdges"])
var[:] = basisGradientU[:]

var = fileOut.createVariable("basisGradientV","d",dimensions=["nCells", "maxEdges", "maxEdges"])
var[:] = basisGradientV[:]

fileOut.close()
