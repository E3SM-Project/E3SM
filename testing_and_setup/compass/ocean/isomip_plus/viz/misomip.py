import numpy
from netCDF4 import Dataset
from shapely.geometry import Polygon, LineString

from progressbar import ProgressBar, Percentage, Bar, ETA

import os
import glob


def compute_misomip_interp_coeffs(folder, expt):

    def getTransectWeights(outFileName, axis):

        pbar = ProgressBar(
            widgets=[
                'compute MPAS and MISOMIP {} transect intersections'.format(
                    axis),
                Percentage(),
                Bar(),
                ETA()],
            maxval=nCells).start()

        if axis == 'x':
            slicePos = xTransect
            outNOther = outNy
            outOtherAxis = y
        else:
            slicePos = yTransect
            outNOther = outNx
            outOtherAxis = x

        sliceCount = numpy.zeros(outNOther, int)
        cellIndices = []
        otherIndices = []
        weights = []
        sliceIndices = []

        for iCell in range(nCells):
            verts = verticesOnCell[iCell, 0:nEdgesOnCell[iCell]]
            verts = numpy.append(verts, verts[0])
            xVert = xVertex[verts]
            yVert = yVertex[verts]

            if axis == 'x':
                sliceAxisVerts = xVert
                otherAxisVerts = yVert
            else:
                sliceAxisVerts = yVert
                otherAxisVerts = xVert

            if(numpy.amax(sliceAxisVerts) < slicePos) or \
                    (numpy.amin(sliceAxisVerts) > slicePos):
                # this polygon doesn't intersect the slice
                continue

            mpasPolygon = Polygon(zip(xVert, yVert))

            indices = numpy.nonzero(
                outOtherAxis < numpy.amin(otherAxisVerts))[0]
            if len(indices) == 0:
                lower = 0
            else:
                lower = indices[-1]

            indices = numpy.nonzero(
                outOtherAxis > numpy.amax(otherAxisVerts))[0]
            if len(indices) == 0:
                upper = outNOther
            else:
                upper = indices[0]

            for otherIndex in range(lower, upper):
                if axis == 'x':
                    vertices = ((slicePos, outOtherAxis[otherIndex]),
                                (slicePos, outOtherAxis[otherIndex + 1]))
                else:
                    vertices = ((outOtherAxis[otherIndex], slicePos),
                                (outOtherAxis[otherIndex + 1], slicePos))

                line = LineString(vertices)
                if not mpasPolygon.intersects(line):
                    continue
                length = mpasPolygon.intersection(line).length
                if length == 0.0:
                    continue

                cellIndices.append(iCell)
                otherIndices.append(otherIndex)

                weights.append(length / outDx)

                sliceIndex = sliceCount[otherIndex]
                sliceCount[otherIndex] += 1
                sliceIndices.append(sliceIndex)
            pbar.update(iCell + 1)

        pbar.finish()

        cellIndices = numpy.array(cellIndices)
        otherIndices = numpy.array(otherIndices)
        sliceIndices = numpy.array(sliceIndices)
        weights = numpy.array(weights)

        # sort the intersections first by otherIndex, then by sliceIndex
        # for efficiency
        nSlices = numpy.amax(sliceCount)
        sortedIndices = numpy.zeros(0, int)
        for sliceIndex in range(nSlices):
            intersectionsInSlice = numpy.nonzero(sliceIndices == sliceIndex)[0]
            indices = numpy.argsort(otherIndices[intersectionsInSlice])
            sortedIndices = numpy.append(
                sortedIndices, intersectionsInSlice[indices])

        outFile = Dataset(outFileName, 'w', format='NETCDF4')
        outFile.createDimension('nIntersections', len(cellIndices))
        outFile.createVariable('cellIndices', 'i4', ('nIntersections',))
        if axis == 'x':
            outFile.createVariable('yIndices', 'i4', ('nIntersections',))
        else:
            outFile.createVariable('xIndices', 'i4', ('nIntersections',))
        outFile.createVariable('sliceIndices', 'i4', ('nIntersections',))
        outFile.createVariable(
            'mpasToMisomipWeights', 'f8', ('nIntersections',))

        outVars = outFile.variables
        outVars['cellIndices'][:] = cellIndices[sortedIndices]
        if axis == 'x':
            outVars['yIndices'][:] = otherIndices[sortedIndices]
        else:
            outVars['xIndices'][:] = otherIndices[sortedIndices]

        outVars['sliceIndices'][:] = sliceIndices[sortedIndices]
        outVars['mpasToMisomipWeights'][:] = weights[sortedIndices]

        outFile.close()

    try:
        os.makedirs('{}/misomip'.format(folder))
    except OSError:
        pass

    meshFileName = '%s/init.nc' % folder
    interpWeightsFileName = '%s/misomip/horiz_map.nc' % folder
    xTransectFileName = '%s/misomip/x_trans_map.nc' % folder
    yTransectFileName = '%s/misomip/y_trans_map.nc' % folder

    outNx, outNy, outNz, x, y, z, xTransect, yTransect, outDx, outDz = \
        _get_out_grid(corners=True)

    inFile = Dataset(meshFileName, 'r')

    nCells = len(inFile.dimensions['nCells'])
    inVars = inFile.variables
    nEdgesOnCell = inVars['nEdgesOnCell'][:]
    verticesOnCell = inVars['verticesOnCell'][:, :] - 1
    xVertex = inVars['xVertex'][:]
    yVertex = inVars['yVertex'][:]

    inFile.close()
    if(not os.path.exists(interpWeightsFileName)):

        cellIndices = []
        xIndices = []
        yIndices = []
        weights = []
        sliceIndices = []

        sliceCount = numpy.zeros((outNy, outNx), int)

        pbar = ProgressBar(
            widgets=[
                'compute MPAS and MISOMIP intersections',
                Percentage(),
                Bar(),
                ETA()],
            maxval=nCells).start()

        for iCell in range(nCells):
            verts = verticesOnCell[iCell, 0:nEdgesOnCell[iCell]]
            verts = numpy.append(verts, verts[0])
            xVert = xVertex[verts]
            yVert = yVertex[verts]
            mpasPolygon = Polygon(zip(xVert, yVert))

            # find the out indices that bound the MPAS polygon
            indices = numpy.nonzero(x < numpy.amin(xVert))[0]
            if len(indices) == 0:
                xl = 0
            else:
                xl = indices[-1]

            indices = numpy.nonzero(x > numpy.amax(xVert))[0]
            if len(indices) == 0:
                xu = outNx
            else:
                xu = indices[0]

            indices = numpy.nonzero(y < numpy.amin(yVert))[0]
            if len(indices) == 0:
                yl = 0
            else:
                yl = indices[-1]

            indices = numpy.nonzero(y > numpy.amax(yVert))[0]
            if len(indices) == 0:
                yu = outNy
            else:
                yu = indices[0]

            for yIndex in range(yl, yu):
                for xIndex in range(xl, xu):
                    vertices = ((x[xIndex], y[yIndex]),
                                (x[xIndex + 1], y[yIndex]),
                                (x[xIndex + 1], y[yIndex + 1]),
                                (x[xIndex], y[yIndex + 1]),
                                (x[xIndex], y[yIndex]))
                    outPoly = Polygon(vertices)
                    if not mpasPolygon.intersects(outPoly):
                        continue
                    intersectionArea = mpasPolygon.intersection(outPoly).area
                    if intersectionArea == 0.:
                        continue

                    cellIndices.append(iCell)
                    xIndices.append(xIndex)
                    yIndices.append(yIndex)

                    weights.append(intersectionArea / outDx**2)

                    sliceIndex = sliceCount[yIndex, xIndex]
                    sliceCount[yIndex, xIndex] += 1
                    sliceIndices.append(sliceIndex)
            pbar.update(iCell + 1)

        pbar.finish()

        cellIndices = numpy.array(cellIndices)
        xIndices = numpy.array(xIndices)
        yIndices = numpy.array(yIndices)
        sliceIndices = numpy.array(sliceIndices)
        weights = numpy.array(weights)

        # sort the intersections first by xIndex, then by yIndex, then by
        # sliceIndex for efficiency
        nSlices = numpy.amax(sliceCount)
        sortedIndices = numpy.zeros(0, int)
        xyIndices = xIndices + outNx * yIndices
        for sliceIndex in range(nSlices):
            intersectionsInSlice = numpy.nonzero(sliceIndices == sliceIndex)[0]
            indices = numpy.argsort(xyIndices[intersectionsInSlice])
            sortedIndices = numpy.append(
                sortedIndices, intersectionsInSlice[indices])

        outFile = Dataset(interpWeightsFileName, 'w', format='NETCDF4')
        outFile.createDimension('nIntersections', len(cellIndices))
        outFile.createVariable('cellIndices', 'i4', ('nIntersections',))
        outFile.createVariable('xIndices', 'i4', ('nIntersections',))
        outFile.createVariable('yIndices', 'i4', ('nIntersections',))
        outFile.createVariable('sliceIndices', 'i4', ('nIntersections',))
        outFile.createVariable(
            'mpasToMisomipWeights', 'f8', ('nIntersections',))

        outVars = outFile.variables
        outVars['cellIndices'][:] = cellIndices[sortedIndices]
        outVars['xIndices'][:] = xIndices[sortedIndices]
        outVars['yIndices'][:] = yIndices[sortedIndices]
        outVars['sliceIndices'][:] = sliceIndices[sortedIndices]
        outVars['mpasToMisomipWeights'][:] = weights[sortedIndices]

        outFile.close()

    if(not os.path.exists(xTransectFileName)):
        getTransectWeights(xTransectFileName, axis='x')

    if(not os.path.exists(yTransectFileName)):
        getTransectWeights(yTransectFileName, axis='y')


def interp_misomip(folder, expt):

    def interpHoriz(field, inMask=None, outFraction=None):
        if inMask is not None:
            field = field * inMask
        outField = numpy.zeros((outNy, outNx))
        for sliceIndex in range(xyNSlices):
            mask = xySliceIndices == sliceIndex
            cellsSlice = xyCellIndices[mask]
            fieldSlice = field[cellsSlice]
            outField[xyYIndices[mask], xyXIndices[mask]
                     ] += (xyMpasToMisomipWeights[mask] * fieldSlice)
        if outFraction is not None:
            mask = outFraction > normalizationThreshold
            outField[mask] /= outFraction[mask]
            outField[numpy.logical_not(mask)] = 0.
        return outField

    def interpHorizOcean(field):
        return interpHoriz(field, cellOceanMask, xyOceanFraction)

    def interpHorizCavity(field):
        return interpHoriz(field, inCavityFraction, outCavityFraction)

    def interpXZTransect(field, normalize=True):
        outField = numpy.zeros((outNz, outNx))

        for sliceIndex in range(xzNSlices):
            mask = xzSliceIndices == sliceIndex
            cellsSlice = xzCellIndices[mask]
            xIndices = xzXIndices[mask]
            weights = xzMpasToMisomipWeights[mask]
            for index in range(len(cellsSlice)):
                iCell = cellsSlice[index]
                xIndex = xIndices[index]
                weight = weights[index]
                layerThick = layerThickness[iCell, :]
                layerThick[maxLevelCell[iCell] + 1:] = 0.
                zInterface = numpy.zeros(2 * nVertLevels)
                fieldColumn = numpy.zeros(2 * nVertLevels)
                fieldColumn[0::2] = field[iCell, :]
                fieldColumn[1::2] = field[iCell, :]
                zInterface[0] = ssh[iCell]
                zInterface[1::2] = ssh[iCell] - numpy.cumsum(layerThick[:])
                zInterface[2::2] = ssh[iCell] - numpy.cumsum(layerThick[:-1])
                outField[:,
                         xIndex] += weight * numpy.interp(z,
                                                          zInterface[::-1],
                                                          fieldColumn[::-1],
                                                          left=0.,
                                                          right=0.)

        if normalize:
            outField[xzOceanMask] /= xzOceanFraction[xzOceanMask]
            outField[numpy.logical_not(xzOceanMask)] = 0.

        return outField

    def interpYZTransect(field, normalize=True):
        outField = numpy.zeros((outNz, outNy))

        for sliceIndex in range(yzNSlices):
            mask = yzSliceIndices == sliceIndex
            cellsSlice = yzCellIndices[mask]
            yIndices = yzYIndices[mask]
            weights = yzMpasToMisomipWeights[mask]
            for index in range(len(cellsSlice)):
                iCell = cellsSlice[index]
                yIndex = yIndices[index]
                weight = weights[index]
                layerThick = layerThickness[iCell, :]
                layerThick[maxLevelCell[iCell] + 1:] = 0.
                zInterface = numpy.zeros(2 * nVertLevels)
                fieldColumn = numpy.zeros(2 * nVertLevels)
                fieldColumn[0::2] = field[iCell, :]
                fieldColumn[1::2] = field[iCell, :]
                zInterface[0] = ssh[iCell]
                zInterface[1::2] = ssh[iCell] - numpy.cumsum(layerThick[:])
                zInterface[2::2] = ssh[iCell] - numpy.cumsum(layerThick[:-1])
                outField[:,
                         yIndex] += weight * numpy.interp(z,
                                                          zInterface[::-1],
                                                          fieldColumn[::-1],
                                                          left=0.,
                                                          right=0.)

        if normalize:
            outField[yzOceanMask] /= yzOceanFraction[yzOceanMask]
            outField[numpy.logical_not(yzOceanMask)] = 0.

        return outField

    def writeMetric(varName, metric):
        vars[varName][tIndex] = metric

    def writeVar(varName, varField, varMask=None):
        if varMask is None:
            maskedVar = varField
        else:
            maskedVar = numpy.ma.masked_array(varField,
                                              mask=numpy.logical_not(varMask))
        vars[varName][tIndex, :, :] = maskedVar

    inFileNames = sorted(glob.glob('%s/timeSeriesStatsMonthly.*.nc' % folder))

    initFile = Dataset('{}/init.nc'.format(folder), 'r')
    outFileName = '%s/Ocean%s_COM_MPAS-Ocean.nc' % (folder, expt)

    rho_fw = 1000.
    secPerDay = 24 * 60 * 60

    normalizationThreshold = 0.001

    outNx, outNy, outNz, x, y, z, xTransect, yTransect, outDx, outDz = \
        _get_out_grid(corners=False)

    inFile = Dataset('%s/misomip/horiz_map.nc' % folder, 'r')
    inVars = inFile.variables
    xyCellIndices = inVars['cellIndices'][:]
    xyXIndices = inVars['xIndices'][:]
    xyYIndices = inVars['yIndices'][:]
    xySliceIndices = inVars['sliceIndices'][:]
    xyMpasToMisomipWeights = inVars['mpasToMisomipWeights'][:]
    inFile.close()
    xyNSlices = numpy.amax(xySliceIndices) + 1

    inFile = Dataset('%s/misomip/x_trans_map.nc' % folder, 'r')
    inVars = inFile.variables
    yzCellIndices = inVars['cellIndices'][:]
    yzYIndices = inVars['yIndices'][:]
    yzSliceIndices = inVars['sliceIndices'][:]
    yzMpasToMisomipWeights = inVars['mpasToMisomipWeights'][:]
    inFile.close()
    yzNSlices = numpy.amax(yzSliceIndices) + 1

    inFile = Dataset('%s/misomip/y_trans_map.nc' % folder, 'r')
    inVars = inFile.variables
    xzCellIndices = inVars['cellIndices'][:]
    xzXIndices = inVars['xIndices'][:]
    xzSliceIndices = inVars['sliceIndices'][:]
    xzMpasToMisomipWeights = inVars['mpasToMisomipWeights'][:]
    inFile.close()
    xzNSlices = numpy.amax(xzSliceIndices) + 1

    dynamicTopo = False

    initFile = Dataset('%s/init.nc' % folder, 'r')
    osfFile = Dataset('%s/overturningStreamfunction.nc' % folder, 'r')
    bsfFile = Dataset('%s/barotropicStreamfunction.nc' % folder, 'r')
    continueOutput = os.path.exists(outFileName)
    if(continueOutput):
        outFile = Dataset(outFileName, 'r+', format='NETCDF4')
        vars = outFile.variables
    else:
        outFile = Dataset(outFileName, 'w', format='NETCDF4')

        outFile.createDimension('nTime', None)
        outFile.createDimension('nx', outNx)
        outFile.createDimension('ny', outNy)
        outFile.createDimension('nz', outNz)

        for varName in ['x', 'y', 'z']:
            var = outFile.createVariable(varName, 'f4', ('n%s' % varName,))
            var.units = 'm'
            var.description = '%s location of cell centers' % varName
        vars = outFile.variables
        print(x.shape, vars['x'])
        vars['x'][:] = x
        vars['y'][:] = y
        vars['z'][:] = z

        var = outFile.createVariable('time', 'f4', ('nTime',))
        var.units = 's'
        var.description = 'time since start of simulation'

        var = outFile.createVariable('meanMeltRate', 'f4', ('nTime'))
        var.units = 'm/s'
        var.description = 'mean melt rate averaged over area of floating ' \
                          'ice, positive for melting'

        var = outFile.createVariable('totalMeltFlux', 'f4', ('nTime'))
        var.units = 'kg/s'
        var.description = 'total flux of melt water summed over area of ' \
                          'floating ice, positive for melting'

        var = outFile.createVariable('totalOceanVolume', 'f4', ('nTime'))
        var.units = 'm^3'
        var.description = 'total volume of ocean'

        var = outFile.createVariable('meanTemperature', 'f4', ('nTime'))
        var.units = 'deg C'
        var.description = 'the potential temperature averaged over the ' \
                          'ocean volume'

        var = outFile.createVariable('meanSalinity', 'f4', ('nTime'))
        var.units = 'PSU'
        var.description = 'the salinity averaged over the ocean volume'

        if dynamicTopo:
            outFile.createVariable('iceDraft', 'f4', ('nTime', 'ny', 'nx',))
            outFile.createVariable('bathymetry', 'f4', ('nTime', 'ny', 'nx',))
        else:
            outFile.createVariable('iceDraft', 'f4', ('ny', 'nx',))
            outFile.createVariable('bathymetry', 'f4', ('ny', 'nx',))

        var = vars['iceDraft']
        var.units = 'm'
        var.description = 'elevation of the ice-ocean interface'

        var = vars['bathymetry']
        var.units = 'm'
        var.description = 'elevation of the bathymetry'

        var = outFile.createVariable('meltRate', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'm/s'
        var.description = 'melt rate, positive for melting'

        var = outFile.createVariable(
            'frictionVelocity', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'm/s'
        var.description = 'friction velocity u* used in melt calculations'

        var = outFile.createVariable(
            'thermalDriving', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'deg C'
        var.description = 'thermal driving used in the melt calculation'

        var = outFile.createVariable(
            'halineDriving', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'PSU'
        var.description = 'haline driving used in the melt calculation'

        var = outFile.createVariable(
            'uBoundaryLayer', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'm/s'
        var.description = 'x-velocity in the boundary layer used to compute u*'

        var = outFile.createVariable(
            'vBoundaryLayer', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'm/s'
        var.description = 'y-velocity in the boundary layer used to compute u*'

        var = outFile.createVariable(
            'barotropicStreamfunction', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'm^3/s'
        var.description = 'barotropic streamfunction'

        var = outFile.createVariable(
            'overturningStreamfunction', 'f4', ('nTime', 'nz', 'nx',))
        var.units = 'm^3/s'
        var.description = 'overturning (meridional) streamfunction'

        var = outFile.createVariable(
            'bottomTemperature', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'deg C'
        var.description = 'temperature in the bottom grid cell of each ' \
                          'ocean column'

        var = outFile.createVariable(
            'bottomSalinity', 'f4', ('nTime', 'ny', 'nx',))
        var.units = 'PSU'
        var.description = 'salinity in the bottom grid cell of each ocean ' \
                          'column'

        var = outFile.createVariable(
            'temperatureXZ', 'f4', ('nTime', 'nz', 'nx',))
        var.units = 'deg C'
        var.description = 'temperature slice in x-z plane through the center' \
                          ' of the domain (y = 40 km)'

        var = outFile.createVariable('salinityXZ', 'f4',
                                     ('nTime', 'nz', 'nx',))
        var.units = 'PSU'
        var.description = 'salinity slice in x-z plane through the center of' \
                          ' the domain (y = 40 km)'

        var = outFile.createVariable(
            'temperatureYZ', 'f4', ('nTime', 'nz', 'ny',))
        var.units = 'deg C'
        var.description = 'temperature slice in y-z plane through x = 500 km'

        var = outFile.createVariable('salinityYZ', 'f4',
                                     ('nTime', 'nz', 'ny',))
        var.units = 'PSU'
        var.description = 'salinity slice in y-z plane through x = 500 km'

    nCells = len(initFile.dimensions['nCells'])
    nVertLevels = len(initFile.dimensions['nVertLevels'])
    nTimeIn = len(inFileNames)
    nTimeIn = min(nTimeIn, len(bsfFile.dimensions['Time']))
    nTimeIn = min(nTimeIn, len(osfFile.dimensions['Time']))

    if(continueOutput):
        nTimeOut = max(0, len(outFile.dimensions['nTime']) - 2)
    else:
        nTimeOut = 0

    areaCell = initFile.variables['areaCell'][:]
    bathymetry = -initFile.variables['bottomDepth'][:]
    maxLevelCell = initFile.variables['maxLevelCell'][:] - 1

    cellMask = numpy.zeros((nCells, nVertLevels))
    for iCell in range(nCells):
        k = maxLevelCell[iCell]
        if (k >= 0):
            cellMask[iCell, 0:k + 1] = 1.0

    cellOceanMask = cellMask[:, 0]
    xyOceanFraction = interpHoriz(cellOceanMask)
    xyOceanMask = xyOceanFraction > normalizationThreshold

    if not dynamicTopo and (nTimeOut == 0):
        vars['iceDraft'][:, :] = interpHorizOcean(
            initFile.variables['ssh'][0, :])
        vars['bathymetry'][:, :] = interpHorizOcean(bathymetry)

    initFile.close()

    daysInMonth = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    daysBeforeMonth = [0] + list(numpy.cumsum(daysInMonth))

    pbar = ProgressBar(
        widgets=[
            'interpolate from MPAS to MISOMIP',
            Percentage(),
            Bar(),
            ETA()],
        maxval=nTimeIn).start()
    for tIndex in range(nTimeOut, nTimeIn):
        inFile = Dataset(inFileNames[tIndex], 'r')
        inVars = inFile.variables
        days = 365 * int(tIndex / 12) + daysBeforeMonth[numpy.mod(tIndex, 12)]
        vars['time'][tIndex] = secPerDay * days

        freshwaterFlux = inVars['timeMonthly_avg_landIceFreshwaterFlux'][0, :]
        inCavityFraction = inVars['timeMonthly_avg_landIceFraction'][0, :]
        outCavityFraction = interpHoriz(inCavityFraction)
        outCavityMask = outCavityFraction > normalizationThreshold
        meltRate = freshwaterFlux / rho_fw

        if not numpy.all(inCavityFraction == 0.):
            writeMetric('meanMeltRate', numpy.sum(meltRate * areaCell)
                        / numpy.sum(inCavityFraction * areaCell))

        writeMetric('totalMeltFlux', numpy.sum(freshwaterFlux * areaCell))

        ssh = inVars['timeMonthly_avg_ssh'][0, :]
        columnThickness = ssh - bathymetry

        writeMetric('totalOceanVolume', numpy.sum(columnThickness * areaCell))

        thermalDriving = (inVars['timeMonthly_avg_landIceBoundaryLayerTracers'
                                 '_landIceBoundaryLayerTemperature'][0, :]
                          - inVars['timeMonthly_avg_landIceInterfaceTracers'
                                   '_landIceInterfaceTemperature'][0, :])
        halineDriving = (inVars['timeMonthly_avg_landIceBoundaryLayerTracers'
                                '_landIceBoundaryLayerSalinity'][0, :]
                         - inVars['timeMonthly_avg_landIceInterfaceTracers'
                                  '_landIceInterfaceSalinity'][0, :])
        frictionVelocity = \
            inVars['timeMonthly_avg_landIceFrictionVelocity'][0, :]

        bsfCell = 1e6 * bsfFile.variables['bsfCell'][tIndex, :]

        # meltRate is already multiplied by inCavityFraction, so no in masking
        writeVar(
            'meltRate',
            interpHoriz(
                meltRate,
                inMask=None,
                outFraction=outCavityFraction),
            outCavityMask)
        writeVar(
            'thermalDriving',
            interpHorizCavity(thermalDriving),
            outCavityMask)
        writeVar(
            'halineDriving',
            interpHorizCavity(halineDriving),
            outCavityMask)
        writeVar('frictionVelocity', interpHorizCavity(frictionVelocity),
                 outCavityMask)
        writeVar('barotropicStreamfunction', interpHorizOcean(bsfCell),
                 xyOceanMask)

        temperature = \
            inVars['timeMonthly_avg_activeTracers_temperature'][0, :, :]
        salinity = inVars['timeMonthly_avg_activeTracers_salinity'][0, :, :]
        layerThickness = inVars['timeMonthly_avg_layerThickness'][0, :, :]

        indices = numpy.arange(nCells)
        bottomTemperature = temperature[indices, maxLevelCell]
        bottomSalinity = salinity[indices, maxLevelCell]

        writeVar('bottomTemperature', interpHorizOcean(bottomTemperature),
                 xyOceanMask)
        writeVar('bottomSalinity', interpHorizOcean(bottomSalinity),
                 xyOceanMask)

        writeMetric(
            'meanTemperature',
            numpy.sum(
                cellMask *
                layerThickness *
                temperature) /
            numpy.sum(
                cellMask *
                layerThickness))
        writeMetric(
            'meanSalinity',
            numpy.sum(
                cellMask *
                layerThickness *
                salinity) /
            numpy.sum(
                cellMask *
                layerThickness))

        uTop = inVars['timeMonthly_avg_velocityX'][0, :, 0]
        vTop = inVars['timeMonthly_avg_velocityY'][0, :, 0]
        writeVar('uBoundaryLayer', interpHorizCavity(uTop), outCavityMask)
        writeVar('vBoundaryLayer', interpHorizCavity(vTop), outCavityMask)

        osf = 1e6 * osfFile.variables['osf'][tIndex, :, :]
        osfX = osfFile.variables['x'][:]
        osfZ = osfFile.variables['z'][:]
        assert(numpy.all(osfX == x))
        # the first and last OSF z values have been tweaked...
        assert(numpy.all(osfZ[1:-1] == z[1:-1]))

        writeVar('overturningStreamfunction', osf)

        xzOceanFraction = interpXZTransect(cellMask, normalize=False)
        xzOceanMask = xzOceanFraction > 0.001

        yzOceanFraction = interpYZTransect(cellMask, normalize=False)
        yzOceanMask = yzOceanFraction > 0.001

        writeVar('temperatureXZ', interpXZTransect(temperature), xzOceanMask)
        writeVar('salinityXZ', interpXZTransect(salinity), xzOceanMask)

        writeVar('temperatureYZ', interpYZTransect(temperature), yzOceanMask)
        writeVar('salinityYZ', interpYZTransect(salinity), yzOceanMask)
        pbar.update(tIndex + 1)

    pbar.finish()

    outFile.close()
    osfFile.close()
    bsfFile.close()


def _get_out_grid(corners):

    outDx = 2e3
    outX0 = 320e3
    outY0 = 0.
    xTransect = 520.01e3  # shifted slightly to avoid exact touching
    yTransect = 40.01e3
    outDz = 5.

    outNx = 240
    outNy = 40
    outNz = 144

    if corners:
        # x, y and z of corners of grid cells
        x = outX0 + outDx * (numpy.arange(outNx + 1))
        y = outY0 + outDx * (numpy.arange(outNy + 1))
        z = -outDz * (numpy.arange(outNz + 1))
    else:
        x = outX0 + outDx*(numpy.arange(outNx)+0.5)
        y = outY0 + outDx*(numpy.arange(outNy)+0.5)
        z = -outDz*(numpy.arange(outNz)+0.5)

    return outNx, outNy, outNz, x, y, z, xTransect, yTransect, outDx, outDz
