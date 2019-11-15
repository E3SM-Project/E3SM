import xarray
import numpy
import scipy.sparse
import scipy.sparse.linalg
import progressbar
import os

from viz.io import write_netcdf, file_complete


def compute_barotropic_streamfunction(dsMesh, ds, folder):
    '''
    compute the barotropic streamfunction for the given mesh and monthly-mean
    data set
    '''

    bsfFileName = '{}/barotropicStreamfunction.nc'.format(folder)
    if file_complete(ds, bsfFileName):
        return

    bsfVertex = _compute_barotropic_streamfunction_vertex(dsMesh, ds)
    bsfCell = _compute_barotropic_streamfunction_cell(dsMesh, bsfVertex)
    dsBSF = xarray.Dataset()
    dsBSF['xtime_startMonthly'] = ds.xtime_startMonthly
    dsBSF['xtime_endMonthly'] = ds.xtime_endMonthly
    dsBSF['bsfVertex'] = bsfVertex
    dsBSF.bsfVertex.attrs['units'] = 'Sv'
    dsBSF.bsfVertex.attrs['description'] = 'barotropic streamfunction ' \
        'on vertices'
    dsBSF['bsfCell'] = bsfCell
    dsBSF.bsfCell.attrs['units'] = 'Sv'
    dsBSF.bsfCell.attrs['description'] = 'barotropic streamfunction ' \
        'on cells'
    dsBSF = dsBSF.transpose('Time', 'nCells', 'nVertices')
    write_netcdf(dsBSF, bsfFileName)


def compute_overturning_streamfunction(dsMesh, ds, folder, dx=2e3, dz=5.):
    '''
    compute the overturning streamfunction for the given mesh and monthly-mean
    data set.

    dx and dz are the resolutions of the OSF in meters
    '''

    osfFileName = '{}/overturningStreamfunction.nc'.format(folder)

    if file_complete(ds, osfFileName):
        return

    xMin = 320e3 + 0.5*dx
    xMax = 800e3 - 0.5*dx
    nx = int((xMax - xMin)/dx + 1)
    x = numpy.linspace(xMin, xMax, nx)

    zMin = -720.0 + 0.5*dz
    zMax = 0.0 - 0.5*dz
    nz = int((zMax - zMin)/dz + 1)
    z = numpy.linspace(zMax, zMin, nz)

    try:
        os.makedirs('{}/cache'.format(folder))
    except OSError:
        pass

    mpasTransportFileName = '{}/cache/osf_mpas_transport.nc'.format(folder)
    _compute_horizontal_transport_mpas(ds, dsMesh, mpasTransportFileName)
    ds = xarray.open_dataset(mpasTransportFileName)

    zlevelTransportFileName = '{}/cache/osf_zlevel_transport.nc'.format(folder)
    _interpolate_horizontal_transport_zlevel(ds, z, zlevelTransportFileName)
    ds = xarray.open_dataset(zlevelTransportFileName)

    cumsumTransportFileName = '{}/cache/osf_cumsum_transport.nc'.format(folder)
    _vertical_cumsum_horizontal_transport(ds, cumsumTransportFileName)
    ds = xarray.open_dataset(cumsumTransportFileName)

    cacheFileName = '{}/cache/osf_vert_slice.nc'.format(folder)
    _horizontally_bin_overturning_streamfunction(ds, dsMesh, x, osfFileName,
                                                 cacheFileName)


def _compute_barotorpic_transport(dsMesh, ds):
    '''
    Compute the barotropic transport for inner edges (not on the domain
    boundary) for the given mesh and monthly-mean data set
    '''

    cellsOnEdge = dsMesh.cellsOnEdge - 1
    innerEdges = numpy.logical_and(cellsOnEdge[:, 0] >= 0,
                                   cellsOnEdge[:, 1] >= 0)

    # convert from boolean mask to indices
    innerEdges = numpy.nonzero(innerEdges.values)[0]

    cell0 = cellsOnEdge[innerEdges, 0]
    cell1 = cellsOnEdge[innerEdges, 1]

    layerThickness = ds.timeMonthly_avg_layerThickness.chunk(
        chunks={'Time': 1})
    normalVelocity = ds.timeMonthly_avg_normalVelocity[:, innerEdges, :].chunk(
        chunks={'Time': 1})

    layerThicknessEdge = 0.5*(layerThickness[:, cell0, :] +
                              layerThickness[:, cell1, :])
    transport = dsMesh.dvEdge[innerEdges] * \
        (layerThicknessEdge * normalVelocity).sum(dim='nVertLevels')

    return innerEdges, transport


def _compute_barotropic_streamfunction_vertex(dsMesh, ds):
    innerEdges, transport = _compute_barotorpic_transport(dsMesh, ds)

    nVertices = dsMesh.sizes['nVertices']
    nTime = ds.sizes['Time']

    cellsOnVertex = dsMesh.cellsOnVertex - 1
    verticesOnEdge = dsMesh.verticesOnEdge - 1
    boundaryVertices = numpy.logical_or(cellsOnVertex[:, 0] == -1,
                                        cellsOnVertex[:, 1] == -1)
    boundaryVertices = numpy.logical_or(boundaryVertices,
                                        cellsOnVertex[:, 2] == -1)

    # convert from boolean mask to indices
    boundaryVertices = numpy.nonzero(boundaryVertices.values)[0]

    nBoundaryVertices = len(boundaryVertices)
    nInnerEdges = len(innerEdges)

    indices = numpy.zeros((2, 2*nInnerEdges+nBoundaryVertices), dtype=int)
    data = numpy.zeros(2*nInnerEdges+nBoundaryVertices, dtype=float)

    # The difference between the streamfunction at vertices on an inner edge
    # should be equal to the transport
    v0 = verticesOnEdge[innerEdges, 0].values
    v1 = verticesOnEdge[innerEdges, 1].values

    ind = numpy.arange(nInnerEdges)
    indices[0, 2*ind] = ind
    indices[1, 2*ind] = v1
    data[2*ind] = 1.

    indices[0, 2*ind+1] = ind
    indices[1, 2*ind+1] = v0
    data[2*ind+1] = -1.

    # the streamfunction should be zero at all boundary vertices
    ind = numpy.arange(nBoundaryVertices)
    indices[0, 2*nInnerEdges + ind] = nInnerEdges + ind
    indices[1, 2*nInnerEdges + ind] = boundaryVertices
    data[2*nInnerEdges + ind] = 1.

    bsfVertex = xarray.DataArray(numpy.zeros((nTime, nVertices)),
                                 dims=('Time', 'nVertices'))

    widgets = ['barotropic streamfunction: ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=nTime).start()

    for tIndex in range(nTime):
        rhs = numpy.zeros(nInnerEdges+nBoundaryVertices, dtype=float)

        # convert to Sv
        ind = numpy.arange(nInnerEdges)
        rhs[ind] = 1e-6*transport.isel(Time=tIndex)

        ind = numpy.arange(nBoundaryVertices)
        rhs[nInnerEdges + ind] = 0.

        M = scipy.sparse.csr_matrix((data, indices),
                                    shape=(nInnerEdges+nBoundaryVertices,
                                           nVertices))

        solution = scipy.sparse.linalg.lsqr(M, rhs)

        bsfVertex[tIndex, :] = -solution[0]
        bar.update(tIndex+1)
    bar.finish()

    return bsfVertex


def _compute_barotropic_streamfunction_cell(dsMesh, bsfVertex):
    '''
    Interpolate the barotropic streamfunction from vertices to cells
    '''
    nEdgesOnCell = dsMesh.nEdgesOnCell
    edgesOnCell = dsMesh.edgesOnCell - 1
    verticesOnCell = dsMesh.verticesOnCell - 1
    areaEdge = dsMesh.dcEdge*dsMesh.dvEdge
    prevEdgesOnCell = edgesOnCell.copy(deep=True)
    prevEdgesOnCell[:, 1:] = edgesOnCell[:, 0:-1]
    prevEdgesOnCell[:, 0] = edgesOnCell[:, nEdgesOnCell-1]

    mask = verticesOnCell >= 0
    areaVert = mask*0.5*(areaEdge[edgesOnCell] + areaEdge[prevEdgesOnCell])

    bsfCell = ((areaVert * bsfVertex[:, verticesOnCell]).sum(dim='maxEdges') /
               areaVert.sum(dim='maxEdges'))

    return bsfCell


def _compute_horizontal_transport_mpas(ds, dsMesh, outFileName):
    '''
    compute the horizontal transport through edges on the native MPAS grid.
    '''

    if file_complete(ds, outFileName):
        return

    nVertLevels = dsMesh.sizes['nVertLevels']
    cellsOnEdge = dsMesh.cellsOnEdge - 1
    maxLevelCell = dsMesh.maxLevelCell - 1

    cell0 = cellsOnEdge[:, 0]
    cell1 = cellsOnEdge[:, 1]

    internalEdgeIndices = xarray.DataArray(
        numpy.nonzero(numpy.logical_and(cell0.values >= 0,
                                        cell1.values >= 0))[0],
        dims=('nInternalEdges',))

    cell0 = cell0[internalEdgeIndices]
    cell1 = cell1[internalEdgeIndices]

    bottomDepth = dsMesh.bottomDepth

    maxLevelEdgeTop = maxLevelCell[cell0]
    mask = numpy.logical_or(cell0 == -1,
                            maxLevelCell[cell1] < maxLevelEdgeTop)
    maxLevelEdgeTop[mask] = maxLevelCell[cell1][mask]

    nVertLevels = dsMesh.sizes['nVertLevels']

    vertIndex = \
        xarray.DataArray.from_dict({'dims': ('nVertLevels',),
                                    'data': numpy.arange(nVertLevels)})

    ds = ds.chunk({'Time': 1})

    chunks = {'nInternalEdges': 1024}
    maxLevelEdgeTop = maxLevelEdgeTop.chunk(chunks)
    dvEdge = dsMesh.dvEdge[internalEdgeIndices].chunk(chunks)
    bottomDepthEdge = 0.5*(bottomDepth[cell0] +
                           bottomDepth[cell1]).chunk(chunks)

    chunks = {'Time': 1, 'nInternalEdges': 1024}

    normalVelocity = ds.timeMonthly_avg_normalVelocity.isel(
        nEdges=internalEdgeIndices).chunk(chunks)
    layerThickness = ds.timeMonthly_avg_layerThickness.chunk()

    layerThicknessEdge = 0.5*(layerThickness.isel(nCells=cell0) +
                              layerThickness.isel(nCells=cell1)).chunk(chunks)

    layerThicknessEdge = layerThicknessEdge.where(
        vertIndex <= maxLevelEdgeTop, other=0.)

    thicknessSum = layerThicknessEdge.sum(dim='nVertLevels')

    thicknessCumSum = layerThicknessEdge.cumsum(dim='nVertLevels')

    zSurface = thicknessSum - bottomDepthEdge

    zInterfaceEdge = -thicknessCumSum + zSurface

    zInterfaceEdge = xarray.concat(
        [zSurface.expand_dims(dim='nVertLevelsP1', axis=2),
         zInterfaceEdge.rename({'nVertLevels': 'nVertLevelsP1'})],
        dim='nVertLevelsP1')

    transportPerDepth = dvEdge*normalVelocity

    dsOut = xarray.Dataset()
    dsOut['xtime_startMonthly'] = ds.xtime_startMonthly
    dsOut['xtime_endMonthly'] = ds.xtime_endMonthly
    dsOut['zInterfaceEdge'] = zInterfaceEdge
    dsOut['layerThicknessEdge'] = layerThicknessEdge
    dsOut['transportPerDepth'] = transportPerDepth
    dsOut['transportVertSum'] = \
        (transportPerDepth*layerThicknessEdge).sum(dim='nVertLevels')

    dsOut = dsOut.transpose('Time', 'nInternalEdges', 'nVertLevels',
                            'nVertLevelsP1')

    print('compute and caching transport on MPAS grid:')
    write_netcdf(dsOut, outFileName, progress=True)


def _interpolate_horizontal_transport_zlevel(ds, z, outFileName):
    '''
    interpolate the horizontal transport through edges onto a z-level grid.
    '''

    if file_complete(ds, outFileName):
        return

    ds = ds.chunk({'Time': 1, 'nInternalEdges': None, 'nVertLevels': 1,
                   'nVertLevelsP1': 1})

    nz = len(z)
    z = xarray.DataArray.from_dict({'dims': ('nz',), 'data': z})

    # make sure we don't miss anything
    z[0] = max(z[0].values, ds.zInterfaceEdge.max())
    z[-1] = min(z[-1].values, ds.zInterfaceEdge.min())

    z0 = z[0:-1].rename({'nz': 'nzM1'})
    z1 = z[1:].rename({'nz': 'nzM1'})

    nTime = ds.sizes['Time']
    nInternalEdges = ds.sizes['nInternalEdges']
    nVertLevels = ds.sizes['nVertLevels']

    widgets = ['interpolating tansport on z-level grid: ',
               progressbar.Percentage(), ' ', progressbar.Bar(), ' ',
               progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets, maxval=nTime).start()

    fileNames = []

    for tIndex in range(nTime):
        fileName = outFileName.replace('.nc', '_{}.nc'.format(tIndex))
        fileNames.append(fileName)
        if os.path.exists(fileName):
            continue

        outTransport = xarray.DataArray(
            numpy.zeros((nInternalEdges, nz-1)),
            dims=('nInternalEdges', 'nzM1'))

        dzSum = xarray.DataArray(
            numpy.zeros((nInternalEdges, nz-1)),
            dims=('nInternalEdges', 'nzM1'))

        dsIn = ds.isel(Time=tIndex)
        for inZIndex in range(nVertLevels):
            zTop = dsIn.zInterfaceEdge.isel(nVertLevelsP1=inZIndex)
            zBot = dsIn.zInterfaceEdge.isel(nVertLevelsP1=inZIndex+1)
            inTransportPerDepth = \
                dsIn.transportPerDepth.isel(nVertLevels=inZIndex)

            zt = numpy.minimum(zTop, z0)
            zb = numpy.maximum(zBot, z1)
            dz = numpy.maximum(zt - zb, 0.)

            outTransport = outTransport + dz*inTransportPerDepth

            dzSum = dzSum + dz

        outTransport.compute()
        dzSum.compute()

        dsOut = xarray.Dataset()
        dsOut['mask'] = dzSum > 0
        dsOut['transport'] = outTransport
        dsOut['transportVertSum'] = outTransport.sum('nzM1')
        dsOut['transportVertSumCheck'] = \
            dsIn.transportVertSum - dsOut.transportVertSum

        dsOut = dsOut.transpose('nzM1', 'nInternalEdges')

        write_netcdf(dsOut, fileName, progress=False)

        assert(numpy.abs(dsOut.transportVertSumCheck).max().values < 1e-9)

        bar.update(tIndex + 1)

    bar.finish()

    dsOut = xarray.open_mfdataset(fileNames, concat_dim='Time')

    dsOut['xtime_startMonthly'] = ds.xtime_startMonthly
    dsOut['xtime_endMonthly'] = ds.xtime_endMonthly
    dsOut['z'] = z

    dsOut = dsOut.transpose('Time', 'nzM1', 'nz', 'nInternalEdges')

    print('caching transport on z-level grid:')
    write_netcdf(dsOut, outFileName, progress=True)


def _vertical_cumsum_horizontal_transport(ds, outFileName):
    '''
    compute the cumsum in the vertical of the horizontal transport
    '''

    if file_complete(ds, outFileName):
        return

    chunks = {'Time': 1, 'nInternalEdges': 32768}
    ds = ds.chunk(chunks)

    nTime = ds.sizes['Time']
    nInternalEdges = ds.sizes['nInternalEdges']
    nz = ds.sizes['nz']

    transport = ds.transport.rename({'nzM1': 'nz'})

    transportSumTop = xarray.DataArray(
        numpy.zeros((nTime, nInternalEdges, 1)),
        dims=('Time', 'nInternalEdges', 'nz')).chunk(chunks)

    # zeros on top and then the cumsum for the rest
    transportSum = xarray.concat([transportSumTop,
                                  transport.cumsum(dim='nz')], dim='nz')

    # mask out locations on the output where no input-grid layers overlap
    # with either the output layer above or the one below
    mask = ds.mask.rename({'nzM1': 'nz'})
    maskTop = mask.isel(nz=0)
    maskBot = mask.isel(nz=nz-2)

    outMask = xarray.concat([maskTop,
                             numpy.logical_or(mask[:, 0:-1, :],
                                              mask[:, 1:, :]),
                             maskBot], dim='nz')

    dsOut = xarray.Dataset()
    dsOut['xtime_startMonthly'] = ds.xtime_startMonthly
    dsOut['xtime_endMonthly'] = ds.xtime_endMonthly
    dsOut['z'] = ds.z
    dsOut['transportSum'] = transportSum
    dsOut['mask'] = outMask

    dsOut = dsOut.transpose('Time', 'nz', 'nInternalEdges')

    print('compute and caching vertical transport sum on z-level grid:')
    write_netcdf(dsOut, outFileName, progress=True)


def _compute_region_boundary_edges(dsMesh, cellMask):
    '''
    Given a mask of cells in a region, find the indices and signs (indicating
    fluxes into the masked region) of non-boundary edges
    '''

    cellsOnEdge = dsMesh.cellsOnEdge.values - 1
    cellsOnEdgeMask = cellMask.values[cellsOnEdge]

    # first, we only want interior edges, not boundary edges
    edgeMask = numpy.logical_and(cellsOnEdge[:, 0] >= 0,
                                 cellsOnEdge[:, 1] >= 0)

    cellMask0 = cellsOnEdgeMask[:, 0][edgeMask]
    cellMask1 = cellsOnEdgeMask[:, 1][edgeMask]

    # second, we only want edges where one side is in the region and the other
    # is not
    edgeMask = cellMask0 != cellMask1

    # convert from boolean mask to indices
    edgeIndices = numpy.nonzero(edgeMask)[0]

    if(len(edgeIndices) == 0):
        return numpy.array([]), numpy.array([])
    else:
        # according to the mesh spec, normals point from cell 0 to cell 1 on a
        # given edge:
        # https://mpas-dev.github.io/files/documents/MPAS-MeshSpec.pdf
        edgeSign = numpy.ones(len(edgeIndices), int)

        # if the first cell on the edge is in the region, we need to reverse
        # the edge direction to point into the region
        mask = cellMask0[edgeIndices]
        edgeSign[mask] = -1

        return edgeIndices, edgeSign


def _horizontally_bin_overturning_streamfunction(ds, dsMesh, x, osfFileName,
                                                 cacheFileName):
    '''
    bin and sum the vertically cumsummed horizontal transport on the z-level
    grid to get the OSF.
    '''

    chunks = {'Time': 1}
    ds = ds.chunk(chunks)

    nTime = ds.sizes['Time']
    nz = ds.sizes['nz']
    nx = len(x)
    x = xarray.DataArray.from_dict({'dims': ('nx',), 'data': x})

    widgets = ['bin overturning streamfunction: ',
               progressbar.Percentage(), ' ', progressbar.Bar(), ' ',
               progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets, maxval=nx).start()

    # sum up the transport into the region bounded by each x value
    # on the output grid

    fileNames = []

    for xIndex in range(nx):
        fileName = cacheFileName.replace('.nc', '_{}.nc'.format(xIndex))
        fileNames.append(fileName)
        if file_complete(ds, fileName):
            continue

        cellMask = dsMesh.xCell >= x[xIndex]
        edgeIndices, edgeSigns = _compute_region_boundary_edges(dsMesh,
                                                                cellMask)

        if len(edgeIndices) == 0:

            localOSF = numpy.nan*xarray.DataArray(numpy.ones((nTime, nz)),
                                                  dims=('Time', 'nz'))
        else:
            # convert to Sv
            transportSum = 1e-6 * \
                edgeSigns*ds.transportSum.isel(nInternalEdges=edgeIndices)

            localOSF = transportSum.sum(dim='nInternalEdges')

            localMask = ds.mask.isel(nInternalEdges=edgeIndices).sum(
                dim='nInternalEdges') > 0
            localOSF = localOSF.where(localMask)

        dsOSF = xarray.Dataset()
        dsOSF['osf'] = localOSF
        write_netcdf(dsOSF, fileName, progress=False)

        bar.update(xIndex+1)

    bar.finish()

    dsOSF = xarray.open_mfdataset(fileNames, concat_dim='nx')
    dsOSF['xtime_startMonthly'] = ds.xtime_startMonthly
    dsOSF['xtime_endMonthly'] = ds.xtime_endMonthly
    dsOSF['x'] = x
    dsOSF['z'] = ds.z
    dsOSF.osf.attrs['units'] = 'Sv'
    dsOSF.osf.attrs['description'] = 'overturning streamfunction '

    dsOSF = dsOSF.transpose('Time', 'nz', 'nx')

    write_netcdf(dsOSF, osfFileName)
