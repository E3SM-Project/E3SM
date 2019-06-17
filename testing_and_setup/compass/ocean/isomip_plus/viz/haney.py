import xarray
import numpy
import progressbar
import os

from viz.io import write_netcdf, file_complete


def compute_haney_number(dsMesh, ds, folder):
    '''
    compute the Haney number rx1 for each edge, and interpolate it to cells.
    '''

    haneyFileName = '{}/haney.nc'.format(folder)
    if file_complete(ds, haneyFileName):
        return

    haneyEdge, haneyCell = _compute_haney_number(dsMesh, ds)
    dsHaney = xarray.Dataset()
    dsHaney['xtime_startMonthly'] = ds.xtime_startMonthly
    dsHaney['xtime_endMonthly'] = ds.xtime_endMonthly
    dsHaney['haneyEdge'] = haneyEdge
    dsHaney.haneyEdge.attrs['units'] = 'unitless'
    dsHaney.haneyEdge.attrs['description'] = 'Haney number on edges'
    dsHaney['haneyCell'] = haneyCell
    dsHaney.haneyCell.attrs['units'] = 'unitless'
    dsHaney.haneyCell.attrs['description'] = 'Haney number on cells'
    dsHaney = dsHaney.transpose('Time', 'nCells', 'nEdges', 'nVertLevels')
    write_netcdf(dsHaney, haneyFileName)


def _compute_haney_number(dsMesh, ds):

    nEdges = dsMesh.sizes['nEdges']
    nCells = dsMesh.sizes['nCells']
    nVertLevels = dsMesh.sizes['nVertLevels']
    nTime = ds.sizes['Time']

    cellsOnEdge = dsMesh.cellsOnEdge - 1
    maxLevelCell = dsMesh.maxLevelCell - 1
    edgesOnCell = dsMesh.edgesOnCell - 1

    internalMask = numpy.logical_and(cellsOnEdge[:, 0] >= 0,
                                     cellsOnEdge[:, 1] >= 1)

    cell0 = cellsOnEdge[:, 0]
    cell1 = cellsOnEdge[:, 1]

    maxLevelEdge = numpy.minimum(maxLevelCell[cell0], maxLevelCell[cell1])
    vertIndex = \
        xarray.DataArray.from_dict({'dims': ('nVertLevels',),
                                    'data': numpy.arange(nVertLevels)})

    edgeMask = vertIndex <= maxLevelEdge

    cell0 = cell0[internalMask]
    cell1 = cell1[internalMask]

    haneyEdge = xarray.DataArray(numpy.zeros((nTime, nEdges, nVertLevels)),
                                 dims=('Time', 'nEdges', 'nVertLevels'))

    haneyCell = xarray.DataArray(numpy.zeros((nTime, nCells, nVertLevels)),
                                 dims=('Time', 'nCells', 'nVertLevels'))

    widgets = ['Haney number: ', progressbar.Percentage(), ' ',
               progressbar.Bar(), ' ', progressbar.ETA()]
    bar = progressbar.ProgressBar(widgets=widgets,
                                  maxval=nTime).start()

    for tIndex in range(nTime):

        zMid = numpy.zeros((nCells, nVertLevels+1))
        layerThickness = ds.timeMonthly_avg_layerThickness.isel(
            Time=tIndex).values
        ssh = ds.timeMonthly_avg_ssh.isel(Time=tIndex).values
        bottomDepth = dsMesh.bottomDepth.values
        zBot = -bottomDepth
        for zIndex in range(nVertLevels-1, -1, -1):
            zMid[:, zIndex+1] = zBot + 0.5*layerThickness[:, zIndex]
            zBot += layerThickness[:, zIndex]
        zMid[:, 0] = ssh

        dzVert1 = zMid[cell0, 0:-1] - zMid[cell0, 1:]
        dzVert2 = zMid[cell1, 0:-1] - zMid[cell1, 1:]
        dzEdge = zMid[cell1, :] - zMid[cell0, :]

        dzVert1[:, 0] *= 2
        dzVert2[:, 0] *= 2

        rx1 = numpy.zeros((nEdges, nVertLevels))

        rx1[internalMask, :] = (numpy.abs(dzEdge[:, 0:-1] + dzEdge[:, 1:]) / 
	                     (dzVert1 + dzVert2))

        haneyEdge[tIndex, :, :] = rx1
        haneyEdge[tIndex, :, :] = haneyEdge[tIndex, :, :].where(edgeMask)
        haneyCell[tIndex, :, :] = haneyEdge[tIndex, edgesOnCell, :].max(
            dim='maxEdges')
        bar.update(tIndex+1)
    bar.finish()

    return haneyEdge, haneyCell


