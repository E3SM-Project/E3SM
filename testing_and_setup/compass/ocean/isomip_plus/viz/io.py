import numpy
import netCDF4
from datetime import datetime
import sys
from dask.diagnostics import ProgressBar
import os
import xarray


def write_netcdf(ds, fileName, fillValues=netCDF4.default_fillvals,
                 progress=False):
    '''Write an xarray Dataset with NetCDF4 fill values where needed'''
    encodingDict = {}
    variableNames = list(ds.data_vars.keys()) + list(ds.coords.keys())
    for variableName in variableNames:
        isNumeric = numpy.issubdtype(ds[variableName].dtype, numpy.number)
        if isNumeric and numpy.any(numpy.isnan(ds[variableName])):
            dtype = ds[variableName].dtype
            for fillType in fillValues:
                if dtype == numpy.dtype(fillType):
                    encodingDict[variableName] = \
                        {'_FillValue': fillValues[fillType]}
                    break
        else:
            encodingDict[variableName] = {'_FillValue': None}

    update_history(ds)

    if progress:
        delayed_obj = ds.to_netcdf(fileName, encoding=encodingDict,
                                   compute=False)
        with ProgressBar():
            delayed_obj.compute()
    else:
        ds.to_netcdf(fileName, encoding=encodingDict)


def update_history(ds):
    '''Add or append history to attributes of a data set'''

    thiscommand = datetime.now().strftime("%a %b %d %H:%M:%S %Y") + ": " + \
        " ".join(sys.argv[:])
    if 'history' in ds.attrs:
        newhist = '\n'.join([thiscommand, ds.attrs['history']])
    else:
        newhist = thiscommand
    ds.attrs['history'] = newhist


def file_complete(ds, fileName):
    '''
    Find out if the file already has the same number of time slices as the
    monthly-mean data set
    '''
    complete = False
    if os.path.exists(fileName):
        with xarray.open_dataset(fileName) as dsCompare:
            if ds.sizes['Time'] == dsCompare.sizes['Time']:
                complete = True

    return complete
