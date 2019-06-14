import numpy
import genutil
import cdutil


def corr(model, obs, axis='xy'):
    corr = -numpy.infty
    try:
        corr = float(genutil.statistics.correlation(
            model, obs, axis=axis, weights='generate'))
    except Exception as err:
        print(err)

    return corr


def max_cdms(variable):
    return float(variable.max())


def mean(variable, axis='xy'):
    return cdutil.averager(variable, axis=axis, weights='generate')


def min_cdms(variable):
    return float(variable.min())


def rmse(model, obs, axis='xy'):
    rmse = -numpy.infty
    try:
        rmse = float(genutil.statistics.rms(
            model, obs, axis=axis, weights='generate'))
    except Exception as err:
        print(err)
    return rmse


def std(variable, axis='xy'):
    std = -numpy.infty
    try:
        std = float(genutil.statistics.std(
            variable, axis=axis, weights='generate'))
    except Exception as err:
        print(err)

    return std
