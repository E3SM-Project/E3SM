from __future__ import print_function

from e3sm_diags.plot.cartopy.zonal_mean_2d_plot import plot as base_plot


def plot(reference, test, diff, metrics_dict, parameter):
    return base_plot(reference, test, diff, metrics_dict, parameter)
