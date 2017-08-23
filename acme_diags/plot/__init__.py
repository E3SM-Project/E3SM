"""Manages plotting, provides a single interface
for different plots with different backends."""
from __future__ import print_function

import importlib
from acme_diags.driver.utils import get_set_name

def _get_plot_fcn(backend, set_num):
    """Get the actual plot() function based on the backend and set_num."""
    try:
        if backend in ['matplotlib', 'mpl']:
            backend = 'cartopy'

        set_num = get_set_name(set_num) 
        mod_str = 'acme_diags.plot.{}.{}_plot'.format(backend, set_num)
        module = importlib.import_module(mod_str)
        return module.plot

    except ImportError as e:
        print(e)
        print('Plotting for set {} with {} is not supported'.format(set_num, backend))

def plot(set_num, ref, test, diff, metrics_dict, parameter):
    """Based on set_num and parameter.backend,
    call the correct plotting function."""
    if hasattr(parameter, 'plot'):
        parameter.plot(ref, test, diff, metrics_dict, parameter)
    else:
        if parameter.backend not in ['vcs', 'cartopy', 'mpl', 'matplotlib']:
            raise RuntimeError('Invalid backend, choose either "vcs" or "matplotlib"/"mpl"/"cartopy"')

        plot_fcn = _get_plot_fcn(parameter.backend, set_num)
        plot_fcn(ref, test, diff, metrics_dict, parameter)
