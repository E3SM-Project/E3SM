"""Manages plotting, provides a single interface
for different plots with different backends."""
from __future__ import print_function, absolute_import

import os
import sys
import importlib
import traceback
import numpy
from matplotlib.colors import LinearSegmentedColormap
from vcs.colors import matplotlib2vcs
import acme_diags


def _get_plot_fcn(backend, set_name):
    """Get the actual plot() function based on the backend and set_name."""
    try:
        if backend in ['matplotlib', 'mpl']:
            backend = 'cartopy'

        mod_str = 'acme_diags.plot.{}.{}_plot'.format(backend, set_name)
        module = importlib.import_module(mod_str)
        return module.plot

    except ModuleNotFoundError:
        print(
            'Plotting for set {} with {} is not supported'.format(
                set_name, backend))
        traceback.print_exc()

def plot(set_name, ref, test, diff, metrics_dict, parameter):
    """Based on set_name and parameter.backend,
    call the correct plotting function."""
    if hasattr(parameter, 'plot'):
        parameter.plot(ref, test, diff, metrics_dict, parameter)
    else:
        if parameter.backend not in ['vcs', 'cartopy', 'mpl', 'matplotlib']:
            raise RuntimeError(
                'Invalid backend, choose either "vcs" or "matplotlib"/"mpl"/"cartopy"')

        plot_fcn = _get_plot_fcn(parameter.backend, set_name)
        if plot_fcn:
            try:
                plot_fcn(ref, test, diff, metrics_dict, parameter)
            except Exception as e:
                print('Error while plotting {} with backend {}'.format(set_name, parameter.backend))
                traceback.print_exc()
                if parameter.debug:
                    sys.exit()

def get_colormap(colormap, parameters):
    """Get the colormap (string, list for vcs, or mpl colormap obj), which can be
    loaded from a local file in the cwd, installed file, or a predefined mpl/vcs one."""
    colormap = str(
        colormap)  # unicode don't seem to work well with string.endswith()
    if not colormap.endswith('.rgb'):  # predefined vcs/mpl colormap
        return colormap

    installed_colormap = os.path.join(acme_diags.INSTALL_PATH, 'colormaps', colormap)

    if os.path.exists(colormap):
        # colormap is an .rgb in the current directory
        pass
    elif not os.path.exists(colormap) and os.path.exists(installed_colormap):
        # use the colormap from /plot/colormaps
        colormap = installed_colormap
    elif not os.path.exists(colormap) and not os.path.exists(installed_colormap):
        pth = os.path.join(acme_diags.INSTALL_PATH, 'colormaps')
        msg = "File {} isn't in the current working directory or installed in {}"
        raise IOError(msg.format(colormap, pth))

    rgb_arr = numpy.loadtxt(colormap)
    rgb_arr = rgb_arr / 255.0

    if parameters.backend in ['cartopy', 'mpl', 'matplotlib']:
        cmap = LinearSegmentedColormap.from_list(name=colormap, colors=rgb_arr)
        return cmap

    elif parameters.backend in ['vcs']:
        n_levels = 240
        cmap = LinearSegmentedColormap.from_list(name=colormap, colors=rgb_arr, N=n_levels)
        vcs_cmap = matplotlib2vcs(cmap, vcs_name=colormap)

        return vcs_cmap, list(range(n_levels))

    else:
        raise RuntimeError('Invalid backend: {}'.format(parameters.backend))
