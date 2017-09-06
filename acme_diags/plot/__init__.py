"""Manages plotting, provides a single interface
for different plots with different backends."""
from __future__ import print_function, absolute_import

import os
import sys
import importlib
import numpy
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


def get_colormap(colormap, parameters):
    """Get the colormap (string, list for vcs, or mpl colormap obj), which can be
    loaded from a local file in the cwd, installed file, or a predefined mpl/vcs one."""
    if not colormap.endswith('.rgb'):  # predefined vcs/mpl colormap
        return colormap

    installed_colormap = os.path.join(sys.prefix, 'share', 'acme_diags', 'colormaps', colormap)

    if os.path.exists(colormap):
        # colormap is an .rgb in the current directory
        pass
    elif not os.path.exists(colormap) and os.path.exists(installed_colormap):
        # use the colormap from /plot/colormaps
        colormap = installed_colormap
    elif not os.path.exists(colormap) and not os.path.exists(installed_colormap):
        pth = os.path.join(sys.prefix, 'share', 'acme_diags', 'colormaps')
        msg = "File {} isn't in the current working directory or installed in {}"
        raise IOError(msg.format(colormap, pth))

    rgb_arr = numpy.loadtxt(colormap)
    rgb_arr = rgb_arr / 255.0

    if parameters.backend in ['cartopy', 'mpl', 'matplotlib']:
        from matplotlib.colors import LinearSegmentedColormap

        cmap = LinearSegmentedColormap.from_list(name=colormap, colors=rgb_arr)
        return cmap

    elif parameters.backend in ['vcs']:
        from vcs import createcolormap

        cmap = createcolormap()
        cmap_range = range(len(rgb_arr))
        for i in cmap_range:
            r, g, b = rgb_arr[i]
            cmap.setcolorcell(i, r, g, b, 100)
        return cmap, cmap_range

    else:
        raise RuntimeError('Invalid backend: {}'.format(parameters.backend))
