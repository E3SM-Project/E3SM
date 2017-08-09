"""Manages plotting, provides a single interface
for different plots with different backends."""

def _get_plot_fcn(backend, set_num):
    """Get the actual plot() function based on the backend and set_num."""
    set_num = str(set_num)
    try:
        import importlib
        if backend in ['matplotlib', 'mpl']:
            backend = 'cartopy'

        if set_num in ['zonal_mean_xy', '3']:
            mod_str = 'acme_diags.plot.{}.zonal_mean_xy_plot'.format(backend)
        elif set_num in ['zonal_mean_2d', '4']:
            mod_str = 'acme_diags.plot.{}.zonal_mean_2d_plot'.format(backend)
        elif set_num in ['lat_lon', '5']:
            mod_str = 'acme_diags.plot.{}.lat_lon_plot'.format(backend)
        elif set_num in ['polar', '7']:
            mod_str = 'acme_diags.plot.{}.polar_plot'.format(backend)
        elif set_num in ['cosp_histogram', '13']:
            mod_str = 'acme_diags.plot.{}.cosp_histogram_plot'.format(backend)
        else:
            raise ImportError('There is no module {} for the {} backend'.format(set_num, backend))

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
