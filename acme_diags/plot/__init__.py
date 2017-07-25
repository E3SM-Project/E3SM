"""Manages plotting, provides a single interface
for different plots with different backends."""

def _get_plot_fcn(backend, set_num):
    """Get the actual plot() function based on the backend and set_num."""
    try:
        import importlib
        if backend == 'matplotlib' or backend == 'mpl':
            backend = 'cartopy'
        if set_num == 'zonal_mean_xy':
            set_num = '3'
        if set_num == 'zonal_mean_2d':
            set_num = '4'
        if set_num == 'lat_lon':
            set_num = '5'
        if set_num == 'polar':
            set_num = '7'

        mod_str = 'acme_diags.plot.{}.set{}'.format(backend, set_num)
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
