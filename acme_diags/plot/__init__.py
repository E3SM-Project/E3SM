"""Manages plotting, provides a single interface
for different plots with different backends."""

def _get_plot_fcn(backend, set_num):
    """Get the actual plot() function based on the backend and set_num."""
    try:
        import importlib
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
        if not 1 <= int(set_num) <= 15:
            raise RuntimeError('set_num must be between 1 and 15')
        if not parameter.backend == 'vcs' and not parameter.backend == 'cartopy':
            raise RuntimeError('Invalid backend, choose either "vcs" or "cartopy"')

        plot_fcn = _get_plot_fcn(parameter.backend, set_num)
        plot_fcn(ref, test, diff, metrics_dict, parameter)
