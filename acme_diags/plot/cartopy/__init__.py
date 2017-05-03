def get_plot_set(set_num, ref, test, diff, metrics_dict, parameter):
    if str(set_num) == '5':
        from . import set5cartopy as setplot
    if str(set_num) == '7':
        from . import set7cartopy as setplot
    else:
        raise RuntimeError('Plot set {} not yet supported for cartopy'.format(set_num))
    setplot.plot(ref, test, diff, metrics_dict, parameter)