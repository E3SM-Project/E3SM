import vcs
from acme_diags.plot import get_colormap

def plot_min_max_mean(tc_dict, canvas, metrics_dict, ref_test_or_diff):
    """Based on ref_test_or_diff, plot the metrics on that plot."""
    var_min = '%.2f' % metrics_dict[ref_test_or_diff]['min']
    var_max = '%.2f' % metrics_dict[ref_test_or_diff]['max']
    var_mean = '%.2f' % metrics_dict[ref_test_or_diff]['mean']

    if ref_test_or_diff == 'ref':  # Remove this when vcdat is done
        ref_test_or_diff = 'reference'
 
    # can be either 'reference', 'test' or 'diff'
    plot = ref_test_or_diff
    min_label = managetextcombined(tc_dict, plot + '_min_label', plot + '_min_label', canvas)
    max_label = managetextcombined(tc_dict, plot + '_max_label', plot + '_max_label', canvas)
    mean_label = managetextcombined(tc_dict, plot + '_mean_label', plot + '_mean_label', canvas)

    min_value = managetextcombined(tc_dict, plot + '_min_value', plot + '_min_value', canvas)
    max_value = managetextcombined(tc_dict, plot + '_max_value', plot + '_max_value', canvas)
    mean_value = managetextcombined(tc_dict, plot + '_mean_value', plot + '_mean_value', canvas)

    min_value.string = var_min
    max_value.string = var_max
    mean_value.string = var_mean
    canvas.plot(min_value)
    canvas.plot(min_label)
    canvas.plot(max_value)
    canvas.plot(max_label)
    canvas.plot(mean_value)
    canvas.plot(mean_label)

def plot_rmse_and_corr(tc_dict, canvas, metrics_dict):
    """Plot the RMSE and CORR metrics."""

    rmse_str = '%.2f' % metrics_dict['misc']['rmse']
    corr_str = '%.2f' % metrics_dict['misc']['corr']

    rmse_label = managetextcombined(tc_dict, 'diff_plot_comment1_title', 'diff_plot_comment1_title', canvas)
    corr_label = managetextcombined(tc_dict, 'diff_plot_comment2_title', 'diff_plot_comment2_title', canvas)
    rmse_label.string = 'RMSE'
    corr_label.string = 'CORR'

    rmse_value = managetextcombined(tc_dict, 'diff_plot_comment1_value', 'diff_plot_comment1_value', canvas)
    corr_value = managetextcombined(tc_dict, 'diff_plot_comment2_value', 'diff_plot_comment2_value', canvas)

    rmse_value.string = rmse_str
    corr_value.string = corr_str

    canvas.plot(rmse_label)
    canvas.plot(corr_label)
    canvas.plot(rmse_value)
    canvas.plot(corr_value)

def get_color_range(gm):
    """For some graphic method colors, it needs a range of (16, 240)"""
    if gm.colormap in ['AMIP',
        'NCAR',
        'bl_to_darkred',
        'bl_to_drkorang',
        'blends',
        'blue2darkorange',
        'blue2darkred',
        'blue2green',
        'blue2grey',
        'blue2orange',
        'blue2orange2red',
        'blue_to_grey',
        'blue_to_grn',
        'blue_to_orange',
        'blue_to_orgred',
        'brown2blue',
        'brown_to_blue',
        'categorical',
        'classic',
        'green2magenta',
        'grn_to_magenta',
        'inferno',
        'lightblue2darkblue',
        'ltbl_to_drkbl',
        'rainbow',
        'rainbow_no_grn',
        'rainbownogreen',
        'sequential',
        'white2blue',
        'white2green',
        'white2magenta',
        'white2red',
        'white2yellow',
        'white_to_blue',
        'white_to_green',
        'white_to_magenta',
        'white_to_red',
        'white_to_yellow']:
        return range(16, 240)
    else:
        return range(256)

def set_colormap_of_graphics_method(canvas, parameter_colormap, method, parameters):
    if parameter_colormap is not '':
        cmap = get_colormap(parameter_colormap, parameters)
        if isinstance(cmap, str):
            if cmap in vcs.listelements('colormap'):
                method.colormap = vcs.getcolormap(cmap)
            else:
                method.colormap = vcs.matplotlib2vcs(cmap)
            cmap_range = get_color_range(method)
        else:
            # cmap is of type vcs.colormap.Cp
            cmap, cmap_range = cmap
            method.colormap = cmap
        colors = vcs.getcolors(method.levels, colors=cmap_range, white=[100, 100, 100], split=0)
        method.fillareacolors = colors

def set_levels_of_graphics_method(method, levels, data, data2=None):
    if levels != []:
        method.levels = levels

    if method.levels == [[1.0000000200408773e+20, 1.0000000200408773e+20]]:
        if data2 is None:
            method.levels = vcs.mkscale(data.min(), data.max())
        else:
            method.levels = vcs.mkscale(min(data.min(), data2.min()), max(data.max(), data2.max()))

def set_units(ref_or_test, units):
    if units != '':
        ref_or_test.units = units

def add_cyclic(var):
    lon = var.getLongitude()
    return var(longitude=(lon[0],lon[0]+360.0,'coe'))

def managetextcombined(tc_dict, tt_name, to_name, vcs_canvas):
    """Caches textcombined objects"""
    new_name = "%s:::%s" % (tt_name, to_name)
    mytc = tc_dict.get(new_name, None)
    if mytc is None:
        mytc = vcs_canvas.createtextcombined(Tt_source=tt_name, To_source=to_name)
        tc_dict[new_name] = mytc
    return mytc
