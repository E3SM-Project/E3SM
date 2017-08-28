from __future__ import print_function

import os
import sys
import numpy
import cdutil
import vcs
import cdms2
import genutil.statistics
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean
from acme_diags.driver.utils import get_output_dir, _chown

textcombined_objs = {}

def managetextcombined(tt_name, to_name, vcs_canvas):
    """Caches textcombined objects"""
    new_name = "%s:::%s" % (tt_name, to_name)
    mytc = textcombined_objs.get(new_name, None)
    if mytc is None:
        mytc = vcs_canvas.createtextcombined(Tt_source=tt_name, To_source=to_name)
        textcombined_objs[new_name] = mytc
    return mytc

def plot_min_max_mean(canvas, metrics_dict, ref_test_or_diff):
    """canvas is a vcs.Canvas, metrics_dict is a dict and
    ref_test_or_diff is a string."""
    var_min = '%.2f' % metrics_dict[ref_test_or_diff]['min']
    var_max = '%.2f' % metrics_dict[ref_test_or_diff]['max']
    var_mean = '%.2f' % metrics_dict[ref_test_or_diff]['mean']

    if ref_test_or_diff == 'ref':  # Remove this when vcdat is done
        ref_test_or_diff = 'reference'
 
    # can be either 'reference', 'test' or 'diff'
    plot = ref_test_or_diff
    min_label = managetextcombined(plot + '_min_label', plot + '_min_label', canvas)
    max_label = managetextcombined(plot + '_max_label', plot + '_max_label', canvas)
    mean_label = managetextcombined(plot + '_mean_label', plot + '_mean_label', canvas)

    min_value = managetextcombined(plot + '_min_value', plot + '_min_value', canvas)
    max_value = managetextcombined(plot + '_max_value', plot + '_max_value', canvas)
    mean_value = managetextcombined(plot + '_mean_value', plot + '_mean_value', canvas)

    min_value.string = var_min
    max_value.string = var_max
    mean_value.string = var_mean
    canvas.plot(min_value)
    canvas.plot(min_label)
    canvas.plot(max_value)
    canvas.plot(max_label)
    canvas.plot(mean_value)
    canvas.plot(mean_label)

def plot_rmse_and_corr(canvas, metrics_dict):
    """canvas is a vcs.Canvas, metrics_dict is a dict"""

    rmse_str = '%.2f' % metrics_dict['misc']['rmse']
    corr_str = '%.2f' % metrics_dict['misc']['corr']

    rmse_label = managetextcombined('diff_plot_comment1_title', 'diff_plot_comment1_title', canvas)
    corr_label = managetextcombined('diff_plot_comment2_title', 'diff_plot_comment2_title', canvas)
    rmse_label.string = 'RMSE'
    corr_label.string = 'CORR'

    rmse_value = managetextcombined('diff_plot_comment1_value', 'diff_plot_comment1_value', canvas)
    corr_value = managetextcombined('diff_plot_comment2_value', 'diff_plot_comment2_value', canvas)

    rmse_value.string = rmse_str
    corr_value.string = corr_str

    canvas.plot(rmse_label)
    canvas.plot(corr_label)
    canvas.plot(rmse_value)
    canvas.plot(corr_value)

def set_colormap_of_graphics_method(canvas, parameter_colormap, method):
    if parameter_colormap is not '':
        if parameter_colormap in vcs.listelements('colormap'):
            method.colormap = vcs.getcolormap(parameter_colormap)
        else:
            method.colormap = vcs.matplotlib2vcs(parameter_colormap)
        colors = vcs.getcolors(method.levels, colors=range(6, 240), white=[100, 100, 100])
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

def plot(reference, test, diff, metrics_dict, parameter):
    vcs_canvas = vcs.init(bg=True, geometry=(parameter.canvas_size_w, parameter.canvas_size_h))
    case_id = parameter.case_id


    file_path = os.path.join(sys.prefix, 'share', 'acme_diags', 'lat_lon')
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_5.json'))
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_5_new.json'))

    template_test = vcs_canvas.gettemplate('plotset5_0_x_0')
    template_ref = vcs_canvas.gettemplate('plotset5_0_x_1')
    template_diff = vcs_canvas.gettemplate('plotset5_0_x_2')

    set_units(test, parameter.test_units)
    set_units(reference, parameter.reference_units)
    set_units(diff, parameter.diff_units)

    test.long_name = parameter.test_title
    reference.long_name = parameter.reference_title
    diff.long_name = parameter.diff_title

    test.id = parameter.test_name
    reference.id = parameter.reference_name
    diff.id = parameter.diff_name

    # model and observation graph
    plot_min_max_mean(vcs_canvas, metrics_dict, 'test')
    plot_min_max_mean(vcs_canvas, metrics_dict, 'ref')
    plot_min_max_mean(vcs_canvas, metrics_dict, 'diff')

    reference_isofill = vcs.getisofill('reference_isofill')
    reference_isofill.missing = 'grey'
    test_isofill = vcs.getisofill('test_isofill')
    test_isofill.missing = 'grey'
    diff_isofill = vcs.getisofill('diff_isofill')
    diff_isofill.missing = 'grey'


    set_levels_of_graphics_method(reference_isofill, parameter.contour_levels, reference, test)
    set_levels_of_graphics_method(test_isofill, parameter.contour_levels, test, reference)
    set_levels_of_graphics_method(diff_isofill, parameter.diff_levels, diff)

    if parameter.arrows:
        reference_isofill.ext_1 = True
        reference_isofill.ext_2 = True
        test_isofill.ext_1 = True
        test_isofill.ext_2 = True
        diff_isofill.ext_1 = True
        diff_isofill.ext_2 = True

    set_colormap_of_graphics_method(vcs_canvas, parameter.reference_colormap, reference_isofill)
    set_colormap_of_graphics_method(vcs_canvas, parameter.test_colormap, test_isofill)
    set_colormap_of_graphics_method(vcs_canvas, parameter.diff_colormap, diff_isofill)

    vcs_canvas.plot(add_cyclic(test), template_test, test_isofill)
    vcs_canvas.plot(add_cyclic(reference), template_ref, reference_isofill)
    vcs_canvas.plot(add_cyclic(diff), template_diff, diff_isofill)

    plot_rmse_and_corr(vcs_canvas, metrics_dict)
    vcs_canvas.portrait()  # for some reason, this needs to be before a call to vcs_canvas.plot()

    # Plotting the main title
    main_title = managetextcombined('main_title', 'main_title', vcs_canvas)
    main_title.string = parameter.main_title
    vcs_canvas.plot(main_title)

    if not parameter.logo:
        vcs_canvas.drawlogooff()

    fnm = os.path.join(get_output_dir(parameter.current_set, parameter), parameter.output_file)
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        if f == 'png':
            vcs_canvas.png(fnm)
            _chown(fnm + '.png', parameter.user)
        elif f == 'pdf':
            vcs_canvas.pdf(fnm)
            _chown(fnm + '.pdf', parameter.user)
        elif f == 'svg':
            vcs_canvas.svg(fnm)
            _chown(fnm + '.svg', parameter.user)
        print('Plot saved in: ' + fnm + '.' + f)
    vcs_canvas.clear()
    
