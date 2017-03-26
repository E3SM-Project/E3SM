import os
import sys
import numpy
import cdutil
import vcs
import cdms2
import genutil.statistics
from acme_diags.metrics import rmse, corr, min_cdms, max_cdms, mean


def plot_min_max_mean(canvas, metrics_dict, ref_test_or_diff):
    """ canvas is a vcs.Canvas, metrics_dict is a dict and 
    ref_test_or_diff is a string """
    var_min = '%.2f' % metrics_dict[ref_test_or_diff]['min']
    var_max = '%.2f' % metrics_dict[ref_test_or_diff]['max']
    var_mean = '%.2f' % metrics_dict[ref_test_or_diff]['mean']

    if ref_test_or_diff == 'ref':  # Remove this when vcdat is done
        ref_test_or_diff = 'reference'
 
    # can be either 'reference', 'test' or 'diff'
    plot = ref_test_or_diff
    min_label = canvas.createtextcombined(Tt_source = plot + '_min_label',
                                          To_source = plot + '_min_label')
    max_label = canvas.createtextcombined(Tt_source = plot + '_max_label',
                                          To_source = plot + '_max_label')
    mean_label = canvas.createtextcombined(Tt_source = plot + '_mean_label',
                                           To_source = plot + '_mean_label')

    min_value = canvas.createtextcombined(Tt_source = plot + '_min_value',
                                          To_source = plot + '_min_value')
    max_value = canvas.createtextcombined(Tt_source = plot + '_max_value',
                                          To_source = plot + '_max_value')
    mean_value = canvas.createtextcombined(Tt_source = plot + '_mean_value',
                                           To_source = plot + '_mean_value')

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
    """ canvas is a vcs.Canvas, metrics_dict is a dict """

    rmse_str = '%.2f' % metrics_dict['misc']['rmse']
    corr_str = '%.2f' % metrics_dict['misc']['corr']

    rmse_label = canvas.createtextcombined(Tt_source = 'diff_plot_comment1_title',
                                           To_source = 'diff_plot_comment1_title')
    corr_label = canvas.createtextcombined(Tt_source = 'diff_plot_comment2_title',
                                           To_source = 'diff_plot_comment2_title')
    rmse_label.string = 'RMSE'
    corr_label.string = 'CORR'

    rmse_value = canvas.createtextcombined(Tt_source = 'diff_plot_comment1_value',
                                           To_source = 'diff_plot_comment1_value')
    corr_value = canvas.createtextcombined(Tt_source = 'diff_plot_comment2_value',
                                           To_source = 'diff_plot_comment2_value')

    rmse_value.string = rmse_str
    corr_value.string = corr_str

    canvas.plot(rmse_label)
    canvas.plot(corr_label)
    canvas.plot(rmse_value)
    canvas.plot(corr_value)

def set_colormap_of_graphics_method(canvas, parameter_colormap, method):
    if parameter_colormap is not '':
        method.colormap = vcs.getcolormap(parameter_colormap)
        colors = vcs.getcolors(method.levels, colors=range(6, 240))
        method.fillareacolors = colors

def set_levels_of_graphics_method(method, levels, data):
    if levels != []:
        method.levels = levels

    if method.levels == [[1.0000000200408773e+20, 1.0000000200408773e+20]]:
        method.levels = vcs.mkscale(data.min(), data.max())

def set_units(ref_or_test, units):
    if units != '':
        ref_or_test.units = units

def plot(reference, test, diff, metrics_dict, parameter):

    case_id = parameter.case_id
    if not os.path.exists(case_id):
        os.makedirs(case_id)

    # Plotting
    vcs_canvas = vcs.init(bg=True, geometry=(parameter.canvas_size_w, parameter.canvas_size_h))
    if not parameter.logo:
        vcs_canvas.drawlogooff()

    file_path = os.path.abspath(os.path.dirname(__file__))
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


    set_levels_of_graphics_method(reference_isofill, parameter.reference_levels, reference)
    set_levels_of_graphics_method(test_isofill, parameter.test_levels, test)
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

    vcs_canvas.plot(test, template_test, test_isofill)
    vcs_canvas.plot(reference, template_ref, reference_isofill)
    vcs_canvas.plot(diff, template_diff, diff_isofill)

    plot_rmse_and_corr(vcs_canvas, metrics_dict)

    # Plotting the main title
    main_title = vcs_canvas.createtextcombined(Tt_source = 'main_title',
                                               To_source = 'main_title')
    main_title.string = parameter.main_title
    vcs_canvas.plot(main_title)

    vcs_canvas.png(case_id + '/' + parameter.output_file)
    print 'Plot saved in: ' + case_id + '/' + parameter.output_file
    
    # Save files being plotted
    # Set cdms preferences - no compression, no shuffling, no complaining
    cdms2.setNetcdfDeflateFlag(1)
    cdms2.setNetcdfDeflateLevelFlag(0) ; # 1-9, min to max - Comes at heavy IO (read/write time cost)
    cdms2.setNetcdfShuffleFlag(0)
    cdms2.setCompressionWarnings(0) ; # Turn off warning messages
    # Save test file
    file_test = cdms2.open(case_id + '/' + parameter.output_file + '_test.nc','w+')
    test.id = parameter.variables
    file_test.write(test) 
    file_test.close()

    # Save reference file
    file_ref = cdms2.open(case_id + '/' + parameter.output_file + '_ref.nc','w+')
    reference.id = parameter.variables
    file_ref.write(reference) 
    file_ref.close()

    # Save difference file
    file_diff = cdms2.open(case_id + '/' + parameter.output_file + '_diff.nc','w+')
    diff.id = parameter.variables+'(test - reference)'
    file_diff.write(diff) 
    file_diff.close()


