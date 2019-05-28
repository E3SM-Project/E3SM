from __future__ import print_function

import os
import sys
import numpy
import vcs
import acme_diags
import acme_diags.plot.vcs as utils
from acme_diags.driver.utils.general import get_output_dir

textcombined_objs = {}


def rotate_180(data):
    """Rotate the data 180 degrees."""
    return data[::-1]


def log_yaxis(data, method):
    """Make the y axis logrithmic."""
    axis = data.getAxis(-2)
    new_axis = axis.clone()
    new_axis[:] = numpy.ma.log10(new_axis[:])
    data.setAxis(-2, new_axis)

    normal_labels = vcs.mkscale(axis[0], axis[-1])
    new_labels = {}
    for l in normal_labels:
        new_labels[numpy.log10(l)] = int(l)
    method.yticlabels1 = new_labels

    # needed until cdms2 supports in-place operations
    return data


def plot(reference, test, diff, metrics_dict, parameter):
    vcs_canvas = vcs.init(bg=True, geometry=(
        parameter.canvas_size_w, parameter.canvas_size_h))
    parameter.case_id

    reference = rotate_180(reference)
    test = rotate_180(test)
    diff = rotate_180(diff)

    file_path = os.path.join(acme_diags.INSTALL_PATH, 'zonal_mean_2d')
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_4.json'))
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_4_new.json'))

    template_test = vcs_canvas.gettemplate('plotset4_0_x_0')
    template_ref = vcs_canvas.gettemplate('plotset4_0_x_1')
    template_diff = vcs_canvas.gettemplate('plotset4_0_x_2')

    # Turn off the units of the axes in the plots.
    template_test.xunits.priority = 0
    template_test.yunits.priority = 0
    template_ref.xunits.priority = 0
    template_ref.yunits.priority = 0
    template_diff.xunits.priority = 0
    template_diff.yunits.priority = 0

    utils.set_units(test, parameter.test_units)
    utils.set_units(reference, parameter.reference_units)
    utils.set_units(diff, parameter.diff_units)

    test.long_name = parameter.test_title
    reference.long_name = parameter.reference_title
    diff.long_name = parameter.diff_title

    test.id = parameter.test_name_yrs
    reference.id = parameter.ref_name_yrs
    diff.id = parameter.diff_name

    # model and observation graph
    utils.plot_min_max_mean(
        textcombined_objs, vcs_canvas, metrics_dict, 'test')
    utils.plot_min_max_mean(textcombined_objs, vcs_canvas, metrics_dict, 'ref')
    utils.plot_min_max_mean(
        textcombined_objs, vcs_canvas, metrics_dict, 'diff')

    reference_isofill = vcs.getisofill('reference_isofill')
    reference_isofill.missing = 'grey'
    reference = log_yaxis(reference, reference_isofill)

    test_isofill = vcs.getisofill('test_isofill')
    test_isofill.missing = 'grey'
    test = log_yaxis(test, test_isofill)

    diff_isofill = vcs.getisofill('diff_isofill')
    diff_isofill.missing = 'grey'
    diff = log_yaxis(diff, diff_isofill)

    utils.set_levels_of_graphics_method(
        reference_isofill, parameter.contour_levels, reference, test)
    utils.set_levels_of_graphics_method(
        test_isofill, parameter.contour_levels, test, reference)
    utils.set_levels_of_graphics_method(
        diff_isofill, parameter.diff_levels, diff)

    if parameter.arrows:
        reference_isofill.ext_1 = True
        reference_isofill.ext_2 = True
        test_isofill.ext_1 = True
        test_isofill.ext_2 = True
        diff_isofill.ext_1 = True
        diff_isofill.ext_2 = True

    utils.set_colormap_of_graphics_method(
        vcs_canvas, parameter.reference_colormap, reference_isofill, parameter)
    utils.set_colormap_of_graphics_method(
        vcs_canvas, parameter.test_colormap, test_isofill, parameter)
    utils.set_colormap_of_graphics_method(
        vcs_canvas, parameter.diff_colormap, diff_isofill, parameter)

    vcs_canvas.plot(test, template_test, test_isofill)
    vcs_canvas.plot(reference, template_ref, reference_isofill)
    vcs_canvas.plot(diff, template_diff, diff_isofill)

    utils.plot_rmse_and_corr(textcombined_objs, vcs_canvas, metrics_dict)

    # Plotting the main title
    main_title = utils.managetextcombined(
        textcombined_objs, 'main_title', 'main_title', vcs_canvas)
    main_title.string = parameter.main_title
    # for some reason, this needs to be before a call to vcs_canvas.plot()
    vcs_canvas.portrait()
    vcs_canvas.plot(main_title)

    if not parameter.logo:
        vcs_canvas.drawlogooff()

    fnm = os.path.join(get_output_dir(parameter.current_set,
                                      parameter), parameter.output_file)
    for f in parameter.output_format:
        f = f.lower().split('.')[-1]
        if f == 'png':
            vcs_canvas.png(fnm)
        elif f == 'pdf':
            vcs_canvas.pdf(fnm)
        elif f == 'svg':
            vcs_canvas.svg(fnm)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(get_output_dir(parameter.current_set, parameter,
            ignore_container=True), parameter.output_file)
        print('Plot saved in: ' + fnm + '.' + f)

    vcs_canvas.clear()
