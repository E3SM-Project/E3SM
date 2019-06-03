from __future__ import print_function

import os
import sys
import vcs
import acme_diags
import acme_diags.plot.vcs as utils
from acme_diags.driver.utils.general import get_output_dir

textcombined_objs = {}


def plot(reference, test, diff, metrics_dict, parameter):
    vcs_canvas = vcs.init(bg=True, geometry=(
        parameter.canvas_size_w, parameter.canvas_size_h))
    parameter.case_id

    file_path = os.path.join(acme_diags.INSTALL_PATH, 'polar')
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_7.json'))
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_7_new.json'))

    template_test = vcs_canvas.gettemplate('plotset7_0_x_0')
    template_ref = vcs_canvas.gettemplate('plotset7_0_x_1')
    template_diff = vcs_canvas.gettemplate('plotset7_0_x_2')

    # Turn off the units of the axes in the plots.
    template_test.xunits.priority = 0
    template_test.yunits.priority = 0
    template_ref.xunits.priority = 0
    template_ref.yunits.priority = 0
    template_diff.xunits.priority = 0
    template_diff.yunits.priority = 0

    template_test.title.x = 0.01
    template_test.dataname.x = 0.01
    template_test.dataname.y = template_test.title.y - 0.02
    template_test.data.y1 -= 0.025
    template_test.data.y2 -= 0.02

    template_ref.title.x = 0.01
    template_ref.dataname.x = 0.01
    template_ref.dataname.y = template_ref.title.y - 0.02
    template_ref.data.y1 -= 0.025
    template_ref.data.y2 -= 0.025

    template_diff.title.x = 0.01
    template_diff.dataname.x = 0.01
    template_diff.dataname.y = template_diff.title.y - 0.02
    template_diff.data.y1 -= 0.025
    template_diff.data.y2 -= 0.025

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
    reference_isofill.xticlabels1 = {
        0: "0", 30: "30$^\circ$E", 60: "60$^\circ$E", 90: "90$^\circ$E", 120: "120$^\circ$E",
        150: "150$^\circ$E", 180: "180$^\circ$W", 210: "150$^\circ$W", 240: "120$^\circ$W",
        270: "90$^\circ$W", 300: "60$^\circ$W", 330: "30$^\circ$W"
    }
    
    test_isofill = vcs.getisofill('test_isofill')
    test_isofill.missing = 'grey'
    test_isofill.xticlabels1 = {
        0: "0", 30: "30$^\circ$E", 60: "60$^\circ$E", 90: "90$^\circ$E", 120: "120$^\circ$E",
        150: "150$^\circ$E", 180: "180$^\circ$W", 210: "150$^\circ$W", 240: "120$^\circ$W",
        270: "90$^\circ$W", 300: "60$^\circ$W", 330: "30$^\circ$W"
    }

    diff_isofill = vcs.getisofill('diff_isofill')
    diff_isofill.missing = 'grey'
    diff_isofill.xticlabels1 = {
        0: "0", 30: "30$^\circ$E", 60: "60$^\circ$E", 90: "90$^\circ$E", 120: "120$^\circ$E",
        150: "150$^\circ$E", 180: "180$^\circ$W", 210: "150$^\circ$W", 240: "120$^\circ$W",
        270: "90$^\circ$W", 300: "60$^\circ$W", 330: "30$^\circ$W"
    }

    if parameter.var_region.lower().find('polar') != -1:
        reference_isofill.projection = 'polar'
        test_isofill.projection = 'polar'
        diff_isofill.projection = 'polar'
        if parameter.var_region.find('S') != -1:
            lat_y1 = -90
            lat_y2 = -55
        elif parameter.var_region.find('N') != -1:
            lat_y1 = 90
            lat_y2 = 50

        # this should extracted from selected domain
        reference_isofill.datawc_y1 = lat_y1
        reference_isofill.datawc_y2 = lat_y2
        test_isofill.datawc_y1 = lat_y1  # this should extracted from selected domain
        test_isofill.datawc_y2 = lat_y2
        diff_isofill.datawc_y1 = lat_y1  # this should extracted from selected domain
        diff_isofill.datawc_y2 = lat_y2

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

    vcs_canvas.plot(utils.add_cyclic(test), template_test, test_isofill)
    vcs_canvas.plot(utils.add_cyclic(reference),
                    template_ref, reference_isofill)
    vcs_canvas.plot(utils.add_cyclic(diff), template_diff, diff_isofill)

    utils.plot_rmse_and_corr(textcombined_objs, vcs_canvas, metrics_dict)

    # Plotting the main title
    main_title = utils.managetextcombined(
        textcombined_objs, 'main_title', 'main_title', vcs_canvas)
    main_title.string = parameter.main_title
    main_title.y = [0.985]
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
