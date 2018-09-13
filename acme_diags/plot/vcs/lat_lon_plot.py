from __future__ import print_function

import os
import sys
import vcs
import acme_diags
import acme_diags.plot.vcs as utils
from acme_diags.driver.utils import get_output_dir

textcombined_objs = {}


def plot(reference, test, diff, metrics_dict, parameter):
    vcs_canvas = vcs.init(bg=True, geometry=(
        parameter.canvas_size_w, parameter.canvas_size_h))
    parameter.case_id

    file_path = os.path.join(acme_diags.INSTALL_PATH, 'lat_lon')
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_5.json'))
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_5_new.json'))

    template_test = vcs_canvas.gettemplate('plotset5_0_x_0')
    template_ref = vcs_canvas.gettemplate('plotset5_0_x_1')
    template_diff = vcs_canvas.gettemplate('plotset5_0_x_2')

    utils.set_units(test, parameter.test_units)
    utils.set_units(reference, parameter.reference_units)
    utils.set_units(diff, parameter.diff_units)

    test.long_name = parameter.test_title
    reference.long_name = parameter.reference_title
    diff.long_name = parameter.diff_title

    test.id = parameter.test_name_yrs
    reference.id = parameter.reference_name
    diff.id = parameter.diff_name

    # model and observation graph
    utils.plot_min_max_mean(
        textcombined_objs, vcs_canvas, metrics_dict, 'test')
    utils.plot_min_max_mean(textcombined_objs, vcs_canvas, metrics_dict, 'ref')
    utils.plot_min_max_mean(
        textcombined_objs, vcs_canvas, metrics_dict, 'diff')

    reference_isofill = vcs.getisofill('reference_isofill')
    reference_isofill.missing = 'grey'
    reference_isofill.yticlabels1 = {-90:"90$^\circ$S", -60:"60$^\circ$S", -30:"30$^\circ$S", 0:"0$^\circ$", 30:"30$^\circ$N", 60:"60$^\circ$N", 90:"90$^\circ$N"}
    test_isofill = vcs.getisofill('test_isofill')
    test_isofill.missing = 'grey'
    test_isofill.yticlabels1 = {-90:"90$^\circ$S", -60:"60$^\circ$S", -30:"30$^\circ$S", 0:"0$^\circ$", 30:"30$^\circ$N", 60:"60$^\circ$N", 90:"90$^\circ$N"}
    diff_isofill = vcs.getisofill('diff_isofill')
    diff_isofill.missing = 'grey'
    diff_isofill.yticlabels1 = {-90:"90$^\circ$S", -60:"60$^\circ$S", -30:"30$^\circ$S", 0:"0$^\circ$", 30:"30$^\circ$N", 60:"60$^\circ$N", 90:"90$^\circ$N"}

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
    # for some reason, this needs to be before a call to vcs_canvas.plot()
    vcs_canvas.portrait()

    # Plotting the main title
    main_title = utils.managetextcombined(
        textcombined_objs, 'main_title', 'main_title', vcs_canvas)
    main_title.string = parameter.main_title
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
        print('Plot saved in: ' + fnm + '.' + f)
    vcs_canvas.clear()
