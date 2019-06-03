from __future__ import print_function

import os
import sys
import vcs
import acme_diags
import acme_diags.plot.vcs as utils
from acme_diags.driver.utils.general import get_output_dir

textcombined_objs = {}


def managetextcombined(tt_name, to_name, vcs_canvas):
    """Caches textcombined objects"""
    new_name = "%s:::%s" % (tt_name, to_name)
    mytc = textcombined_objs.get(new_name, None)
    if mytc is None:
        mytc = vcs_canvas.createtextcombined(
            Tt_source=tt_name, To_source=to_name)
        textcombined_objs[new_name] = mytc
    return mytc


def plot(ref, test, diff, metrics_dict, parameters):
    vcs_canvas = vcs.init(bg=True, geometry=(
        parameters.canvas_size_w, parameters.canvas_size_h))

    # Line options, see here: https://uvcdat.llnl.gov/documentation/vcs/vcs-10.html
    # Other options not in the above link:
    # https://uvcdat.llnl.gov/docs/vcs/graphics/unified1D.html
    ref_plot_linetype = 0
    ref_plot_color = 215  # 6 to 239
    ref_plot_width = 3  # 1 to 100
    ref_plot_marker = 1
    ref_plot_markersize = 1
    ref_plot_markercolor = 215

    test_plot_linetype = 0
    test_plot_color = 1
    test_plot_width = 3
    test_plot_marker = 1
    test_plot_markersize = 1
    test_plot_markercolor = 1

    diff_plot_linetype = 0
    diff_plot_color = 1
    diff_plot_width = 3
    diff_plot_marker = 1
    diff_plot_markersize = 1
    diff_plot_markercolor = 1

    file_path = os.path.join(acme_diags.INSTALL_PATH, 'zonal_mean_xy')
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_3.json'))

    utils.set_units(test, parameters.test_units)
    utils.set_units(ref, parameters.reference_units)
    utils.set_units(diff, parameters.diff_units)

    if hasattr(test, 'long_name'):
        test.long_name = parameters.test_title if parameters.test_title is not '' else test.long_name
    if hasattr(ref, 'long_name'):
        ref.long_name = parameters.reference_title if parameters.reference_title is not '' else ref.long_name
    if hasattr(diff, 'long_name'):
        diff.long_name = parameters.diff_title if parameters.diff_title is not '' else diff.long_name

    test.id = str(
        parameters.test_name_yrs) if parameters.test_name_yrs is not '' else test.id
    ref.id = str(
        parameters.ref_name_yrs) if parameters.ref_name_yrs is not '' else ref.id
    diff.id = str(
        parameters.diff_name) if parameters.diff_name is not '' else diff.id

    # use vcs_canvas.show('colormap') to view all colormaps
    # 6 to 239 are purple to red in rainbow order
    vcs_canvas.setcolormap('rainbow')

    ref_test_template = vcs.gettemplate('ref_test_template')

    # Turn off the units of the axes in the plots.
    ref_test_template.xunits.priority = 0
    ref_test_template.yunits.priority = 0

    ref_test_yaxis_title = managetextcombined(
        'ref_test_yaxis_title', 'ref_test_yaxis_title', vcs_canvas)
    ref_test_yaxis_title.angle = 270
    ref_test_yaxis_title.halign = 'center'
    ref_test_yaxis_title.y = (
        ref_test_template.data.y1 + ref_test_template.data.y2) / 2
    ref_test_yaxis_title.x = ref_test_template.data.x1 - 0.08
    ref_test_yaxis_title.string = test.long_name + ' (' + test.units + ')'
    vcs_canvas.plot(ref_test_yaxis_title)

    ref_test_template.legend.priority = 0
    ref_test_template.title.priority = 0

    # the actual box around the plot
    ref_test_template.box1.x1 = 0.1223
    ref_test_template.box1.x2 = 0.96
    ref_test_template.box1.y1 = 0.55
    ref_test_template.box1.y2 = 0.90

    # data (the lines) need to be offset accordingly
    ref_test_template.data.x1 = 0.1223
    ref_test_template.data.x2 = 0.96
    ref_test_template.data.y1 = 0.55
    ref_test_template.data.y2 = 0.90

    ref_test_template.units.textorientation = 'defright'
    ref_test_template.units.x = 0.96
    ref_test_template.units.y = 0.91

    # labels on xaxis
    ref_test_template.xlabel1.y = (0.55) - 0.02  # no xlabel1.x attribute

    # actual ticks on xaxis
    ref_test_template.xtic1.y1 = (0.55 - 0.005) + 0.01
    ref_test_template.xtic1.y2 = (0.55 - 0.005)

    # name of xaxis
    # ref_test_template.xname.y += 0.29

    # labels on yaxis
    ref_test_template.ylabel1.x = 0.11  # no ylabel1.y attribute

    # actual ticks on yaxis
    ref_test_template.ytic1.x1 = (0.1223 - 0.006) + 0.01
    ref_test_template.ytic1.x2 = (0.1223 - 0.006)

    diff_template = vcs.gettemplate('diff_template')

    # Turn off the units of the axes in the plots.
    diff_template.xunits.priority = 0
    diff_template.yunits.priority = 0

    diff_yaxis_title = managetextcombined(
        'diff_yaxis_title', 'diff_yaxis_title', vcs_canvas)
    diff_yaxis_title.angle = 270
    diff_yaxis_title.halign = 'center'
    diff_yaxis_title.y = (diff_template.data.y1 + diff_template.data.y2) / 2
    diff_yaxis_title.x = diff_template.data.x1 - 0.08
    diff_yaxis_title.string = test.long_name + ' (' + test.units + ')'
    vcs_canvas.plot(diff_yaxis_title)

    diff_template.units.textorientation = 'defright'
    diff_template.units.x += 0.01
    diff_template.legend.priority = 0

    diff_template.ytic1.x1 = (0.1223 - 0.006) + 0.01
    diff_template.ytic1.x2 = (0.1223 - 0.006)
    diff_template.ylabel1.x = 0.11  # no ylabel1.y attribute
    diff_template.units.textorientation = 'defright'
    diff_template.units.x = 0.96

    ref_line = vcs_canvas.getxvsy('ref_plot')
    ref_line.datawc_y1 = min(ref.min(), test.min())
    ref_line.datawc_y2 = max(ref.max(), test.max())
    ref_line.datawc_x1 = -90
    ref_line.datawc_x2 = 90
    ref_line.xticlabels1 = {-90: "90S", -60: "60S",
                            -30: "30S", 0: "Eq", 30: "30N",
                            60: "60N", 90: "90N"}

    test_line = vcs_canvas.getxvsy('test_plot')
    test_line.datawc_y1 = min(ref.min(), test.min())
    test_line.datawc_y2 = max(ref.max(), test.max())
    test_line.datawc_x1 = -90
    test_line.datawc_x2 = 90

    diff_line = vcs_canvas.getxvsy('diff_plot')
    diff_line.datawc_y1 = diff.min()
    diff_line.datawc_y2 = diff.max()
    diff_line.datawc_x1 = -90
    diff_line.datawc_x2 = 90
    diff_line.xticlabels1 = {-90: "90S", -60: "60S",
                             -30: "30S", 0: "Eq", 30: "30N",
                             60: "60N", 90: "90N"}

    ref_line.linetype = ref_plot_linetype
    ref_line.linecolor = ref_plot_color
    ref_line.linewidth = ref_plot_width
    ref_line.marker = ref_plot_marker
    ref_line.markersize = ref_plot_markersize
    ref_line.markercolor = ref_plot_markercolor

    test_line.linetype = test_plot_linetype
    test_line.linecolor = test_plot_color
    test_line.linewidth = test_plot_width
    test_line.marker = test_plot_marker
    test_line.markersize = test_plot_markersize
    test_line.markercolor = test_plot_markercolor

    diff_line.linetype = diff_plot_linetype
    diff_line.linecolor = diff_plot_color
    diff_line.linewidth = diff_plot_width
    diff_line.marker = diff_plot_marker
    diff_line.markersize = diff_plot_markersize
    diff_line.markercolor = diff_plot_markercolor

    blank_template = vcs_canvas.gettemplate('blank_template')
    blank_template.legend.priority = 0
    blank_template.data.priority = 1

    vcs_canvas.plot(ref, ref_line, ref_test_template)
    vcs_canvas.plot(test, test_line, blank_template)
    vcs_canvas.plot(diff, diff_line, diff_template)

    # Plot the main title
    main_title = managetextcombined('main_title', 'main_title', vcs_canvas)
    main_title.string = parameters.main_title
    # for some reason, this needs to be before a call to vcs_canvas.plot()
    vcs_canvas.portrait()
    vcs_canvas.plot(main_title)

    test_title = managetextcombined('test_title', 'test_title', vcs_canvas)
    test_title.string = "Test: " + str(parameters.test_name)
    test_title.color = 1
    test_title.x = ref_test_template.data.x1 - 0.05
    test_title.y = ref_test_template.data.y2 + 0.045
    test_title.height = 12
    vcs_canvas.plot(test_title)

    ref_title = managetextcombined('ref_title', 'ref_title', vcs_canvas)
    ref_title.string = "Reference: " + str(parameters.reference_name)
    ref_title.color = 215
    ref_title.x = ref_test_template.data.x1 - 0.05
    ref_title.y = ref_test_template.data.y2 + 0.025
    ref_title.height = 12
    vcs_canvas.plot(ref_title)

    if not parameters.logo:
        vcs_canvas.drawlogooff()

    fnm = os.path.join(get_output_dir(parameters.current_set,
                                      parameters), parameters.output_file)
    for f in parameters.output_format:
        f = f.lower().split('.')[-1]
        if f == 'png':
            vcs_canvas.png(fnm)
        elif f == 'pdf':
            vcs_canvas.pdf(fnm)
        elif f == 'svg':
            vcs_canvas.svg(fnm)
        # Get the filename that the user has passed in and display that.
        # When running in a container, the paths are modified.
        fnm = os.path.join(get_output_dir(parameters.current_set, parameters,
            ignore_container=True), parameters.output_file)
        print('Plot saved in: ' + fnm + '.' + f)
    vcs_canvas.clear()
