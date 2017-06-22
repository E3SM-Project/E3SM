import os
import sys
import vcs
from acme_diags.driver.utils import get_output_dir

vcs_canvas = vcs.init()
textcombined_objs = {}

def set_units(ref_or_test, units):
    if units != '':
        ref_or_test.units = units

def managetextcombined(tt_name, to_name):
    """Caches textcombined objects"""
    new_name = "%s:::%s" % (tt_name, to_name)
    mytc = textcombined_objs.get(new_name, None)
    if mytc is None:
        mytc = vcs_canvas.createtextcombined(Tt_source=tt_name, To_source=to_name)
        textcombined_objs[new_name] = mytc
    return mytc

def plot(ref, test, diff, metrics_dict, parameters):
    # Line options, see here: https://uvcdat.llnl.gov/documentation/vcs/vcs-10.html
    # Other options not in the above link: https://uvcdat.llnl.gov/docs/vcs/graphics/unified1D.html
    ref_plot_linetype = 0
    ref_plot_color = 1  # 6 to 239
    ref_plot_width = 3  # 1 to 100
    ref_plot_marker = 1
    ref_plot_markersize = 1
    ref_plot_markercolor = 1

    test_plot_linetype = 0
    test_plot_color = 215
    test_plot_width = 3
    test_plot_marker = 1
    test_plot_markersize = 1
    test_plot_markercolor = 215

    diff_plot_linetype = 0
    diff_plot_color = 1
    diff_plot_width = 3
    diff_plot_marker = 1
    diff_plot_markersize = 1
    diff_plot_markercolor = 1

    file_path = os.path.join(sys.prefix, 'share', 'acme_diags', 'set3')
    vcs_canvas.scriptrun(os.path.join(file_path, 'plot_set_3.json'))

    set_units(test, parameters.test_units)
    set_units(ref, parameters.reference_units)
    set_units(diff, parameters.diff_units)

    test.long_name = parameters.test_title
    ref.long_name = parameters.reference_title
    diff.long_name = parameters.diff_title

    test.id = str(parameters.test_name)
    ref.id = str(parameters.reference_name)
    diff.id = str(parameters.diff_name)

    # use vcs_canvas.show('colormap') to view all colormaps
    vcs_canvas.setcolormap('rainbow')  # 6 to 239 are purple to red in rainbow order

    ref_test_template = vcs.gettemplate('ref_test_template')
    # make all of the elements listed have priority = 0
    ref_test_template.blank(["mean", "max", "min", "zvalue", "dataname", "crtime", "ytic2", "xtic2"])
    
    # the actual box around the plot
    ref_test_template.box1.x1 = 0.123
    ref_test_template.box1.x2 = 0.86
    ref_test_template.box1.y1 = 0.55
    ref_test_template.box1.y2 = 0.90

    # data (the lines) need to be offset accordingly
    ref_test_template.data.x1 = 0.123
    ref_test_template.data.x2 = 0.86
    ref_test_template.data.y1 = 0.55
    ref_test_template.data.y2 = 0.90

    ref_test_template.legend.x1 = 0.88
    ref_test_template.legend.x2 = 0.98
    ref_test_template.legend.y1 = 0.86
    ref_test_template.legend.y2 = 0.88
    ref_test_template.legend.textorientation = 'defright'

    ref_test_template.title.x = 0.5
    ref_test_template.title.textorientation = 'defcenter'

    ref_test_template.units.textorientation = 'defright'
    ref_test_template.units.x = 0.86
    ref_test_template.units.y = 0.91

    # labels on xaxis
    ref_test_template.xlabel1.y = (0.55) - 0.02  # no xlabel1.x attribute

    # actual ticks on xaxis
    ref_test_template.xtic1.y1 = (0.55 - 0.005) + 0.01
    ref_test_template.xtic1.y2 = (0.55 - 0.005)

    # name of xaxis
    # ref_test_template.xname.y += 0.29

    # labels on yaxis
    ref_test_template.ylabel1.x = 0.1108  # no ylabel1.y attribute

    # actual ticks on yaxis
    ref_test_template.ytic1.x1 = (0.123 - 0.005) + 0.01
    ref_test_template.ytic1.x2 = (0.123 - 0.005)

    # name of yaxis
    ref_test_template.yname.priority = 1
    ref_test_template.yname.x = ref_test_template.xname.x
    ref_test_template.yname.y = ref_test_template.xname.y
    
    #ref_test_template.xname.y += 0.20
    # ref_test_template.yname.x += 0.05
    # ref_test_template.yname.y += 0.17


    diff_template = vcs.gettemplate('diff_template')
    diff_template.units.textorientation = 'defright'
    diff_template.units.x += 0.01
    diff_template.legend.priority = 0

    '''
    diff_template.box1.y1 -= 0.47
    diff_template.box1.y2 -= 0.47

    diff_template.data.y1 -= 0.47
    diff_template.data.y2 -= 0.47

    diff_template.legend.y1 -= 0.47
    diff_template.legend.y2 -= 0.47

    diff_template.title.y -= 0.47
    diff_template.units.y -= 0.47

    diff_template.xlabel1.y -= 0.47

    diff_template.xtic1.y1 -= 0.47
    diff_template.xtic1.y2 -= 0.47

    diff_template.xname.y -= 0.47
    diff_template.yname.y -= 0.47
    '''

    ref_line = vcs_canvas.getxvsy('ref_plot')
    ref_line.datawc_y1 = min(ref.min(), test.min())
    ref_line.datawc_y2 = max(ref.max(), test.max())
    ref_line.datawc_x1 = -90
    ref_line.datawc_x2 = 90
    ref_line.xticlabels1 = {-90: "90S", -80: "80S", -60: "60S",
                        -40: "40S", -20: "20S", 0: "Eq",
                        20: "20N", 40: "40N", 60: "60N",
                        80: "80N", 90: "90N"}


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
    diff_line.xticlabels1 = {-90: "90S", -80: "80S", -60: "60S",
                        -40: "40S", -20: "20S", 0: "Eq",
                        20: "20N", 40: "40N", 60: "60N",
                        80: "80N", 90: "90N"}



    #ref_line.line = ref_plot_linetype
    ref_line.linetype = ref_plot_linetype
    ref_line.linecolor = ref_plot_color
    ref_line.linewidth = ref_plot_width
    ref_line.marker = ref_plot_marker
    ref_line.markersize = ref_plot_markersize
    ref_line.markercolor = ref_plot_markercolor

    #test_line.line = test_plot_linetype
    test_line.linetype = test_plot_linetype
    test_line.linecolor = test_plot_color
    test_line.linewidth = test_plot_width
    test_line.marker = test_plot_marker
    test_line.markersize = test_plot_markersize
    test_line.markercolor = test_plot_markercolor
    # test_line.smooth = 1

    #diff_line.line = diff_plot_linetype
    diff_line.linetype = diff_plot_linetype
    diff_line.linecolor = diff_plot_color
    diff_line.linewidth = diff_plot_width
    diff_line.marker = diff_plot_marker
    diff_line.markersize = diff_plot_markersize
    diff_line.markercolor = diff_plot_markercolor

    blank_template = vcs_canvas.gettemplate('blank_template')
    #blank_template.blank()
    #blank_template.legend.priority = 1
    #blank_template.legend.y1 -= 0.05
    #blank_template.legend.y2 -= 0.05

    vcs_canvas.plot(ref, ref_line, ref_test_template)
    vcs_canvas.plot(test, test_line, blank_template)
    vcs_canvas.plot(diff, diff_line, diff_template)

    # Plot the main title
    main_title = managetextcombined('main_title', 'main_title')
    main_title.string = parameters.main_title
    #main_title.height = 18
    #main_title.halign = 'center'
    #main_title.x = 0.5
    #main_title.y = 0.97
    vcs_canvas.portrait()  # for some reason, this needs to be before a call to vcs_canvas.plot()
    vcs_canvas.plot(main_title)

    #ref_test_template.script('plot_set_3.json')
    #blank_template.script('plot_set_3.json')
    #diff_template.script('plot_set_3.json')
    #ref_line.script('plot_set_3.json')
    #test_line.script('plot_set_3.json')
    #diff_line.script('plot_set_3.json')
    #main_title.script('plot_set_3.json')

    vcs_canvas.bgX = parameters.canvas_size_w
    vcs_canvas.bgY = parameters.canvas_size_h
    if not parameters.logo:
        vcs_canvas.drawlogooff()

    fnm = os.path.join(get_output_dir('3', parameters), parameters.output_file)
    for f in parameters.output_format:
        f = f.lower().split('.')[-1]
        if f == 'png':
            vcs_canvas.png(fnm)
        elif f == 'pdf':
            vcs_canvas.pdf(fnm)
        elif f == 'svg':
            vcs_canvas.svg(fnm)

        print('Plot saved in: ' + fnm + '.' + f)
    vcs_canvas.clear()
