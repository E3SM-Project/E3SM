#!/usr/bin/env python
import os
import cdms2
import cdutil
import vcs
import genutil.statistics
import numpy
import acme_diags.acme_parser

def compute_rmse(model, obs):
    rmse = -numpy.infty
    try:
        rmse = float(genutil.statistics.rms(model, obs, axis='xy', weights='generate'))
    except Exception, err:
        print err
    return rmse

def compute_corr(model, obs):
    corr = -numpy.infty
    try:
        corr = float(genutil.statistics.correlation(model, obs, axis='xy', weights='generate'))
    except Exception, err:
        print err
    return corr

def plot_min_max_mean(canvas, variable, ref_test_or_diff):
    """canvas is a vcs.Canvas, variable is a
    cdms2.tvariable.TransientVariable and ref_test_or_diff is a string"""
    var_min = '%.2f' % float(variable.min())
    var_max = '%.2f' % float(variable.max())

    # Doesn't work: var_mean = str(round(float(variable.mean()), 2))
    var_mean = '%.2f' % cdutil.averager(variable, axis='xy', weights='generate')

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

def plot_rmse_and_corr(canvas, model, obs):
    """canvas is a vcs.Canvas, model and obs are
    a cdms2.tvariable.TransientVariable"""

    rmse = '%.2f' % compute_rmse(obs, model)
    corr = '%.2f' % compute_corr(obs, model)

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

    rmse_value.string = rmse
    corr_value.string = corr

    canvas.plot(rmse_label)
    canvas.plot(corr_label)
    canvas.plot(rmse_value)
    canvas.plot(corr_value)

def set_colormap_of_graphics_method(canvas, parameter_colormap, method):
    if parameter_colormap is not '':
        method.colormap = canvas.getcolormap(parameter_colormap)
        _fix_levels_on_new_colormap(method)

def _fix_levels_on_new_colormap(method):
    colors = vcs.getcolors(method.levels, colors=range(6, 240))
    method.fillareacolors = colors

def set_levels_of_graphics_method(method, levels):
    if levels != []:
        method.levels = levels


def set_units(ref_or_test, units):
    if units != '':
        ref_or_test.units = units

parser = acme_diags.acme_parser.ACMEParser()
parameter = parser.get_parameter()

reference_data_path = parameter.reference_data_path
reference_data_set = parameter.reference_data_set # observation
test_data_path = parameter.test_data_path
test_data_set = parameter.test_data_set # model

case_id = parameter.case_id
if not os.path.exists(case_id):
    os.makedirs(case_id)
var = parameter.variables
season = parameter.season

# Below should be read from metadata
mod_name = '1850_alpha6_01 (yrs0070-0099)'
obs_name = 'GPCP (yrs1979-2009)'

f_obs = cdms2.open(reference_data_path + reference_data_set)
f_mod = cdms2.open(test_data_path + test_data_set)

obs_pr = f_obs(var, longitude=(-180, 540))
mod_pr = (f_mod('PRECC', longitude=(-180, 540)) + f_mod('PRECL', longitude=(-180, 540)))*3600.0*24.0*1000.0

set_units(obs_pr, parameter.reference_units)
set_units(mod_pr, parameter.test_units)

# For plotting, original grid is plotted for model observation, differece plot is regridded to coaser grid. Need if statement to evaluate grid size. aminusb_2ax from uvcmetrics takes care of this,which also considers complex corner cases.
axes1 = mod_pr.getAxisList()
axes2 = obs_pr.getAxisList()
if len(axes1[1]) <= len(axes2[1]): # use nlat to decide data resolution, higher number means higher data resolution. For the difference plot, regrid toward lower resolution
    model_grid = mod_pr.getGrid()
    mod_pr_reg = mod_pr
    obs_pr_reg = obs_pr.regrid(model_grid, regridTool=parameter.regrid_tool, regridMethod=parameter.regrid_method)
else:
    obs_grid = obs_pr.getGrid()
    obs_pr_reg = obs_pr
    mod_pr_reg = mod_pr.regrid(obs_grid, regridTool=parameter.regrid_tool, regridMethod=parameter.regrid_method)
dif_pr = mod_pr_reg - obs_pr_reg

# Plotting
vcs_canvas = vcs.init(bg=True, geometry=(parameter.canvas_size_w, parameter.canvas_size_h))
if not parameter.logo:
    vcs_canvas.drawlogooff()

vcs_canvas.scriptrun('plot_set_5.json')
vcs_canvas.scriptrun('plot_set_5_new.json')
template_test = vcs_canvas.gettemplate('plotset5_0_x_0')
template_ref = vcs_canvas.gettemplate('plotset5_0_x_1')
template_diff = vcs_canvas.gettemplate('plotset5_0_x_2')

mod_pr.long_name = parameter.test_title
obs_pr.long_name = parameter.reference_title
dif_pr.long_name = parameter.diff_title

mod_pr.id = parameter.test_name
obs_pr.id = parameter.reference_name
dif_pr.id = parameter.diff_name

# model and observation graph
plot_min_max_mean(vcs_canvas, mod_pr, 'test')
plot_min_max_mean(vcs_canvas, obs_pr, 'reference')
plot_min_max_mean(vcs_canvas, dif_pr, 'diff')

reference_isofill = vcs.getisofill('reference_isofill')
test_isofill = vcs.getisofill('test_isofill')
diff_isofill = vcs.getisofill('diff_isofill')

set_levels_of_graphics_method(reference_isofill, parameter.reference_levels)
set_levels_of_graphics_method(test_isofill, parameter.test_levels)
set_levels_of_graphics_method(diff_isofill, parameter.diff_levels)

set_colormap_of_graphics_method(vcs_canvas, parameter.reference_colormap, reference_isofill)
set_colormap_of_graphics_method(vcs_canvas, parameter.test_colormap, test_isofill)
set_colormap_of_graphics_method(vcs_canvas, parameter.diff_colormap, diff_isofill)



vcs_canvas.plot(mod_pr, template_test, reference_isofill)
vcs_canvas.plot(obs_pr, template_ref, test_isofill)
vcs_canvas.plot(dif_pr, template_diff, diff_isofill)

plot_rmse_and_corr(vcs_canvas, mod_pr_reg, obs_pr_reg)

# Plotting the main title
main_title = vcs_canvas.createtextcombined(Tt_source = 'main_title',
                                           To_source = 'main_title')
main_title.string = parameter.main_title
vcs_canvas.plot(main_title)

vcs_canvas.png(case_id + '/' + parameter.output_file)
