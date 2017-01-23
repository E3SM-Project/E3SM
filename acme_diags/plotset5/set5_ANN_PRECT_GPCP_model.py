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

def plot_min_max_mean(canvas, template, variable):
    """canvas is a vcs.Canvas, template is a vcs.Template,
    variable is a cdms2.tvariable.TransientVariable"""
    var_min = '%.2f' % float(variable.min())
    var_max = '%.2f' % float(variable.max())

    # Doesn't work: var_mean = str(round(float(variable.mean()), 2))
    var_mean = '%.2f' % cdutil.averager(variable, axis='xy', weights='generate')

    # Todo: Just have a textorientation object in the script titled 'min_max_mean'
    # which will have the height and everything
    # EX: text_orientation = canvas.gettextorientation("min_max_mean")
    text_orientation = canvas.gettextorientation(template.mean.textorientation)
    height = text_orientation.height

    # Draw the "Max", "Mean", and "Min" labels
    plot_text(canvas, 'Min', template.min.x, template.min.y, height, "left")
    plot_text(canvas, 'Max', template.max.x, template.max.y, height, "left")
    plot_text(canvas, 'Mean', template.mean.x, template.mean.y, height, "left")

    # Draw the actual mean, max, min labels
    plot_text(canvas, var_min, template.min.x+0.12, template.min.y, height, "right")
    plot_text(canvas, var_max, template.max.x+0.12, template.max.y, height, "right")
    plot_text(canvas, var_mean, template.mean.x+0.12, template.mean.y, height, "right")

def plot_rmse_and_corr(canvas, model, obs):
    """canvas is a vcs.Canvas, template is a vcs.Template,
    model and obs are a cdms2.tvariable.TransientVariable"""

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

def plot_text(canvas, label_string, x, y, height, align):
    label = vcs.createtextcombined()
    label.x = x
    label.y = y
    label.string = label_string
    label.height = height
    label.halign = align
    canvas.plot(label)

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
mod_pr.units = 'mm/day'

# For plotting, original grid is plotted for model observation, differece plot is regridded to coaser grid. Need if statement to evaluate grid size. aminusb_2ax from uvcmetrics takes care of this,which also considers complex corner cases.
axes1 = mod_pr.getAxisList()
axes2 = obs_pr.getAxisList()
if len(axes1[1]) <= len(axes2[1]): # use nlat to decide data resolution, higher number means higher data resolution. For the difference plot, regrid toward lower resolution
    model_grid = mod_pr.getGrid()
    mod_pr_reg = mod_pr
    obs_pr_reg = obs_pr.regrid(model_grid, regridTool='esmf', regridMethod='linear')
else:
    obs_grid = obs_pr.getGrid()
    obs_pr_reg = obs_pr
    mod_pr_reg = mod_pr.regrid(obs_grid, regridTool='esmf', regridMethod='linear')
dif_pr = mod_pr_reg - obs_pr_reg

# Plotting
x = vcs.init(bg=True, geometry=(1212,1628))

x.scriptrun('plot_set_5.json')
x.scriptrun('plot_set_5_new.json')
template_0 = x.gettemplate('plotset5_0_x_0')
template_1 = x.gettemplate('plotset5_0_x_1')
template_2 = x.gettemplate('plotset5_0_x_2')

mod_pr.long_name = 'model'
obs_pr.long_name = 'observation'
dif_pr.long_name = 'model-observation'

template_0.title.priority = 1
template_1.title.priority = 1
template_2.title.priority = 1

template_0.units.priority = 1
template_1.units.priority = 1
template_2.units.priority = 1

mod_pr.id = mod_name
obs_pr.id = obs_name
template_0.dataname.priority = 1
template_1.dataname.priority = 1

# model and observation graph
plot_min_max_mean(x, template_0, mod_pr)
plot_min_max_mean(x, template_1, obs_pr)
x.plot(mod_pr, template_0, vcs.getisofill('reference_isofill'))
x.plot(obs_pr, template_1, vcs.getisofill('test_isofill'))

# Create main title for the 3 plots
plot_text(x, ' '.join([var, season]), 0.42, 0.98, 18, "left")

plot_min_max_mean(x, template_2, dif_pr)
#x.plot(dif_pr, template_2, isofill)
x.plot(dif_pr, template_2, vcs.getisofill('diff_plot'))

plot_rmse_and_corr(x, mod_pr_reg, obs_pr_reg)
x.png(case_id + '/' + parameter.output_file)
