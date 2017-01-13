import cdms2
import cdutil
import vcs
import MV2
import EzTemplate
from metrics.computation.reductions import aminusb
import genutil.statistics
import numpy
import cdutil

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
    #print 'type(variable.mean)'
    #print type(variable.mean)
    var_min = str(round(float(variable.min()), 2))
    var_max = str(round(float(variable.max()), 2))

    #var_mean = str(round(float(variable.mean()), 2))
    var_mean = str(round(cdutil.averager(variable, axis='xy', weights='generate'), 2))

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

def plot_rmse_and_corr(canvas, template, model, obs):
    """canvas is a vcs.Canvas, template is a vcs.Template,
    model and obs are a cdms2.tvariable.TransientVariable"""

    #rmse = str(round(compute_rmse(obs, model), 2))
    rmse = str(compute_rmse(obs, model))
    #corr = str(round(compute_corr(obs, model), 2))
    corr = str(compute_corr(obs, model))

    text_orientation = canvas.gettextorientation(template.mean.textorientation)
    height = text_orientation.height

    #plot_text(canvas, "RMSE", template.comment1.x, template.comment1.y, height, "left")
    #plot_text(canvas, "CORR", template.comment2.x, template.comment2.y, height, "left")

    plot_text(canvas, rmse, template.comment1.x+0.12, template.comment1.y, height, "right")
    plot_text(canvas, corr, template.comment2.x+0.12, template.comment2.y, height, "right")


def plot_text(canvas, label_string, x, y, height, align):
    label = vcs.createtextcombined()
    label.x = x
    label.y = y
    label.string = label_string
    label.height = height
    label.halign = align
    canvas.plot(label)

#reference_data_path='/space1/test_data/obs_for_diagnostics/'  # observation
#test_data_path='/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'  # model
reference_data_path='../'
test_data_path='../'

#Pre defined variables
var='PRECT'
season='ANN'
#Below should be read from metadata
mod_name='1850_alpha6_01 (yrs0070-0099)'
obs_name='GPCP (yrs1979-2009)'

#Read in data
#reference_data_set='GPCP_v2.2_ANN_climo.nc'  # observation
reference_data_set='GPCP_ANN_climo.nc'  # observation
test_data_set='20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model

f_obs=cdms2.open(reference_data_path + reference_data_set)
f_mod=cdms2.open(test_data_path + test_data_set)

obs_pr=f_obs('PRECT',longitude=(-180, 540))
#obs_pr=f_obs('PRECT')
mod_pr=(f_mod('PRECC',longitude=(-180, 540))+f_mod('PRECL', longitude=(-180, 540)))*3600.0*24.0*1000.0
#mod_pr=(f_mod('PRECC')+f_mod('PRECL'))*3600.0*24.0*1000.0
mod_pr.units='mm/day'

#For plotting, original grid is plotted for model observation, differece plot is regridded to coaser grid. Need if statement to evaluate grid size. aminusb_2ax from uvcmetrics takes care of this,which also considers complex corner cases.
axes1=mod_pr.getAxisList()
axes2=obs_pr.getAxisList()
if len(axes1[1])<=len(axes2[1]): #use nlat to decide data resolution, higher number means higher data resolution. For the difference plot, regrid toward lower resolution
    model_grid=mod_pr.getGrid()
    mod_pr_reg=mod_pr
    obs_pr_reg=obs_pr.regrid(model_grid,regridTool='esmf',regridMethod='linear')
else:
    obs_grid=obs_pr.getGrid()
    obs_pr_reg=obs_pr
    mod_pr_reg=mod_pr.regrid(obs_grid,regridTool='esmf',regridMethod='linear')
dif_pr=mod_pr_reg-obs_pr_reg

#Plotting
x = vcs.init(bg=True, geometry=(1212,1628))
x.portrait()

x.scriptrun('plot_set_5.json')
template_0 = x.gettemplate('plotset5_0_x_0')
template_1 = x.gettemplate('plotset5_0_x_1')
template_2 = x.gettemplate('plotset5_0_x_2')

#It turns out the long_name attribute of the mv appears as title in .json.
mod_pr.long_name='model'
obs_pr.long_name='observation'
dif_pr.long_name='model-observation'

template_0.title.priority=1
template_1.title.priority=1
template_2.title.priority=1

template_0.units.priority=1
template_1.units.priority=1
template_2.units.priority=1

mod_pr.id=mod_name
obs_pr.id=obs_name
template_0.dataname.priority=1
template_1.dataname.priority=1

# model and observation graph
isofill = x.createisofill()
isofill.datawc_x1=0
isofill.datawc_x2=360
isofill.datawc_y1=-90
isofill.datawc_y2=90
isofill.levels=[0, 0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]
# NOTE: because of the obs and model data files used,
# there is no 360 degree value, so we use 358 as 0.
# Same for 0 where we use 2 instead.
isofill.xticlabels1 = {0: '0', 30: '30E', 60: '60E', 90: '90E',
                       120: '120E', 150: '150E', 180: '180W', 210: '150W',
                       240: '120W', 270: '90W', 300: '60W', 330: '30W', 360: '0'}
isofill.yticlabels1 = {-90: '90S', -80: '80S', -60: '60S', -40: '40S',
                       -20:'20S', 0: 'Eq', 20: '20N', 40: '40N', 60: '60N',
                       80: '80N', 90: '90N'}

#ext_1 and ext_2 are arrows
isofill.ext_1 = True
isofill.ext_2 = True

plot_min_max_mean(x, template_0, mod_pr)
plot_min_max_mean(x, template_1, obs_pr)
x.plot(mod_pr, template_0, isofill)
x.plot(obs_pr, template_1, isofill)

#Create main title for the 3 plots
plot_text(x, ' '.join([var, season]), 0.42, 0.98, 20, "left")
#For plotting, original grid is plotted for model observation, differece plot is regridded to coaser grid. Need if statement to evaluate grid size. aminusb_2ax from uvcmetrics takes care of this,which also considers complex corner cases.
axes1=mod_pr.getAxisList()
axes2=obs_pr.getAxisList()
if len(axes1[1])<=len(axes2[1]): #use nlat to decide data resolution, higher number means higher data resolution. For the difference plot, regrid toward lower resolution
    model_grid=mod_pr.getGrid()
    mod_pr_reg=mod_pr
    obs_pr_reg=obs_pr.regrid(model_grid,regridTool='esmf',regridMethod='linear')
else:
    obs_grid=obs_pr.getGrid()
    obs_pr_reg=obs_pr
    mod_pr_reg=mod_pr.regrid(obs_grid,regridTool='esmf',regridMethod='linear')
dif_pr=mod_pr_reg-obs_pr_reg

# difference graph
isofill = x.createisofill()
isofill.datawc_x1=0
isofill.datawc_x2=360
isofill.datawc_y1=-90
isofill.datawc_y2=90

isofill.levels=[-6, -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6]
# After you set arrows, need to enable arrows again
###isofill.ext_1 = True
###isofill.ext_2 = True
isofill.colormap = x.getcolormap('bl_to_darkred')
isofill.ext_1 = True
isofill.ext_2 = True

plot_min_max_mean(x, template_2, dif_pr)
x.plot(dif_pr, template_2, isofill)

plot_rmse_and_corr(x, template_2, mod_pr_reg, obs_pr_reg)
x.png('test.png')
