import cdms2,cdutil
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
        weights = cdutil.area_weights(model)
        rmse = float(genutil.statistics.rms(model, obs, axis='xy', weights=weights))
    except Exception, err:
        print err
    return rmse

def compute_corr(model, obs):
    corr = -numpy.infty
    try:
        weights = cdutil.area_weights(model)
        corr = float(genutil.statistics.correlation(model, obs, axis='xy', weights=weights))
    except Exception, err:
        print err
    return corr

def plot_min_max_mean(canvas, template, variable):
    """canvas is a vcs.Canvas, template is a vcs.Template,
    variable is a cdms2.tvariable.TransientVariable"""

    #var_min = str(round(float(variable.min()), 2))
    #var_max = str(round(float(variable.max()), 2))
    #var_mean = str(round(float(variable.mean), 2))
    #The reason why I think '%.2f' is better than round is because round somehow ignore zeros at the end of the number.
    var_min =  '%.2f' %variable.min
    var_max =  '%.2f' %variable.max
    var_mean = '%.2f' %variable.mean


    
    print variable.id,var_min,var_max,var_mean

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
    #corr = str(round(compute_corr(obs, model), 2))
    rmse ='%.2f' %compute_rmse(obs, model)
    corr ='%.2f' %compute_corr(obs, model)

    text_orientation = canvas.gettextorientation(template.mean.textorientation)
    height = text_orientation.height

    plot_text(canvas, "RMSE", template.comment1.x, template.comment1.y, height, "left")
    plot_text(canvas, "CORR", template.comment2.x, template.comment2.y, height, "left")

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
reference_data_path='./'
test_data_path='./'

#Read in data
reference_data_set='figure-set5_Global_ANN_T_plot--obs.nc'  # observation
test_data_set='figure-set5_Global_ANN_T_plot--model.nc'  # model

f_obs=cdms2.open(reference_data_path + reference_data_set)
f_mod=cdms2.open(test_data_path + test_data_set)

#obs_pr=f_obs('rv_T_ANN_NCEP')
obs_pr=f_obs('rv_T_ANN_NCEP', longitude=(-180, 540))
mod_pr=f_mod('rv_T_ANN_cam_output', longitude=(-180, 540))

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
#In uvcmetrics, the substraction of data is handled, 'min,max' are calculated correctly for the difference map, however the mean is not. we want to calculate all these explicitly
dif_pr=mod_pr_reg-obs_pr_reg
#calculate metrics and pass in as mv attribute,failed
obs_pr.mean=round(cdutil.averager(obs_pr, axis='xy', weights='generate'),2) #area-weighting
mod_pr.mean=round(cdutil.averager(mod_pr, axis='xy', weights='generate'),2) #area-weighting
dif_pr.mean=round(cdutil.averager(dif_pr, axis='xy', weights='generate'),2) #area-weighting

obs_pr.max=round(mod_pr.max(),2)
mod_pr.max=round(mod_pr.max(),2)
dif_pr.max=round(dif_pr.max(),2)

obs_pr.min=round(mod_pr.min(),2)
mod_pr.min=round(mod_pr.min(),2)
dif_pr.min=round(dif_pr.min(),2)



#Plotting
x = vcs.init(bg=True, geometry=(1212,1628))
x.portrait()

x.scriptrun('plot_set_5.json')
template_0 = x.gettemplate('plotset5_0_x_0')
template_1 = x.gettemplate('plotset5_0_x_1')
template_2 = x.gettemplate('plotset5_0_x_2')

plot_rmse_and_corr(x, template_2, mod_pr_reg, obs_pr_reg)

#It turns out the long_name attribute of the mv appears as title in .json.
mod_pr.long_name='model'
obs_pr.long_name='observation'
dif_pr.long_name='model-observation'
template_0.title.priority=1
template_1.title.priority=1
template_2.title.priority=1

# model and observation graph
isofill = x.createisofill()
isofill.datawc_x1=0
isofill.datawc_x2=360
isofill.datawc_y1=-90
isofill.datawc_y2=90

isofill.levels=[212, 216, 220, 224, 228, 232, 236, 240, 244, 248, 252, 256, 260, 264]#[x+4 for x in range(208, 264, 4)]
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
#isofill.ext_1 = True
#isofill.ext_2 = True
plot_min_max_mean(x, template_0, mod_pr)
plot_min_max_mean(x, template_1, obs_pr)
x.plot(mod_pr, template_0, isofill)
x.plot(obs_pr, template_1, isofill)

isofill = x.createisofill()
isofill.colormap = x.getcolormap('lightblue2darkblue')
if isofill.colormap in ['AMIP',
            'NCAR',
            'bl_to_darkred',
            'bl_to_drkorang',
            'blends',
            'blue2darkorange',
            'blue2darkred',
            'blue2green',
            'blue2grey',
            'blue2orange',
            'blue2orange2red',
            'blue_to_grey',
            'blue_to_grn',
            'blue_to_orange',
            'blue_to_orgred',
            'brown2blue',
            'brown_to_blue',
            'categorical',
            'classic',
            'green2magenta',
            'grn_to_magenta',
            'inferno',
            'lightblue2darkblue',
            'ltbl_to_drkbl',
            'rainbow',
            'rainbow_no_grn',
            'rainbownogreen',
            'sequential',
            'white2blue',
            'white2green',
            'white2magenta',
            'white2red',
            'white2yellow',
            'white_to_blue',
            'white_to_green',
            'white_to_magenta',
            'white_to_red',
            'white_to_yellow']:
        # Restore old behaviour of default colors for old colormaps
        vcs.utils.defaultColorsRange = range(239,15,-1)
else:
    vcs.utils.defaultColorsRange = range(256)

isofill.datawc_x1=0
isofill.datawc_x2=360
isofill.datawc_y1=-90
isofill.datawc_y2=90

# difference graph
###isofill.levels=[-18, -16, -14, -12, -10, -8, -6, -4, -2, 0, 2]
# After you set arrows, need to enable arrows again
#isofill.ext_1 = True
#isofill.ext_2 = True
#isofill.fillareacolors=vcs.getcolors(range(16,240))
#x.setcolormap('bl_to_darkred')
plot_min_max_mean(x, template_2, dif_pr)
x.plot(dif_pr, template_2, isofill)

x.png('test.png')
