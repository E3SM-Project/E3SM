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

<<<<<<< HEAD:acme_diags/uvcmetric_test/plotset5/plotset5_uvcmetrics.py
reference_data_path='./'
test_data_path='./'
=======
#reference_data_path='/space1/test_data/obs_for_diagnostics/'  # observation
#test_data_path='/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'  # model
reference_data_path='../'
test_data_path='../'

#Pre defined variables
var='PRECT'
season='ANN'
>>>>>>> 2f74ec1ef67ef52dcb7447cf568062b5492cf4c2:acme_diags/jill_test/set5_ANN_PRECT_GPCP_model_tmp.py

#Read in data
reference_data_set='figure-set5_Global_ANN_T_plot--obs.nc'  # observation
test_data_set='figure-set5_Global_ANN_T_plot--model.nc'  # model

f_obs=cdms2.open(reference_data_path + reference_data_set)
f_mod=cdms2.open(test_data_path + test_data_set)

<<<<<<< HEAD:acme_diags/uvcmetric_test/plotset5/plotset5_uvcmetrics.py
obs_pr=f_obs('rv_T_ANN_NCEP')
mod_pr=f_mod('rv_T_ANN_cam_output')
=======
obs_pr=f_obs('PRECT')
mod_pr=(f_mod('PRECC')+f_mod('PRECL'))*3600.0*24.0*1000.0 #this is derived variable convert m/s to mm/d
mod_pr.units='mm/d'
>>>>>>> 2f74ec1ef67ef52dcb7447cf568062b5492cf4c2:acme_diags/jill_test/set5_ANN_PRECT_GPCP_model_tmp.py


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
<<<<<<< HEAD:acme_diags/uvcmetric_test/plotset5/plotset5_uvcmetrics.py

#CORR and RMSE need to be calculated after reduction to ensure same array shapes.
rmse= 'RMSE:'+'%.2f' %round(compute_rmse(obs_pr_reg, mod_pr_reg),2)
corr= 'CORR:'+'%.2f' %round(compute_corr(obs_pr_reg, mod_pr_reg),2)
=======

#calculate metrics and pass in as mv attribute,failed
#obs_pr.mean=round(cdutil.averager(obs_pr, axis='xy', weights='generate'),2) #area-weighting
#mod_pr.mean=round(cdutil.averager(mod_pr, axis='xy', weights='generate'),2) #area-weighting
#dif_pr.mean=round(cdutil.averager(dif_pr, axis='xy', weights='generate'),2) #area-weighting
#
#obs_pr.max=round(mod_pr.max(),2)
#mod_pr.max=round(mod_pr.max(),2)
#dif_pr.max=round(dif_pr.max(),2)
#
#obs_pr.min=round(mod_pr.min(),2)
#mod_pr.min=round(mod_pr.min(),2)
#dif_pr.min=round(dif_pr.min(),2)


#CORR and RMSE need to be calculated after reduction to ensure same array shapes.
rmse= 'RMSE:'+'%.2f' %compute_rmse(obs_pr_reg, mod_pr_reg)
corr= 'CORR:'+'%.2f' %compute_corr(obs_pr_reg, mod_pr_reg)
print rmse, corr
>>>>>>> 2f74ec1ef67ef52dcb7447cf568062b5492cf4c2:acme_diags/jill_test/set5_ANN_PRECT_GPCP_model_tmp.py

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

#Or we can change title attribute, adding 3 more lines though...
#mod_pr.title='model'
#obs_pr.title='observation'
#dif_pr.title='model-observation'
#
#del mod_pr.long_name
#del obs_pr.long_name
#del dif_pr.long_name
template_0.title.priority=1
template_1.title.priority=1
template_2.title.priority=1

<<<<<<< HEAD:acme_diags/uvcmetric_test/plotset5/plotset5_uvcmetrics.py
# model and observation graph
=======
template_0.units.priority=1
template_1.units.priority=1
template_2.units.priority=1

#The dataname part should pass in basic model info: model_version...and basic obs info:obsname_version_years
#Those information will be provided by metadata.
template_0.dataname.priority=1
template_1.dataname.priority=1

mod_pr.dataname='test'


>>>>>>> 2f74ec1ef67ef52dcb7447cf568062b5492cf4c2:acme_diags/jill_test/set5_ANN_PRECT_GPCP_model_tmp.py
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
isofill.ext_1 = True
isofill.ext_2 = True

<<<<<<< HEAD:acme_diags/uvcmetric_test/plotset5/plotset5_uvcmetrics.py
x.plot(mod_pr, template_0, isofill)
=======
#x.plot(obs_pr, template_1, isofill)
>>>>>>> 2f74ec1ef67ef52dcb7447cf568062b5492cf4c2:acme_diags/jill_test/set5_ANN_PRECT_GPCP_model_tmp.py
x.plot(obs_pr, template_1, isofill)
isofill = x.createisofill()
isofill.datawc_x1=0
isofill.datawc_x2=360
isofill.datawc_y1=-90
isofill.datawc_y2=90

# difference graph
#isofill.levels=[-1e20, -6, -5, -4, -3, -2, -1, -0.5, 0, 0.5, 1, 2, 3, 4, 5, 6, 1e20]
# After you set arrows, need to enable arrows again
isofill.ext_1 = True
isofill.ext_2 = True
isofill.colormap = x.getcolormap('ltbl_to_drkbl')
x.plot(dif_pr, template_2, isofill, comment1=rmse, comment2=corr)

#Create main title for the 3 plots 
main_title=vcs.createtext()
main_title.string=' '.join([var,season])
main_title.x=0.42
main_title.y=0.98
x.plot(main_title)

x.png('test.png')
