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

reference_data_path='/space1/test_data/obs_for_diagnostics/'  # observation
test_data_path='/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'  # model
#reference_data_path='../'
#test_data_path='../'

#Read in data
reference_data_set='GPCP_v2.2_ANN_climo.nc'  # observation
test_data_set='20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'  # model

f_obs=cdms2.open(reference_data_path + reference_data_set)
f_mod=cdms2.open(test_data_path + test_data_set)

obs_pr=f_obs('PRECT')
mod_pr=(f_mod('PRECC')+f_mod('PRECL'))*3600.0*24.0*1000.0

print obs_pr.shape, mod_pr.shape

mean_obs=cdutil.averager(obs_pr, axis='xy', weights='generate') #area-weighting
mean_mod=cdutil.averager(mod_pr, axis='xy', weights='generate') #area-weighting

max_obs=obs_pr.max()
min_obs=obs_pr.min()

max_mod=mod_pr.max()
min_mod=mod_pr.min()
print min_obs, max_obs

max_obs=obs_pr.max()
min_obs=obs_pr.min()
print min_obs, max_obs

#For plotting, original grid is plotted for model observation, differece plot is regridded to coaser grid. Need if statement to evaluate grid size
#model_grid=mod_pr.getGrid()
#mod_pr_reg=mod_pr
#obs_pr_reg=obs_pr.regrid(model_grid,regridTool='esmf',regridMethod='linear')

obs_grid=obs_pr.getGrid()
obs_pr_reg=obs_pr
mod_pr_reg=mod_pr.regrid(obs_grid,regridTool='esmf',regridMethod='linear')
diff_pr_reg=mod_pr_reg-obs_pr_reg

#CORR and RMSE need to be calculated after reduction to ensure same array shapes.
print round(compute_rmse(obs_pr_reg, mod_pr_reg),2)
print round(compute_corr(obs_pr_reg, mod_pr_reg),2)

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
diff_pr_reg.long_name='model-observation'
template_0.title.priority=1
template_1.title.priority=1
template_2.title.priority=1

isofill = x.createisofill()
isofill.levels=[0,0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]
#x.setcolormap('rainbow')
x.plot(mod_pr, template_0, isofill)

isofill = x.createisofill()
isofill.levels=[0,0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]

x.plot(obs_pr, template_1, isofill)
#x.setcolormap('bl_to_darkred')
#
isofill = x.createisofill()
isofill.levels=[-10,-8, -6, -4, -3, -2, -1,-0.5, 0, 0.5, 1, 2, 3, 4, 6, 8,10]
x.plot(diff_pr_reg, template_2, isofill)

x.png('test.png')



