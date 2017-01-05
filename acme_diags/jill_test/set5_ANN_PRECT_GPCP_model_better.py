import cdms2,cdutil
import vcs
import MV2
import EzTemplate
from metrics.computation.reductions import set_mean
#from regrid import Regridder


reference_data_path='/space1/test_data/obs_for_diagnostics/'  # observation
test_data_path='/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'  # model
#reference_data_path='./'
#test_data_path='./'

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

model_grid=mod_pr.getGrid()
mod_pr_reg=mod_pr
obs_pr_reg=obs_pr.regrid(model_grid,regridTool='esmf',regridMethod='linear')
max_obs=obs_pr_reg.max()
min_obs=obs_pr_reg.min()
print min_obs, max_obs

x = vcs.init(bg=True, geometry=(1212,1628))
x.portrait()

x.scriptrun('plot_set_5.json')
template_0 = x.gettemplate('plotset5_0_x_0')
template_1 = x.gettemplate('plotset5_0_x_1')
template_2 = x.gettemplate('plotset5_0_x_2')
#template_0.blank(["title","mean","min","max","dataname","crdate","crtime","units"]) ## Turn off additional information
isofill = x.createisofill()
isofill.levels=[0,0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]
#x.setcolormap('rainbow')
x.plot(mod_pr_reg, template_0, isofill)
isofill = x.createisofill()
isofill.levels=[0,0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]
#x.setcolormap('rainbow')
x.plot(obs_pr_reg, template_1, isofill)
#x.setcolormap('bl_to_darkred')
isofill = x.createisofill()
isofill.levels=[-10,-8, -6, -4, -3, -2, -1,-0.5, 0, 0.5, 1, 2, 3, 4, 6, 8,10]
x.plot(mod_pr_reg-obs_pr_reg, template_2, isofill)

x.png('test.png')




##x.show('isoline')
##x.createisoline([0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17])
#
#png_filename='mod.png'
#x.png(png_filename)
#
#y=vcs.init()
#y.scriptrun('plot_set_5.json')
#y.plot(obs_pr_reg)
#
#png_filename='obs.png'
#y.png(png_filename)
#
#
#z=vcs.init()
#z.scriptrun('plot_set_5.json')
#z.plot(diff_pr_reg)
#
#png_filename='diff.png'
#z.png(png_filename)
