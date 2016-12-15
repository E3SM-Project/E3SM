import cdms2
import vcs
import MV2
import EzTemplate

reference_data_path='./uvcmetric_test/plotset5/'
test_data_path='./uvcmetric_test/plotset5/'

reference_data_set='figure-set5_Global_ANN_T_plot--obs.nc'  # observation
test_data_set='figure-set5_Global_ANN_T_plot--model.nc'  # model

f_obs=cdms2.open(reference_data_path + reference_data_set)  # vars are: ['bounds_lon', 'bounds_lat', 'rv_T_ANN_NCEP']
f_mod=cdms2.open(test_data_path + test_data_set)  # vars are: ['bounds_lon', 'bounds_lat', 'rv_T_ANN_cam_output']

obs_pr=f_obs('rv_T_ANN_NCEP')
mod_pr=f_mod('rv_T_ANN_cam_output')
print obs_pr.shape, mod_pr.shape


mod_pr_reg=mod_pr#.regrid(grid_normal,regridTool='esmf',regridMethod='conserve')
obs_pr_reg=obs_pr#.regrid(grid_normal,regridTool='esmf',regridMethod='conserve')

diff_pr_reg=mod_pr_reg-obs_pr_reg

#Plot
x=vcs.init()
x.scriptrun('plot_set_5.json')
x.plot(mod_pr_reg)
png_filename='mod.png'
x.png(png_filename)
x.clear()

x.plot(obs_pr_reg)
png_filename='obs.png'
x.png(png_filename)
x.clear()

x.plot(diff_pr_reg)
png_filename='diff.png'
x.png(png_filename)
