import cdms2
import vcs
import MV2
import EzTemplate

# Data stuff
reference_data_path='./uvcmetric_test/plotset5/'
test_data_path='./uvcmetric_test/plotset5/'

reference_data_set='figure-set5_Global_ANN_T_plot--obs.nc'  # observation
test_data_set='figure-set5_Global_ANN_T_plot--model.nc'  # model

f_obs=cdms2.open(reference_data_path + reference_data_set)  # vars are: ['bounds_lon', 'bounds_lat', 'rv_T_ANN_NCEP']
f_mod=cdms2.open(test_data_path + test_data_set)  # vars are: ['bounds_lon', 'bounds_lat', 'rv_T_ANN_cam_output']

obs_pr=f_obs('rv_T_ANN_NCEP')
mod_pr=f_mod('rv_T_ANN_cam_output')

mod_pr_reg=mod_pr#.regrid(grid_normal,regridTool='esmf',regridMethod='conserve')
obs_pr_reg=obs_pr#.regrid(grid_normal,regridTool='esmf',regridMethod='conserve')

diff_pr_reg=mod_pr_reg-obs_pr_reg

# Graphing stuff

# How to create, store and load a vcs script
#
# Create a vcs object with a name
#template = x.createtemplate('plot_set_5_0')
#
# Modify the vcs object if needed
#template.data.x1 = 0.12300000000000001
#template.data.x2 = 0.86
#template.data.y1 = 0.7016666666666667
#template.data.y2 = 0.95
#
#template.legend.x1 = 0.88
#template.legend.x2 = 0.9193066666666667
#template.legend.y1 = 0.7265
#template.legend.y2 = 0.9251666666666666
#
# Store the vcs object as a .json file
#x.scriptobject(template, 'plot_set_5_new')
#
# Later you can load the template object like:
#x.scriptrun('plot_set_5_new.json')
#new_template = x.gettemplate('plot_set_5_0')

x = vcs.init(bg=True, geometry=(1212,1628))
x.portrait()
isofill = x.createisofill()

x.scriptrun('plot_set_5.json')
template_0 = x.gettemplate('plotset5_0_x_0')
template_1 = x.gettemplate('plotset5_0_x_1')
template_2 = x.gettemplate('plotset5_0_x_2')
x.plot(mod_pr_reg, template_0, isofill)
x.plot(obs_pr_reg, template_1, isofill)
x.plot(diff_pr_reg, template_2, isofill)

x.png('out.png')
