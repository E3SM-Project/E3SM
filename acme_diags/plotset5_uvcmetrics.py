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

'''
x = vcs.init()
isofill = x.createisofill()
x.plot(diff_pr_reg, isofill)
x.png('temp.png')
quit()
'''

# Plotting stuff
three_row_multi_template = EzTemplate.Multi(rows=3, columns=1)

parent_canvas = vcs.init(bg=True, geometry=(1212,1628))  # canvas which hold 3 children canveses
#parent_canvas.setcolormap('bl_to_darkred')  # NCAR colors
parent_canvas.portrait()
isofill = parent_canvas.createisofill()
#isofill.levels=[0,0.2, 0.5, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 17]

three_row_template = three_row_multi_template.get(legend='none')
parent_canvas.plot(mod_pr_reg, three_row_template, isofill)

three_row_template = three_row_multi_template.get(legend='none')
parent_canvas.plot(obs_pr_reg, three_row_template, isofill)

three_row_template = three_row_multi_template.get(legend='none')
parent_canvas.plot(diff_pr_reg, three_row_template, isofill)

png_filename='diff.png'
parent_canvas.png(png_filename)
