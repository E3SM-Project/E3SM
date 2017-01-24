#!/usr/bin/env python
import cdms2
import acme_diags.acme_parser
import plot_set_5

parser = acme_diags.acme_parser.ACMEParser()
parameter = parser.get_parameter()

reference_data_path = parameter.reference_data_path
reference_data_set = parameter.reference_data_set # observation
test_data_path = parameter.test_data_path
test_data_set = parameter.test_data_set # model

var = parameter.variables
season = parameter.season

f_obs = cdms2.open(reference_data_path + reference_data_set)
f_mod = cdms2.open(test_data_path + test_data_set)

obs_pr = f_obs(var, longitude=(-180, 540))
mod_pr = (f_mod('PRECC', longitude=(-180, 540)) + f_mod('PRECL', longitude=(-180, 540)))*3600.0*24.0*1000.0


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

plot_set_5.plot(obs_pr, mod_pr, obs_pr_reg, mod_pr_reg, parameter)
