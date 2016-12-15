import cdms2
import vcs
import MV2
import EzTemplate

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

#Regrid both model and obs by creating a common grid nlat=36, nlat = 72
nlat=72
nlon=144
lat=mod_pr.getLatitude()
lon=mod_pr.getLongitude()
model_grid=mod_pr.getGrid()
deltalat=(lat[-1]-lat[0])/(nlat-1)
deltalon=(lon[-1]-lon[0])/(nlon-1)
grid_normal=cdms2.createUniformGrid(lat[0], nlat, deltalat, lon[0], nlon, deltalon, order='yx')

mod_pr_reg=mod_pr.regrid(grid_normal,regridTool='esmf',regridMethod='conserve')
obs_pr_reg=obs_pr.regrid(grid_normal,regridTool='esmf',regridMethod='conserve')

diff_pr_reg=mod_pr_reg-obs_pr_reg

#Plot
x=vcs.init()
#x.scriptrun('plot_set_5.json')
x.plot(diff_pr_reg)

png_filename='test.png'
x.png(png_filename)
