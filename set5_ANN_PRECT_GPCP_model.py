import cdms2
import vcs
import MV2
import EzTemplate

obs_path='/space1/test_data/obs_for_diagnostics/'
mod_path='/space/golaz1/ACME_simulations/20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01/pp/clim_rgr/0070-0099/'

#Read in data
obs_name='GPCP_v2.2_ANN_climo.nc'
mod_name='20160520.A_WCYCL1850.ne30_oEC.edison.alpha6_01_ANN_climo.nc'

f_obs=cdms2.open(obs_path+obs_name)
f_mod=cdms2.open(mod_path+mod_name)

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
x.setcolormap('rainbow')
aa1=x.createboxfill()
aa1.level_1=-15
aa1.level_2=15
aa1.ext_1='y'
aa1.ext_2='y'
x.plot(mod_pr_reg,aa1)

#t=M.get(legend='local')
#aa2=x.createboxfill()
##aa2.level_1=-15
##aa2.level_2=15
#aa2.ext_1='y'
#aa2.ext_2='y'
##x.plot(obs_pr_reg,t,aa2)
#x.plot(obs_pr_reg,aa2)
#
#aa3=x.createboxfill()
#aa3.level_1=-5
#aa3.level_2=5
#aa3.ext_1='y'
#aa3.ext_2='y'
#t=M.get(legend='local')
#x.plot(diff_pr_reg,t,aa3)
png_filename='test.png'#"".join([figure_directory, 'WC_Diag3_', model_case, '_', season, '.png'])
x.png(png_filename)
x.clear()
