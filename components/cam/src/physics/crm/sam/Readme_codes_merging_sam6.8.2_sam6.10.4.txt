
Here we merge CRM in SPCAM5 (https://svn-ccsm-models.cgd.ucar.edu/cam1/branches/spcam1_5_00_cam5_2_09_pnnl) 
from the version of sam6.8.2 (sam_clubb trunk revision r763) to sam6.10.4 (the pnnl branch of sam_clubb revision tag r1130: 
 http://carson.math.uwm.edu/repos/sam_repos/branches/sam_clubb_r1061_pnnl)

steps to do this: 
1. compare sam_clubb r763 with the pnnl branch of sam_CLUBB r1130
2. compare sam_clubb r763 with crm in SPCAM5 
3. compare sam_clubb r1130 with crm in SPCAM5 

copy r763, r1130 to the src directory (models/atm/cam/src/physics/)

July 1st, 2013: 
advect_mom.F90: no change from spcam5_2_09
advect_all_scalars.F90: not in r763, so copy it directly from r1130. DONE
./ADV_MPDATA/advect_scalar.F90: remove statistical part
            /advect_scalar2D.F90: no change from r1130
            /advect_scalar3D.F90: no change from r1130
            /advection.F90: no change from r1130
./ADV_UM5/advect_scalar.F90: remove statistical part
         /advect_scalar2D.F90: no change from r1130
         /advect_scalar3D.F90: no change from r1130
The above three files listed under ./crm are removed. 

boudaries.F90: copy "use grid, only: dompi" from r1130. So now boudaries.F90
	are identifical for spcam5_2_09 and r1130. 
buoyancy.F90: add betu, betd part from r1130 

clubb_sgs.F90: Incorporate changes from r1130 

clubbvars.F90: incorporate changes from r1130 to spcam5_2_09
clubb_silhs_vars.F90: directly copy it from r1130. This is not enabled in MMF. 

comparess3D.F90: the same as spcam5_2_09. No change. 
coriolis.F90: update dvdt formula from r1130. 

crm_module.F90: DONE 
crmsurface.F90: the same as spcam5_2_09. surface.F90 in r763 and r1130 are
	different, but these differences are not relevant to SPCAM. 
crmtracers.F90: the same as spcam5_2_09. No change from r763 to r1130. 
damping.F90: No change. Note: the damping of t and micro_filed is removed from
	r1130. Need to check with Marat to see whether we should incoroprate
	this change to SPCAM5 as well. 
diagnose.F90: incorporate changes from r763 to r1130 to spcam5_2_09. 

create two new subdirectories: SGS_TKE; SGS_CLUBBkvhkvm for subgrid treatment 
./SGS_TKE/diffuse_mom.F90: the same as that in spcam5_2_09
./SGS_TKE/diffuse_mom2D.F90: the same as r1130. No clubb-related codes, as
	CLUBB-related code is added in a separate directory (SGS_CLUBBkvhkvm)
./SGS_TKE/diffuse_mom3d.F90: the same as r1130. No clubb-related codes, as
        CLUBB-related code is added in a separate directory (SGS_CLUBBkvhkvm) 
./SGS_TKE/diffuse_scalar.F90: the same as spcam5_2_09 (except tkh from sgs)
./SGS_TKE/diffuse_scalar2D.F90: the same as r1130;
./SGS_TKE/diffuse_scalar3D.F90: the same as r1130;
./SGS_TKE/shear_prod2D.F90: no change from r1130 or spcam5_2_09
./SGS_TKE/shear_prod3D.F90: no change from r1130 or spcam5_2_09
./SGS_TKE/tke_full.F90: adopted the one from r1130, but add changes from
	spcam5_2_09 in terms of *_crm subroutine. In r763, tke is only updated 
	if .not.doscalar when dosmagor is true, but in r1130, no such
	restriction.  
./SGS_TKE/sgs.F90: 
	i) sgs_setparm: comment out reading namelist
       ii) no change in sgs_init. This is now called in crm_module.F90, after
	micro_init, and grdf_x, grdf_y, grdf_z are calcluated in
        sgs_init. These were calcluated in crm_module.F90 in SPCAM5. 
       iii) sgs_statistics: this may need to be removed

./SGS_CLUBBkvhkvm/sgs.F90: add docam_sfc_fluxes flag
./SGS_CLUBBkvhkvm/tke_full.F90: the same as the one from ./SGS_TKE/
./SGS_CLUBBkvhkvm/diffuse_mom.F90: remove statistics
./SGS_CLUBBkvhkvm/diffuse_mom2D.F90: add docam_sfc_fluxes
./SGS_CLUBBkvhkvm/diffuse_mom2D_xy.F90: remove CLUBB-related.
./SGS_CLUBBkvhkvm/diffuse_mom2D_z.F90: add docam_sfc_fluxes
./SGS_CLUBBkvhkvm/diffuse_mom3D.F90: add docam_sfc_fluxes
./SGS_CLUBBkvhkvm/diffuse_mom3D_xy.F90: remove clubb-related
./SGS_CLUBBkvhkvm/diffuse_mom3D_z.F90: add docam_sfc_fluxes
./SGS_CLUBBkvhkvm/diffuse_scalar.F90: incorporate changes from spcam5_2_09
./SGS_CLUBBkvhkvm/diffuse_scalar_xy.F90: incorporate changes from spcam5_2_09
./SGS_CLUBBkvhkvm/diffuse_scalar_z.F90: incorporate chagnes from spcam5_2_09
./SGS_CLUBBkvhkvm/fluxes_scalar_z.F90: incorporate changes from spcam5_2_09

domain.F90: no change from spcam5_2_09
ftt.F: no change from spcam5_2_09
forcing.F90: no change from spcam5_2_09
gammaff.c:  no change from spcam5_2_09 (seems not included in r1130 or r763)
grid.F90: Identifical to the one from r1130. There are large difference
	between r1130 and r763. Need to double check whether there is any
	potential issues. 
ice_fall.F90: no change from spcam5_2_09 
init.F90: add qtostor to the one from spcam5_2_09 
kurant.F90: adopt one from r1130
params.F90: adopt one from r1130, but add CRM-related codes. This is quite
	different from r763 and spcam5_2_09. Need to double check to see
	whether there is any poential issues
periodic.F90: adopt the one from r1130, and change CLUBB to CLUBB_CRM
precip_fall.F90: No change from spcam5_2_09 
press_grad.F90: the same as spcam5_2_09, but adopte changes from r763 to
	r1130 ( a fix by P. Bloss). 
press_rhs.F90: the same as the one from r1130 
pressure.F90: the same as the one from spcam5_2_09, but add "use params, only:
	dowallx, dowally, docolumn". Probably need to check with Marat to see
	whether we need update this. Pressure-related subroutines have littles
	change from r763 to r1130. 
random.F90: no changes from either spcam5_2_09 or r1130 
sat.F90: the same as spcam5_2_09 (quite different from r1130. But no change
	from r763 to r1130). 

NO SETDATA.F90 in spcam5, but sgs_init is called in setdata. 
	so sgs_inti is called in crm_module.F90  

setparm.F90: adopt from r1130, and add MMF-related from spcam5_2_09 
        Things to note: sgs_setparm; forz and fcor are not caclcualted here in
        r1130 any more (they are calculated in setgrid.F), but this is still
        kept here.
setperturb.F90: Tke is now treated by calling setperturb_sgs. Otherwise, it is
	the same as spcam5_2_09.  

stat_clubb.F90: Copy it from r1130. NO CHANGE YET. NEED TO BE CHANGED

stepout.F90: No change from spcam5_2_09. It is not used in spcam5. so we may
	remove it in the future. 
task_init.F90: No change from spcam5_2_09
task_util_NOMPI.F90: No change from spacm5_2_09
tke_full.F90: deleted, as this has been added into ./SGS_TKE/
utils.F90: Incorporate changes from r1130. 
vars.F90: Incorporate changes from r1130. fcory(ny) is changed to fcory(0:ny). So the calculation of fcory in
crm_module is changed as well.

./MICRO_SAM1MOM/cloud.F90: the same as spcam5_2_09
./MICRO_SAM1MOM/micro_params.F90: the same as spcam5_2_09
./MICRO_SAM1MOM/microphysics.F90: adopt changes from r1130. 
	s_ar is removed from micro_precip_fall 
./MICRO_SAM1MOM/precip_init.F90: the same as spcam5_2_09
./MICRO_SAM1MOM/precip_proc.F90: the same as spcam5_2_09
./MICRO_SAM1MOM/precip_proc_clubb.F90: adopt from r1130
./MICRO_M2005/microphysics.F90: incorporates changes from r1130
./MICRO_M2005/module_mp_graupel.F90: incorporate change from r1130. Those
	changes are quite minor, except a scaling factor is applied to contact
	freezing nucleaiton rate and homogeneous freezing of cloud droplets. 

./CLUBB/: create a new CLUBB directory for the latest CLUBB used in MMF
