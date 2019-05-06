module runtime_opts

!-----------------------------------------------------------------------
!
! Provide driver level routine for making calls to the namelist readers
! for the infrastructure and the dycore and physics parameterizations.
!
!-----------------------------------------------------------------------

use shr_kind_mod,    only: r8=>shr_kind_r8

implicit none
private
save

public :: read_namelist

!=======================================================================
contains
!=======================================================================

!subroutine read_namelist(nlfilename, single_column, scmlat, scmlon)
subroutine read_namelist(nlfilename_in, single_column_in, scmlat_in, scmlon_in)

   use cam_initfiles,       only: cam_initfiles_readnl
   use constituents,        only: cnst_readnl

   use phys_grid,           only: phys_grid_readnl

   use chem_surfvals,       only: chem_surfvals_readnl
   use check_energy,        only: check_energy_readnl
   use radiation,           only: radiation_readnl
   use carma_flags_mod,     only: carma_readnl
   use co2_cycle,           only: co2_cycle_readnl
   !use scamMod,             only: scam_readnl

   use spmd_utils,          only: spmd_utils_readnl
   use cam_history,         only: history_readnl
   use physconst,           only: physconst_readnl
   use physics_buffer,      only: pbuf_readnl
   use phys_control,        only: phys_ctl_readnl
   use wv_saturation,       only: wv_sat_readnl
   use ref_pres,            only: ref_pres_readnl
   use cam3_aero_data,      only: cam3_aero_data_readnl
   use cam3_ozone_data,     only: cam3_ozone_data_readnl
   use dadadj_cam,          only: dadadj_readnl
   use macrop_driver,       only: macrop_driver_readnl
   use microp_driver,       only: microp_driver_readnl
   use microp_aero,         only: microp_aero_readnl
   use subcol,              only: subcol_readnl
   use cloud_fraction,      only: cldfrc_readnl
   use cldfrc2m,            only: cldfrc2m_readnl
   use rk_stratiform,       only: rk_stratiform_readnl
   use unicon_cam,          only: unicon_cam_readnl
   use zm_conv_intr,        only: zm_conv_readnl
   use hk_conv,             only: hkconv_readnl
   use uwshcu,              only: uwshcu_readnl
   use pkg_cld_sediment,    only: cld_sediment_readnl
   use gw_drag,             only: gw_drag_readnl
   use qbo,                 only: qbo_readnl
   use iondrag,             only: iondrag_readnl
   use phys_debug_util,     only: phys_debug_readnl
   use conv_water,          only: conv_water_readnl
   use rad_constituents,    only: rad_cnst_readnl
   use radiation_data,      only: rad_data_readnl
   use modal_aer_opt,       only: modal_aer_opt_readnl
   use clubb_intr,          only: clubb_readnl
   use chemistry,           only: chem_readnl
   use prescribed_volcaero, only: prescribed_volcaero_readnl
   use prescribed_strataero,only: prescribed_strataero_readnl
   use aerodep_flx,         only: aerodep_flx_readnl
   use solar_data,          only: solar_data_readnl
   use tropopause,          only: tropopause_readnl
   use aoa_tracers,         only: aoa_tracers_readnl
   use prescribed_ozone,    only: prescribed_ozone_readnl
   use prescribed_aero,     only: prescribed_aero_readnl
   use prescribed_ghg,      only: prescribed_ghg_readnl
   use aircraft_emit,       only: aircraft_emit_readnl
   use cospsimulator_intr,  only: cospsimulator_intr_readnl
   use vertical_diffusion,  only: vd_readnl
   use rayleigh_friction,   only: rayleigh_friction_readnl

   use cam_diagnostics,     only: diag_readnl
   use radheat,             only: radheat_readnl
#if ( defined OFFLINE_DYN )
   use metdata,             only: metdata_readnl
#endif
   !use offline_driver,      only: offline_driver_readnl
   !use inic_analytic_utils, only: analytic_ic_readnl
   !use rate_diags,          only: rate_diags_readnl
   use tracers,             only: tracers_readnl

   use dyn_comp,            only: dyn_readnl
   !use ionosphere_interface,only: ionosphere_readnl
   use qneg_module,         only: qneg_readnl

   use time_manager,        only: dtime
   !---------------------------Arguments-----------------------------------

   character(len=*), intent(in) :: nlfilename_in
   logical,          intent(in) :: single_column_in
   real(r8),         intent(in) :: scmlat_in
   real(r8),         intent(in) :: scmlon_in

   !---------------------------Local variables-----------------------------
   character(len=*), parameter ::  subname = "read_namelist"

   !-----------------------------------------------------------------------

   ! Call subroutines for modules to read their own namelist.
   ! In some cases namelist default values may depend on settings from
   ! other modules, so there may be an order dependence in the following
   ! calls.
   ! ***N.B.*** In particular, physconst_readnl should be called before
   !            the other readnl methods in case that method is used to set
   !            physical constants, some of which are set at runtime
   !            by the physconst_readnl method.
   ! Modules that read their own namelist are responsible for making sure
   ! all processes receive the values.

   call spmd_utils_readnl(nlfilename_in)
   call phys_grid_readnl(nlfilename_in)
   call physconst_readnl(nlfilename_in)
!++bee 13 Oct 2015, need to fix the pbuf_global_allocate functionality, then
!                   can uncomment the pbuf_readnl line
!   call pbuf_readnl(nlfilename)
   call cam_initfiles_readnl(nlfilename_in)
   call cnst_readnl(nlfilename_in)
   call history_readnl(nlfilename_in,dtime)
   call chem_surfvals_readnl(nlfilename_in)
   call phys_ctl_readnl(nlfilename_in)
   call wv_sat_readnl(nlfilename_in)
   call ref_pres_readnl(nlfilename_in)
   call cam3_aero_data_readnl(nlfilename_in)
   call cam3_ozone_data_readnl(nlfilename_in)
   call dadadj_readnl(nlfilename_in)
   call macrop_driver_readnl(nlfilename_in)
   call microp_driver_readnl(nlfilename_in)
   call microp_aero_readnl(nlfilename_in)
   call clubb_readnl(nlfilename_in)
   call subcol_readnl(nlfilename_in)
   call cldfrc_readnl(nlfilename_in)
   call cldfrc2m_readnl(nlfilename_in)
   call unicon_cam_readnl(nlfilename_in)
   !call zm_conv_readnl(nlfilename_in)   ! mdb:  comment out since it's not used with q3d
   call rk_stratiform_readnl(nlfilename_in)
   call hkconv_readnl(nlfilename_in)
   call uwshcu_readnl(nlfilename_in)
   call cld_sediment_readnl(nlfilename_in)
   call gw_drag_readnl(nlfilename_in)
   call qbo_readnl(nlfilename_in)
   call iondrag_readnl(nlfilename_in)
   call phys_debug_readnl(nlfilename_in)
   call conv_water_readnl(nlfilename_in)
   call radiation_readnl(nlfilename_in)
   call rad_cnst_readnl(nlfilename_in)
   call rad_data_readnl(nlfilename_in)
   call modal_aer_opt_readnl(nlfilename_in)
   call chem_readnl(nlfilename_in)
   call prescribed_volcaero_readnl(nlfilename_in)
   call prescribed_strataero_readnl(nlfilename_in)
   call solar_data_readnl(nlfilename_in)
   call carma_readnl(nlfilename_in)
   call tropopause_readnl(nlfilename_in)
   call aoa_tracers_readnl(nlfilename_in)
   call tracers_readnl(nlfilename_in)
   call aerodep_flx_readnl(nlfilename_in)
   call prescribed_ozone_readnl(nlfilename_in)
   call prescribed_aero_readnl(nlfilename_in)
   call prescribed_ghg_readnl(nlfilename_in)
   call co2_cycle_readnl(nlfilename_in)
   call aircraft_emit_readnl(nlfilename_in)
   call cospsimulator_intr_readnl(nlfilename_in)
   call diag_readnl(nlfilename_in)
   call check_energy_readnl(nlfilename_in)
   call radheat_readnl(nlfilename_in)
   call vd_readnl(nlfilename_in)
   call rayleigh_friction_readnl(nlfilename_in)
#if ( defined OFFLINE_DYN )
   call metdata_readnl(nlfilename_in)
#endif
!### for q3d
   !call offline_driver_readnl(nlfilename_in)
   !call analytic_ic_readnl(nlfilename_in)
   !call rate_diags_readnl(nlfilename_in)
   !call scam_readnl(nlfilename_in, single_column_in, scmlat_in, scmlon_in)

   call dyn_readnl(nlfilename_in)
   !call ionosphere_readnl(nlfilename_in)
   call qneg_readnl(nlfilename_in)

end subroutine read_namelist


!=======================================================================

end module runtime_opts
