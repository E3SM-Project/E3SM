module phys_control
!-----------------------------------------------------------------------
! Purpose:
!
! Provides a control interface to CAM physics packages
!
! Revision history:
! 2006-05-01  D. B. Coleman,  Creation of module
! 2009-02-13  Eaton           Replace *_{default,set}opts methods with module namelist.
!                             Add vars to indicate physics version and chemistry type.
!-----------------------------------------------------------------------

use spmd_utils,    only: masterproc
use cam_logfile,   only: iulog
use cam_abortutils,    only: endrun
use shr_kind_mod,  only: r8 => shr_kind_r8

implicit none
private
save

public :: &
   phys_ctl_readnl,   &! read namelist from file
   phys_getopts,      &! generic query method
   phys_deepconv_pbl, &! return true if deep convection is allowed in the PBL
   phys_do_flux_avg,  &! return true to average surface fluxes
   cam_physpkg_is,    &! query for the name of the physics package
   cam_chempkg_is,    &! query for the name of the chemistry package
   waccmx_is

! Private module data

character(len=16), parameter :: unset_str = 'UNSET'
integer,           parameter :: unset_int = huge(1)

! Namelist variables:
character(len=16) :: cam_physpkg          = unset_str  ! CAM physics package [cam3 | cam4 | cam5 |
                                                       !   ideal | adiabatic].
character(len=32) :: cam_chempkg          = unset_str  ! CAM chemistry package [waccm_mozart | 
                                                       !  waccm_ghg | trop_mozart | trop_ghg | 
                                                       !  trop_bam | trop_mam3 | trop_mam4 | 
                                                       !  trop_mam4_resus | trop_mam4_resus_soag |
                                                       !  trop_mam4_mom |
                                                       !  trop_mam4_resus_mom | trop_mam7 |
                                                       !  trop_mam9 |
                                                       !  linoz_mam3 | linoz_mam4_resus |
                                                       !  linoz_mam4_resus_mom |
                                                       !  linoz_mam4_resus_mom_soag |
                                                       !  super_fast_llnl | super_fast_llnl_mam3 | 
                                                       !  waccm_mozart_mam3 | none
character(len=16) :: waccmx_opt           = unset_str  ! WACCMX run option [ionosphere | neutral | off
character(len=16) :: deep_scheme          = unset_str  ! deep convection package
character(len=16) :: shallow_scheme       = unset_str  ! shallow convection package
character(len=16) :: eddy_scheme          = unset_str  ! vertical diffusion package
character(len=16) :: microp_scheme        = unset_str  ! microphysics package
character(len=16) :: macrop_scheme        = unset_str  ! macrophysics package
character(len=16) :: radiation_scheme     = unset_str  ! radiation package
integer           :: srf_flux_avg         = unset_int  ! 1 => smooth surface fluxes, 0 otherwise
integer           :: conv_water_in_rad    = unset_int  ! 0==> No; 1==> Yes-Arithmetic average;
                                                       ! 2==> Yes-Average in emissivity.

logical           :: use_subcol_microp    = .false.    ! if .true. then use sub-columns in microphysics

logical           :: atm_dep_flux         = .true.     ! true => deposition fluxes will be provided
                                                       ! to the coupler
logical           :: history_amwg         = .true.     ! output the variables used by the AMWG diag package
logical           :: history_vdiag        = .false.    ! output the variables used by the AMWG variability diag package
logical           :: history_aerosol      = .false.    ! output the MAM aerosol variables and tendencies
logical           :: history_aero_optics  = .false.    ! output the aerosol
logical           :: history_eddy         = .false.    ! output the eddy variables
logical           :: history_budget       = .false.    ! output tendencies and state variables for CAM4
                                                       ! temperature, water vapor, cloud ice and cloud
                                                       ! liquid budgets.
logical           :: ssalt_tuning         = .false.    ! sea salt tuning flag for progseasalts_intr.F90
logical           :: resus_fix            = .false.    ! to address resuspension bug fix in wetdep.F90 
logical           :: convproc_do_aer      = .false.    ! to apply unified convective transport/removal for aerosols
logical           :: convproc_do_gas      = .false.    ! to apply unified convective transport/removal for trace gases  
                                                       ! *** the unified conv. trans/removal currently does not do gases ***
!  convproc_method_activate - 1=apply abdulrazzak-ghan to entrained aerosols for lowest nlayers
!                             2=do secondary activation with prescribed supersat
integer           :: convproc_method_activate = 2      ! controls activation in the unified convective transport/removal
integer           :: mam_amicphys_optaa   = 0          ! <= 0 -- use old microphysics code (separate calls to gasaerexch, 
                                                       !                                    newnuc, and coag routines) 
                                                       !  > 0 -- use new microphysics code (single call to amicphys routine)
real(r8)          :: n_so4_monolayers_pcage = huge(1.0_r8) ! number of so4(+nh4) monolayers needed to "age" a carbon particle
real(r8)          :: micro_mg_accre_enhan_fac = huge(1.0_r8) !!Accretion enhancement factor
logical           :: liqcf_fix            = .false.    ! liq cld fraction fix calc.                     
logical           :: regen_fix            = .false.    ! aerosol regeneration bug fix for ndrop.F90 
logical           :: demott_ice_nuc       = .false.    ! use DeMott ice nucleation treatment in microphysics 
integer           :: history_budget_histfile_num = 1   ! output history file number for budget fields
logical           :: history_waccm        = .true.     ! output variables of interest for WACCM runs
logical           :: history_clubb        = .true.     ! output default CLUBB-related variables
logical           :: do_clubb_sgs
real(r8)          :: prc_coef1            = huge(1.0_r8)
real(r8)          :: prc_exp              = huge(1.0_r8)
real(r8)          :: prc_exp1             = huge(1.0_r8)
real(r8)          :: cld_sed              = huge(1.0_r8)
logical           :: do_tms
logical           :: micro_do_icesupersat
logical           :: state_debug_checks   = .false.    ! Extra checks for validity of physics_state objects
                                                       ! in physics_update.
logical, public, protected :: use_mass_borrower    = .false.     ! switch on tracer borrower, instead of using the QNEG3 clipping
logical, public, protected :: use_qqflx_fixer      = .false.     ! switch on water vapor fixer to compensate changes in qflx
logical, public, protected :: print_fixer_message  = .false.     ! switch on error message printout in log file

! Macro/micro-physics co-substeps
integer           :: cld_macmic_num_steps = 1

logical :: prog_modal_aero ! determines whether prognostic modal aerosols are present in the run.

!BSINGH -  Bugfix flags (Must be removed once the bug fix is accepted for master merge)
logical :: fix_g1_err_ndrop = .false. !BSINGH - default is false
! Option to use heterogeneous freezing
logical, public, protected :: use_hetfrz_classnuc = .false.

! Which gravity wave sources are used?
! Orographic
logical, public, protected :: use_gw_oro = .true.
! Frontogenesis
logical, public, protected :: use_gw_front = .false.
! Convective
logical, public, protected :: use_gw_convect = .false.

! Switches that turn on/off individual parameterizations.
!
! Comment by Hui Wan (PNNL, 2014-12):
! This set of switches were implemeted in a very simplistic way
! for a short-term time-step convergence test performed 
! with the "standard" CAM5 as of 2014. 
! The purpose was to identify which moist processes 
! were responsible for the poor convergence of the full model. 
! We did not make any attempt to test details of MAM
! or the non-standard model configurations/components such as 
! WACCM, CLUBB, CARMA. It is unlikely that the switches will work
! for those configurations. 

logical :: l_tracer_aero   = .true.
logical :: l_vdiff         = .true.
logical :: l_rayleigh      = .true.
logical :: l_gw_drag       = .true.
logical :: l_ac_energy_chk = .true.
logical :: l_bc_energy_fix = .true.
logical :: l_dry_adj       = .true.
logical :: l_st_mac        = .true.
logical :: l_st_mic        = .true.
logical :: l_rad           = .true.


!======================================================================= 
contains
!======================================================================= 

subroutine phys_ctl_readnl(nlfile)

   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand
   use cam_control_mod, only: cam_ctrl_set_physics_type

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'phys_ctl_readnl'

   namelist /phys_ctl_nl/ cam_physpkg, cam_chempkg, waccmx_opt, deep_scheme, shallow_scheme, &
      eddy_scheme, microp_scheme,  macrop_scheme, radiation_scheme, srf_flux_avg, &
      use_subcol_microp, atm_dep_flux, history_amwg, history_vdiag, history_aerosol, history_aero_optics, &
      history_eddy, history_budget,  history_budget_histfile_num, history_waccm, &
      conv_water_in_rad, history_clubb, do_clubb_sgs, do_tms, state_debug_checks, &
      use_mass_borrower, & 
      use_qqflx_fixer, & 
      print_fixer_message, & 
      use_hetfrz_classnuc, use_gw_oro, use_gw_front, use_gw_convect, &
      cld_macmic_num_steps, micro_do_icesupersat, &
      fix_g1_err_ndrop, ssalt_tuning, resus_fix, convproc_do_aer, &
      convproc_do_gas, convproc_method_activate, liqcf_fix, regen_fix, demott_ice_nuc, &
      mam_amicphys_optaa, n_so4_monolayers_pcage,micro_mg_accre_enhan_fac, &
      l_tracer_aero, l_vdiff, l_rayleigh, l_gw_drag, l_ac_energy_chk, &
      l_bc_energy_fix, l_dry_adj, l_st_mac, l_st_mic, l_rad, prc_coef1,prc_exp,prc_exp1,cld_sed
   !-----------------------------------------------------------------------------

   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'phys_ctl_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, phys_ctl_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
      call freeunit(unitn)
   end if

#ifdef SPMD
   ! Broadcast namelist variables
   call mpibcast(deep_scheme,      len(deep_scheme)      , mpichar, 0, mpicom)
   call mpibcast(cam_physpkg,      len(cam_physpkg)      , mpichar, 0, mpicom)
   call mpibcast(cam_chempkg,      len(cam_chempkg)      , mpichar, 0, mpicom)
   call mpibcast(waccmx_opt,       len(waccmx_opt)       , mpichar, 0, mpicom)
   call mpibcast(shallow_scheme,   len(shallow_scheme)   , mpichar, 0, mpicom)
   call mpibcast(eddy_scheme,      len(eddy_scheme)      , mpichar, 0, mpicom)
   call mpibcast(microp_scheme,    len(microp_scheme)    , mpichar, 0, mpicom)
   call mpibcast(radiation_scheme, len(radiation_scheme) , mpichar, 0, mpicom)
   call mpibcast(macrop_scheme,    len(macrop_scheme)    , mpichar, 0, mpicom)
   call mpibcast(srf_flux_avg,                    1 , mpiint,  0, mpicom)
   call mpibcast(use_subcol_microp,               1 , mpilog,  0, mpicom)
   call mpibcast(atm_dep_flux,                    1 , mpilog,  0, mpicom)
   call mpibcast(history_amwg,                    1 , mpilog,  0, mpicom)
   call mpibcast(history_vdiag,                   1 , mpilog,  0, mpicom)
   call mpibcast(history_eddy,                    1 , mpilog,  0, mpicom)
   call mpibcast(history_aerosol,                 1 , mpilog,  0, mpicom)
   call mpibcast(history_aero_optics,             1 , mpilog,  0, mpicom)
   call mpibcast(history_budget,                  1 , mpilog,  0, mpicom)
   call mpibcast(history_budget_histfile_num,     1 , mpiint,  0, mpicom)
   call mpibcast(history_waccm,                   1 , mpilog,  0, mpicom)
   call mpibcast(history_clubb,                   1 , mpilog,  0, mpicom)
   call mpibcast(do_clubb_sgs,                    1 , mpilog,  0, mpicom)
   call mpibcast(conv_water_in_rad,               1 , mpiint,  0, mpicom)
   call mpibcast(do_tms,                          1 , mpilog,  0, mpicom)
   call mpibcast(use_mass_borrower,               1 , mpilog,  0, mpicom)
   call mpibcast(use_qqflx_fixer,                 1 , mpilog,  0, mpicom)
   call mpibcast(print_fixer_message,             1 , mpilog,  0, mpicom)
   call mpibcast(micro_do_icesupersat,            1 , mpilog,  0, mpicom)
   call mpibcast(state_debug_checks,              1 , mpilog,  0, mpicom)
   call mpibcast(use_hetfrz_classnuc,             1 , mpilog,  0, mpicom)
   call mpibcast(use_gw_oro,                      1 , mpilog,  0, mpicom)
   call mpibcast(use_gw_front,                    1 , mpilog,  0, mpicom)
   call mpibcast(use_gw_convect,                  1 , mpilog,  0, mpicom)
   call mpibcast(fix_g1_err_ndrop,                1 , mpilog,  0, mpicom)
   call mpibcast(ssalt_tuning,                    1 , mpilog,  0, mpicom)
   call mpibcast(resus_fix,                       1 , mpilog,  0, mpicom)
   call mpibcast(convproc_do_aer,                 1 , mpiint,  0, mpicom)
   call mpibcast(convproc_do_gas,                 1 , mpilog,  0, mpicom)
   call mpibcast(convproc_method_activate,        1 , mpilog,  0, mpicom)
   call mpibcast(mam_amicphys_optaa,              1 , mpilog,  0, mpicom)
   call mpibcast(n_so4_monolayers_pcage,          1 , mpir8,   0, mpicom)
   call mpibcast(micro_mg_accre_enhan_fac,        1 , mpir8,   0, mpicom)
   call mpibcast(liqcf_fix,                       1 , mpilog,  0, mpicom)
   call mpibcast(regen_fix,                       1 , mpilog,  0, mpicom)
   call mpibcast(demott_ice_nuc,                  1 , mpilog,  0, mpicom)
   call mpibcast(l_tracer_aero,                   1 , mpilog,  0, mpicom)
   call mpibcast(l_vdiff,                         1 , mpilog,  0, mpicom)
   call mpibcast(l_rayleigh,                      1 , mpilog,  0, mpicom)
   call mpibcast(l_gw_drag,                       1 , mpilog,  0, mpicom)
   call mpibcast(l_ac_energy_chk,                 1 , mpilog,  0, mpicom)
   call mpibcast(l_bc_energy_fix,                 1 , mpilog,  0, mpicom)
   call mpibcast(l_dry_adj,                       1 , mpilog,  0, mpicom)
   call mpibcast(l_st_mac,                        1 , mpilog,  0, mpicom)
   call mpibcast(l_st_mic,                        1 , mpilog,  0, mpicom)
   call mpibcast(l_rad,                           1 , mpilog,  0, mpicom)
   call mpibcast(cld_macmic_num_steps,            1 , mpiint,  0, mpicom)
   call mpibcast(prc_coef1,                       1 , mpir8,   0, mpicom)
   call mpibcast(prc_exp,                         1 , mpir8,   0, mpicom)
   call mpibcast(prc_exp1,                        1 , mpir8,   0, mpicom)
   call mpibcast(cld_sed,                         1 , mpir8,   0, mpicom)
#endif

   call cam_ctrl_set_physics_type(cam_physpkg)

   ! Error checking:

   ! Defaults for PBL and microphysics are set in build-namelist.  Check here that
   ! values have been set to guard against problems with hand edited namelists.

   ! WACCM-X run option set in build-namelist. Check for valid values
   if (.not. (waccmx_opt == 'ionosphere' .or. waccmx_opt == 'neutral' .or. waccmx_opt == 'off')) then
      write(iulog,*)'waccm: illegal value of waccmx_opt:', waccmx_opt
      call endrun('waccm: illegal value of waccmx_opt')
   endif
   if (.not. (shallow_scheme .eq. 'Hack' .or. shallow_scheme .eq. 'UW' .or. &
              shallow_scheme .eq. 'CLUBB_SGS' .or. shallow_scheme == 'UNICON' &
       .or. shallow_scheme.eq.'off')) then
      write(iulog,*)'phys_setopts: illegal value of shallow_scheme:', shallow_scheme
      call endrun('phys_setopts: illegal value of shallow_scheme')
   endif
   if (.not. (eddy_scheme .eq. 'HB' .or. eddy_scheme .eq. 'HBR' .or. eddy_scheme .eq. 'diag_TKE' .or. &
              eddy_scheme .eq. 'CLUBB_SGS') ) then
      write(iulog,*)'phys_setopts: illegal value of eddy_scheme:', eddy_scheme
      call endrun('phys_setopts: illegal value of eddy_scheme')
   endif
   if ((microp_scheme /= 'MG' .and. microp_scheme /= 'RK')) then
      write(iulog,*)'phys_setopts: illegal value of microp_scheme:', microp_scheme
      call endrun('phys_setopts: illegal value of microp_scheme')
   endif

   ! Check compatibility of eddy & shallow schemes
   if (( shallow_scheme .eq. 'UW' ) .and. ( eddy_scheme .ne. 'diag_TKE' )) then
      write(iulog,*)'Do you really want to run UW shallow scheme without diagnostic TKE eddy scheme? Quiting'
      call endrun('shallow convection and eddy scheme may be incompatible')
   endif

   if (( shallow_scheme .eq. 'Hack' ) .and. ( ( eddy_scheme .ne. 'HB' ) .and. ( eddy_scheme .ne. 'HBR' ))) then
      write(iulog,*)'Do you really want to run Hack shallow scheme with a non-standard eddy scheme? Quiting.'
      call endrun('shallow convection and eddy scheme may be incompatible')
   endif

   ! Check compatibility of PBL and Microphysics schemes
   if (( eddy_scheme .eq. 'diag_TKE' ) .and. ( microp_scheme .eq. 'RK' )) then
      write(iulog,*)'UW PBL is not compatible with RK microphysics.  Quiting'
      call endrun('PBL and Microphysics schemes incompatible')
   endif
   
   ! Add a check to make sure CLUBB and MG are used together
   if ( do_clubb_sgs .and. ( microp_scheme .ne. 'MG')) then
      write(iulog,*)'CLUBB is only compatible with MG microphysics.  Quiting'
      call endrun('CLUBB and microphysics schemes incompatible')
   endif

   ! Check that eddy_scheme, macrop_scheme, shallow_scheme are all set to CLUBB_SGS if do_clubb_sgs is true
   if (do_clubb_sgs) then
      if (eddy_scheme .ne. 'CLUBB_SGS' .or. macrop_scheme .ne. 'CLUBB_SGS' .or. shallow_scheme .ne. 'CLUBB_SGS') then
         write(iulog,*)'eddy_scheme, macrop_scheme and shallow_scheme must all be CLUBB_SGS.  Quiting'
         call endrun('CLUBB and eddy, macrop or shallow schemes incompatible')
      endif
   endif

   ! Macro/micro co-substepping support.
   if (cld_macmic_num_steps > 1) then
      if (microp_scheme /= "MG" .or. (macrop_scheme /= "park" .and. macrop_scheme /= "CLUBB_SGS")) then
         call endrun ("Setting cld_macmic_num_steps > 1 is only &
              &supported with Park or CLUBB macrophysics and MG microphysics.")
      end if
   end if

   ! prog_modal_aero determines whether prognostic modal aerosols are present in the run.
   prog_modal_aero = (     cam_chempkg_is('trop_mam3') &
                      .or. cam_chempkg_is('trop_mam4') &
                      .or. cam_chempkg_is('trop_mam4_resus') &
                      .or. cam_chempkg_is('trop_mam4_resus_soag') &
                      .or. cam_chempkg_is('trop_mam4_mom') &
                      .or. cam_chempkg_is('trop_mam4_resus_mom') &
                      .or. cam_chempkg_is('trop_mam7') &
                      .or. cam_chempkg_is('trop_mam9') &
                      .or. cam_chempkg_is('linoz_mam3') &
                      .or. cam_chempkg_is('linoz_mam4_resus') &
                      .or. cam_chempkg_is('linoz_mam4_resus_soag') &
                      .or. cam_chempkg_is('linoz_mam4_resus_mom') &
                      .or. cam_chempkg_is('linoz_mam4_resus_mom_soag') &
                      .or. cam_chempkg_is('super_fast_llnl_mam3') &
                      .or. cam_chempkg_is('trop_mozart_mam3') &
                      .or. cam_chempkg_is('trop_strat_mam3') &
                      .or. cam_chempkg_is('trop_strat_mam7') &
                      .or. cam_chempkg_is('waccm_mozart_mam3'))
end subroutine phys_ctl_readnl

!===============================================================================

logical function cam_physpkg_is(name)

   ! query for the name of the physics package

   character(len=*) :: name
   
   cam_physpkg_is = (trim(name) == trim(cam_physpkg))
end function cam_physpkg_is

!===============================================================================

logical function cam_chempkg_is(name)

   ! query for the name of the chemics package

   character(len=*) :: name
   
   cam_chempkg_is = (trim(name) == trim(cam_chempkg))
end function cam_chempkg_is

!===============================================================================

logical function waccmx_is(name)

   ! query for the name of the waccmx run option

   character(len=*) :: name
   
   waccmx_is = (trim(name) == trim(waccmx_opt))
end function waccmx_is

!===============================================================================

subroutine phys_getopts(deep_scheme_out, shallow_scheme_out, eddy_scheme_out, microp_scheme_out, &
                        radiation_scheme_out, use_subcol_microp_out, atm_dep_flux_out, &
                        history_amwg_out, history_vdiag_out, history_aerosol_out, history_aero_optics_out, history_eddy_out, &
                        history_budget_out, history_budget_histfile_num_out, history_waccm_out, &
                        history_clubb_out, conv_water_in_rad_out, cam_chempkg_out, prog_modal_aero_out, macrop_scheme_out, &
                        do_clubb_sgs_out, do_tms_out, state_debug_checks_out, &
                        use_mass_borrower_out, & 
                        use_qqflx_fixer_out, & 
                        print_fixer_message_out, & 
                        cld_macmic_num_steps_out, micro_do_icesupersat_out, &
                        fix_g1_err_ndrop_out, ssalt_tuning_out,resus_fix_out,convproc_do_aer_out,  &
                        convproc_do_gas_out, convproc_method_activate_out, mam_amicphys_optaa_out, n_so4_monolayers_pcage_out, &
                        micro_mg_accre_enhan_fac_out, liqcf_fix_out, regen_fix_out,demott_ice_nuc_out      &
                       ,l_tracer_aero_out, l_vdiff_out, l_rayleigh_out, l_gw_drag_out, l_ac_energy_chk_out  &
                       ,l_bc_energy_fix_out, l_dry_adj_out, l_st_mac_out, l_st_mic_out, l_rad_out  &
                       ,prc_coef1_out,prc_exp_out,prc_exp1_out, cld_sed_out)

!-----------------------------------------------------------------------
! Purpose: Return runtime settings
!          deep_scheme_out   : deep convection scheme
!          shallow_scheme_out: shallow convection scheme
!          eddy_scheme_out   : vertical diffusion scheme
!          microp_scheme_out : microphysics scheme
!          radiation_scheme_out : radiation_scheme
!-----------------------------------------------------------------------

   character(len=16), intent(out), optional :: deep_scheme_out
   character(len=16), intent(out), optional :: shallow_scheme_out
   character(len=16), intent(out), optional :: eddy_scheme_out
   character(len=16), intent(out), optional :: microp_scheme_out
   character(len=16), intent(out), optional :: radiation_scheme_out
   character(len=16), intent(out), optional :: macrop_scheme_out
   logical,           intent(out), optional :: use_subcol_microp_out
   logical,           intent(out), optional :: atm_dep_flux_out
   logical,           intent(out), optional :: history_amwg_out
   logical,           intent(out), optional :: history_vdiag_out
   logical,           intent(out), optional :: history_eddy_out
   logical,           intent(out), optional :: history_aerosol_out
   logical,           intent(out), optional :: history_aero_optics_out
   logical,           intent(out), optional :: history_budget_out
   integer,           intent(out), optional :: history_budget_histfile_num_out
   logical,           intent(out), optional :: history_waccm_out
   logical,           intent(out), optional :: history_clubb_out
   logical,           intent(out), optional :: do_clubb_sgs_out
   logical,           intent(out), optional :: micro_do_icesupersat_out
   integer,           intent(out), optional :: conv_water_in_rad_out
   character(len=32), intent(out), optional :: cam_chempkg_out
   logical,           intent(out), optional :: prog_modal_aero_out
   logical,           intent(out), optional :: do_tms_out
   logical,           intent(out), optional :: use_mass_borrower_out
   logical,           intent(out), optional :: use_qqflx_fixer_out
   logical,           intent(out), optional :: print_fixer_message_out
   logical,           intent(out), optional :: state_debug_checks_out
   logical,           intent(out), optional :: fix_g1_err_ndrop_out!BSINGH - bugfix for ndrop.F90
   logical,           intent(out), optional :: ssalt_tuning_out    
   logical,           intent(out), optional :: resus_fix_out       
   logical,           intent(out), optional :: convproc_do_aer_out 
   logical,           intent(out), optional :: convproc_do_gas_out 
   integer,           intent(out), optional :: convproc_method_activate_out 
   integer,           intent(out), optional :: mam_amicphys_optaa_out
   real(r8),          intent(out), optional :: n_so4_monolayers_pcage_out
   real(r8),          intent(out), optional :: micro_mg_accre_enhan_fac_out
   logical,           intent(out), optional :: liqcf_fix_out       
   logical,           intent(out), optional :: regen_fix_out       
   logical,           intent(out), optional :: demott_ice_nuc_out  


   logical,           intent(out), optional :: l_tracer_aero_out
   logical,           intent(out), optional :: l_vdiff_out
   logical,           intent(out), optional :: l_rayleigh_out
   logical,           intent(out), optional :: l_gw_drag_out
   logical,           intent(out), optional :: l_ac_energy_chk_out
   logical,           intent(out), optional :: l_bc_energy_fix_out
   logical,           intent(out), optional :: l_dry_adj_out
   logical,           intent(out), optional :: l_st_mac_out
   logical,           intent(out), optional :: l_st_mic_out
   logical,           intent(out), optional :: l_rad_out
   integer,           intent(out), optional :: cld_macmic_num_steps_out
   real(r8),          intent(out), optional :: prc_coef1_out
   real(r8),          intent(out), optional :: prc_exp_out
   real(r8),          intent(out), optional :: prc_exp1_out
   real(r8),          intent(out), optional :: cld_sed_out

   if ( present(deep_scheme_out         ) ) deep_scheme_out          = deep_scheme
   if ( present(shallow_scheme_out      ) ) shallow_scheme_out       = shallow_scheme
   if ( present(eddy_scheme_out         ) ) eddy_scheme_out          = eddy_scheme
   if ( present(microp_scheme_out       ) ) microp_scheme_out        = microp_scheme
   if ( present(radiation_scheme_out    ) ) radiation_scheme_out     = radiation_scheme

   if ( present(use_subcol_microp_out   ) ) use_subcol_microp_out    = use_subcol_microp
   if ( present(macrop_scheme_out       ) ) macrop_scheme_out        = macrop_scheme
   if ( present(atm_dep_flux_out        ) ) atm_dep_flux_out         = atm_dep_flux
   if ( present(history_aerosol_out     ) ) history_aerosol_out      = history_aerosol
   if ( present(history_aero_optics_out ) ) history_aero_optics_out  = history_aero_optics
   if ( present(history_budget_out      ) ) history_budget_out       = history_budget
   if ( present(history_amwg_out        ) ) history_amwg_out         = history_amwg
   if ( present(history_vdiag_out       ) ) history_vdiag_out        = history_vdiag
   if ( present(history_eddy_out        ) ) history_eddy_out         = history_eddy
   if ( present(history_budget_histfile_num_out ) ) history_budget_histfile_num_out = history_budget_histfile_num
   if ( present(history_waccm_out       ) ) history_waccm_out        = history_waccm
   if ( present(history_clubb_out       ) ) history_clubb_out        = history_clubb
   if ( present(do_clubb_sgs_out        ) ) do_clubb_sgs_out         = do_clubb_sgs
   if ( present(micro_do_icesupersat_out )) micro_do_icesupersat_out = micro_do_icesupersat
   if ( present(conv_water_in_rad_out   ) ) conv_water_in_rad_out    = conv_water_in_rad
   if ( present(cam_chempkg_out         ) ) cam_chempkg_out          = cam_chempkg
   if ( present(prog_modal_aero_out     ) ) prog_modal_aero_out      = prog_modal_aero
   if ( present(do_tms_out              ) ) do_tms_out               = do_tms
   if ( present(use_mass_borrower_out   ) ) use_mass_borrower_out    = use_mass_borrower
   if ( present(use_qqflx_fixer_out     ) ) use_qqflx_fixer_out      = use_qqflx_fixer
   if ( present(print_fixer_message_out ) ) print_fixer_message_out  = print_fixer_message
   if ( present(state_debug_checks_out  ) ) state_debug_checks_out   = state_debug_checks
   if ( present(fix_g1_err_ndrop_out    ) ) fix_g1_err_ndrop_out     = fix_g1_err_ndrop
   if ( present(ssalt_tuning_out        ) ) ssalt_tuning_out         = ssalt_tuning   
   if ( present(resus_fix_out           ) ) resus_fix_out            = resus_fix      
   if ( present(convproc_do_aer_out     ) ) convproc_do_aer_out      = convproc_do_aer
   if ( present(convproc_do_gas_out     ) ) convproc_do_gas_out      = convproc_do_gas
   if ( present(convproc_method_activate_out ) ) convproc_method_activate_out = convproc_method_activate
   if ( present(mam_amicphys_optaa_out  ) ) mam_amicphys_optaa_out  = mam_amicphys_optaa
   if ( present(n_so4_monolayers_pcage_out  ) ) n_so4_monolayers_pcage_out = n_so4_monolayers_pcage
   if ( present(micro_mg_accre_enhan_fac_out)) micro_mg_accre_enhan_fac_out = micro_mg_accre_enhan_fac
   if ( present(liqcf_fix_out           ) ) liqcf_fix_out            = liqcf_fix      
   if ( present(regen_fix_out           ) ) regen_fix_out            = regen_fix      
   if ( present(demott_ice_nuc_out      ) ) demott_ice_nuc_out       = demott_ice_nuc 
   if ( present(l_tracer_aero_out       ) ) l_tracer_aero_out     = l_tracer_aero
   if ( present(l_vdiff_out             ) ) l_vdiff_out           = l_vdiff
   if ( present(l_rayleigh_out          ) ) l_rayleigh_out        = l_rayleigh
   if ( present(l_gw_drag_out           ) ) l_gw_drag_out         = l_gw_drag
   if ( present(l_ac_energy_chk_out     ) ) l_ac_energy_chk_out   = l_ac_energy_chk
   if ( present(l_bc_energy_fix_out     ) ) l_bc_energy_fix_out   = l_bc_energy_fix
   if ( present(l_dry_adj_out           ) ) l_dry_adj_out         = l_dry_adj
   if ( present(l_st_mac_out            ) ) l_st_mac_out          = l_st_mac
   if ( present(l_st_mic_out            ) ) l_st_mic_out          = l_st_mic
   if ( present(l_rad_out               ) ) l_rad_out             = l_rad
   if ( present(cld_macmic_num_steps_out) ) cld_macmic_num_steps_out = cld_macmic_num_steps
   if ( present(prc_coef1_out           ) ) prc_coef1_out            = prc_coef1
   if ( present(prc_exp_out             ) ) prc_exp_out              = prc_exp
   if ( present(prc_exp1_out            ) ) prc_exp1_out             = prc_exp1 
   if ( present(cld_sed_out             ) ) cld_sed_out              = cld_sed

end subroutine phys_getopts

!===============================================================================

function phys_deepconv_pbl()

  logical phys_deepconv_pbl

   ! Don't allow deep convection in PBL if running UW PBL scheme
   if ( (eddy_scheme .eq. 'diag_TKE' ) .or. (shallow_scheme .eq. 'UW' ) ) then
      phys_deepconv_pbl = .true.
   else
      phys_deepconv_pbl = .false.
   endif

   return

end function phys_deepconv_pbl

!===============================================================================

function phys_do_flux_avg()

   logical :: phys_do_flux_avg
   !----------------------------------------------------------------------

   phys_do_flux_avg = .false.
   if (srf_flux_avg == 1) phys_do_flux_avg = .true.

end function phys_do_flux_avg

!===============================================================================
end module phys_control
