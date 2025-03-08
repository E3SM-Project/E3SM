module physpkg
  !-----------------------------------------------------------------------
  ! Purpose:
  !
  ! Provides the interface to CAM physics package
  !
  ! Revision history:
  ! Aug  2005,  E. B. Kluzek,  Creation of module from physpkg subroutine
  ! 2005-10-17  B. Eaton       Add contents of inti.F90 to phys_init().  Add
  !                            initialization of grid info in phys_state.
  ! Nov 2010    A. Gettelman   Put micro/macro physics into separate routines
  ! July 2015   B. Singh       Added code for unified convective transport
  !-----------------------------------------------------------------------


  use shr_kind_mod,     only: i8 => shr_kind_i8, r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use physconst,        only: latvap, latice, rh2o
  use physics_types,    only: physics_state, physics_tend, physics_state_set_grid, &
       physics_ptend, physics_tend_init,    &
       physics_type_alloc, physics_ptend_dealloc,&
       physics_state_alloc, physics_state_dealloc, physics_tend_alloc, physics_tend_dealloc
  use conditional_diag, only: cnd_diag_t, cnd_diag_info
  use conditional_diag_main, only: cnd_diag_checkpoint
  use physics_update_mod,  only: physics_update, physics_update_init, hist_vars, nvars_prtrb_hist, get_var
  use phys_grid,        only: get_ncols_p, print_cost_p, update_cost_p, phys_proc_cost
  use phys_gmean,       only: gmean_mass
  use ppgrid,           only: begchunk, endchunk, pcols, pver, pverp, psubcols
  use constituents,     only: pcnst, cnst_name, cnst_get_ind, &
                              setup_moist_indices, icldice, icldliq, irain, isnow
  use camsrfexch,       only: cam_out_t, cam_in_t

  use cam_control_mod,  only: ideal_phys, adiabatic
  use phys_control,     only: phys_do_flux_avg, phys_getopts, waccmx_is
  use zm_conv,          only: do_zmconv_dcape_ull => trigdcape_ull, &
                              do_zmconv_dcape_only => trig_dcape_only
  use iop_data_mod,     only: single_column
  use flux_avg,         only: flux_avg_init
  use infnan,           only: posinf, assignment(=)
#ifdef SPMD
  use mpishorthand
#endif
  use perf_mod
  use cam_logfile,     only: iulog
  use camsrfexch,      only: cam_export

  use modal_aero_calcsize,    only: modal_aero_calcsize_init, &
                                    modal_aero_calcsize_reg
  use modal_aero_wateruptake, only: modal_aero_wateruptake_init, &
                                    modal_aero_wateruptake_reg
  use chem_mods,    only : gas_pcnst
  use mo_tracname,  only : solsym
  use mo_chm_diags, only : aer_species
  use mo_gas_phase_chemdr, only : gas_ac_name, gas_ac_name_2D

  use check_energy,    only: nstep_ignore_diagn1, nstep_ignore_diagn2, &
                             check_energy_gmean, check_energy_gmean_additional_diagn, &
                             print_additional_diagn

  implicit none
  private

  !  Physics buffer index
  integer ::  teout_idx          = 0  

  integer ::  tini_idx           = 0 
  integer ::  qini_idx           = 0 
  integer ::  cldliqini_idx      = 0 
  integer ::  cldiceini_idx      = 0 
  integer ::  static_ener_ac_idx = 0
  integer ::  water_vap_ac_idx   = 0

  integer ::  prec_str_idx       = 0
  integer ::  snow_str_idx       = 0
  integer ::  prec_sed_idx       = 0
  integer ::  snow_sed_idx       = 0
  integer ::  prec_pcw_idx       = 0
  integer ::  snow_pcw_idx       = 0
  integer ::  prec_dp_idx        = 0
  integer ::  snow_dp_idx        = 0
  integer ::  prec_sh_idx        = 0
  integer ::  snow_sh_idx        = 0
  integer ::  rice2_idx          = 0
  integer ::  gas_ac_idx         = 0
  integer :: species_class(pcnst)  = -1 !BSINGH: Moved from modal_aero_data.F90 as it is being used in second call to zm deep convection scheme (convect_deep_tend_2)

  save

  ! Public methods
  public phys_register ! was initindx  - register physics methods
  public phys_init   ! Public initialization method
  public phys_run1   ! First phase of the public run method
  public phys_run2   ! Second phase of the public run method
  public phys_final  ! Public finalization method
  !
  ! Private module data
  !
  ! Physics package options
  character(len=16) :: shallow_scheme
  character(len=16) :: macrop_scheme
  character(len=16) :: microp_scheme 
  integer           :: cld_macmic_num_steps    ! Number of macro/micro substeps
  logical           :: do_clubb_sgs
  logical           :: do_shoc_sgs
  logical           :: use_subcol_microp   ! if true, use subcolumns in microphysics
  logical           :: state_debug_checks  ! Debug physics_state.
  logical           :: clim_modal_aero     ! climate controled by prognostic or prescribed modal aerosols
  logical           :: prog_modal_aero     ! Prognostic modal aerosols present
  logical           :: micro_do_icesupersat
  logical           :: pergro_test_active= .false.
  logical           :: pergro_mods = .false.
  logical           :: is_cmip6_volc !true if cmip6 style volcanic file is read otherwise false
  logical :: history_gaschmbudget ! output gas chemistry tracer concentrations and tendencies
  logical :: history_gaschmbudget_2D ! output 2D gas chemistry tracer concentrations and tendencies
  logical :: history_gaschmbudget_2D_levels ! output 2D gas chemistry tracer concentrations and tendencies within certain layers
  logical :: history_chemdyg_summary

  !======================================================================= 
contains

subroutine phys_register
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: Register constituents and physics buffer fields.
    ! 
    ! Author:    CSM Contact: M. Vertenstein, Aug. 1997
    !            B.A. Boville, Oct 2001
    !            A. Gettelman, Nov 2010 - put micro/macro physics into separate routines
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer,     only: pbuf_init_time
    use physics_buffer,     only: pbuf_add_field, dtype_r8, pbuf_register_subcol
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use spmd_utils,         only: masterproc
    use constituents,       only: pcnst, cnst_add, cnst_chk_dim, cnst_name

    use cam_control_mod,    only: moist_physics
    use chemistry,          only: chem_register
    use cloud_fraction,     only: cldfrc_register
    use stratiform,         only: stratiform_register
    use microp_driver,      only: microp_driver_register
    use microp_aero,        only: microp_aero_register
    use macrop_driver,      only: macrop_driver_register
    use clubb_intr,         only: clubb_register_cam
    use shoc_intr,          only: shoc_register_e3sm
    use conv_water,         only: conv_water_register
    use physconst,          only: mwdry, cpair, mwh2o, cpwv
    use tracers,            only: tracers_register
    use check_energy,       only: check_energy_register
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_register
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_register
    use ghg_data,           only: ghg_data_register
    use vertical_diffusion, only: vd_register
    use convect_deep,       only: convect_deep_register
    use convect_shallow,    only: convect_shallow_register
    use radiation,          only: radiation_register
    use co2_cycle,          only: co2_register
    use co2_diagnostics,    only: co2_diags_register
    use gw_drag,            only: gw_register
    use flux_avg,           only: flux_avg_register
    use iondrag,            only: iondrag_register
    use ionosphere,         only: ionos_register
    use string_utils,       only: to_lower
    use prescribed_ozone,   only: prescribed_ozone_register
    use prescribed_volcaero,only: prescribed_volcaero_register
    use prescribed_aero,    only: prescribed_aero_register
    use read_spa_data,      only: read_spa_data_register
    use prescribed_ghg,     only: prescribed_ghg_register
    use sslt_rebin,         only: sslt_rebin_register
    use aoa_tracers,        only: aoa_tracers_register
    use aircraft_emit,      only: aircraft_emit_register
    use cam_diagnostics,    only: diag_register
    use cloud_diagnostics,  only: cloud_diagnostics_register
    use cospsimulator_intr, only: cospsimulator_intr_register
    use rad_constituents,   only: rad_cnst_get_info ! Added to query if it is a modal aero sim or not
    use subcol,             only: subcol_register
    use subcol_utils,       only: is_subcol_on
    use output_aerocom_aie, only: output_aerocom_aie_register, do_aerocom_ind3
    use mo_chm_diags,       only: chm_diags_inti_ac

    !---------------------------Local variables-----------------------------
    !
    integer  :: m        ! loop index
    integer  :: mm       ! constituent index 
    !-----------------------------------------------------------------------

    integer :: nmodes
    character(len=16) :: spc_name

    call phys_getopts(shallow_scheme_out       = shallow_scheme, &
                      macrop_scheme_out        = macrop_scheme,   &
                      microp_scheme_out        = microp_scheme,   &
                      cld_macmic_num_steps_out = cld_macmic_num_steps, &
                      do_clubb_sgs_out         = do_clubb_sgs,     &
		      do_shoc_sgs_out          = do_shoc_sgs, &
                      do_aerocom_ind3_out      = do_aerocom_ind3,  &
                      use_subcol_microp_out    = use_subcol_microp, &
                      state_debug_checks_out   = state_debug_checks, &
                      micro_do_icesupersat_out = micro_do_icesupersat, &
                      pergro_test_active_out   = pergro_test_active, &
                      pergro_mods_out          = pergro_mods, &
                      history_gaschmbudget_out = history_gaschmbudget, &
                   history_gaschmbudget_2D_out = history_gaschmbudget_2D, &
            history_gaschmbudget_2D_levels_out = history_gaschmbudget_2D_levels, &
                   history_chemdyg_summary_out = history_chemdyg_summary )

    ! Initialize dyn_time_lvls
    call pbuf_init_time()

    ! Register the subcol scheme
    call subcol_register()

    ! Register water vapor.
    ! ***** N.B. ***** This must be the first call to cnst_add so that
    !                  water vapor is constituent 1.
    if (moist_physics) then
       call cnst_add('Q', mwh2o, cpwv, 1.E-12_r8, mm, &
            longname='Specific humidity', readiv=.true., is_convtran1=.true.)
    else
       call cnst_add('Q', mwh2o, cpwv, 0.0_r8, mm, &
            longname='Specific humidity', readiv=.false., is_convtran1=.true.)
    end if

    ! Fields for physics package diagnostics
    call pbuf_add_field('TINI',      'physpkg', dtype_r8, (/pcols,pver/), tini_idx)
    call pbuf_add_field('QINI',      'physpkg', dtype_r8, (/pcols,pver/), qini_idx)
    call pbuf_add_field('CLDLIQINI', 'physpkg', dtype_r8, (/pcols,pver/), cldliqini_idx)
    call pbuf_add_field('CLDICEINI', 'physpkg', dtype_r8, (/pcols,pver/), cldiceini_idx)
    call pbuf_add_field('static_ener_ac', 'global', dtype_r8, (/pcols/), static_ener_ac_idx)
    call pbuf_add_field('water_vap_ac',   'global', dtype_r8, (/pcols/), water_vap_ac_idx)


    ! check energy package
    call check_energy_register

    ! If using an ideal/adiabatic physics option, the CAM physics parameterizations 
    ! aren't called.
    if (moist_physics) then

       ! register fluxes for saving across time
       if (phys_do_flux_avg()) call flux_avg_register()

       call cldfrc_register()

       ! cloud water
       if( microp_scheme == 'RK' ) then
          call stratiform_register()
       elseif( microp_scheme == 'MG' .or. microp_scheme == 'P3' ) then
          if (.not. do_clubb_sgs .and. .not. do_shoc_sgs) call macrop_driver_register()
          call microp_aero_register()
          call microp_driver_register()
       end if
       
       ! Register CLUBB_SGS here
       if (do_clubb_sgs) call clubb_register_cam()
       
       ! Register SHOC_SGS here
       if (do_shoc_sgs) call shoc_register_e3sm()
       

       call pbuf_add_field('PREC_STR',  'physpkg',dtype_r8,(/pcols/),prec_str_idx)
       call pbuf_add_field('SNOW_STR',  'physpkg',dtype_r8,(/pcols/),snow_str_idx)
       call pbuf_add_field('PREC_PCW',  'physpkg',dtype_r8,(/pcols/),prec_pcw_idx)
       call pbuf_add_field('SNOW_PCW',  'physpkg',dtype_r8,(/pcols/),snow_pcw_idx)
       call pbuf_add_field('PREC_SED',  'physpkg',dtype_r8,(/pcols/),prec_sed_idx)
       call pbuf_add_field('SNOW_SED',  'physpkg',dtype_r8,(/pcols/),snow_sed_idx)
       if (is_subcol_on()) then
         call pbuf_register_subcol('PREC_STR', 'phys_register', prec_str_idx)
         call pbuf_register_subcol('SNOW_STR', 'phys_register', snow_str_idx)
         call pbuf_register_subcol('PREC_PCW', 'phys_register', prec_pcw_idx)
         call pbuf_register_subcol('SNOW_PCW', 'phys_register', snow_pcw_idx)
         call pbuf_register_subcol('PREC_SED', 'phys_register', prec_sed_idx)
         call pbuf_register_subcol('SNOW_SED', 'phys_register', snow_sed_idx)
       end if

    ! Who should add FRACIS? 
    ! -- It does not seem that aero_intr should add it since FRACIS is used in convection
    !     even if there are no prognostic aerosols ... so do it here for now 
       call pbuf_add_field('FRACIS','physpkg',dtype_r8,(/pcols,pver,pcnst/),m)

       call conv_water_register()
       
       ! Determine whether its a 'modal' aerosol simulation  or not
       call rad_cnst_get_info(0, nmodes=nmodes)
       clim_modal_aero = (nmodes > 0)

       if (clim_modal_aero) then
          call modal_aero_calcsize_reg()
          call modal_aero_wateruptake_reg()
       endif

       ! register chemical constituents including aerosols ...
       call chem_register(species_class)

       ! NB: has to be after chem_register to use tracer names
       ! Fields for gas chemistry tracers
       if (history_gaschmbudget .or. history_gaschmbudget_2D .or. history_gaschmbudget_2D_levels &
            .or. history_chemdyg_summary) then
         call chm_diags_inti_ac() ! to get aer_species
         do m = 1,gas_pcnst
            if (.not. any( aer_species == m )) then
              spc_name = trim(solsym(m))
              if (history_gaschmbudget .or. history_gaschmbudget_2D_levels .or. history_chemdyg_summary) then
                gas_ac_name(m) = 'ac_'//spc_name
                call pbuf_add_field(gas_ac_name(m), 'global', dtype_r8, (/pcols,pver/), gas_ac_idx)
              end if

              if (history_gaschmbudget_2D .or. history_chemdyg_summary) then
                gas_ac_name_2D(m) = 'ac_2D_'//spc_name
                call pbuf_add_field(gas_ac_name_2D(m), 'global', dtype_r8, (/pcols/), gas_ac_idx)
              end if
            end if
         enddo
       end if

       ! co2 constituents
       call co2_register()
       call co2_diags_register()

       call gw_register()

       ! register data model ozone with pbuf
       if (cam3_ozone_data_on) then
          call cam3_ozone_data_register()
       end if
       call prescribed_volcaero_register()
       call prescribed_ozone_register()
       call prescribed_aero_register()
       call read_spa_data_register()
       call prescribed_ghg_register()
       call sslt_rebin_register

       ! CAM3 prescribed aerosols
       if (cam3_aero_data_on) then
          call cam3_aero_data_register()
       end if

       ! register various data model gasses with pbuf
       call ghg_data_register()

       ! Register iondrag variables with pbuf
       call iondrag_register()

       ! Register ionosphere variables with pbuf if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call ionos_register()
       endif

       call aircraft_emit_register()

       ! deep convection
       call convect_deep_register

       !  shallow convection
       call convect_shallow_register

       ! radiation
       call radiation_register
       call cloud_diagnostics_register

       ! COSP
       call cospsimulator_intr_register

       ! vertical diffusion
       if (.not. do_clubb_sgs .and. .not. do_shoc_sgs) call vd_register()

       if (do_aerocom_ind3) call output_aerocom_aie_register()
    
    end if

    ! Register diagnostics PBUF
    call diag_register()

    ! Register age of air tracers
    call aoa_tracers_register()

    ! Register test tracers
    ! ***** N.B. ***** This is the last call to register constituents because
    !                  the test tracers fill the remaining available slots up
    !                  to constituent number PCNST -- regardless of what PCNST is set to.
    call tracers_register()

    ! All tracers registered, check that the dimensions are correct
    call cnst_chk_dim()

    ! ***NOTE*** No registering constituents after the call to cnst_chk_dim.

end subroutine phys_register



  !======================================================================= 

subroutine phys_inidat( cam_out, pbuf2d )
    use cam_abortutils, only : endrun

    use physics_buffer, only : pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_set_field, dyn_time_lvls


    use cam_initfiles,       only: initial_file_get_id, topo_file_get_id
    use cam_grid_support,    only: cam_grid_check, cam_grid_id
    use cam_grid_support,    only: cam_grid_get_dim_names
    use pio,                 only: file_desc_t
    use ncdio_atm,           only: infld
    use dycore,              only: dycore_is
    use polar_avg,           only: polar_average
    use short_lived_species, only: initialize_short_lived_species
    use comsrf,              only: sgh, sgh30
    use cam_control_mod,     only: aqua_planet

    type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    integer :: lchnk, m, n, i, k, ncol
    type(file_desc_t), pointer :: fh_ini, fh_topo
    character(len=8) :: fieldname
    real(r8), pointer :: cldptr(:,:,:,:), convptr_3d(:,:,:,:)
    real(r8), pointer :: tptr(:,:), tptr3d(:,:,:), tptr3d_2(:,:,:)
    real(r8), pointer :: qpert(:,:)

    character*11 :: subname='phys_inidat' ! subroutine name
    integer :: tpert_idx, qpert_idx, pblh_idx, vmag_gust_idx

    logical :: found=.false., found2=.false.
    integer :: ierr
    character(len=8) :: dim1name, dim2name
    integer                   :: grid_id  ! grid ID for data mapping
    nullify(tptr,tptr3d,tptr3d_2,cldptr,convptr_3d)

    fh_ini=>initial_file_get_id()

    !   dynamics variables are handled in dyn_init - here we read variables needed for physics 
    !   but not dynamics

    grid_id = cam_grid_id('physgrid')
    if (.not. cam_grid_check(grid_id)) then
      call endrun(trim(subname)//': Internal error, no "physgrid" grid')
    end if
    call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

    if(aqua_planet) then
       sgh = 0._r8
       sgh30 = 0._r8
       if (masterproc) write(iulog,*) 'AQUA_PLANET simulation, sgh, sgh30 initialized to 0.'
    else
       if (masterproc) write(iulog,*) 'NOT AN AQUA_PLANET simulation, initialize &
                                      &sgh, sgh30, land m using data from file.'
       fh_topo=>topo_file_get_id()
       call t_startf('phys_inidat_infld')
       call infld('SGH', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
            sgh, found, gridname='physgrid')
       call t_stopf('phys_inidat_infld')
       if(.not. found) call endrun('ERROR: SGH not found on topo file')

       call t_startf('phys_inidat_infld')
       call infld('SGH30', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
            sgh30, found, gridname='physgrid')
       call t_stopf('phys_inidat_infld')
       
       if(.not. found) then
          if (masterproc) write(iulog,*) 'Warning: Error reading SGH30 from topo file.'
          if (masterproc) write(iulog,*) 'The field SGH30 will be filled using data from SGH.'
          sgh30 = sgh
       end if
    end if

    allocate(tptr(1:pcols,begchunk:endchunk))

    call t_startf('phys_inidat_infld')
    call infld('PBLH', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, gridname='physgrid')
    call t_stopf('phys_inidat_infld')
    
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'PBLH initialized to 0.'
    end if
    pblh_idx = pbuf_get_index('pblh')

    call pbuf_set_field(pbuf2d, pblh_idx, tptr)

    call infld('TPERT', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, gridname='physgrid')
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'TPERT initialized to 0.'
    end if
    tpert_idx = pbuf_get_index( 'tpert')
    call pbuf_set_field(pbuf2d, tpert_idx, tptr)


    call infld('vmag_gust', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, gridname='physgrid')
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'vmag_gust initialized to 1.'
    end if
    vmag_gust_idx = pbuf_get_index( 'vmag_gust')
    call pbuf_set_field(pbuf2d, vmag_gust_idx, tptr)


    fieldname='QPERT'  
    qpert_idx = pbuf_get_index( 'qpert',ierr)
    if (qpert_idx > 0) then
       call infld(fieldname, fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
            tptr, found, gridname='physgrid')
       if(.not. found) then
          tptr=0_r8
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
       end if

       allocate(tptr3d_2(pcols,pcnst,begchunk:endchunk))
       tptr3d_2 = 0_r8
       tptr3d_2(:,1,:) = tptr(:,:)

       call pbuf_set_field(pbuf2d, qpert_idx, tptr3d_2)
       deallocate(tptr3d_2)
    end if

    fieldname='CUSH'
    m = pbuf_get_index('cush')
    call infld(fieldname, fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, &
         tptr, found, gridname='physgrid')
    if(.not.found) then
       if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 1000.'
       tptr=1000._r8
    end if
    do n=1,dyn_time_lvls
       call pbuf_set_field(pbuf2d, m, tptr, start=(/1,n/), kount=(/pcols,1/))
    end do
    deallocate(tptr)

    !
    ! 3-D fields
    !

    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname='CLOUD'
    m = pbuf_get_index('CLD')
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if(found) then
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname='QCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, gridname='physgrid')
       if(.not. found) then
          call infld('Q',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, gridname='physgrid')
          if (found) then
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with Q'
             if(dycore_is('LR')) call polar_average(pver, tptr3d) 	
          else
             call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
          end if
       end if
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    fieldname = 'ICCWAT'
    m = pbuf_get_index(fieldname, ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
          tptr3d, found, gridname='physgrid')
       if(found) then
          do n = 1, dyn_time_lvls
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          call infld('CLDICE',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
             tptr3d, found, gridname='physgrid')
          if(found) then
             do n = 1, dyn_time_lvls
                call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
             end do
          else
             call pbuf_set_field(pbuf2d, m, 0._r8)
          end if
          if (masterproc) then
             if (found) then
                write(iulog,*) trim(fieldname), ' initialized with CLDICE'
             else
                write(iulog,*) trim(fieldname), ' initialized to 0.0'
             end if
          end if
       end if
    end if

    fieldname = 'LCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, gridname='physgrid')
       if(found) then
          do n = 1, dyn_time_lvls
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          allocate(tptr3d_2(pcols,pver,begchunk:endchunk))     
          call infld('CLDICE',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, gridname='physgrid')
          call infld('CLDLIQ',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d_2, found2, gridname='physgrid')
          if(found .and. found2) then
             tptr3d(:,:,:)=tptr3d(:,:,:)+tptr3d_2(:,:,:)
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE + CLDLIQ'
          else if (found) then ! Data already loaded in tptr3d
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDICE only'
          else if (found2) then
             tptr3d(:,:,:)=tptr3d_2(:,:,:)
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with CLDLIQ only'
          end if

          if (found .or. found2) then
             do n = 1, dyn_time_lvls
                call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
             end do
             if(dycore_is('LR')) call polar_average(pver, tptr3d) 	
          else
             call pbuf_set_field(pbuf2d, m, 0._r8)
             if (masterproc)  write(iulog,*) trim(fieldname), ' initialized to 0.0'
          end if
          deallocate(tptr3d_2)
       end if
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname = 'TCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, gridname='physgrid')
       if(.not.found) then
          call infld('T', fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, gridname='physgrid')
          if(dycore_is('LR')) call polar_average(pver, tptr3d) 	
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized with T'
       end if
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pverp,begchunk:endchunk))

    fieldname = 'TKE'
    m = pbuf_get_index( 'tke')
    call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0.01_r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.01'
    end if


    fieldname = 'KVM'
    m = pbuf_get_index('kvm')
    call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if


    fieldname = 'KVH'
    m = pbuf_get_index('kvh')
    call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname = 'CONCLD'
    m = pbuf_get_index('CONCLD')
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, gridname='physgrid')
    if(found) then
       do n = 1, dyn_time_lvls
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    deallocate (tptr3d)

    call initialize_short_lived_species(fh_ini, pbuf2d)
end subroutine phys_inidat


subroutine phys_init( phys_state, phys_tend, pbuf2d, cam_out )

    !----------------------------------------------------------------------- 
    ! 
    ! Initialization of physics package.
    ! 
    !-----------------------------------------------------------------------

    use physics_buffer,     only: physics_buffer_desc, pbuf_initialize, pbuf_get_index
    use physics_buffer,     only: pbuf_get_chunk
    use time_manager,       only: get_nstep, is_first_step, is_first_restart_step
    use physconst,          only: rair, cpair, gravit, stebol, tmelt, &
                                  latvap, latice, rh2o, rhoh2o, pstd, zvir, &
                                  karman, rhodair, physconst_init 
    use ref_pres,           only: pref_edge, pref_mid
    use cloud_rad_props,    only: cloud_rad_props_init
    use cam_control_mod,    only: nsrest  ! restart flag
    use check_energy,       only: check_energy_init
    use chemistry,          only: chem_init
    use prescribed_ozone,   only: prescribed_ozone_init
    use prescribed_ghg,     only: prescribed_ghg_init
    use prescribed_aero,    only: prescribed_aero_init
    use read_spa_data,      only: read_spa_data_init
    use seasalt_model,      only: init_ocean_data, has_mam_mom
    use aerodep_flx,        only: aerodep_flx_init
    use aircraft_emit,      only: aircraft_emit_init
    use prescribed_volcaero,only: prescribed_volcaero_init
    use cloud_fraction,     only: cldfrc_init
    use cldfrc2m,           only: cldfrc2m_init
    use co2_cycle,          only: co2_init, co2_transport
    use convect_deep,       only: convect_deep_init
    use convect_shallow,    only: convect_shallow_init
    use cam_diagnostics,    only: diag_init
    use gw_drag,            only: gw_init
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_init
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_init
    use radheat,            only: radheat_init
    use radiation,          only: radiation_init
    use cloud_diagnostics,  only: cloud_diagnostics_init
    use co2_diagnostics,    only: co2_diags_init
    use stratiform,         only: stratiform_init
    use wv_saturation,      only: wv_sat_init
    use microp_driver,      only: microp_driver_init
    use microp_aero,        only: microp_aero_init
    use macrop_driver,      only: macrop_driver_init
    use conv_water,         only: conv_water_init
    use tracers,            only: tracers_init
    use aoa_tracers,        only: aoa_tracers_init
    use rayleigh_friction,  only: rayleigh_friction_init
    use pbl_utils,          only: pbl_utils_init
    use vertical_diffusion, only: vertical_diffusion_init
    use dycore,             only: dycore_is
    use phys_debug_util,    only: phys_debug_init
    use rad_constituents,   only: rad_cnst_init
    use aer_rad_props,      only: aer_rad_props_init
    use subcol,             only: subcol_init
    use qbo,                only: qbo_init
    use iondrag,            only: iondrag_init
#if ( defined OFFLINE_DYN )
    use metdata,            only: metdata_phys_init
#endif
    use ionosphere,	    only: ionos_init  ! Initialization of ionosphere module (WACCM-X)
    use majorsp_diffusion,  only: mspd_init   ! Initialization of major species diffusion module (WACCM-X)
    use clubb_intr,         only: clubb_ini_cam
    use shoc_intr,          only: shoc_init_e3sm
    use sslt_rebin,         only: sslt_rebin_init
    use tropopause,         only: tropopause_init
    use solar_data,         only: solar_data_init
    use rad_solar_var,      only: rad_solar_var_init
    use nudging,            only: Nudge_Model,nudging_init
    use output_aerocom_aie, only: output_aerocom_aie_init, do_aerocom_ind3
    use misc_diagnostics,   only: dcape_diags_init
    use conditional_diag_output_utils, only: cnd_diag_output_init
    use phys_grid_ctem,     only: phys_grid_ctem_init

    ! Input/output arguments
    type(physics_state), pointer       :: phys_state(:)
    type(physics_tend ), pointer       :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(physics_buffer_desc), pointer :: pbuf1d(:)

    type(cam_out_t),intent(inout)      :: cam_out(begchunk:endchunk)

    ! local variables
    integer :: lchnk
    real(r8) :: dp1 = huge(1.0_r8) !set in namelist, assigned in cloud_fraction.F90
    character(len=16)  :: deep_scheme      ! Default set in phys_control.F90

    !-----------------------------------------------------------------------

    call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, pcols)

    do lchnk = begchunk, endchunk
       call physics_state_set_grid(lchnk, phys_state(lchnk))
    end do

    !-------------------------------------------------------------------------------------------
    ! Initialize any variables in physconst which are not temporally and/or spatially constant
    !------------------------------------------------------------------------------------------- 
    call physconst_init()

    ! Initialize debugging a physics column
    call phys_debug_init()

    call pbuf_initialize(pbuf2d)

    !initialize physics update interface routine
    call physics_update_init()

    ! Initialize subcol scheme
    call subcol_init(pbuf2d)

    call setup_moist_indices()

    ! diag_init makes addfld calls for dynamics fields that are output from
    ! the physics decomposition
    call diag_init()

    call check_energy_init(phys_state)

    call tracers_init()

    ! age of air tracers
    call aoa_tracers_init()

    teout_idx = pbuf_get_index( 'TEOUT')

    ! "addfld" and "add_default" calls for CondiDiag
    call cnd_diag_output_init( pver, cnd_diag_info )

    ! For adiabatic or ideal physics don't need to initialize any of the
    ! parameterizations below:
    if (adiabatic .or. ideal_phys) return

    if (nsrest .eq. 0) then
       call t_startf ('phys_inidat')
       call phys_inidat(cam_out, pbuf2d) 
       call t_stopf ('phys_inidat')
    end if
    
    ! wv_saturation is relatively independent of everything else and
    ! low level, so init it early. Must at least do this before radiation.
    call wv_sat_init

    ! CAM3 prescribed aerosols
    if (cam3_aero_data_on) call cam3_aero_data_init(phys_state)

    ! Initialize rad constituents and their properties
    call t_startf ('rad_cnst_init')
    call rad_cnst_init()
    call t_stopf ('rad_cnst_init')
    
    call aer_rad_props_init()
    call cloud_rad_props_init()

    ! solar irradiance data modules
    call t_startf ('solar_data_init')
    call solar_data_init()
    call t_stopf ('solar_data_init')

    ! Prognostic chemistry.
    call t_startf ('chem_init')
    call chem_init(phys_state,pbuf2d, species_class)
    call t_stopf ('chem_init')

    ! Prescribed tracers
    call prescribed_ozone_init()
    call prescribed_ghg_init()
    call prescribed_aero_init()
    call read_spa_data_init()
    call aerodep_flx_init()
    call aircraft_emit_init(phys_state,pbuf2d)
    !when is_cmip6_volc is true ,cmip6 style volcanic file is read
    !Initialized to .false. here but it gets its values from prescribed_volcaero_init
    is_cmip6_volc = .false. 
    call t_startf ('prescribed_volcaero_init')
    call prescribed_volcaero_init(is_cmip6_volc)
    call t_stopf ('prescribed_volcaero_init')

    ! Initialize ocean data
    if (has_mam_mom) then
       call t_startf ('init_ocean_data')
       call init_ocean_data()
       call t_stopf ('init_ocean_data')
    end if

    ! co2 cycle            
    if (co2_transport()) then
       call co2_init()
    end if
    call co2_diags_init(phys_state)

    ! CAM3 prescribed ozone
    if (cam3_ozone_data_on) call cam3_ozone_data_init(phys_state)

    call gw_init(pbuf2d)

    call rayleigh_friction_init()

    call pbl_utils_init(gravit, karman, cpair, rair, zvir)
    if (.not. do_clubb_sgs .and. .not. do_shoc_sgs) call vertical_diffusion_init(pbuf2d)

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call mspd_init ()
       ! Initialization of ionosphere module if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call ionos_init()
       endif
    endif

    call tsinti(tmelt, latvap, rair, stebol, latice)

    call t_startf ('radiation_init')
    call radiation_init(phys_state,pbuf2d)
    call t_stopf ('radiation_init')

    call rad_solar_var_init()

    call cloud_diagnostics_init()

    call radheat_init(pref_mid)

    call convect_shallow_init(pref_edge, pbuf2d)

    call cldfrc_init(dp1)! for passing dp1 on to clubb
    call cldfrc2m_init()

    call convect_deep_init(pref_edge, pbuf2d)

    ! If this is the first timestep of an initial run,
    ! set some old-timestep values for dCAPE diagnostics.
    ! In restart/branch/hybrid runs, these values are read from
    ! the restart file.

    call phys_getopts( deep_scheme_out   = deep_scheme )

    if (is_first_step() .and. deep_scheme .eq. 'ZM') then
       do lchnk = begchunk, endchunk
          pbuf1d => pbuf_get_chunk(pbuf2d, lchnk)
          call dcape_diags_init( pbuf1d, pver )
       end do
    end if
    !----

    if( microp_scheme == 'RK' ) then
       call stratiform_init()
    elseif( microp_scheme == 'MG' .or. microp_scheme == 'P3' ) then 
       if (.not. do_clubb_sgs .and. .not. do_shoc_sgs) call macrop_driver_init(pbuf2d)
       call microp_aero_init()
       call microp_driver_init(pbuf2d)
       call conv_water_init(pbuf2d)
    end if

    ! initiate CLUBB within CAM
    if (do_clubb_sgs) call clubb_ini_cam(pbuf2d,dp1)
    
    ! initiate SHOC within E3SM
    if (do_shoc_sgs) call shoc_init_e3sm(pbuf2d,dp1)

    call qbo_init

    call iondrag_init(pref_mid)

#if ( defined OFFLINE_DYN )
    call metdata_phys_init()
#endif
    call sslt_rebin_init()
    
    call t_startf ('tropopause_init')
    call tropopause_init()
    call t_stopf ('tropopause_init')

    if(do_aerocom_ind3) call output_aerocom_aie_init()

    prec_dp_idx  = pbuf_get_index('PREC_DP')
    snow_dp_idx  = pbuf_get_index('SNOW_DP')
    prec_sh_idx  = pbuf_get_index('PREC_SH')
    snow_sh_idx  = pbuf_get_index('SNOW_SH')

    call phys_getopts(prog_modal_aero_out=prog_modal_aero)

    if (clim_modal_aero) then

       ! If climate calculations are affected by prescribed modal aerosols, the
       ! the initialization routine for the dry mode radius calculation is called
       ! here.  For prognostic MAM the initialization is called from
       ! modal_aero_initialize
       if (.not. prog_modal_aero) then
          call modal_aero_calcsize_init(pbuf2d, species_class)
       endif

       call modal_aero_wateruptake_init(pbuf2d)

    end if

    ! Initialize Nudging Parameters
    !--------------------------------
    if(Nudge_Model) call nudging_init

    ! Initialize Transformed Eularian Mean (TEM) diagnostics
    call phys_grid_ctem_init()

    
   !BSINGH -  addfld and adddefault calls for perturb growth testing    
    if(pergro_test_active)call add_fld_default_calls()

    if(is_first_step().or.is_first_restart_step())then
      nstep_ignore_diagn1 = get_nstep()
      nstep_ignore_diagn2 = nstep_ignore_diagn1 + 1
    endif


end subroutine phys_init

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out, phys_diag)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! First part of atmospheric physics package before updating of surface models
    ! 
    !-----------------------------------------------------------------------
    use time_manager,   only: get_nstep
    use cam_diagnostics,only: diag_allocate, diag_physvar_ic

    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, pbuf_allocate
#if (defined E3SM_SCM_REPLAY )
    use cam_history,    only: outfld
#endif
    use comsrf,         only: fsns, fsnt, flns, sgh, sgh30, flnt, fsds
    use cam_abortutils,     only: endrun
#if ( defined OFFLINE_DYN )
     use metdata,       only: get_met_srf1
#endif

    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
    type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
    type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
    type(cnd_diag_t),    intent(inout), dimension(begchunk:endchunk) :: phys_diag
    !-----------------------------------------------------------------------
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! indices
    integer :: ncol                              ! number of columns
    integer :: nstep                             ! current timestep number
#if (! defined SPMD)
    integer :: mpicom = 0
#endif
    integer(i8) :: beg_proc_cnt                  ! start time for the process
    integer(i8) :: end_proc_cnt                  ! stop time for the process
    integer(i8) :: beg_chnk_cnt                  ! start time for a chunk
    integer(i8) :: end_chnk_cnt                  ! stop time for a chunk
    integer(i8) :: sysclock_rate                 ! system clock rate
    integer(i8) :: sysclock_max                  ! system clock max value
    real(r8)    :: chunk_cost                    ! measured cost per chunk
    type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

    logical     :: print_additional_diagn_local

    call t_startf ('physpkg_st1')
    nstep = get_nstep()

    print_additional_diagn_local = .true.
    if(nstep == nstep_ignore_diagn1 .or. nstep == nstep_ignore_diagn2) print_additional_diagn_local = .false.
    if(.not. print_additional_diagn) print_additional_diagn_local = .false.

#if ( defined OFFLINE_DYN )
    !
    ! if offline mode set SNOWH and TS for micro-phys
    !
    call get_met_srf1( cam_in )
#endif

    ! The following initialization depends on the import state (cam_in)
    ! being initialized.  This isn't true when cam_init is called, so need
    ! to postpone this initialization to here.
    if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in,  pbuf2d)

    ! Compute total energy of input state and previous output state
    call t_startf ('chk_en_gmean')
    call check_energy_gmean(phys_state, pbuf2d, ztodt, nstep)
    call t_stopf ('chk_en_gmean')

    
    if ( adiabatic .or. ideal_phys )then
       call t_stopf ('physpkg_st1')

       call t_startf ('bc_physics')
       call phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend,  pbuf2d)
       call t_stopf ('bc_physics')
    else
       call pbuf_allocate(pbuf2d, 'physpkg')
       call diag_allocate()

       !-----------------------------------------------------------------------
       ! Advance time information
       !-----------------------------------------------------------------------

       call t_startf('phys_timestep_init')
       call phys_timestep_init( phys_state, cam_out, pbuf2d)
       call t_stopf('phys_timestep_init')

       call t_stopf ('physpkg_st1')

#ifdef TRACER_CHECK
       call gmean_mass ('before tphysbc DRY', phys_state)
#endif


       !-----------------------------------------------------------------------
       ! Tendency physics before flux coupler invocation
       !-----------------------------------------------------------------------
       !

#if (defined BFB_CAM_SCAM_IOP )
       do c=begchunk, endchunk
          call outfld('Tg',cam_in(c)%ts,pcols   ,c     )
       end do
#endif
       call t_barrierf('sync_bc_physics', mpicom)
       call t_startf ('bc_physics')
       !call t_adj_detailf(+1)

       call system_clock(count=beg_proc_cnt)
       
!$OMP PARALLEL DO SCHEDULE(STATIC,1) &
!$OMP PRIVATE (c, beg_chnk_cnt, phys_buffer_chunk, end_chnk_cnt, sysclock_rate, sysclock_max, chunk_cost)
       do c=begchunk, endchunk

          call system_clock(count=beg_chnk_cnt)

          !
          ! Output physics terms to IC file
          !
          phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

          call t_startf ('diag_physvar_ic')
          call diag_physvar_ic ( c,  phys_buffer_chunk, cam_out(c), cam_in(c) )
          call t_stopf ('diag_physvar_ic')

          call tphysbc (ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), phys_state(c),        &
                       phys_tend(c), phys_buffer_chunk, phys_diag(c), fsds(1,c),                 &
                       sgh(1,c), sgh30(1,c), cam_out(c), cam_in(c) )

          call system_clock(count=end_chnk_cnt, count_rate=sysclock_rate, count_max=sysclock_max)
          if ( end_chnk_cnt < beg_chnk_cnt ) end_chnk_cnt = end_chnk_cnt + sysclock_max
          chunk_cost = real( (end_chnk_cnt-beg_chnk_cnt), r8)/real(sysclock_rate, r8)
          call update_cost_p(c, chunk_cost)

       end do

       call system_clock(count=end_proc_cnt, count_rate=sysclock_rate, count_max=sysclock_max)
       if ( end_proc_cnt < beg_proc_cnt ) end_proc_cnt = end_proc_cnt + sysclock_max
       phys_proc_cost = phys_proc_cost + real( (end_proc_cnt-beg_proc_cnt), r8)/real(sysclock_rate, r8)

       !call t_adj_detailf(-1)
       call t_stopf ('bc_physics')

#ifdef TRACER_CHECK
       call gmean_mass ('between DRY', phys_state)
#endif

       if( print_additional_diagn_local ) then
          call t_startf ('chk_en_gmean_additional')
          call check_energy_gmean_additional_diagn(phys_state, pbuf2d, ztodt, nstep)
          call t_stopf ('chk_en_gmean_additional')
       endif

    end if !not adiabatic or ideal

end subroutine phys_run1

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend,  pbuf2d)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Physics for adiabatic or idealized physics case.
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only : physics_buffer_desc, pbuf_set_field, pbuf_get_chunk, pbuf_old_tim_idx
    use time_manager,     only: get_nstep
    use cam_diagnostics,  only: diag_phys_writeout
    use check_energy,     only: check_energy_fix, check_energy_chng
    use dycore,           only: dycore_is

    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    !-----------------------------------------------------------------------
    !---------------------------Local workspace-----------------------------
    !
    integer             :: c               ! indices
    integer             :: nstep           ! current timestep number
    type(physics_ptend) :: ptend(begchunk:endchunk) ! indivdual parameterization tendencies
    real(r8)            :: flx_heat(pcols) ! effective sensible heat flux
    real(r8)            :: zero(pcols)     ! array of zeros

    integer(i8)         :: beg_proc_cnt    ! start time for the process
    integer(i8)         :: end_proc_cnt    ! stop time for the process
    integer(i8)         :: beg_chnk_cnt    ! start time for a chunk
    integer(i8)         :: end_chnk_cnt    ! stop time for a chunk
    integer(i8)         :: sysclock_rate   ! system clock rate
    integer(i8)         :: sysclock_max    ! system clock max value
    real(r8)            :: chunk_cost      ! measured cost per chunk

    character(len=128)  :: ideal_phys_option

    ! physics buffer field for total energy
    real(r8), pointer, dimension(:) :: teout
    logical, SAVE :: first_exec_of_phys_run1_adiabatic_or_ideal  = .TRUE.
    !-----------------------------------------------------------------------

    nstep = get_nstep()
    zero  = 0._r8

    ! Associate pointers with physics buffer fields
    if (first_exec_of_phys_run1_adiabatic_or_ideal) then
       first_exec_of_phys_run1_adiabatic_or_ideal  = .FALSE.
    endif

    call system_clock(count=beg_proc_cnt)

    call phys_getopts(ideal_phys_option_out=ideal_phys_option)

!$OMP PARALLEL DO SCHEDULE(STATIC,1) &
!$OMP PRIVATE (c, beg_chnk_cnt, flx_heat, end_chnk_cnt, sysclock_rate, sysclock_max, chunk_cost)
    do c=begchunk, endchunk

       call system_clock(count=beg_chnk_cnt)

       ! Initialize the physics tendencies to zero.
       call physics_tend_init(phys_tend(c))

       ! Dump dynamics variables to history buffers
       call diag_phys_writeout(phys_state(c))

       if (dycore_is('LR') .or. dycore_is('SE') ) then
          call check_energy_fix(phys_state(c), ptend(c), nstep, flx_heat)
          call physics_update(phys_state(c), ptend(c), ztodt, phys_tend(c))
          call check_energy_chng(phys_state(c), phys_tend(c), "chkengyfix", nstep, ztodt, &
               zero, zero, zero, flx_heat)
          call physics_ptend_dealloc(ptend(c))
       end if

       if ( ideal_phys )then
          call t_startf('tphysidl')
          call tphysidl(ztodt, phys_state(c), phys_tend(c), trim(ideal_phys_option))
          call t_stopf('tphysidl')
       end if

       ! Save total enery after physics for energy conservation checks
       call pbuf_set_field(pbuf_get_chunk(pbuf2d, c), teout_idx, phys_state(c)%te_cur)

       call system_clock(count=end_chnk_cnt, count_rate=sysclock_rate, count_max=sysclock_max)
       if ( end_chnk_cnt < beg_chnk_cnt ) end_chnk_cnt = end_chnk_cnt + sysclock_max
       chunk_cost = real( (end_chnk_cnt-beg_chnk_cnt), r8)/real(sysclock_rate, r8)
       call update_cost_p(c, chunk_cost)

    end do

    call system_clock(count=end_proc_cnt, count_rate=sysclock_rate, count_max=sysclock_max)
    if ( end_proc_cnt < beg_proc_cnt ) end_proc_cnt = end_proc_cnt + sysclock_max
    phys_proc_cost = phys_proc_cost + real( (end_proc_cnt-beg_proc_cnt), r8)/real(sysclock_rate, r8)

end subroutine phys_run1_adiabatic_or_ideal

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run2(phys_state, ztodt, phys_tend, pbuf2d,  cam_out, &
       cam_in, phys_diag )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Second part of atmospheric physics package after updating of surface models
    ! 
    ! Modified by Kai Zhang 2017-03: add IEFLX fixer treatment 
    !-----------------------------------------------------------------------
    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, pbuf_deallocate, pbuf_update_tim_idx
    use mo_lightning,   only: lightning_no_prod


    use cam_diagnostics,only: diag_deallocate, diag_surf
    use comsrf,         only: trefmxav, trefmnav, sgh, sgh30, fsds
    use physconst,      only: stebol, latvap
#if ( defined OFFLINE_DYN )
    use metdata,        only: get_met_srf2
#endif
    use time_manager,   only: get_nstep, is_first_step, is_end_curr_month, &
                              is_first_restart_step, is_last_step
    use check_energy,   only: ieflx_gmean, check_ieflx_fix 
    use phys_control,   only: ieflx_opt
    use co2_diagnostics,only: get_total_carbon, print_global_carbon_diags, &
                              co2_diags_store_fields, co2_diags_read_fields
    use co2_cycle,      only: co2_transport
    !
    ! Input arguments
    !
    real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
    !
    ! Input/Output arguments
    !
    type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
    type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
    type(physics_buffer_desc),pointer, dimension(:,:)     :: pbuf2d

    type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
    type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
    type(cnd_diag_t),    intent(inout), dimension(begchunk:endchunk) :: phys_diag

    !
    !-----------------------------------------------------------------------
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! chunk index
    integer :: ncol                              ! number of columns
    integer :: nstep                             ! current timestep number
#if (! defined SPMD)
    integer :: mpicom = 0
#endif
    integer(i8) :: beg_proc_cnt                  ! start time for the process
    integer(i8) :: end_proc_cnt                  ! stop time for the process
    integer(i8) :: beg_chnk_cnt                  ! start time for a chunk
    integer(i8) :: end_chnk_cnt                  ! stop time for a chunk
    integer(i8) :: sysclock_rate                 ! system clock rate
    integer(i8) :: sysclock_max                  ! system clock max value
    real(r8)    :: chunk_cost                    ! measured cost per chunk
    type(physics_buffer_desc),pointer, dimension(:)     :: phys_buffer_chunk
    !
    ! If exit condition just return
    !

    if ( adiabatic .or. ideal_phys ) return
    !-----------------------------------------------------------------------
    ! Tendency physics after coupler 
    ! Not necessary at terminal timestep.
    !-----------------------------------------------------------------------
    !
#if ( defined OFFLINE_DYN )
    !
    ! if offline mode set SHFLX QFLX TAUX TAUY for vert diffusion
    !
    call get_met_srf2( cam_in )
#endif
    ! Set lightning production of NO
    call t_startf ('lightning_no_prod')
    call lightning_no_prod( phys_state, pbuf2d,  cam_in )
    call t_stopf ('lightning_no_prod')

    call t_barrierf('sync_ac_physics', mpicom)
    call t_startf ('ac_physics')
    !call t_adj_detailf(+1)

    nstep = get_nstep()


    !! calculate the global mean ieflx 

    if(ieflx_opt>0) then
       call ieflx_gmean(phys_state, phys_tend, pbuf2d, cam_in, cam_out, nstep)
    end if

    ! Get carbon conservation fields from pbuf if restarting
    if ( co2_transport() ) then
       if ( is_first_restart_step() ) then
          call co2_diags_read_fields(phys_state, pbuf2d)
       end if
    end if

    call system_clock(count=beg_proc_cnt)

!$OMP PARALLEL DO SCHEDULE(STATIC,1) &
!$OMP PRIVATE (c, beg_chnk_cnt, ncol, phys_buffer_chunk, end_chnk_cnt, sysclock_rate, sysclock_max, chunk_cost)
    do c=begchunk,endchunk

       call system_clock(count=beg_chnk_cnt)

       ncol = get_ncols_p(c)
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

       !! 
       !! add the implied internal energy flux to sensible heat flux
       !! 

       if(ieflx_opt>0) then
          call check_ieflx_fix(c, ncol, nstep, cam_in(c)%shf)
       end if

       !
       ! surface diagnostics for history files
       !
       call t_startf('diag_surf')
       call diag_surf(cam_in(c), cam_out(c), phys_state(c)%ps,trefmxav(1,c), trefmnav(1,c))
       call t_stopf('diag_surf')

       call tphysac(ztodt, cam_in(c),  &
            sgh(1,c), sgh30(1,c), cam_out(c),                              &
            phys_state(c), phys_tend(c), phys_buffer_chunk, phys_diag(c),  &
            fsds(1,c))

       call system_clock(count=end_chnk_cnt, count_rate=sysclock_rate, count_max=sysclock_max)
       if ( end_chnk_cnt < beg_chnk_cnt ) end_chnk_cnt = end_chnk_cnt + sysclock_max
       chunk_cost = real( (end_chnk_cnt-beg_chnk_cnt), r8)/real(sysclock_rate, r8)
       call update_cost_p(c, chunk_cost)

    end do                    ! Chunk loop

    call system_clock(count=end_proc_cnt, count_rate=sysclock_rate, count_max=sysclock_max)
    if ( end_proc_cnt < beg_proc_cnt ) end_proc_cnt = end_proc_cnt + sysclock_max
    phys_proc_cost = phys_proc_cost + real( (end_proc_cnt-beg_proc_cnt), r8)/real(sysclock_rate, r8)

    !call t_adj_detailf(-1)
    call t_stopf('ac_physics')

#ifdef TRACER_CHECK
    call gmean_mass ('after tphysac FV:WET)', phys_state)
#endif

    !
    ! Check for carbon conservation
    !
    if ( co2_transport() ) then
       do c = begchunk, endchunk
          call get_total_carbon(phys_state(c), 'wet')
       end do
       call print_global_carbon_diags(phys_state, ztodt, nstep)
       do c = begchunk, endchunk
          ncol = get_ncols_p(c)
          phys_state(c)%tc_prev(:ncol) = phys_state(c)%tc_curr(:ncol)
          if ( is_first_step() ) then
             phys_state(c)%tc_init(:ncol) = phys_state(c)%tc_curr(:ncol)
          end if
          if ( is_end_curr_month() ) then
             phys_state(c)%tc_mnst(:ncol) = phys_state(c)%tc_curr(:ncol)
          end if
       end do
       call co2_diags_store_fields(phys_state, pbuf2d)
    end if

    call t_startf ('physpkg_st2')
    call pbuf_deallocate(pbuf2d, 'physpkg')

    call pbuf_update_tim_idx()
    call diag_deallocate()
    call t_stopf ('physpkg_st2')

end subroutine phys_run2

  !
  !----------------------------------------------------------------------- 
  !

subroutine phys_final( phys_state, phys_tend, pbuf2d, phys_diag )
    use physics_buffer, only : physics_buffer_desc, pbuf_deallocate
    use chemistry, only : chem_final
    use wv_saturation, only : wv_sat_final
    use radiation, only: radiation_final
    use phys_grid_ctem, only : phys_grid_ctem_final
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Finalization of physics package
    ! 
    !-----------------------------------------------------------------------
    ! Input/output arguments
    type(physics_state), pointer :: phys_state(:)
    type(physics_tend ), pointer :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)
    type(cnd_diag_t),          pointer :: phys_diag(:)

    if(associated(pbuf2d)) then
       call pbuf_deallocate(pbuf2d,'global')
       deallocate(pbuf2d)
    end if
    deallocate(phys_state)
    deallocate(phys_tend)
    deallocate(phys_diag)

    call t_startf ('chem_final')
    call chem_final
    call t_stopf ('chem_final')

    call t_startf ('wv_sat_final')
    call wv_sat_final
    call t_stopf ('wv_sat_final')

    call t_startf ('radiation_final')
    call radiation_final()
    call t_stopf ('radiation_final')

    call t_startf ('print_cost_p')
    call print_cost_p
    call t_stopf ('print_cost_p')

    call t_startf ('phys_grid_ctem_final')
    call phys_grid_ctem_final()
    call t_stopf ('phys_grid_ctem_final')

end subroutine phys_final


subroutine tphysac (ztodt,   cam_in,  &
       sgh,     sgh30,                                     &
       cam_out,  state,   tend,    pbuf,   diag,           &
       fsds    )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Tendency physics after coupling to land, sea, and ice models.
    ! Computes the following:
    !   o Radon surface flux and decay (optional)
    !   o Vertical diffusion and planetary boundary layer
    !   o Multiple gravity wave drag
    ! 
    ! Method: 
    ! <Describe the algorithm(s) used in the routine.> 
    ! <Also include any applicable external references.> 
    ! 
    ! Author: CCM1, CMS Contact: J. Truesdale
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer, only: physics_buffer_desc, pbuf_set_field, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use chemistry,          only: chem_is_active, chem_timestep_tend, chem_emissions
    use cam_diagnostics,    only: diag_phys_tend_writeout
    use gw_drag,            only: gw_tend
    use vertical_diffusion, only: vertical_diffusion_tend
    use rayleigh_friction,  only: rayleigh_friction_tend
    use physics_types,      only: physics_state, physics_tend, physics_ptend,    &
         physics_dme_adjust, set_dry_to_wet, physics_state_check
    use majorsp_diffusion,  only: mspd_intr  ! WACCM-X major diffusion
    use ionosphere,         only: ionos_intr ! WACCM-X ionosphere
    use tracers,            only: tracers_timestep_tend
    use aoa_tracers,        only: aoa_tracers_timestep_tend
    use physconst,          only: rhoh2o, latvap,latice, rga, gravit
    use aero_model,         only: aero_model_drydep
    use check_energy,       only: check_energy_chng, check_water, & 
                                  check_prect, check_qflx , &
                                  check_tracers_data, check_tracers_init, &
                                  check_tracers_chng, check_tracers_fini
    use time_manager,       only: get_nstep, is_first_step, is_end_curr_month
    use cam_abortutils,         only: endrun
    use dycore,             only: dycore_is
    use cam_control_mod,    only: aqua_planet 
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_set
    use charge_neutrality,  only: charge_fix
    use qbo,                only: qbo_relax
    use iondrag,            only: iondrag_calc, do_waccm_ions
    use clubb_intr,         only: clubb_surface
    use cflx,               only: cflx_tend
    use perf_mod
    use flux_avg,           only: flux_avg_run
    use nudging,            only: Nudge_Model,Nudge_ON,nudging_timestep_tend
    use phys_control,       only: use_qqflx_fixer
    use co2_cycle,          only: co2_cycle_set_ptend, co2_transport
    use co2_diagnostics,    only: get_carbon_sfc_fluxes, get_carbon_air_fluxes

    implicit none

    !
    ! Arguments
    !
    real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
    real(r8), intent(in) :: fsds(pcols)            ! down solar flux
    real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
    real(r8), intent(in) :: sgh30(pcols)           ! Std. deviation of 30s orography for tms

    type(cam_in_t),      intent(inout) :: cam_in
    type(cam_out_t),     intent(inout) :: cam_out
    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(cnd_diag_t),    intent(inout) :: diag

    !
    !---------------------------Local workspace-----------------------------
    !
    type(check_tracers_data):: tracerint           ! tracer mass integrals and cummulative boundary fluxes
    type(physics_ptend)     :: ptend               ! indivdual parameterization tendencies

    integer  :: nstep                              ! current timestep number
    real(r8) :: zero(pcols)                        ! array of zeros

    integer :: lchnk                                ! chunk identifier
    integer :: ncol                                 ! number of atmospheric columns
    integer i,k,m                 ! Longitude, level indices
    integer :: yr, mon, day, tod       ! components of a date

    logical :: labort                            ! abort flag

    real(r8) tvm(pcols,pver)           ! virtual temperature
    real(r8) prect(pcols)              ! total precipitation
    real(r8) surfric(pcols)            ! surface friction velocity
    real(r8) obklen(pcols)             ! Obukhov length
    real(r8) :: fh2o(pcols)            ! h2o flux to balance source from methane chemistry
    real(r8) :: tmp_q     (pcols,pver) ! tmp space
    real(r8) :: tmp_cldliq(pcols,pver) ! tmp space
    real(r8) :: tmp_cldice(pcols,pver) ! tmp space
    real(r8) :: tmp_t     (pcols,pver) ! tmp space
    real(r8) :: ftem      (pcols,pver) ! tmp space
    real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
    real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor

    ! physics buffer fields for total energy and mass adjustment
    integer itim_old, ifld

    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: cld
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore
    real(r8), pointer, dimension(:,:) :: ast     ! relative humidity cloud fraction 

    logical :: do_clubb_sgs 

    !DCAPE-ULL: physics buffer fields to compute tendencies for dcape
    real(r8), pointer, dimension(:,:) :: t_star   ! temperature
    real(r8), pointer, dimension(:,:) :: q_star   ! moisture

    character(len=16)  :: deep_scheme             ! Default set in phys_control.F90

    ! Debug physics_state.
    logical :: state_debug_checks

    logical :: l_tracer_aero
    logical :: l_vdiff
    logical :: l_rayleigh
    logical :: l_gw_drag
    logical :: l_ac_energy_chk

    real(r8) :: cc,mm,qq,pp,adjust

    ! Numerical schemes for process coupling
    integer :: cflx_cpl_opt  ! When to apply surface tracer fluxes  (not including water vapor).
                             ! The default for aerosols is to do this 
                             ! after tphysac:clubb_surface and before aerosol dry removal.
                             ! For chemical gases, different versions of EAM 
                             ! might use different process ordering.
    !
    !-----------------------------------------------------------------------
    !
    lchnk = state%lchnk
    ncol  = state%ncol

    nstep = get_nstep()
    
    call phys_getopts( do_clubb_sgs_out       = do_clubb_sgs, &
                       do_shoc_sgs_out        = do_shoc_sgs, &
                       state_debug_checks_out = state_debug_checks &
                      ,deep_scheme_out        = deep_scheme        &
                      ,cflx_cpl_opt_out       = cflx_cpl_opt       &
                      ,l_tracer_aero_out      = l_tracer_aero      &
                      ,l_vdiff_out            = l_vdiff            &
                      ,l_rayleigh_out         = l_rayleigh         &
                      ,l_gw_drag_out          = l_gw_drag          &
                      ,l_ac_energy_chk_out    = l_ac_energy_chk    &
                     )

    ! Adjust the surface fluxes to reduce instabilities in near sfc layer
    if (phys_do_flux_avg()) then 
       call flux_avg_run(state, cam_in,  pbuf, nstep, ztodt)
    endif

    ! Validate the physics state.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysac")

    ! CondiDiag checkpoint before the after-coupling parameterization suite
    call cnd_diag_checkpoint( diag, 'MCTCPL', state, pbuf, cam_in, cam_out )

    call t_startf('tphysac_init')
    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()


    ifld = pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, start=(/1,1,itim_old/),kount=(/pcols,pver,1/))

    ifld = pbuf_get_index('AST')
    call pbuf_get_field(pbuf, ifld, ast, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    !
    ! accumulate fluxes into net flux array for spectral dycores
    ! jrm Include latent heat of fusion for snow
    !
    do i=1,ncol
       tend%flx_net(i) = tend%flx_net(i) + cam_in%shf(i) + (cam_out%precc(i) &
            + cam_out%precl(i))*latvap*rhoh2o &
            + (cam_out%precsc(i) + cam_out%precsl(i))*latice*rhoh2o
    end do

if (l_tracer_aero) then

    ! emissions of aerosols and gas-phase chemistry constituents at surface
    call chem_emissions( state, cam_in )

end if ! l_tracer_aero

    call cnd_diag_checkpoint( diag, 'CHEMEMIS', state, pbuf, cam_in, cam_out )

    ! get nstep and zero array for energy checker
    zero = 0._r8
    nstep = get_nstep()
    call check_tracers_init(state, tracerint)

    if(print_additional_diagn) then
       !save TE before physstep
       call check_energy_chng(state, tend, "assign_te_before_physstep", nstep, ztodt, zero, zero, zero, zero)
       call additional_diagn_before_step_part1(state, cam_in)
    endif

!!== KZ_WCON

    !just output
    call check_qflx(state, tend, "PHYAC01", nstep, ztodt, cam_in%cflx(:,1))

    if(.not.use_qqflx_fixer) then 
       ! Check if latent heat flux exceeds the total moisture content of the
       ! lowest model layer, thereby creating negative moisture.
       call qneg4('TPHYSAC '       ,lchnk               ,ncol  ,ztodt ,               &
            state%q(1,pver,1),state%rpdel(1,pver) ,cam_in%shf ,         &
            cam_in%lhf , cam_in%cflx )
    end if 

    if(print_additional_diagn) then
       call additional_diagn_before_step_part2(state, cam_in)
    endif

    !just output
    call check_qflx(state, tend, "PHYAC02", nstep, ztodt, cam_in%cflx(:,1))

!!== KZ_WCON

    call t_stopf('tphysac_init')
    call cnd_diag_checkpoint( diag, 'PACINI', state, pbuf, cam_in, cam_out )

if (l_tracer_aero) then
    !===================================================
    ! Source/sink terms for advected tracers.
    !===================================================
    call t_startf('adv_tracer_src_snk')
    ! Test tracers

    call tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "tracers_timestep_tend", nstep, ztodt,   &
         cam_in%cflx)

    call aoa_tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "aoa_tracers_timestep_tend", nstep, ztodt,   &
         cam_in%cflx)
    
    ! add tendency from aircraft emissions
    call co2_cycle_set_ptend(state, pbuf, ptend)
    call physics_update(state, ptend, ztodt, tend)
    call get_carbon_air_fluxes(state, pbuf, ztodt)

    ! Chemistry calculation
    if (chem_is_active()) then
       call chem_timestep_tend(state, ptend, cam_in, cam_out, ztodt, &
            pbuf,  fh2o, fsds)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chem", nstep, ztodt, fh2o, zero, zero, zero)
       call check_tracers_chng(state, tracerint, "chem_timestep_tend", nstep, ztodt, &
            cam_in%cflx)
    end if
    call t_stopf('adv_tracer_src_snk')

end if ! l_tracer_aero

    call cnd_diag_checkpoint( diag, 'CHEM', state, pbuf, cam_in, cam_out )

    !===================================================
    ! Vertical diffusion/pbl calculation
    ! Call vertical diffusion code (pbl, free atmosphere and molecular)
    !===================================================

    ! If CLUBB (or SHOC) is called, do not call vertical diffusion, but obukov length and
    !   surface friction velocity still need to be computed.  In addition, 
    !   surface fluxes need to be updated here for constituents 
    if (do_clubb_sgs .or. do_shoc_sgs) then

       ! If CLUBB is called, do not call vertical diffusion, but still
       ! calculate surface friction velocity (ustar) and Obukhov length
       call clubb_surface ( state, cam_in, pbuf, surfric, obklen)

       ! Diagnose tracer mixing ratio tendencies from surface fluxes, 
       ! then update the mixing ratios. (If cflx_cpl_opt==2, these are done in 
       ! tphysbc after deep convection before the cloud mac-mic subcycles
       ! so that the emission-induced updates of state can be closer to turbulent transport.)
       ! Note that the two subroutine calls below do not touch water vapor. 
       ! They also have no effects on tracers for which cam_in%cflx(:,m) 
       ! is zero at this point.

       if (cflx_cpl_opt==1) then
          call cflx_tend( state, cam_in, ztodt, ptend)       
          call physics_update(state, ptend, ztodt, tend)
       end if

       call cnd_diag_checkpoint( diag, 'CFLXAPP', state, pbuf, cam_in, cam_out )

    else
    if (l_vdiff) then

       call t_startf('vertical_diffusion_tend')
       call vertical_diffusion_tend (ztodt ,state ,cam_in%wsx, cam_in%wsy,   &
            cam_in%shf     ,cam_in%cflx     ,surfric  ,obklen   ,ptend    ,ast    ,&
            cam_in%ocnfrac  , cam_in%landfrac ,        &
            sgh30    ,pbuf )

    !------------------------------------------
    ! Call major diffusion for extended model
    !------------------------------------------
    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call mspd_intr (ztodt    ,state    ,ptend)
    endif

       call physics_update(state, ptend, ztodt, tend)
       call t_stopf ('vertical_diffusion_tend')
       call cnd_diag_checkpoint( diag, 'VDIFF', state, pbuf, cam_in, cam_out )
    
    end if ! l_vdiff
    endif

    ! collect surface carbon fluxes
    call get_carbon_sfc_fluxes(state, cam_in, ztodt)


if (l_rayleigh) then
    !===================================================
    ! Rayleigh friction calculation
    !===================================================
    call t_startf('rayleigh_friction')
    call rayleigh_friction_tend( ztodt, state, ptend)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('rayleigh_friction')

    if (do_clubb_sgs .or. do_shoc_sgs) then
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, zero, zero, zero, zero)
    else
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, cam_in%cflx(:,1), zero, &
           zero, cam_in%shf)
    endif
    
    call check_tracers_chng(state, tracerint, "vdiff", nstep, ztodt, cam_in%cflx)

end if ! l_rayleigh

    call cnd_diag_checkpoint( diag, 'RAYLEIGH', state, pbuf, cam_in, cam_out )

if (l_tracer_aero) then

    !  aerosol dry deposition processes
    call t_startf('aero_drydep')
    call aero_model_drydep( state, pbuf, obklen, surfric, cam_in, ztodt, cam_out, ptend )
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('aero_drydep')

    !---------------------------------------------------------------------------------
    !	... enforce charge neutrality
    !---------------------------------------------------------------------------------
    call charge_fix( ncol, state%q(:,:,:) )

end if ! l_tracer_aero

    call cnd_diag_checkpoint( diag, 'AERDRYRM', state, pbuf, cam_in, cam_out )

if (l_gw_drag) then
    !===================================================
    ! Gravity wave drag
    !===================================================
    call t_startf('gw_tend')

    call gw_tend(state, sgh, pbuf, ztodt, ptend, cam_in)

    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, zero, zero, zero)
    call t_stopf('gw_tend')
    call cnd_diag_checkpoint( diag, 'GWDRAG', state, pbuf, cam_in, cam_out )

    ! QBO relaxation
    call qbo_relax(state, pbuf, ptend)
    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "qborelax", nstep, ztodt, zero, zero, zero, zero)

    ! Ion drag calculation
    call t_startf ( 'iondrag' )

    if ( do_waccm_ions ) then
       call iondrag_calc( lchnk, ncol, state, ptend, pbuf,  ztodt )
    else
       call iondrag_calc( lchnk, ncol, state, ptend)
    endif
    !----------------------------------------------------------------------------
    ! Call ionosphere routines for extended model if mode is set to ionosphere
    !----------------------------------------------------------------------------
    if( waccmx_is('ionosphere') ) then
       call ionos_intr(state, ptend, pbuf, ztodt)
    endif

    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "iondrag", nstep, ztodt, zero, zero, zero, zero)
    call t_stopf  ( 'iondrag' )

    call cnd_diag_checkpoint( diag, 'IONDRAG', state, pbuf, cam_in, cam_out )

end if ! l_gw_drag

    !===================================================
    ! Update Nudging values, if needed
    !===================================================
    if((Nudge_Model).and.(Nudge_ON)) then
      call nudging_timestep_tend(state,ptend)
      call physics_update(state,ptend,ztodt,tend)
    endif
    call cnd_diag_checkpoint( diag, 'NDG', state, pbuf, cam_in, cam_out )

if (l_ac_energy_chk) then
    !-------------- Energy budget checks vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    call pbuf_set_field(pbuf, teout_idx, state%te_cur, (/1,itim_old/),(/pcols,1/))       

    tmp_t(:ncol,:pver) = state%t(:ncol,:pver)

    ! store dse after tphysac in buffer
    do k = 1,pver
       dtcore(:ncol,k) = state%t(:ncol,k)
    end do

    !
    ! FV: convert dry-type mixing ratios to moist here because physics_dme_adjust
    !     assumes moist. This is done in p_d_coupling for other dynamics. Bundy, Feb 2004.


    if ( dycore_is('LR') .or. dycore_is('SE')) call set_dry_to_wet(state)    ! Physics had dry, dynamics wants moist


    ! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
    tmp_q     (:ncol,:pver) = state%q(:ncol,:pver,1)
    tmp_cldliq(:ncol,:pver) = state%q(:ncol,:pver,icldliq)
    tmp_cldice(:ncol,:pver) = state%q(:ncol,:pver,icldice)
    call physics_dme_adjust(state, tend, qini, ztodt)
!!!   REMOVE THIS CALL, SINCE ONLY Q IS BEING ADJUSTED. WON'T BALANCE ENERGY. TE IS SAVED BEFORE THIS
!!!   call check_energy_chng(state, tend, "drymass", nstep, ztodt, zero, zero, zero, zero)

    !-------------- Energy budget checks ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
end if ! l_ac_energy_chk

    call cnd_diag_checkpoint( diag, 'DRYWET', state, pbuf, cam_in, cam_out )

    if (aqua_planet) then
       labort = .false.
       do i=1,ncol
          if (cam_in%ocnfrac(i) /= 1._r8) labort = .true.
       end do
       if (labort) then
          call endrun ('TPHYSAC error:  grid contains non-ocean point')
       endif
    endif

    call diag_phys_tend_writeout (state, pbuf,  tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, &
         tmp_t, qini, cldliqini, cldiceini)

    ! DCAPE-ULL: record current state of T and q for computing dynamical tendencies
    !            the calculation follows the same format as in diag_phys_tend_writeout
    if (deep_scheme .eq. 'ZM' .and. (do_zmconv_dcape_ull .or. do_zmconv_dcape_only)) then
      ifld = pbuf_get_index('T_STAR')
      call pbuf_get_field(pbuf, ifld, t_star, (/1,1/),(/pcols,pver/))
      ifld = pbuf_get_index('Q_STAR')
      call pbuf_get_field(pbuf, ifld, q_star, (/1,1/),(/pcols,pver/))
      do k=1,pver
         q_star(:ncol,k) = state%q(:ncol,k,1)
         t_star(:ncol,k) = state%t(:ncol,k)
      enddo
    endif

    call clybry_fam_set( ncol, lchnk, map2chm, state%q, pbuf )

    static_ener_ac_idx = pbuf_get_index('static_ener_ac')
    call pbuf_get_field(pbuf, static_ener_ac_idx, static_ener_ac_2d )
    water_vap_ac_idx   = pbuf_get_index('water_vap_ac')
    call pbuf_get_field(pbuf, water_vap_ac_idx, water_vap_ac_2d )

    !Integrate column static energy
    ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
    do k=2,pver
       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    static_ener_ac_2d(:ncol) = ftem(:ncol,1)

    !Integrate water vapor
    ftem(:ncol,:) = state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
    do k=2,pver
       ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
    end do
    water_vap_ac_2d(:ncol) = ftem(:ncol,1)

    call check_tracers_fini(tracerint)
    call cnd_diag_checkpoint( diag, 'PACEND', state, pbuf, cam_in, cam_out )

end subroutine tphysac



subroutine additional_diagn_before_step_part1(state, cam_in)

    type(physics_state), intent(inout) :: state
    type(cam_in_t), intent(in)         :: cam_in

    integer ncol                               ! number of atmospheric columns

    ncol  = state%ncol

    !this assumes a very recent call to check_energy
    state%te_before_physstep(:ncol)=state%te_cur(:ncol)
    !save water mass here
    state%tw_before(:ncol) = state%tw_cur(:ncol)

    !save values as received from the coupler, before limiters on cflx
    state%cflx_raw(:ncol)  = cam_in%cflx(:ncol,1)
    state%shf_raw(:ncol) = cam_in%shf(:ncol)

end subroutine additional_diagn_before_step_part1



subroutine additional_diagn_before_step_part2(state, cam_in)

    type(physics_state), intent(inout) :: state
    type(cam_in_t), intent(in)         :: cam_in

    integer ncol                               ! number of atmospheric columns

    ncol  = state%ncol

    !save diffs after limiters on cflx
    state%cflx_diff(:ncol) = state%cflx_raw(:ncol) - cam_in%cflx(:ncol,1)
    state%shf_diff(:ncol)  = state%shf_raw(:ncol)  - cam_in%shf(:ncol)

end subroutine additional_diagn_before_step_part2



subroutine additional_diagn_after_step(state, cam_in, cam_out, fsns, fsnt, flns, flnt, ztodt)

    type(physics_state), intent(inout) :: state
    type(cam_in_t), intent(in)         :: cam_in
    type(cam_out_t), intent(in)        :: cam_out

    real(r8), intent(in) :: fsns(pcols)                   ! Surface solar absorbed flux
    real(r8), intent(in) :: fsnt(pcols)                   ! Net column abs solar flux at model top
    real(r8), intent(in) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
    real(r8), intent(in) :: flnt(pcols)                   ! Net outgoing lw flux at model top

    real(r8), intent(in) :: ztodt 

    integer ncol                               ! number of atmospheric columns

    ncol  = state%ncol

    !if including effects from qneg4
    !state%deltaw_flux(:ncol) = state%cflx_raw(:ncol)
    !if excluding qneg4 from diagnostics
    state%deltaw_flux(:ncol) = cam_in%cflx(:ncol,1)

    state%deltaw_flux(:ncol) = state%deltaw_flux(:ncol) - &
    1000.0*(cam_out%precc(:ncol)+cam_out%precl(:ncol))

    state%deltaw_flux(:ncol) = state%deltaw_flux(:ncol) * ztodt
    state%deltaw_step(:ncol) = state%tw_cur(:ncol) - state%tw_before(:ncol)

    !!!! compute net energy budget here
    !!!! expected TE2 - TE1 = restom - ressurf
    !    for TE2-TE1 we will use te_cur - te_before_physstep
    !    for restom - ressurf we will use fluxes to/from atm following AMWG energy scripts

    state%delta_te(1:ncol)=state%te_cur(1:ncol)-state%te_before_physstep(:ncol)

    !if including effects from qneg4
    !state%rr(1:ncol) = state%shf_raw(1:ncol)
    !if excluding qneg4 from diagnostisc -- do it consistently with cflx_raw
    state%rr(1:ncol) = cam_in%shf(1:ncol)

    state%rr(1:ncol) = state%rr(1:ncol)                    &
      + ( fsnt(1:ncol) - flnt(1:ncol) )                    &
      - ( fsns(1:ncol) - flns(1:ncol)                      &
      -        (latvap+latice)*cam_in%cflx(:ncol,1)        &
      + 1000.0* latice        *( cam_out%precc(1:ncol)+cam_out%precl(1:ncol) - cam_out%precsc(1:ncol) - cam_out%precsl(1:ncol) ) )

    ! deltaw_flux and deltaw_step are in kg/m2
    ! delta_te is in J/m2
    ! rr is in J/m2/sec = W/m2

end subroutine additional_diagn_after_step



subroutine tphysbc (ztodt,               &
       fsns,    fsnt,    flns,    flnt,    state,   &
       tend,    pbuf,    diag,    fsds,             &
       sgh, sgh30, cam_out, cam_in )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Evaluate and apply physical processes that are calculated BEFORE 
    ! coupling to land, sea, and ice models.  
    !
    ! Processes currently included are: 
    ! dry adjustment, moist convection, stratiform, wet deposition, radiation
    !
    ! Pass surface fields for separate surface flux calculations
    ! Dump appropriate fields to history file.
    ! 
    ! Method: 
    !
    ! Each parameterization should be implemented with this sequence of calls:
    !  1)  Call physics interface
    !  2)  Check energy
    !  3)  Call physics_update
    ! See Interface to Column Physics and Chemistry Packages 
    !   http://www.ccsm.ucar.edu/models/atm-cam/docs/phys-interface/index.html
    ! 
    ! Author: CCM1, CMS Contact: J. Truesdale
    !         modified by A. Gettelman and C. Craig Nov 2010 to separate micro/macro physics
    ! 
    !-----------------------------------------------------------------------

    use physics_buffer,          only : physics_buffer_desc, pbuf_get_field
    use physics_buffer,          only : pbuf_get_index, pbuf_old_tim_idx
    use physics_buffer,          only : col_type_subcol, dyn_time_lvls
    use shr_kind_mod,    only: r8 => shr_kind_r8

    use stratiform,      only: stratiform_tend
    use microp_driver,   only: microp_driver_tend
    use microp_aero,     only: microp_aero_run
    use macrop_driver,   only: macrop_driver_tend
    use physics_types,   only: physics_state, physics_tend, physics_ptend, &
         physics_ptend_init, physics_ptend_sum, physics_state_check, physics_ptend_scale
    use cam_diagnostics, only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, diag_export, diag_state_b4_phys_write
    use cam_history,     only: outfld, fieldname_len
    use physconst,       only: cpair, latvap, gravit, rga
    use constituents,    only: pcnst, qmin, cnst_get_ind
    use convect_deep,    only: convect_deep_tend, convect_deep_tend_2, deep_scheme_does_scav_trans
    use time_manager,    only: get_nstep
    use convect_shallow, only: convect_shallow_tend
    use check_energy,    only: check_energy_chng, check_energy_fix, &
                               check_qflx, check_water, check_prect, & 
                               check_energy_timestep_init, &
                               check_tracers_data, check_tracers_init, &
                               check_tracers_chng, check_tracers_fini
    use dycore,          only: dycore_is
    use aero_model,      only: aero_model_wetdep
    use froude,          only: calc_uovern
    use radiation,       only: radiation_tend
    use cloud_diagnostics, only: cloud_diagnostics_calc
    use perf_mod
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_adj
    use clubb_intr,      only: clubb_tend_cam
    use shoc_intr,       only: shoc_tend_e3sm
    use sslt_rebin,      only: sslt_rebin_adv
    use tropopause,      only: tropopause_output, tropopause_e90_3d_output
    use output_aerocom_aie, only: do_aerocom_ind3, cloud_top_aerocom
    use cam_abortutils,      only: endrun
    use subcol,          only: subcol_gen, subcol_ptend_avg
    use subcol_utils,    only: subcol_ptend_copy, is_subcol_on
    use phys_control,    only: use_qqflx_fixer, use_mass_borrower
    use nudging,         only: Nudge_Model,Nudge_Loc_PhysOut,nudging_calc_tend
    use debug_info,      only: get_debug_chunk, get_debug_macmiciter
    use lnd_infodata,    only: precip_downscaling_method
    use cflx,            only: cflx_tend

    implicit none

    !
    ! Arguments
    !
    real(r8), intent(in) :: ztodt                            ! 2 delta t (model time increment)
    real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
    real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
    real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
    real(r8), intent(inout) :: fsds(pcols)                   ! Surface solar down flux
    real(r8), intent(in) :: sgh(pcols)                       ! Std. deviation of orography
    real(r8), intent(in) :: sgh30(pcols)                     ! Std. deviation of 30 s orography for tms

    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)
    type(cnd_diag_t),    intent(inout) :: diag

    type(cam_out_t),     intent(inout) :: cam_out
    type(cam_in_t),      intent(in)    :: cam_in


    !
    !---------------------------Local workspace-----------------------------
    !

    type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
    type(physics_state)   :: state_sc         ! state for sub-columns
    type(physics_ptend)   :: ptend_sc         ! ptend for sub-columns
    type(physics_ptend)   :: ptend_aero       ! ptend for microp_aero
    type(physics_ptend)   :: ptend_aero_sc    ! ptend for microp_aero on sub-columns
    type(physics_tend)    :: tend_sc          ! tend for sub-columns

    integer :: nstep                          ! current timestep number

    real(r8) :: net_flx(pcols)

    real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
    real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c

    real(r8) cmfcme(pcols,pver)                ! cmf condensation - evaporation
    real(r8) cmfmc2(pcols,pverp)               ! Moist convection cloud mass flux
    real(r8) dlf(pcols,pver)                   ! Detraining cld H20 from shallow + deep convections
    real(r8) dlf2(pcols,pver)                  ! Detraining cld H20 from shallow convections
    real(r8) pflx(pcols,pverp)                 ! Conv rain flux thru out btm of lev
    real(r8) rtdt                              ! 1./ztodt

    integer lchnk                              ! chunk identifier
    integer ncol                               ! number of atmospheric columns
    integer ierr

    integer  i,k,m,ihist                       ! Longitude, level, constituent indices
    ! for macro/micro co-substepping
    integer :: macmic_it                       ! iteration variables
    real(r8) :: cld_macmic_ztodt               ! modified timestep
    character(len=2) :: char_macmic_it

    ! physics buffer fields to compute tendencies for stratiform package
    integer itim_old, ifld
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction

    character(len=16)  :: deep_scheme      ! Default set in phys_control.F90

    ! physics buffer fields for total energy and mass adjustment
    real(r8), pointer, dimension(:  ) :: teout
    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore

    real(r8), pointer, dimension(:,:,:) :: fracis  ! fraction of transported species that are insoluble

    ! convective precipitation variables
    real(r8),pointer :: prec_dp(:)                ! total precipitation from ZM convection
    real(r8),pointer :: snow_dp(:)                ! snow from ZM convection
    real(r8),pointer :: prec_sh(:)                ! total precipitation from Hack convection
    real(r8),pointer :: snow_sh(:)                ! snow from Hack convection
    real(r8) :: totice_dp(pcols)                  ! total frozen precip from ZM convection

    ! stratiform precipitation variables
    real(r8),pointer :: prec_str(:)    ! sfc flux of precip from stratiform (m/s)
    real(r8),pointer :: snow_str(:)     ! sfc flux of snow from stratiform   (m/s)
    real(r8),pointer :: prec_str_sc(:)  ! sfc flux of precip from stratiform (m/s) -- for subcolumns
    real(r8),pointer :: snow_str_sc(:)  ! sfc flux of snow from stratiform   (m/s) -- for subcolumns
    real(r8),pointer :: prec_pcw(:)     ! total precip from prognostic cloud scheme
    real(r8),pointer :: snow_pcw(:)     ! snow from prognostic cloud scheme
    real(r8),pointer :: prec_sed(:)     ! total precip from cloud sedimentation
    real(r8),pointer :: snow_sed(:)     ! snow from cloud ice sedimentation
    real(r8) :: sh_e_ed_ratio(pcols,pver)       ! shallow conv [ent/(ent+det)] ratio  


    ! Local copies for substepping
    real(r8) :: prec_pcw_macmic(pcols)
    real(r8) :: snow_pcw_macmic(pcols)
    real(r8) :: prec_sed_macmic(pcols)
    real(r8) :: snow_sed_macmic(pcols)

    ! energy checking variables
    real(r8) :: zero(pcols)                    ! array of zeros
    real(r8) :: zero_sc(pcols*psubcols)        ! array of zeros
    real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
    real(r8) :: rice(pcols)                    ! vertical integral of ice not yet in q(ixcldice)
    real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
    real(r8) :: det_s  (pcols)                 ! vertical integral of detrained static energy from ice
    real(r8) :: det_ice(pcols)                 ! vertical integral of detrained ice
    real(r8) :: flx_cnd(pcols)
    real(r8) :: flx_heat(pcols)
    type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes
    real(r8) :: zero_tracers(pcols,pcnst)

    logical   :: lq(pcnst)

    character(len=fieldname_len)   :: varname, vsuffix
    !BSINGH - Following variables are from zm_conv_intr, which are moved here as they are now used
    ! by aero_model_wetdep subroutine. 

    real(r8):: mu(pcols,pver) 
    real(r8):: eu(pcols,pver)
    real(r8):: du(pcols,pver)
    real(r8):: md(pcols,pver)
    real(r8):: ed(pcols,pver)
    real(r8):: dp(pcols,pver)
    
    ! wg layer thickness in mbs (between upper/lower interface).
    real(r8):: dsubcld(pcols)
    
    ! wg layer thickness in mbs between lcl and maxi.    
    integer :: jt(pcols)
    
    ! wg top  level index of deep cumulus convection.
    integer :: maxg(pcols)
    
    ! wg gathered values of maxi.
    integer :: ideep(pcols)
    
    ! w holds position of gathered points vs longitude index
    integer :: lengath

    real(r8)  :: lcldo(pcols,pver)              !Pass old liqclf from macro_driver to micro_driver

    real(r8) :: ftem(pcols,pver)         ! tmp space
    real(r8), pointer, dimension(:) :: static_ener_ac_2d ! Vertically integrated static energy
    real(r8), pointer, dimension(:) :: water_vap_ac_2d   ! Vertically integrated water vapor
    real(r8) :: CIDiff(pcols)            ! Difference in vertically integrated static energy

    !HuiWan (2014/15): added for a short-term time step convergence test ++ 
    logical :: l_bc_energy_fix
    logical :: l_dry_adj
    logical :: l_tracer_aero
    logical :: l_st_mac
    logical :: l_st_mic
    logical :: l_rad
    !HuiWan (2014/15): added for a short-term time step convergence test ==

    ! Numerical schemes for process coupling
    integer :: cflx_cpl_opt  ! When to apply surface tracer fluxes  (not including water vapor).
                             ! The default for aerosols is to do this 
                             ! after tphysac:clubb_surface and before aerosol dry removal.
                             ! For chemical gases, different versions of EAM 
                             ! might use different process ordering.
    !-----------------------------------------------------------------------
    call cnd_diag_checkpoint( diag, 'DYNEND', state, pbuf, cam_in, cam_out )
    !-----------------------------------------------------------------------

    call phys_getopts( microp_scheme_out      = microp_scheme, &
                       macrop_scheme_out      = macrop_scheme, &
                       use_subcol_microp_out  = use_subcol_microp, &
                       deep_scheme_out        = deep_scheme,       &
                       state_debug_checks_out = state_debug_checks &
                      ,cflx_cpl_opt_out       = cflx_cpl_opt       &
                      ,l_bc_energy_fix_out    = l_bc_energy_fix    &
                      ,l_dry_adj_out          = l_dry_adj          &
                      ,l_tracer_aero_out      = l_tracer_aero      &
                      ,l_st_mac_out           = l_st_mac           &
                      ,l_st_mic_out           = l_st_mic           &
                      ,l_rad_out              = l_rad              &
                      )
    
    !-----------------------------------------------------------------------
    call t_startf('bc_init')

    zero = 0._r8
    zero_tracers(:,:) = 0._r8
    zero_sc(:) = 0._r8
    flx_cnd(:) = 0._r8
    flx_heat(:) = 0._r8

    lchnk = state%lchnk
    ncol  = state%ncol

    call get_debug_chunk(lchnk)

    rtdt = 1._r8/ztodt

    nstep = get_nstep()

    if (pergro_test_active) then 
       !call outfld calls
       do ihist = 1 , nvars_prtrb_hist
          vsuffix  = trim(adjustl(hist_vars(ihist)))
          varname  = trim(adjustl(vsuffix))//'_topphysbc' ! form variable name
          call outfld( trim(adjustl(varname)),get_var(state,vsuffix), pcols , lchnk )
       enddo
    endif

    static_ener_ac_idx = pbuf_get_index('static_ener_ac')
    call pbuf_get_field(pbuf, static_ener_ac_idx, static_ener_ac_2d )
    water_vap_ac_idx   = pbuf_get_index('water_vap_ac')
    call pbuf_get_field(pbuf, water_vap_ac_idx, water_vap_ac_2d )

    ! Integrate and compute the difference
    ! CIDiff = difference of column integrated values
    if( nstep == 0 ) then
       CIDiff(:ncol) = 0.0_r8
       call outfld('DTENDTH', CIDiff, pcols, lchnk )
       call outfld('DTENDTQ', CIDiff, pcols, lchnk )
    else
       ! MSE first
       ftem(:ncol,:) = (state%s(:ncol,:) + latvap*state%q(:ncol,:,1)) * state%pdel(:ncol,:)*rga
       do k=2,pver
          ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
       end do
       CIDiff(:ncol) = (ftem(:ncol,1) - static_ener_ac_2d(:ncol))*rtdt

       call outfld('DTENDTH', CIDiff, pcols, lchnk )
       ! Water vapor second
       ftem(:ncol,:) = state%q(:ncol,:,1)*state%pdel(:ncol,:)*rga
       do k=2,pver
          ftem(:ncol,1) = ftem(:ncol,1) + ftem(:ncol,k)
       end do
       CIDiff(:ncol) = (ftem(:ncol,1) - water_vap_ac_2d(:ncol))*rtdt

       call outfld('DTENDTQ', CIDiff, pcols, lchnk )
    end if

    ! Associate pointers with physics buffer fields
    itim_old = pbuf_old_tim_idx()
    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim_old/),(/pcols,pver,1/))

    call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim_old/), (/pcols,1/))

    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld   =  pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

    ifld    = pbuf_get_index('FRACIS')
    call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/)  )
    fracis (:ncol,:,1:pcnst) = 1._r8

    ! Set physics tendencies to 0
    tend %dTdt(:ncol,:pver)  = 0._r8
    tend %dudt(:ncol,:pver)  = 0._r8
    tend %dvdt(:ncol,:pver)  = 0._r8

!!== KZ_WCON
    call check_qflx (state, tend, "PHYBC01", nstep, ztodt, cam_in%cflx(:,1))
    call check_water(state, tend, "PHYBC01", nstep, ztodt)


    if(use_mass_borrower) then 

      !! printout diagnostic information
      !!.................................................................
       call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
            1, pcnst, qmin  ,state%q, .False.)

      !! tracer borrower for mass conservation 
      !!.................................................................

       do m = 1, pcnst 
          call massborrow("PHYBC01",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
       end do

!!      call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
!!           1, pcnst, qmin  ,state%q, .True. )
    else

      !! original fixer to make sure tracers are all positive
      !!.................................................................

      call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
           1, pcnst, qmin  ,state%q, .True. )

    end if 

!!== KZ_WCON

    ! Validate state coming from the dynamics.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysbc (dycore?)")

    call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

!!== KZ_WCON

    if(use_mass_borrower) then

       !! if use_mass_borrower = True, only printout diagnostic information
       !!.................................................................

       call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
            1, pcnst, qmin  ,state%q, .False. )

       !! tracer borrower for mass conservation 
       !!.................................................................
       do m = 1, pcnst
          call massborrow("PHYBC02",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
       end do

!!       call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
!!            1, pcnst, qmin  ,state%q, .True. )

    else

       !! original fixer to make sure tracers are all positive
       !!.................................................................

       call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
            1, pcnst, qmin  ,state%q, .True. )

    end if


    call check_water(state, tend, "PHYBC02", nstep, ztodt)
!!== KZ_WCON

    ! Validate output of clybry_fam_adj.
    if (state_debug_checks) &
         call physics_state_check(state, name="clybry_fam_adj")

    !
    ! Dump out "before physics" state
    !
    call diag_state_b4_phys_write (state)

    ! compute mass integrals of input tracers state
    call check_tracers_init(state, tracerint)

    call t_stopf('bc_init')

    !===================================================
    ! Global mean total energy fixer
    !===================================================
if (l_bc_energy_fix) then

    call t_startf('energy_fixer')

    tini(:ncol,:pver) = state%t(:ncol,:pver)
    if (dycore_is('LR') .or. dycore_is('SE'))  then
       call check_energy_fix(state, ptend, nstep, flx_heat)
       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
    end if
    ! Save state for convective tendency calculations.
    call diag_conv_tend_ini(state, pbuf)

    qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
    cldliqini(:ncol,:pver) = state%q(:ncol,:pver,icldliq)
    cldiceini(:ncol,:pver) = state%q(:ncol,:pver,icldice)


    call outfld('TEOUT', teout       , pcols, lchnk   )
    call outfld('TEINP', state%te_ini, pcols, lchnk   )
    call outfld('TEFIX', state%te_cur, pcols, lchnk   )

    ! set and output the dse change due to dynpkg
    if( nstep > dyn_time_lvls-1 ) then
       do k = 1,pver
          dtcore(:ncol,k) = (tini(:ncol,k) - dtcore(:ncol,k))/(ztodt) + tend%dTdt(:ncol,k)
       end do
       call outfld( 'DTCORE', dtcore, pcols, lchnk )
    end if

    call t_stopf('energy_fixer')

end if

    call cnd_diag_checkpoint( diag, 'PBCINI', state, pbuf, cam_in, cam_out )

    !
    !===================================================
    ! Dry adjustment
    ! This code block is not a good example of interfacing a parameterization
    !===================================================
if (l_dry_adj) then

    call t_startf('dry_adjustment')

    ! Copy state info for input to dadadj
    ! This is a kludge, so that dadadj does not have to be correctly reformulated in dry static energy

    lq(:) = .FALSE.
    lq(1) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)
    ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
    ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

    call dadadj (lchnk, ncol, state%pmid,  state%pint,  state%pdel,  &
         ptend%s, ptend%q(1,1,1))
    ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
    ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt
    call physics_update(state, ptend, ztodt, tend)

    call t_stopf('dry_adjustment')

end if
    call cnd_diag_checkpoint( diag, 'DRYADJ', state, pbuf, cam_in, cam_out )

    !
    !===================================================
    ! Moist convection
    !===================================================
    call t_startf('moist_convection')
    !
    ! Since the PBL doesn't pass constituent perturbations, they
    ! are zeroed here for input to the moist convection routine
    !
    call t_startf ('convect_deep_tend')
    call convect_deep_tend(  &
         cmfmc,      cmfcme,             &
         dlf,        pflx,    zdu,       &
         rliq,       rice, &
         ztodt,   &
         state,   ptend, cam_in%landfrac, pbuf, mu, eu, du, md, ed, dp,   &
         dsubcld, jt, maxg, ideep, lengath) 
    call t_stopf('convect_deep_tend')

    call physics_update(state, ptend, ztodt, tend)

    call pbuf_get_field(pbuf, prec_dp_idx, prec_dp )
    call pbuf_get_field(pbuf, snow_dp_idx, snow_dp )
    call pbuf_get_field(pbuf, prec_sh_idx, prec_sh )
    call pbuf_get_field(pbuf, snow_sh_idx, snow_sh )
    call pbuf_get_field(pbuf, prec_str_idx, prec_str)
    call pbuf_get_field(pbuf, snow_str_idx, snow_str)
    call pbuf_get_field(pbuf, prec_sed_idx, prec_sed)
    call pbuf_get_field(pbuf, snow_sed_idx, snow_sed)
    call pbuf_get_field(pbuf, prec_pcw_idx, prec_pcw )
    call pbuf_get_field(pbuf, snow_pcw_idx, snow_pcw )

    if (use_subcol_microp) then
      call pbuf_get_field(pbuf, prec_str_idx, prec_str_sc, col_type=col_type_subcol)
      call pbuf_get_field(pbuf, snow_str_idx, snow_str_sc, col_type=col_type_subcol)
    end if

    ! Check energy integrals, including "reserved liquid"
    flx_cnd(:ncol) = prec_dp(:ncol) + rliq(:ncol)
    totice_dp(:ncol) = snow_dp(:ncol) + rice(:ncol)
    call check_energy_chng(state, tend, "convect_deep", nstep, ztodt, zero, flx_cnd, totice_dp, zero)

    call cnd_diag_checkpoint( diag, 'DEEPCU', state, pbuf, cam_in, cam_out )

    !
    ! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
    !
    call t_startf ('convect_shallow_tend')

    call convect_shallow_tend (ztodt   , cmfmc,  cmfmc2  ,&
         dlf        , dlf2   ,  rliq   , rliq2, & 
         state      , ptend  ,  pbuf   , sh_e_ed_ratio   , sgh, sgh30, cam_in) 
    call t_stopf ('convect_shallow_tend')

    call physics_update(state, ptend, ztodt, tend)

    flx_cnd(:ncol) = prec_sh(:ncol) + rliq2(:ncol)
    call check_energy_chng(state, tend, "convect_shallow", nstep, ztodt, zero, flx_cnd, snow_sh, zero)

    call check_tracers_chng(state, tracerint, "convect_shallow", nstep, ztodt, zero_tracers)

    call cnd_diag_checkpoint( diag, 'SHCU', state, pbuf, cam_in, cam_out )

    call t_stopf('moist_convection')

if (l_tracer_aero) then

    ! Rebin the 4-bin version of sea salt into bins for coarse and accumulation
    ! modes that correspond to the available optics data.  This is only necessary
    ! for CAM-RT.  But it's done here so that the microphysics code which is called
    ! from the stratiform interface has access to the same aerosols as the radiation
    ! code.
    call sslt_rebin_adv(pbuf,  state)
    
end if


    !========================================================================================
    ! Stratiform cloud macro and microphysics, turbulence, aerosol activation-resuspension
    !========================================================================================
    if( microp_scheme == 'RK' ) then

     if (l_st_mac.or.l_st_mic) then
       !===================================================
       ! Calculate stratiform tendency (sedimentation, detrain, cloud fraction and microphysics )
       !===================================================
       call t_startf('stratiform_tend')

       call stratiform_tend(state, ptend, pbuf, ztodt, &
            cam_in%icefrac, cam_in%landfrac, cam_in%ocnfrac, &
            cam_in%snowhland, & ! sediment
            dlf, dlf2, & ! detrain
            rliq  , & ! check energy after detrain
            cmfmc,   cmfmc2, &
            cam_in%ts,      cam_in%sst,        zdu)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_str(:ncol), snow_str(:ncol), zero)

       call t_stopf('stratiform_tend')
     end if !l_st_mac

    elseif( microp_scheme == 'MG' .or. microp_scheme == 'P3' ) then

       !========================================================================================
       ! Apply surface tracer fluxes to update tracer mixing ratios before turbulent tranport 
       !========================================================================================
       ! Diagnose tracer mixing ratio tendencies from surface fluxes, then update the mixing ratios.
       ! Note that these subroutine calls do not touch water vapor. They also have no effects
       ! on tracers for which cam_in%cflx(:,m) is zero at this point.

      !if ( do_clubb_sgs .and. (cflx_cpl_opt==2) ) then
       if ( cflx_cpl_opt==2 ) then
          call cflx_tend( state, cam_in, ztodt, ptend)
          call physics_update(state, ptend, ztodt, tend)
       end if

       !========================================================================================
       ! Start co-substepping of macrophysics and microphysics
       cld_macmic_ztodt = ztodt/cld_macmic_num_steps

       ! Clear precip fields that should accumulate.
       prec_sed_macmic = 0._r8
       snow_sed_macmic = 0._r8
       prec_pcw_macmic = 0._r8
       snow_pcw_macmic = 0._r8

       do macmic_it = 1, cld_macmic_num_steps

        call get_debug_macmiciter(macmic_it)
        write(char_macmic_it,'(i2.2)') macmic_it

        if (l_st_mac) then

          if (micro_do_icesupersat) then 

            !===================================================
            ! Aerosol Activation
            !===================================================
            call t_startf('microp_aero_run')
            call microp_aero_run(state, ptend, cld_macmic_ztodt, pbuf, lcldo)
            call t_stopf('microp_aero_run')

            call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)

            call physics_update(state, ptend, ztodt, tend)
            call check_energy_chng(state, tend, "mp_aero_tend", nstep, ztodt, zero, zero, zero, zero)      

          endif
          !===================================================
          ! Calculate macrophysical tendency (sedimentation, detrain, cloud fraction)
          !===================================================

          call t_startf('macrop_tend')

          ! don't call Park macrophysics if CLUBB is called
          if (macrop_scheme .ne. 'CLUBB_SGS' .and. macrop_scheme .ne. 'SHOC_SGS') then

             call macrop_driver_tend( &
                  state,           ptend,          cld_macmic_ztodt, &
                  cam_in%landfrac, cam_in%ocnfrac, cam_in%snowhland, & ! sediment
                  dlf,             dlf2,                             & ! detrain
                  cmfmc,           cmfmc2,                           &
                  cam_in%ts,       cam_in%sst,     zdu,              &
                  pbuf,            det_s,          det_ice,          &
                  lcldo)

             !  Since we "added" the reserved liquid back in this routine, we need 
             !    to account for it in the energy checker
             flx_cnd(:ncol) = -1._r8*rliq(:ncol) 
             flx_heat(:ncol) = det_s(:ncol)

             ! Unfortunately, physics_update does not know what time period
             ! "tend" is supposed to cover, and therefore can't update it
             ! with substeps correctly. For now, work around this by scaling
             ! ptend down by the number of substeps, then applying it for
             ! the full time (ztodt).
             call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)          
             call physics_update(state, ptend, ztodt, tend)
             call check_energy_chng(state, tend, "macrop_tend", nstep, ztodt, &
                  zero, flx_cnd/cld_macmic_num_steps, &
                  det_ice/cld_macmic_num_steps, flx_heat/cld_macmic_num_steps)
       
          else ! Calculate CLUBB macrophysics


!!== KZ_WATCON 

    !! qqflx fixer to avoid negative water vapor in the surface layer
    !! due to strong negative qflx  
    !!.................................................................

    if(use_qqflx_fixer) then
       call qqflx_fixer('TPHYSBC ', lchnk, ncol, cld_macmic_ztodt, &
            state%q(1,1,1), state%rpdel(1,1), cam_in%shf, &
            cam_in%lhf , cam_in%cflx/cld_macmic_num_steps )

    end if
!!== KZ_WATCON 

             ! =====================================================
             !    CLUBB call (PBL, shallow convection, macrophysics)
             ! =====================================================  
           if (do_clubb_sgs) then
             call clubb_tend_cam(state,ptend,pbuf,diag,cld_macmic_ztodt,&
                cmfmc, cam_in, cam_out, sgh30, macmic_it, cld_macmic_num_steps, &
                dlf, det_s, det_ice, lcldo)
	   endif
   
           if (do_shoc_sgs) then
             call shoc_tend_e3sm(state,ptend,pbuf,cld_macmic_ztodt,&
                cmfmc, cam_in, sgh30, macmic_it, cld_macmic_num_steps, & 
                dlf, det_s, det_ice, lcldo)
           endif   

                !  Since we "added" the reserved liquid back in this routine, we need 
                !    to account for it in the energy checker
                flx_cnd(:ncol) = -1._r8*rliq(:ncol) 
                flx_heat(:ncol) = cam_in%shf(:ncol) + det_s(:ncol)

                ! Unfortunately, physics_update does not know what time period
                ! "tend" is supposed to cover, and therefore can't update it
                ! with substeps correctly. For now, work around this by scaling
                ! ptend down by the number of substeps, then applying it for
                ! the full time (ztodt).
                call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)
                !    Update physics tendencies and copy state to state_eq, because that is 
                !      input for microphysics              
                call physics_update(state, ptend, ztodt, tend)
                call check_energy_chng(state, tend, "clubb_tend", nstep, ztodt, &
                     cam_in%cflx(:,1)/cld_macmic_num_steps, flx_cnd/cld_macmic_num_steps, &
                     det_ice/cld_macmic_num_steps, flx_heat/cld_macmic_num_steps)


 
          endif

          call t_stopf('macrop_tend')
        end if ! l_st_mac
        call cnd_diag_checkpoint( diag, 'CLDMAC'//char_macmic_it, state, pbuf, cam_in, cam_out )

          !===================================================
          ! Calculate cloud microphysics 
          !===================================================
        if (l_st_mic) then

          if (is_subcol_on()) then
             ! Allocate sub-column structures. 
             call physics_state_alloc(state_sc, lchnk, psubcols*pcols)
             call physics_tend_alloc(tend_sc, psubcols*pcols)

             ! Generate sub-columns using the requested scheme
             call subcol_gen(state, tend, state_sc, tend_sc, pbuf)

             !Initialize check energy for subcolumns
             call check_energy_timestep_init(state_sc, tend_sc, pbuf, col_type_subcol)
          end if

          if (.not. micro_do_icesupersat) then 

            call t_startf('microp_aero_run')
            call microp_aero_run(state, ptend_aero, cld_macmic_ztodt, pbuf, lcldo)
            call t_stopf('microp_aero_run')

          endif
          call cnd_diag_checkpoint( diag, 'CLDAER'//char_macmic_it, state, pbuf, cam_in, cam_out )


          call t_startf('microp_tend')


          if (use_subcol_microp) then
             call microp_driver_tend(state_sc, ptend_sc, cld_macmic_ztodt, pbuf)

             ! Average the sub-column ptend for use in gridded update - will not contain ptend_aero
             call subcol_ptend_avg(ptend_sc, state_sc%ngrdcol, lchnk, ptend)

             ! Copy ptend_aero field to one dimensioned by sub-columns before summing with ptend
             call subcol_ptend_copy(ptend_aero, state_sc, ptend_aero_sc)
             call physics_ptend_sum(ptend_aero_sc, ptend_sc, state_sc%ncol)
             call physics_ptend_dealloc(ptend_aero_sc)

             ! Have to scale and apply for full timestep to get tend right
             ! (see above note for macrophysics).
             call physics_ptend_scale(ptend_sc, 1._r8/cld_macmic_num_steps, ncol)

             call physics_update (state_sc, ptend_sc, ztodt, tend_sc)
             call check_energy_chng(state_sc, tend_sc, "microp_tend_subcol", &
                  nstep, ztodt, zero_sc, prec_str_sc(:ncol)/cld_macmic_num_steps, &
                  snow_str_sc(:ncol)/cld_macmic_num_steps, zero_sc)

             call physics_state_dealloc(state_sc)
             call physics_tend_dealloc(tend_sc)
             call physics_ptend_dealloc(ptend_sc)
          else
             call microp_driver_tend(state, ptend, cld_macmic_ztodt, pbuf)
          end if
          ! combine aero and micro tendencies for the grid
          if (.not. micro_do_icesupersat) then
             call physics_ptend_sum(ptend_aero, ptend, ncol)
             call physics_ptend_dealloc(ptend_aero)
          endif

          ! Have to scale and apply for full timestep to get tend right
          ! (see above note for macrophysics).
          call physics_ptend_scale(ptend, 1._r8/cld_macmic_num_steps, ncol)

          call physics_update (state, ptend, ztodt, tend)
          call check_energy_chng(state, tend, "microp_tend", nstep, ztodt, &
               zero, prec_str(:ncol)/cld_macmic_num_steps, &
               snow_str(:ncol)/cld_macmic_num_steps, zero)

          call t_stopf('microp_tend')

        else 
        ! If microphysics is off, set surface cloud liquid/ice and rain/snow fluxes to zero

          prec_sed = 0._r8
          snow_sed = 0._r8
          prec_pcw = 0._r8
          snow_pcw = 0._r8

        end if ! l_st_mic

          prec_sed_macmic(:ncol) = prec_sed_macmic(:ncol) + prec_sed(:ncol)
          snow_sed_macmic(:ncol) = snow_sed_macmic(:ncol) + snow_sed(:ncol)
          prec_pcw_macmic(:ncol) = prec_pcw_macmic(:ncol) + prec_pcw(:ncol)
          snow_pcw_macmic(:ncol) = snow_pcw_macmic(:ncol) + snow_pcw(:ncol)

          call cnd_diag_checkpoint( diag, 'CLDMIC'//char_macmic_it, state, pbuf, cam_in, cam_out )

       end do ! end substepping over macrophysics/microphysics

       prec_sed(:ncol) = prec_sed_macmic(:ncol)/cld_macmic_num_steps
       snow_sed(:ncol) = snow_sed_macmic(:ncol)/cld_macmic_num_steps
       prec_pcw(:ncol) = prec_pcw_macmic(:ncol)/cld_macmic_num_steps
       snow_pcw(:ncol) = snow_pcw_macmic(:ncol)/cld_macmic_num_steps
       prec_str(:ncol) = prec_pcw(:ncol) + prec_sed(:ncol)
       snow_str(:ncol) = snow_pcw(:ncol) + snow_sed(:ncol)

     end if !microp_scheme

     call cnd_diag_checkpoint( diag, 'STCLD', state, pbuf, cam_in, cam_out )
     !----------------------------------------------------------------------

   if (l_tracer_aero) then
      if ( .not. deep_scheme_does_scav_trans() ) then

         !======================================================================
         ! Aerosol wet chemistry determines scavenging and transformations.
         ! This is followed by convective transport of all trace species except
         ! water vapor and condensate. Scavenging needs to occur prior to
         ! transport in order to determine interstitial fraction.
         !======================================================================

         call t_startf('tphysbc_aerosols')

         if (do_clubb_sgs .or. do_shoc_sgs) then
            sh_e_ed_ratio = 0.0_r8
         endif

         ! Aerosol wet removal (including aerosol water uptake)
         call aero_model_wetdep( ztodt, dlf, dlf2, cmfmc2, state,  & ! inputs
                sh_e_ed_ratio, mu, md, du, eu, ed, dp, dsubcld,    &
                jt, maxg, ideep, lengath, species_class,           &
                cam_out, pbuf, ptend)                      ! outputs
         call physics_update(state, ptend, ztodt, tend)

         ! deep convective aerosol transport
         call convect_deep_tend_2( state, ptend, ztodt, pbuf, &
                mu, eu, du, md, ed, dp, dsubcld, jt, maxg,    &
                ideep, lengath, species_class )
         call physics_update(state, ptend, ztodt, tend)

         ! check tracer integrals
         call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt,  zero_tracers)

         call t_stopf('tphysbc_aerosols')

      end if
   end if ! l_tracer_aero

    call cnd_diag_checkpoint( diag, 'AERWETRM', state, pbuf, cam_in, cam_out )

    !===================================================
    ! Moist physical parameteriztions complete: 
    ! send dynamical variables, and derived variables to history file
    !===================================================

    call t_startf('bc_history_write')
    call diag_phys_writeout(state, cam_out%psl)
    call diag_conv(state, ztodt, pbuf)

    call t_stopf('bc_history_write')

    !===================================
    ! Update Nudging tendency if needed
    !===================================
    if (Nudge_Model .and. Nudge_Loc_PhysOut) then
       call nudging_calc_tend(state)
    endif

    !===================================================
    ! Write cloud diagnostics on history file
    !===================================================

    call t_startf('bc_cld_diag_history_write')

    call cloud_diagnostics_calc(state, pbuf)

    call t_stopf('bc_cld_diag_history_write')

    call cnd_diag_checkpoint( diag, 'PBCDIAG', state, pbuf, cam_in, cam_out )

if (l_rad) then
    !===================================================
    ! Radiation computations
    !===================================================
    call t_startf('radiation')


    call radiation_tend(state,ptend, pbuf, &
         cam_out, cam_in, &
         cam_in%landfrac, cam_in%icefrac, cam_in%snowhland, &
         fsns,    fsnt, flns,    flnt,  &
         fsds, net_flx,is_cmip6_volc, ztodt)

    ! Set net flux used by spectral dycores
    do i=1,ncol
       tend%flx_net(i) = net_flx(i)
    end do
    call physics_update(state, ptend, ztodt, tend)
    call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, net_flx)

    call t_stopf('radiation')

end if ! l_rad

    call cnd_diag_checkpoint( diag, 'RAD', state, pbuf, cam_in, cam_out )

    if(do_aerocom_ind3) then
       call cloud_top_aerocom(state, pbuf) 
    end if

    if (trim(adjustl(precip_downscaling_method)) == "FNM") then
       !if the land model's precip downscaling method is FNM, compute uovern
       call calc_uovern(state,cam_out)
    end if

    ! Diagnose the location of the tropopause and its location to the history file(s).
    call t_startf('tropopause')
    call tropopause_output(state)
    call tropopause_e90_3d_output(state)
    call t_stopf('tropopause')

    ! Save atmospheric fields to force surface models
    call t_startf('cam_export')
    call cam_export (state,cam_out,pbuf)
    call t_stopf('cam_export')

    ! Write export state to history file
    call t_startf('diag_export')
    call diag_export(cam_out)
    call t_stopf('diag_export')

    call check_tracers_fini(tracerint)

    if(print_additional_diagn) then
       call additional_diagn_after_step(state, cam_in, cam_out, fsns, fsnt, flns, flnt, ztodt)
    endif

    call cnd_diag_checkpoint( diag, 'PBCEND', state, pbuf, cam_in, cam_out )

end subroutine tphysbc

subroutine phys_timestep_init(phys_state, cam_out, pbuf2d)
!-----------------------------------------------------------------------------------
!
! Purpose: The place for parameterizations to call per timestep initializations.
!          Generally this is used to update time interpolated fields from boundary
!          datasets.
!
!-----------------------------------------------------------------------------------
  use shr_kind_mod,        only: r8 => shr_kind_r8
  use chemistry,           only: chem_timestep_init
  use chem_surfvals,       only: chem_surfvals_set
  use physics_types,       only: physics_state
  use physics_buffer,      only: physics_buffer_desc
  use ghg_data,            only: ghg_data_timestep_init
  use cam3_aero_data,      only: cam3_aero_data_on, cam3_aero_data_timestep_init
  use cam3_ozone_data,     only: cam3_ozone_data_on, cam3_ozone_data_timestep_init
  use radiation,           only: radiation_do
  use tracers,             only: tracers_timestep_init
  use aoa_tracers,         only: aoa_tracers_timestep_init
  use vertical_diffusion,  only: vertical_diffusion_ts_init
  use radheat,             only: radheat_timestep_init
  use solar_data,          only: solar_data_advance
  use qbo,                 only: qbo_timestep_init
  use efield,              only: get_efield
  use iondrag,             only: do_waccm_ions
  use perf_mod

  use prescribed_ozone,    only: prescribed_ozone_adv
  use prescribed_ghg,      only: prescribed_ghg_adv
  use prescribed_aero,     only: prescribed_aero_adv
  use read_spa_data,       only: read_spa_data_adv
  use aerodep_flx,         only: aerodep_flx_adv
  use aircraft_emit,       only: aircraft_emit_adv
  use prescribed_volcaero, only: prescribed_volcaero_adv
  use nudging,             only: Nudge_Model,nudging_timestep_init

  use seasalt_model,       only: advance_ocean_data, has_mam_mom
  use phys_grid_ctem,      only: phys_grid_ctem_diags

  implicit none

  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
  
  type(physics_buffer_desc), pointer                 :: pbuf2d(:,:)

  !-----------------------------------------------------------------------------

  ! Chemistry surface values
  call chem_surfvals_set()

  ! Solar irradiance
  call solar_data_advance()

  ! Time interpolate for chemistry.
  call t_startf('chem_timestep_init')
  call chem_timestep_init(phys_state, pbuf2d)
  call t_stopf('chem_timestep_init')

  ! Prescribed tracers
  call prescribed_ozone_adv(phys_state, pbuf2d)
  call prescribed_ghg_adv(phys_state, pbuf2d)
  call prescribed_aero_adv(phys_state, pbuf2d)
  call read_spa_data_adv(phys_state, pbuf2d)
  call aircraft_emit_adv(phys_state, pbuf2d)

  call t_startf('prescribed_volcaero_adv')
  call prescribed_volcaero_adv(phys_state, pbuf2d)
  call t_stopf('prescribed_volcaero_adv')

  if (has_mam_mom) then
     call t_startf('advance_ocean_data')
     call advance_ocean_data(phys_state, pbuf2d)
     call t_stopf('advance_ocean_data')
  end if

  ! prescribed aerosol deposition fluxes
  call aerodep_flx_adv(phys_state, pbuf2d, cam_out)

  ! CAM3 prescribed aerosol masses
  if (cam3_aero_data_on) call cam3_aero_data_timestep_init(pbuf2d,  phys_state)

  ! CAM3 prescribed ozone data
  if (cam3_ozone_data_on) call cam3_ozone_data_timestep_init(pbuf2d,  phys_state)

  ! Time interpolate data models of gasses in pbuf2d
  call ghg_data_timestep_init(pbuf2d,  phys_state)

  ! Upper atmosphere radiative processes
  call radheat_timestep_init(phys_state, pbuf2d)
 
  ! Time interpolate for vertical diffusion upper boundary condition
  call vertical_diffusion_ts_init(pbuf2d, phys_state)

  !----------------------------------------------------------------------
  ! update QBO data for this time step
  !----------------------------------------------------------------------
  call qbo_timestep_init

  if (do_waccm_ions) then
     ! Compute the electric field
     call t_startf ('efield')
     call get_efield
     call t_stopf ('efield')
  endif

  ! Time interpolate for tracers, if appropriate
  call tracers_timestep_init(phys_state)

  ! age of air tracers
  call aoa_tracers_timestep_init(phys_state)

  ! Update Nudging values, if needed
  !----------------------------------
  if(Nudge_Model) call nudging_timestep_init(phys_state)

  ! Update Transformed Eularian Mean (TEM) diagnostics
  call phys_grid_ctem_diags(phys_state)

end subroutine phys_timestep_init


subroutine add_fld_default_calls()
  !BSINGH -  For adding addfld and add defualt calls
  use cam_history,        only: addfld, add_default, fieldname_len

  implicit none

  !Add all existing ptend names for the addfld calls
  character(len=20), parameter :: vlist(28) = (/     'topphysbc           '                       ,&
       'chkenergyfix        ','dadadj              ','zm_convr            ','zm_conv_evap        ',&
       'momtran             ','zm_conv_tend        ','UWSHCU              ','convect_shallow     ',&
       'pcwdetrain_mac      ','macro_park          ','macrop              ','micro_mg            ',&
       'cldwat_mic          ','aero_model_wetdep_ma','convtran2           ','cam_radheat         ',&
       'chemistry           ','vdiff               ','rayleigh_friction   ','aero_model_drydep_ma',&
       'Grav_wave_drag      ','convect_shallow_off ','clubb_ice1          ','clubb_det           ',&
       'clubb_ice4          ','clubb_srf           ','micro_p3            ' /)



  character(len=fieldname_len) :: varname
  character(len=1000)          :: s_lngname,stend_lngname,qv_lngname,qvtend_lngname,t_lngname
  
  integer :: iv, ntot, ihist

  ntot = size(vlist)

  do ihist = 1 , nvars_prtrb_hist
     do iv = 1, ntot   
        
        varname  = trim(adjustl(hist_vars(ihist)))//'_'//trim(adjustl(vlist(iv))) ! form variable name

        call addfld (trim(adjustl(varname)), (/ 'lev' /), 'A', 'prg_test_units', 'pergro_longname',flag_xyfill=.true.)!The units and longname are dummy as it is for a test only
        call add_default (trim(adjustl(varname)), 1, ' ')        
     enddo
  enddo

end subroutine add_fld_default_calls

end module physpkg
