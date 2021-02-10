module physpkg
  !-----------------------------------------------------------------------------
  ! Purpose: Provides the interface to MMF physics package
  !-----------------------------------------------------------------------------
#ifdef SPMD
  use mpishorthand
#endif
  use shr_kind_mod,       only: i8 => shr_kind_i8, r8 => shr_kind_r8
  use shr_sys_mod,        only: shr_sys_irtc
  use shr_sys_mod,        only: shr_sys_flush
  use spmd_utils,         only: masterproc, iam
  use physconst,          only: latvap, latice, rh2o
  use physics_types,      only: physics_state, physics_tend, physics_state_set_grid, &
                                physics_ptend, physics_tend_init, physics_type_alloc, physics_ptend_dealloc,&
                                physics_state_alloc, physics_state_dealloc, physics_tend_alloc, physics_tend_dealloc
  use physics_update_mod, only: physics_update, physics_update_init, hist_vars, nvars_prtrb_hist, get_var
  use phys_grid,          only: get_ncols_p, print_cost_p, update_cost_p
  use phys_gmean,         only: gmean_mass
  use ppgrid,             only: begchunk, endchunk, pcols, pver, pverp
  use constituents,       only: pcnst, cnst_name, cnst_get_ind!, species_class
  use camsrfexch,         only: cam_out_t, cam_in_t
  use phys_control,       only: phys_do_flux_avg, phys_getopts
  use zm_conv,            only: do_zmconv_dcape_ull => trigdcape_ull
  use scamMod,            only: single_column, scm_crm_mode
  use flux_avg,           only: flux_avg_init
  use cam_logfile,        only: iulog
  use camsrfexch,         only: cam_export
  use perf_mod
  use modal_aero_calcsize,    only: modal_aero_calcsize_init, &
                                    modal_aero_calcsize_reg, &
                                    modal_aero_calcsize_sub
  use modal_aero_wateruptake, only: modal_aero_wateruptake_init, &
                                    modal_aero_wateruptake_dr, &
                                    modal_aero_wateruptake_reg
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

  integer :: species_class(pcnst)  = -1 

  save
  !-----------------------------------------------------------------------------
  ! Public methods
  !-----------------------------------------------------------------------------
  public phys_register ! was initindx  - register physics methods
  public phys_init   ! Public initialization method
  public phys_run1   ! First phase of the public run method
  public phys_run2   ! Second phase of the public run method
  public phys_final  ! Public finalization method
  !-----------------------------------------------------------------------------
  ! Private module data
  !-----------------------------------------------------------------------------
  ! Physics package options
  character(len=16) :: shallow_scheme
  character(len=16) :: macrop_scheme
  character(len=16) :: microp_scheme 
  integer           :: cld_macmic_num_steps    ! Number of macro/micro substeps
  logical           :: state_debug_checks  ! Debug physics_state.
  logical           :: clim_modal_aero     ! climate controled by prognostic or prescribed modal aerosols
  logical           :: prog_modal_aero     ! Prognostic modal aerosols present
  logical           :: micro_do_icesupersat
  logical           :: pergro_test_active= .false.
  logical           :: pergro_mods = .false.
  logical           :: is_cmip6_volc       ! true if cmip6 style volcanic file is read otherwise false

contains

!===================================================================================================
!===================================================================================================

subroutine phys_register
  !----------------------------------------------------------------------- 
  ! Purpose: Register constituents and physics buffer fields.
  !-----------------------------------------------------------------------
  use physics_buffer,     only: pbuf_init_time
  use physics_buffer,     only: pbuf_add_field, dtype_r8
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
  use conv_water,         only: conv_water_register
  use physconst,          only: mwdry, cpair, mwh2o, cpwv
  use tracers,            only: tracers_register
  use check_energy,       only: check_energy_register
  use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_register
  use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_register
  use ghg_data,           only: ghg_data_register
  use mmf_vertical_diffusion, only: mmf_vertical_diffusion_register
  use convect_deep,       only: convect_deep_register
  use convect_shallow,    only: convect_shallow_register
  use radiation,          only: radiation_register
  use co2_cycle,          only: co2_register
  use flux_avg,           only: flux_avg_register
  use ionosphere,         only: ionos_register
  use string_utils,       only: to_lower
  use prescribed_ozone,   only: prescribed_ozone_register
  use prescribed_volcaero,only: prescribed_volcaero_register
  use prescribed_aero,    only: prescribed_aero_register
  use prescribed_ghg,     only: prescribed_ghg_register
  use sslt_rebin,         only: sslt_rebin_register
  use aoa_tracers,        only: aoa_tracers_register
  use aircraft_emit,      only: aircraft_emit_register
  use cam_diagnostics,    only: diag_register
  use cloud_diagnostics,  only: cloud_diagnostics_register
  use cospsimulator_intr, only: cospsimulator_intr_register
  use rad_constituents,   only: rad_cnst_get_info ! Added to query if it is a modal aero sim or not
  use output_aerocom_aie, only: output_aerocom_aie_register, do_aerocom_ind3
  use crm_physics,        only: crm_physics_register
  !---------------------------------------------------------------------------
  ! Local variables
  !---------------------------------------------------------------------------
  integer :: m,mm     ! loop and constituent indices
  integer :: dummy    ! for unused output from pbuf_add_field calls
  integer :: nmodes
  character(len=16) :: MMF_microphysics_scheme
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  call phys_getopts(cld_macmic_num_steps_out = cld_macmic_num_steps, &
                    microp_scheme_out        = microp_scheme,   &
                    do_aerocom_ind3_out      = do_aerocom_ind3,  &
                    state_debug_checks_out   = state_debug_checks, &
                    micro_do_icesupersat_out = micro_do_icesupersat, &
                    pergro_test_active_out   = pergro_test_active, &
                    pergro_mods_out          = pergro_mods)
  call phys_getopts( MMF_microphysics_scheme_out = MMF_microphysics_scheme)

  ! Initialize dyn_time_lvls
  call pbuf_init_time()

  ! Register water vapor.
  ! This must be the first call to cnst_add so that water vapor is constituent 1.

  call cnst_add('Q', mwh2o, cpwv, 1.E-12_r8, mm, longname='Specific humidity', readiv=.true., is_convtran1=.true.)
  

  ! Fields for physics package diagnostics
  call pbuf_add_field('TINI',      'physpkg', dtype_r8, (/pcols,pver/), tini_idx)
  call pbuf_add_field('QINI',      'physpkg', dtype_r8, (/pcols,pver/), qini_idx)
  call pbuf_add_field('CLDLIQINI', 'physpkg', dtype_r8, (/pcols,pver/), cldliqini_idx)
  call pbuf_add_field('CLDICEINI', 'physpkg', dtype_r8, (/pcols,pver/), cldiceini_idx)
  call pbuf_add_field('static_ener_ac', 'global', dtype_r8, (/pcols/), static_ener_ac_idx)
  call pbuf_add_field('water_vap_ac',   'global', dtype_r8, (/pcols/), water_vap_ac_idx)

  ! check energy package
  call check_energy_register

  ! register fluxes for saving across time
  if (phys_do_flux_avg()) call flux_avg_register()

  call cldfrc_register()

  ! ! cloud water
  ! if( microp_scheme == 'RK' ) then
  !    call stratiform_register()
  ! elseif ( microp_scheme == 'MG' ) then
  !    call macrop_driver_register()
  !    call microp_aero_register()
  !    call microp_driver_register()
  ! end if
   
     
  ! vmag_gust must be registered even if clubb is not used
  call pbuf_add_field('vmag_gust', 'global', dtype_r8, (/pcols/), dummy)
  call pbuf_add_field('PREC_STR',  'physpkg',dtype_r8,(/pcols/),prec_str_idx)
  call pbuf_add_field('SNOW_STR',  'physpkg',dtype_r8,(/pcols/),snow_str_idx)
  call pbuf_add_field('PREC_PCW',  'physpkg',dtype_r8,(/pcols/),prec_pcw_idx)
  call pbuf_add_field('SNOW_PCW',  'physpkg',dtype_r8,(/pcols/),snow_pcw_idx)
  call pbuf_add_field('PREC_SED',  'physpkg',dtype_r8,(/pcols/),prec_sed_idx)
  call pbuf_add_field('SNOW_SED',  'physpkg',dtype_r8,(/pcols/),snow_sed_idx)

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

  ! co2 constituents
  call co2_register()

  ! register data model ozone with pbuf
  if (cam3_ozone_data_on) call cam3_ozone_data_register()

  call prescribed_volcaero_register()
  call prescribed_ozone_register()
  call prescribed_aero_register()
  call prescribed_ghg_register()
  call sslt_rebin_register

  ! CAM3 prescribed aerosols
  if (cam3_aero_data_on) call cam3_aero_data_register()

  ! register various data model gasses with pbuf
  call ghg_data_register()

  call aircraft_emit_register()

  ! ! deep convection
  ! call convect_deep_register

  ! !  shallow convection
  ! call convect_shallow_register

  call crm_physics_register

  ! radiation
  call radiation_register
  call cloud_diagnostics_register

  ! COSP
  call cospsimulator_intr_register

  ! vertical diffusion
  call mmf_vertical_diffusion_register()

  if (do_aerocom_ind3) call output_aerocom_aie_register()

  ! Register diagnostics PBUF
  call diag_register()

  ! Register age of air tracers
  call aoa_tracers_register()

  ! Register test tracers
  ! This is the last call to register constituents because
  ! the test tracers fill the remaining available slots up
  ! to constituent number PCNST -- regardless of what PCNST is set to.
  call tracers_register()

  ! All tracers registered, check that the dimensions are correct
  call cnst_chk_dim()

  ! ***NOTE*** No registering constituents after the call to cnst_chk_dim.

end subroutine phys_register

!===================================================================================================
!===================================================================================================

subroutine phys_inidat( cam_out, pbuf2d )
  use cam_abortutils, only : endrun
  use physics_buffer, only : pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_set_field, dyn_time_lvls
  use cam_initfiles,       only: initial_file_get_id, topo_file_get_id
  use cam_grid_support,    only: cam_grid_check, cam_grid_id
  use cam_grid_support,    only: cam_grid_get_dim_names
  use pio,                 only: file_desc_t
  use ncdio_atm,           only: infld
  use polar_avg,           only: polar_average
  use short_lived_species, only: initialize_short_lived_species
  use comsrf,              only: landm, sgh, sgh30
  use cam_control_mod,     only: aqua_planet
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  type(cam_out_t),     intent(inout) :: cam_out(begchunk:endchunk)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  integer :: lchnk, m, n, i, k, ncol
  type(file_desc_t), pointer :: fh_ini, fh_topo
  character(len=8) :: fieldname
  real(r8), pointer :: cldptr(:,:,:,:), convptr_3d(:,:,:,:)
  real(r8), pointer :: tptr(:,:), tptr3d(:,:,:), tptr3d_2(:,:,:)
  real(r8), pointer :: qpert(:,:)

  character*11 :: subname='phys_inidat' ! subroutine name
  character(len=8) :: dim1name, dim2name
  integer :: tpert_idx, qpert_idx, pblh_idx, vmag_gust_idx
  logical :: found=.false., found2=.false.
  integer :: ierr
  integer :: ixcldice, ixcldliq
  integer :: grid_id  ! grid ID for data mapping
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------
  nullify(tptr,tptr3d,tptr3d_2,cldptr,convptr_3d)

  fh_ini=>initial_file_get_id()

  grid_id = cam_grid_id('physgrid')
  if (.not. cam_grid_check(grid_id)) then
    call endrun(trim(subname)//': Internal error, no "physgrid" grid')
  end if
  call cam_grid_get_dim_names(grid_id, dim1name, dim2name)

  if(aqua_planet) then
    sgh = 0._r8
    sgh30 = 0._r8
    landm = 0._r8
    if (masterproc) write(iulog,*) 'AQUA_PLANET simulation, sgh, sgh30, landm initialized to 0.'
  else    
    if (masterproc) write(iulog,*) 'NOT AN AQUA_PLANET simulation, initialize sgh, sgh30, land m using data from file.'
    fh_topo=>topo_file_get_id()
    call infld('SGH', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, sgh, found, gridname='physgrid')
    if(.not. found) call endrun('ERROR: SGH not found on topo file')

    call infld('SGH30', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, sgh30, found, gridname='physgrid')
    if(.not. found) then
      if (masterproc) write(iulog,*) 'Warning: Error reading SGH30 from topo file.'
      if (masterproc) write(iulog,*) 'The field SGH30 will be filled using data from SGH.'
      sgh30 = sgh
    end if

    call infld('LANDM_COSLAT', fh_topo, dim1name, dim2name, 1, pcols, begchunk, endchunk, landm, found, gridname='physgrid')
    if(.not.found) call endrun(' ERROR: LANDM_COSLAT not found on topo dataset.')
  end if

  allocate(tptr(1:pcols,begchunk:endchunk))

  call infld('PBLH', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr(:,:), found, gridname='physgrid')
  if(.not. found) then
    tptr(:,:) = 0._r8
    if (masterproc) write(iulog,*) 'PBLH initialized to 0.'
  end if
  pblh_idx = pbuf_get_index('pblh')

  call pbuf_set_field(pbuf2d, pblh_idx, tptr)

  call infld('TPERT', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr(:,:), found, gridname='physgrid')
  if(.not. found) then
    tptr(:,:) = 0._r8
    if (masterproc) write(iulog,*) 'TPERT initialized to 0.'
  end if
  tpert_idx = pbuf_get_index( 'tpert')
  call pbuf_set_field(pbuf2d, tpert_idx, tptr)


  call infld('vmag_gust', fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr(:,:), found, gridname='physgrid')
  if(.not. found) then
    tptr(:,:) = 0._r8
    if (masterproc) write(iulog,*) 'vmag_gust initialized to 1.'
  end if
  vmag_gust_idx = pbuf_get_index( 'vmag_gust')
  call pbuf_set_field(pbuf2d, vmag_gust_idx, tptr)


  fieldname='QPERT'  
  qpert_idx = pbuf_get_index( 'qpert',ierr)
  if (qpert_idx > 0) then
    call infld(fieldname, fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr, found, gridname='physgrid')
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

  ! fieldname='CUSH'
  ! m = pbuf_get_index('cush')
  ! call infld(fieldname, fh_ini, dim1name, dim2name, 1, pcols, begchunk, endchunk, tptr, found, gridname='physgrid')
  ! if(.not.found) then
  !   if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 1000.'
  !   tptr=1000._r8
  ! end if
  ! do n=1,dyn_time_lvls
  !   call pbuf_set_field(pbuf2d, m, tptr, start=(/1,n/), kount=(/pcols,1/))
  ! end do
  ! deallocate(tptr)

  !
  ! 3-D fields
  !

  allocate(tptr3d(pcols,pver,begchunk:endchunk))

  fieldname='CLOUD'
  m = pbuf_get_index('CLD')
  call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
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
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
    if(.not. found) then
      call infld('Q',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
      if (found) then
        if (masterproc) write(iulog,*) trim(fieldname), ' initialized with Q'
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
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
    if(found) then
      do n = 1, dyn_time_lvls
        call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
      end do
    else
      call cnst_get_ind('CLDICE', ixcldice)
      call infld('CLDICE',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
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
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
    if(found) then
      do n = 1, dyn_time_lvls
        call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
      end do
    else
      allocate(tptr3d_2(pcols,pver,begchunk:endchunk))     
      call cnst_get_ind('CLDICE', ixcldice)
      call cnst_get_ind('CLDLIQ', ixcldliq)
      call infld('CLDICE',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
      call infld('CLDLIQ',fh_ini,dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d_2, found2, gridname='physgrid')
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
    call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
    if(.not.found) then
      call infld('T', fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
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
  call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, tptr3d, found, gridname='physgrid')
  if (found) then
    call pbuf_set_field(pbuf2d, m, tptr3d)
  else
    call pbuf_set_field(pbuf2d, m, 0.01_r8)
    if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.01'
  end if


  fieldname = 'KVM'
  m = pbuf_get_index('kvm')
  call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, tptr3d, found, gridname='physgrid')
  if (found) then
    call pbuf_set_field(pbuf2d, m, tptr3d)
  else
    call pbuf_set_field(pbuf2d, m, 0._r8)
    if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
  end if


  fieldname = 'KVH'
  m = pbuf_get_index('kvh')
  call infld(fieldname, fh_ini, dim1name, 'ilev', dim2name, 1, pcols, 1, pverp, begchunk, endchunk, tptr3d, found, gridname='physgrid')
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
  call infld(fieldname, fh_ini, dim1name, 'lev', dim2name, 1, pcols, 1, pver, begchunk, endchunk, tptr3d, found, gridname='physgrid')
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

!===================================================================================================
!===================================================================================================

subroutine phys_init( phys_state, phys_tend, pbuf2d, cam_out )
  !-----------------------------------------------------------------------------
  ! Purpose: Initialization of physics package
  !-----------------------------------------------------------------------------
  use physics_buffer,     only: physics_buffer_desc, pbuf_initialize, pbuf_get_index
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
  use mmf_vertical_diffusion, only: mmf_vertical_diffusion_init
  use phys_debug_util,    only: phys_debug_init
  use rad_constituents,   only: rad_cnst_init
  use aer_rad_props,      only: aer_rad_props_init
  use clubb_intr,         only: clubb_ini_cam
  use sslt_rebin,         only: sslt_rebin_init
  use tropopause,         only: tropopause_init
  use solar_data,         only: solar_data_init
  use rad_solar_var,      only: rad_solar_var_init
  use output_aerocom_aie, only: output_aerocom_aie_init, do_aerocom_ind3
  use dyn_grid,           only: fv_nphys
  use cam_history,        only: addfld, add_default, horiz_only 
  use crm_physics,        only: crm_physics_init 
  !-----------------------------------------------------------------------------
  ! Input/output arguments
  !-----------------------------------------------------------------------------
  type(physics_state), pointer       :: phys_state(:)
  type(physics_tend ), pointer       :: phys_tend(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  type(cam_out_t),intent(inout)      :: cam_out(begchunk:endchunk)
  !-----------------------------------------------------------------------------
  ! local variables
  !-----------------------------------------------------------------------------
  integer :: lchnk
  real(r8) :: dp1 = huge(1.0_r8) !set in namelist, assigned in cloud_fraction.F90
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  call physics_type_alloc(phys_state, phys_tend, begchunk, endchunk, pcols)

  do lchnk = begchunk, endchunk
     call physics_state_set_grid(lchnk, phys_state(lchnk))
  end do

  !-----------------------------------------------------------------------------
  ! Initialize variables in physconst which are not temporally and/or spatially constant
  !----------------------------------------------------------------------------- 
  call physconst_init()

  ! Initialize debugging a physics column
  call phys_debug_init()

  call pbuf_initialize(pbuf2d)

  !initialize physics update interface routine
  call physics_update_init()

  ! diag_init makes addfld calls for dyn fields that are output from the physics decomp
  call diag_init()

  call check_energy_init()

  call tracers_init()

  ! age of air tracers
  call aoa_tracers_init()

  teout_idx = pbuf_get_index( 'TEOUT')

  if (nsrest .eq. 0) then
    call phys_inidat(cam_out, pbuf2d) 
  end if
  
  ! wv_saturation is relatively independent of everything else and
  ! low level, so init it early. Must at least do this before radiation.
  call wv_sat_init

  ! CAM3 prescribed aerosols
  if (cam3_aero_data_on) call cam3_aero_data_init(phys_state)

  ! Initialize rad constituents and their properties
  call rad_cnst_init()
  call aer_rad_props_init()
  call cloud_rad_props_init()

  ! solar irradiance data modules
  call solar_data_init()

  ! Prognostic chemistry.
  call chem_init(phys_state,pbuf2d, species_class)

  ! Prescribed tracers
  call prescribed_ozone_init()
  call prescribed_ghg_init()
  call prescribed_aero_init()
  call aerodep_flx_init()
  call aircraft_emit_init()
  !when is_cmip6_volc is true ,cmip6 style volcanic file is read
  !Initialized to .false. here but it gets its values from prescribed_volcaero_init
  is_cmip6_volc = .false. 
  call prescribed_volcaero_init(is_cmip6_volc)

  ! co2 cycle            
  if (co2_transport()) then
     call co2_init(phys_state, pbuf2d)
  end if

  ! CAM3 prescribed ozone
  if (cam3_ozone_data_on) call cam3_ozone_data_init(phys_state)

  call gw_init()

  call rayleigh_friction_init()

  call pbl_utils_init(gravit, karman, cpair, rair, zvir)
  call mmf_vertical_diffusion_init(pbuf2d)

  call tsinti(tmelt, latvap, rair, stebol, latice)

  call radiation_init(phys_state)

  call rad_solar_var_init()

  call cloud_diagnostics_init()

  call radheat_init(pref_mid)

  if (use_MMF) call crm_physics_init( phys_state, species_class )

  call cldfrc_init(dp1)! for passing dp1 on to clubb
  call cldfrc2m_init()

  ! call convect_deep_init(pref_edge)
  ! if( microp_scheme == 'RK' ) then
  !    call stratiform_init()
  ! elseif( microp_scheme == 'MG' ) then 
  !    call macrop_driver_init(pbuf2d)
  !    call microp_aero_init()
  !    call microp_driver_init(pbuf2d)
     call conv_water_init
  ! end if

  call crm_physics_init(pbuf2d, species_class)

  call sslt_rebin_init()
  call tropopause_init()

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
  
 !BSINGH - addfld and adddefault calls for perturb growth testing    
  if(pergro_test_active)call add_fld_default_calls()

end subroutine phys_init

!===================================================================================================
!===================================================================================================

subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out)
  !-----------------------------------------------------------------------------
  ! Purpose: First part of atmos physics before updating of surface components
  !-----------------------------------------------------------------------------
#if (defined E3SM_SCM_REPLAY )
  use cam_history,    only: outfld
#endif
  use time_manager,     only: get_nstep
  use cam_diagnostics,  only: diag_allocate, diag_physvar_ic
  use check_energy,     only: check_energy_gmean
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_chunk, pbuf_allocate
  use comsrf,           only: fsns, fsnt, flns, sgh, sgh30, flnt, landm, fsds
  use cam_abortutils,   only: endrun

  !-----------------------------------------------------------------------------
  ! Interface arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in) :: ztodt            ! physics time step unless nstep=0
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend

  type(physics_buffer_desc), pointer, dimension(:,:) :: pbuf2d
  type(cam_in_t),                     dimension(begchunk:endchunk) :: cam_in
  type(cam_out_t),                    dimension(begchunk:endchunk) :: cam_out
  !-----------------------------------------------------------------------------
  ! Local Variables
  !-----------------------------------------------------------------------------
  integer :: c                                 ! indices
  integer :: ncol                              ! number of columns
  integer :: nstep                             ! current timestep number
#if (! defined SPMD)
  integer :: mpicom = 0
#endif
  integer(i8) :: beg_count                     ! start time for a chunk
  integer(i8) :: end_count                     ! stop time for a chunk
  integer(i8) :: irtc_rate                     ! irtc clock rate
  real(r8)    :: chunk_cost                    ! measured cost per chunk
  type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call t_startf ('physpkg_st1')
  nstep = get_nstep()

  ! The following initialization depends on the import state (cam_in)
  ! being initialized.  This isn't true when cam_init is called, so need
  ! to postpone this initialization to here.
  if (nstep == 0 .and. phys_do_flux_avg()) call flux_avg_init(cam_in,  pbuf2d)

  ! Compute total energy of input state and previous output state
  call t_startf ('chk_en_gmean')
  call check_energy_gmean(phys_state, pbuf2d, ztodt, nstep)
  call t_stopf ('chk_en_gmean')

  call t_stopf ('physpkg_st1')

  call t_startf ('physpkg_st1')

  call pbuf_allocate(pbuf2d, 'physpkg')
  call diag_allocate()

  !-----------------------------------------------------------------------------
  ! Advance time information
  !-----------------------------------------------------------------------------

  call phys_timestep_init( phys_state, cam_out, pbuf2d)

  call t_stopf ('physpkg_st1')

#ifdef TRACER_CHECK
  call gmean_mass ('before tphysbc DRY', phys_state)
#endif

  !-----------------------------------------------------------------------------
  ! Tendency physics before flux coupler invocation
  !-----------------------------------------------------------------------------

  call t_barrierf('sync_bc_physics', mpicom)
  call t_startf ('bc_physics')

!$OMP PARALLEL DO PRIVATE (C, beg_count, phys_buffer_chunk, end_count, chunk_cost)
  do c=begchunk, endchunk

    beg_count = shr_sys_irtc(irtc_rate)

    ! Output physics terms to IC file
    phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

    call t_startf ('diag_physvar_ic')
    call diag_physvar_ic ( c,  phys_buffer_chunk, cam_out(c), cam_in(c) )
    call t_stopf ('diag_physvar_ic')

    call tphysbc(ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), &
                 phys_state(c), phys_tend(c), phys_buffer_chunk, &
                 fsds(1,c), landm(1,c), sgh(1,c), sgh30(1,c), &
                 cam_out(c), cam_in(c) )

    end_count = shr_sys_irtc(irtc_rate)
    chunk_cost = real( (end_count-beg_count), r8)/real(irtc_rate, r8)
    call update_cost_p(c, chunk_cost)

  end do

  call t_stopf ('bc_physics')

  if(single_column.and.scm_crm_mode) return ! Don't call the rest in CRM mode

#ifdef TRACER_CHECK
  call gmean_mass ('between DRY', phys_state)
#endif

end subroutine phys_run1

!===================================================================================================
!===================================================================================================

subroutine phys_run2(phys_state, ztodt, phys_tend, pbuf2d, cam_out, cam_in )
  !-----------------------------------------------------------------------------
  ! Purpose: Second part of atmos physics after updating of surface models
  !-----------------------------------------------------------------------------
  use physics_buffer,   only: physics_buffer_desc, pbuf_get_chunk, pbuf_deallocate, pbuf_update_tim_idx
  use cam_diagnostics,  only: diag_deallocate, diag_surf
  use comsrf,           only: trefmxav, trefmnav, sgh, sgh30, fsds 
  use physconst,        only: stebol, latvap
  use time_manager,     only: get_nstep
  use check_energy,     only: ieflx_gmean, check_ieflx_fix 
  use phys_control,     only: ieflx_opt
  !-----------------------------------------------------------------------------
  ! Interface arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in) :: ztodt                       ! physics time step unless nstep=0
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(physics_tend ), intent(inout), dimension(begchunk:endchunk) :: phys_tend
  type(physics_buffer_desc),pointer, dimension(:,:)     :: pbuf2d
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out
  type(cam_in_t),      intent(inout), dimension(begchunk:endchunk) :: cam_in
  !-----------------------------------------------------------------------------
  ! Local Variables
  !-----------------------------------------------------------------------------
  integer :: c                                 ! chunk index
  integer :: ncol                              ! number of columns
  integer :: nstep                             ! current timestep number
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
#if (! defined SPMD)
  integer :: mpicom = 0
#endif
  integer(i8) :: beg_count, end_count, irtc_rate ! for measuring chunk cost
  real(r8):: chunk_cost
  type(physics_buffer_desc),pointer, dimension(:)     :: phys_buffer_chunk

  if(single_column.and.scm_crm_mode) return

  call t_barrierf('sync_ac_physics', mpicom)
  call t_startf ('ac_physics')

  nstep = get_nstep()

  !-----------------------------------------------------------------------------
  ! calculate the global mean ieflx 
  !-----------------------------------------------------------------------------
  if(ieflx_opt>0) call ieflx_gmean(phys_state, phys_tend, pbuf2d, cam_in, cam_out, nstep)

!$OMP PARALLEL DO PRIVATE (C, beg_count, NCOL, phys_buffer_chunk, end_count, chunk_cost)
  do c=begchunk,endchunk

    beg_count = shr_sys_irtc(irtc_rate)

    ncol = get_ncols_p(c)
    phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

    ! add the implied internal energy flux to sensible heat flux
    if(ieflx_opt>0) call check_ieflx_fix(c, ncol, nstep, cam_in(c)%shf)

    ! surface diagnostics for history files
    call t_startf('diag_surf')
    call diag_surf(cam_in(c), cam_out(c), phys_state(c)%ps,trefmxav(1,c), trefmnav(1,c))
    call t_stopf('diag_surf')

    call tphysac(ztodt, cam_in(c),  &
        sgh(1,c), sgh30(1,c), cam_out(c),                              &
        phys_state(c), phys_tend(c), phys_buffer_chunk,&
        fsds(1,c))

    end_count = shr_sys_irtc(irtc_rate)
    chunk_cost = real( (end_count-beg_count), r8)/real(irtc_rate, r8)
    call update_cost_p(c, chunk_cost)

  end do ! Chunk loop

  call t_stopf('ac_physics')

#ifdef TRACER_CHECK
  call gmean_mass ('after tphysac FV:WET)', phys_state)
#endif

  call t_startf ('physpkg_st2')
  call pbuf_deallocate(pbuf2d, 'physpkg')

  call pbuf_update_tim_idx()
  call diag_deallocate()
  call t_stopf ('physpkg_st2')

end subroutine phys_run2

!===================================================================================================
!===================================================================================================

subroutine phys_final( phys_state, phys_tend, pbuf2d )
  use physics_buffer, only : physics_buffer_desc, pbuf_deallocate
  use chemistry, only : chem_final
  use wv_saturation, only : wv_sat_final
  !----------------------------------------------------------------------- 
  ! Purpose: Finalization of physics package
  !-----------------------------------------------------------------------
  ! Input/output arguments
  type(physics_state), pointer :: phys_state(:)
  type(physics_tend ), pointer :: phys_tend(:)
  type(physics_buffer_desc), pointer :: pbuf2d(:,:)
  !---------------------------------------------------------------------------
  !---------------------------------------------------------------------------

  if(associated(pbuf2d)) then
     call pbuf_deallocate(pbuf2d,'global')
     deallocate(pbuf2d)
  end if
  deallocate(phys_state)
  deallocate(phys_tend)

  call t_startf ('chem_final')
  call chem_final
  call t_stopf ('chem_final')

  call t_startf ('wv_sat_final')
  call wv_sat_final
  call t_stopf ('wv_sat_final')

  call t_startf ('print_cost_p')
  call print_cost_p
  call t_stopf ('print_cost_p')

end subroutine phys_final

!===================================================================================================
!===================================================================================================

subroutine tphysac (ztodt, cam_in, sgh, sgh30, cam_out, state, tend, pbuf, fsds )
  !----------------------------------------------------------------------------- 
  ! Purpose: 
  !   Tendency physics after coupling to land, sea, and ice models.
  !   Computes the following:
  !     o Radon surface flux and decay (optional)
  !     o Vertical diffusion and planetary boundary layer
  !     o Multiple gravity wave drag
  !-----------------------------------------------------------------------------
  use physics_buffer,     only: physics_buffer_desc, pbuf_set_field, pbuf_get_index, pbuf_get_field, pbuf_old_tim_idx
  use shr_kind_mod,       only: r8 => shr_kind_r8
  use chemistry,          only: chem_is_active, chem_timestep_tend, chem_emissions
  use cam_diagnostics,    only: diag_phys_tend_writeout, diag_conv
  use gw_drag,            only: gw_tend
  use mmf_vertical_diffusion, only: mmf_vertical_diffusion_tend
  use mmf_surface_mod,    only: mmf_surface_ac
  use rayleigh_friction,  only: rayleigh_friction_tend
  use constituents,       only: cnst_get_ind
  use physics_types,      only: physics_state, physics_tend, physics_ptend, &
                                physics_dme_adjust, set_dry_to_wet, physics_state_check
  use tracers,            only: tracers_timestep_tend
  use aoa_tracers,        only: aoa_tracers_timestep_tend
  use physconst,          only: rhoh2o, latvap,latice, rga
  use aero_model,         only: aero_model_drydep
  use check_energy,       only: check_energy_chng, check_prect, &
                                check_tracers_data, check_tracers_init, &
                                check_tracers_chng, check_tracers_fini
  use time_manager,       only: get_nstep
  use cam_abortutils,     only: endrun
  use cam_control_mod,    only: aqua_planet 
  use mo_gas_phase_chemdr,only: map2chm
  use clybry_fam,         only: clybry_fam_set
  use charge_neutrality,  only: charge_fix
  use perf_mod
  use flux_avg,           only: flux_avg_run
  use phys_control,       only: use_qqflx_fixer
  use cam_history,        only: outfld 

  implicit none
  !-----------------------------------------------------------------------------
  ! Arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in) :: ztodt                  ! Two times model timestep (2 delta-t)
  real(r8), intent(in) :: fsds(pcols)            ! down solar flux
  real(r8), intent(in) :: sgh(pcols)             ! Std. deviation of orography for gwd
  real(r8), intent(in) :: sgh30(pcols)           ! Std. deviation of 30s orography for tms

  type(cam_in_t),      intent(inout) :: cam_in
  type(cam_out_t),     intent(inout) :: cam_out
  type(physics_state), intent(inout) :: state
  type(physics_tend ), intent(inout) :: tend
  type(physics_buffer_desc), pointer :: pbuf(:)
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  real(r8), pointer, dimension(:)   :: static_ener_ac_2d  ! Vertically integrated static energy
  real(r8), pointer, dimension(:)   :: water_vap_ac_2d    ! Vertically integrated water vapor
  real(r8), pointer, dimension(:,:) :: tini               !
  real(r8), pointer, dimension(:,:) :: cld                !
  real(r8), pointer, dimension(:,:) :: qini               !
  real(r8), pointer, dimension(:,:) :: cldliqini          !
  real(r8), pointer, dimension(:,:) :: cldiceini          !
  real(r8), pointer, dimension(:,:) :: dtcore             !
  real(r8), pointer, dimension(:,:) :: ast                ! relative humidity cloud fraction 
  type(check_tracers_data):: tracerint      ! tracer mass integrals and cummulative boundary fluxes
  type(physics_ptend) :: ptend              ! indivdual parameterization tendencies
  integer  :: nstep                         ! current timestep number
  real(r8) :: zero(pcols)                   ! array of zeros
  integer  :: lchnk, ncol                   ! chunk identifier
  integer  :: i,k,m                         ! loop iterators
  integer  :: yr, mon, day, tod             ! components of a date
  integer  :: ixcldice, ixcldliq            ! constituent indices for cloud liq and ice
  integer  :: itim_old                      ! pbuf time index
  real(r8) :: obklen(pcols)                 ! Obukhov length for aerosol dry deposition
  real(r8) :: surfric(pcols)                ! surface friction velocity for aerosol dry deposition
  real(r8) :: fh2o(pcols)                   ! h2o flux to balance source from methane chemistry
  real(r8) :: tmp_q     (pcols,pver)        ! tmp variable
  real(r8) :: tmp_cldliq(pcols,pver)        ! tmp variable
  real(r8) :: tmp_cldice(pcols,pver)        ! tmp variable
  real(r8) :: tmp_t     (pcols,pver)        ! tmp variable
  real(r8) :: ftem      (pcols,pver)        ! tmp variable
  logical  :: labort                        ! abort flag
  logical  :: state_debug_checks            ! Debug physics_state
  logical  :: l_tracer_aero
  logical  :: l_vdiff
  logical  :: l_rayleigh
  logical  :: l_gw_drag
  logical  :: l_ac_energy_chk
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  lchnk = state%lchnk
  ncol  = state%ncol
  nstep = get_nstep()
  
  call phys_getopts( state_debug_checks_out = state_debug_checks &
                    ,l_tracer_aero_out      = l_tracer_aero      &
                    ,l_vdiff_out            = l_vdiff            &
                    ,l_rayleigh_out         = l_rayleigh         &
                    ,l_gw_drag_out          = l_gw_drag          &
                    ,l_ac_energy_chk_out    = l_ac_energy_chk    &
                   )

  ! Adjust the surface fluxes to reduce instabilities in near sfc layer
  if (phys_do_flux_avg()) call flux_avg_run(state, cam_in,  pbuf, nstep, ztodt)

  ! Validate the physics state.
  if (state_debug_checks) call physics_state_check(state, name="before tphysac")

  ! Associate pointers with physics buffer fields
  itim_old = pbuf_old_tim_idx()

  call pbuf_get_field(pbuf, pbuf_get_index('DTCORE'), dtcore, start=(/1,1,itim_old/), kount=(/pcols,pver,1/) )

  call pbuf_get_field(pbuf, tini_idx, tini)
  call pbuf_get_field(pbuf, qini_idx, qini)
  call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
  call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

  call pbuf_get_field(pbuf, pbuf_get_index('CLD'),cld,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))
  call pbuf_get_field(pbuf, pbuf_get_index('AST'),ast,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))

  !-----------------------------------------------------------------------------
  ! emissions of aerosols and gas-phase chemistry constituents at surface
  !-----------------------------------------------------------------------------
  if (l_tracer_aero) then
    call chem_emissions( state, cam_in )
  end if ! l_tracer_aero

  !-----------------------------------------------------------------------------
  ! get nstep and zero array for energy checker
  !-----------------------------------------------------------------------------
  zero = 0._r8
  nstep = get_nstep()
  call check_tracers_init(state, tracerint)

  ! Check if LHF exceeds the total moisture content of the lowest layer
  call qneg4('TPHYSAC ', lchnk, ncol, ztodt, &
              state%q(1,pver,1), state%rpdel(1,pver), &
              cam_in%shf, cam_in%lhf, cam_in%cflx )

  !-----------------------------------------------------------------------------
  ! Source/sink terms for advected tracers.
  !-----------------------------------------------------------------------------
  if (l_tracer_aero) then

    call tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      
    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "tracers_timestep_tend", nstep, ztodt, cam_in%cflx)

    call aoa_tracers_timestep_tend(state, ptend, cam_in%cflx, cam_in%landfrac, ztodt)      

    call physics_update(state, ptend, ztodt, tend)
    call check_tracers_chng(state, tracerint, "aoa_tracers_timestep_tend", nstep, ztodt, cam_in%cflx)

    ! Chemistry calculation
    if (chem_is_active()) then
       call chem_timestep_tend(state, ptend, cam_in, cam_out, ztodt, pbuf,  fh2o, fsds)
       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chem", nstep, ztodt, fh2o, zero, zero, zero)
       call check_tracers_chng(state, tracerint, "chem_timestep_tend", nstep, ztodt, cam_in%cflx)
    end if ! chem_is_active

  end if ! l_tracer_aero

  !-----------------------------------------------------------------------------
  ! Vertical diffusion/pbl calculation - (pbl, free atmosphere and molecular)
  !-----------------------------------------------------------------------------

#ifndef MMF_SKIP_FLUX

#ifdef MMF_USE_HB

  call mmf_vertical_diffusion_tend(ztodt, state, &
                                   cam_in%wsx, cam_in%wsy,  &
                                   cam_in%shf, cam_in%cflx, &
                                   surfric, obklen, ptend, ast, &
                                   cam_in%landfrac, sgh30, pbuf )
  call physics_update(state, ptend, ztodt, tend)

#else

  call mmf_surface_ac(state, cam_in, ptend)
  call physics_update(state, ptend, ztodt, tend)

#endif

#endif
  
  obklen(:) = 0.
  surfric(:) = 0.

  !-----------------------------------------------------------------------------
  ! Rayleigh friction calculation
  !-----------------------------------------------------------------------------
  if (l_rayleigh) then

    call rayleigh_friction_tend( ztodt, state, ptend)
    call physics_update(state, ptend, ztodt, tend)

    ! call check_energy_chng(state, tend, "vdiff", nstep, ztodt, cam_in%cflx(:,1), zero, zero, cam_in%shf)
    call check_energy_chng(state, tend, "vdiff", nstep, ztodt, zero, zero, zero, zero)
    
    call check_tracers_chng(state, tracerint, "vdiff", nstep, ztodt, cam_in%cflx)

  end if ! l_rayleigh

  !-----------------------------------------------------------------------------
  !  aerosol dry deposition processes
  !-----------------------------------------------------------------------------
  if (l_tracer_aero) then

    call aero_model_drydep( state, pbuf, obklen, surfric, cam_in, ztodt, cam_out, ptend )
    call physics_update(state, ptend, ztodt, tend)

    ! enforce charge neutrality
    call charge_fix( ncol, state%q(:,:,:) )

  end if ! l_tracer_aero

  !-----------------------------------------------------------------------------
  ! Gravity wave drag
  !-----------------------------------------------------------------------------
  if (l_gw_drag) then

    call gw_tend(state, sgh, pbuf, ztodt, ptend, cam_in)
    call physics_update(state, ptend, ztodt, tend)
    
    ! Check energy integrals
    call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, zero, zero, zero)

  end if ! l_gw_drag

  !-----------------------------------------------------------------------------
  ! Energy budget checks
  !-----------------------------------------------------------------------------
  if (l_ac_energy_chk) then

    call pbuf_set_field(pbuf, teout_idx, state%te_cur, (/1,itim_old/),(/pcols,1/))       

    tmp_t(:ncol,:pver) = state%t(:ncol,:pver)

    ! store dse after tphysac in buffer
    do k = 1,pver
      dtcore(:ncol,k) = state%t(:ncol,k)
    end do

    call set_dry_to_wet(state)    ! Physics had dry, dynamics wants moist

    ! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    tmp_q     (:ncol,:pver) = state%q(:ncol,:pver,1)
    tmp_cldliq(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    tmp_cldice(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
    call physics_dme_adjust(state, tend, qini, ztodt)

  end if ! l_ac_energy_chk

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call diag_phys_tend_writeout (state, pbuf,  tend, ztodt, tmp_q, tmp_cldliq, tmp_cldice, tmp_t, qini, cldliqini, cldiceini)

  call clybry_fam_set( ncol, lchnk, map2chm, state%q, pbuf )

  !-----------------------------------------------------------------------------
  !Integrate column static energy and water vapor
  !-----------------------------------------------------------------------------
  call pbuf_get_field(pbuf, pbuf_get_index('static_ener_ac'), static_ener_ac_2d )
  call pbuf_get_field(pbuf, pbuf_get_index('water_vap_ac'), water_vap_ac_2d )

  static_ener_ac_2d(:) = 0
  water_vap_ac_2d(:)   = 0
  do i = 1,ncol
    do k = 1,pver
      static_ener_ac_2d(i) = static_ener_ac_2d(i) + state%pdel(i,k)*rga*( state%s(i,k) + state%q(i,k,1)*latvap )
      water_vap_ac_2d(i)   = water_vap_ac_2d(i)   + state%pdel(i,k)*rga*  state%q(i,k,1)
    end do
  end do

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
  call check_tracers_fini(tracerint)

end subroutine tphysac

!===================================================================================================
!===================================================================================================

subroutine tphysbc(ztodt, fsns, fsnt, flns, flnt, state, tend, pbuf, fsds, &
                   landm, sgh, sgh30, cam_out, cam_in )
  !----------------------------------------------------------------------------- 
  ! Purpose: Evaluate physics processes BEFORE coupling to sfc components
  !
  ! Processes currently included are: 
  ! dry adjustment, moist convection, stratiform, wet deposition, radiation
  !
  ! Pass surface fields for separate surface flux calculations
  ! Dump appropriate fields to history file.
  ! 
  ! Method: 
  !   Each parameterization should be implemented with this sequence of calls:
  !    1)  Call physics interface
  !    2)  Check energy
  !    3)  Call physics_update
  !   See Interface to Column Physics and Chemistry Packages 
  !     http://www.ccsm.ucar.edu/models/atm-cam/docs/phys-interface/index.html
  !-----------------------------------------------------------------------------
  use physics_buffer,         only: physics_buffer_desc, pbuf_get_field
  use physics_buffer,         only: pbuf_get_index, pbuf_old_tim_idx
  use physics_buffer,         only: dyn_time_lvls
  use shr_kind_mod,           only: r8 => shr_kind_r8
  use physics_types,          only: physics_state, physics_tend, physics_ptend, &
                                    physics_ptend_init, physics_ptend_sum, &
                                    physics_state_check, physics_ptend_scale
  use cam_diagnostics,        only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, &
                                    diag_export, diag_state_b4_phys_write
  use cam_history,            only: outfld, fieldname_len
  use physconst,              only: cpair, latvap, gravit, rga
  use constituents,           only: pcnst, qmin, cnst_get_ind
  use convect_deep,           only: convect_deep_tend_2, deep_scheme_does_scav_trans
  use time_manager,           only: is_first_step, get_nstep
  use check_energy,           only: check_energy_chng, check_energy_fix, check_prect, & 
                                    check_energy_timestep_init, &
                                    check_tracers_data, check_tracers_init, &
                                    check_tracers_chng, check_tracers_fini
  use aero_model,             only: aero_model_wetdep
  use radiation,              only: radiation_tend
  use cloud_diagnostics,      only: cloud_diagnostics_calc
  use perf_mod
  use mo_gas_phase_chemdr,    only: map2chm
  use clybry_fam,             only: clybry_fam_adj
  use clubb_intr,             only: clubb_tend_cam
  use sslt_rebin,             only: sslt_rebin_adv
  use tropopause,             only: tropopause_output
  use output_aerocom_aie,     only: do_aerocom_ind3, cloud_top_aerocom
  use cam_abortutils,         only: endrun
  use phys_control,           only: use_qqflx_fixer, use_mass_borrower
  use crmdims,                only: crm_nz, crm_nx, crm_ny, crm_dx, crm_dy, crm_dt
  use crm_physics,            only: crm_physics_tend, crm_surface_flux_bypass_tend
  use crm_ecpp_output_module, only: crm_ecpp_output_type
  use mmf_surface_mod,        only: mmf_surface_bc

  implicit none
  !-----------------------------------------------------------------------------
  ! Interface Arguments
  !-----------------------------------------------------------------------------
  real(r8), intent(in   ) :: ztodt                         ! 2 delta t (model time increment)
  real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
  real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
  real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
  real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
  real(r8), intent(inout) :: fsds(pcols)                   ! Surface solar down flux
  real(r8), intent(in   ) :: landm(pcols)                  ! land fraction ramp
  real(r8), intent(in   ) :: sgh(pcols)                    ! Std. deviation of orography
  real(r8), intent(in   ) :: sgh30(pcols)                  ! Std. deviation of 30 s orography for tms
  type(physics_state), intent(inout) :: state
  type(physics_tend ), intent(inout) :: tend
  type(physics_buffer_desc), pointer :: pbuf(:)
  type(cam_out_t),     intent(inout) :: cam_out
  type(cam_in_t),      intent(in)    :: cam_in
  !-----------------------------------------------------------------------------
  ! Local variables
  !-----------------------------------------------------------------------------
  type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
  type(physics_ptend)   :: ptend_aero       ! ptend for microp_aero
  type(physics_state)   :: state_alt        ! alt state for CRM input
  integer  :: nstep                         ! current timestep number
  real(r8) :: net_flx(pcols)
  real(r8) :: zdu(pcols,pver)               ! detraining mass flux from deep convection
  real(r8) :: cmfmc(pcols,pverp)            ! Convective mass flux--m sub c
  real(r8) :: cmfcme(pcols,pver)            ! cmf condensation - evaporation
  real(r8) :: cmfmc2(pcols,pverp)           ! Moist convection cloud mass flux
  real(r8) :: dlf(pcols,pver)               ! Detraining cld H20 from shallow + deep convections
  real(r8) :: dlf2(pcols,pver)              ! Detraining cld H20 from shallow convections
  real(r8) :: pflx(pcols,pverp)             ! Conv rain flux thru out btm of lev
  real(r8) :: rtdt                          ! 1./ztodt
  integer  ::lchnk                          ! chunk identifier
  integer  ::ncol                           ! number of atmospheric columns
  integer  ::ierr
  integer  :: i,k,m,ihist                   ! Longitude, level, constituent indices
  integer  :: ixcldice, ixcldliq            ! constituent indices for cloud liquid and ice water.

  ! physics buffer fields to compute tendencies for stratiform package
  integer itim_old, ifld
  real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction
  real(r8), pointer, dimension(:,:) :: cldo       ! old cloud fraction

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

  ! stratiform precipitation variables
  real(r8),pointer :: prec_str(:)    ! sfc flux of precip from stratiform (m/s)
  real(r8),pointer :: snow_str(:)     ! sfc flux of snow from stratiform   (m/s)
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
  real(r8) :: rliq(pcols)                    ! vertical integral of liquid not yet in q(ixcldliq)
  real(r8) :: rliq2(pcols)                   ! vertical integral of liquid from shallow scheme
  real(r8) :: det_s  (pcols)                 ! vertical integral of detrained static energy from ice
  real(r8) :: det_ice(pcols)                 ! vertical integral of detrained ice
  real(r8) :: flx_cnd(pcols)
  real(r8) :: flx_heat(pcols)
  type(check_tracers_data):: tracerint             ! energy integrals and cummulative boundary fluxes
  real(r8) :: zero_tracers(pcols,pcnst)

  logical   :: lq(pcnst)

  character(len=fieldname_len)   :: varname, vsuffix

  !BSINGH - these were moved from zm_conv_intr because they are  used by aero_model_wetdep 
  real(r8):: mu(pcols,pver) 
  real(r8):: eu(pcols,pver)
  real(r8):: du(pcols,pver)
  real(r8):: md(pcols,pver)
  real(r8):: ed(pcols,pver)
  real(r8):: dp(pcols,pver)
  real(r8):: dsubcld(pcols) ! wg layer thickness in mbs (between upper/lower interface).
  integer :: jt(pcols)      ! wg layer thickness in mbs between lcl and maxi.    
  integer :: maxg(pcols)    ! wg top  level index of deep cumulus convection.
  integer :: ideep(pcols)   ! wg gathered values of maxi.
  integer :: lengath        ! w holds position of gathered points vs longitude index

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
  
  character(len=16) :: MMF_microphysics_scheme
  real(r8)          :: crm_run_time              ! length of CRM integration
  real(r8), dimension(pcols,pver) :: mmf_clear_rh ! CRM clear air relative humidity used for aerosol water uptake
  real(r8), dimension(pcols) :: mmf_qchk_prec_dp  ! CRM precipitation diagostic (liq+ice)  used for check_energy_chng
  real(r8), dimension(pcols) :: mmf_qchk_snow_dp  ! CRM precipitation diagostic (ice only) used for check_energy_chng
  real(r8), dimension(pcols) :: mmf_rad_flux      ! CRM radiative flux diagnostic used for check_energy_chng
  type(crm_ecpp_output_type) :: crm_ecpp_output   ! CRM output data for ECPP calculations

  call phys_getopts( MMF_microphysics_scheme_out = MMF_microphysics_scheme)

  call phys_getopts( microp_scheme_out      = microp_scheme, &
                     state_debug_checks_out = state_debug_checks &
                    ,l_bc_energy_fix_out    = l_bc_energy_fix    &
                    ,l_dry_adj_out          = l_dry_adj          &
                    ,l_tracer_aero_out      = l_tracer_aero      &
                    ,l_st_mac_out           = l_st_mac           &
                    ,l_st_mic_out           = l_st_mic           &
                    ,l_rad_out              = l_rad              &
                    )
  
  !-----------------------------------------------------------------------------
  call t_startf('bc_init')

  zero = 0._r8
  zero_tracers(:,:) = 0._r8
  lchnk = state%lchnk
  ncol  = state%ncol
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

  call pbuf_get_field(pbuf, pbuf_get_index('static_ener_ac'), static_ener_ac_2d )
  call pbuf_get_field(pbuf, pbuf_get_index('water_vap_ac'), water_vap_ac_2d )

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

  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------

  if(use_mass_borrower) then 
    ! printout diagnostic information
    call qneg3('TPHYSBCb', lchnk, ncol, pcols, pver, 1, pcnst, qmin, state%q, .False.)
    ! tracer borrower for mass conservation 
     do m = 1, pcnst 
        call massborrow("PHYBC01",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
     end do
  else
    ! original fixer to make sure tracers are all positive
    call qneg3('TPHYSBCb', lchnk, ncol, pcols, pver, 1, pcnst, qmin, state%q, .True. )
  end if 

  ! Validate state coming from the dynamics.
  if (state_debug_checks) call physics_state_check(state, name="before tphysbc (dycore?)")

  call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

  if(use_mass_borrower) then
     ! if use_mass_borrower = True, only printout diagnostic information
     call qneg3('TPHYSBCc',lchnk  ,ncol, pcols, pver, 1, pcnst, qmin, state%q, .False. )
     ! tracer borrower for mass conservation 
     do m = 1, pcnst
        call massborrow("PHYBC02",lchnk,ncol,state%psetcols,m,m,qmin(m),state%q(1,1,m),state%pdel)
     end do
  else
    ! original fixer to make sure tracers are all positive
    call qneg3('TPHYSBCc', lchnk, ncol, pcols, pver, 1, pcnst, qmin, state%q, .True. )
  end if

  ! Validate output of clybry_fam_adj.
  if (state_debug_checks) call physics_state_check(state, name="clybry_fam_adj")

  ! Dump out "before physics" state
  call diag_state_b4_phys_write (state)

  ! compute mass integrals of input tracers state
  call check_tracers_init(state, tracerint)

  call t_stopf('bc_init')

  !-----------------------------------------------------------------------------
  ! Global mean total energy fixer
  !-----------------------------------------------------------------------------
  if (l_bc_energy_fix) then

    call t_startf('energy_fixer')

    tini(:ncol,:pver) = state%t(:ncol,:pver)
    
    call check_energy_fix(state, ptend, nstep, flx_heat)
    call physics_update(state, ptend, ztodt, tend)
    call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
    
    ! Save state for convective tendency calculations.
    call diag_conv_tend_ini(state, pbuf)

    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)

    qini     (:ncol,:pver) = state%q(:ncol,:pver,       1)
    cldliqini(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    cldiceini(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)

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
    
  !-----------------------------------------------------------------------------
  ! Dry adjustment
  !-----------------------------------------------------------------------------
  if (l_dry_adj) then

    call t_startf('dry_adjustment')

    ! Copy state info for input to dadadj
    ! This is a kludge so dadadj doesn't have to be reformulated for DSE

    lq(:) = .FALSE.
    lq(1) = .TRUE.
    call physics_ptend_init(ptend, state%psetcols, 'dadadj', ls=.true., lq=lq)
    ptend%s(:ncol,:pver)   = state%t(:ncol,:pver)
    ptend%q(:ncol,:pver,1) = state%q(:ncol,:pver,1)

    call dadadj (lchnk, ncol, state%pmid, state%pint, state%pdel, ptend%s, ptend%q(1,1,1))

    ptend%s(:ncol,:)   = (ptend%s(:ncol,:)   - state%t(:ncol,:)  )/ztodt * cpair
    ptend%q(:ncol,:,1) = (ptend%q(:ncol,:,1) - state%q(:ncol,:,1))/ztodt

    call physics_update(state, ptend, ztodt, tend)

    call t_stopf('dry_adjustment')

  end if

  !-----------------------------------------------------------------------------
  ! MMF surface flux bypass
  !-----------------------------------------------------------------------------
  ! ! Check if LHF exceeds the total moisture content of the lowest layer
  ! call qneg4('TPHYSBC ', lchnk, ncol, ztodt, &
  !             state%q(1,pver,1), state%rpdel(1,pver), &
  !             cam_in%shf, cam_in%lhf, cam_in%cflx )

#ifndef MMF_SKIP_FLUX

  if (.not. is_first_step()) then
    call mmf_surface_bc(state, cam_in, ptend)
    call physics_update(state, ptend, ztodt, tend)
  end if

#endif

  !-----------------------------------------------------------------------------
  ! Run the CRM 
  !-----------------------------------------------------------------------------
  crm_run_time = ztodt
  call crm_physics_tend(ztodt, state, tend,ptend, pbuf, cam_in, cam_out,    &
                        species_class, crm_ecpp_output, mmf_clear_rh,       &
                        mmf_qchk_prec_dp, mmf_qchk_snow_dp, mmf_rad_flux)

  call physics_update(state, ptend, crm_run_time, tend)

  call check_energy_chng(state, tend, "crm_tend", nstep, crm_run_time,  &
                         zero, mmf_qchk_prec_dp, mmf_qchk_snow_dp, mmf_rad_flux)

  !-----------------------------------------------------------------------------
  ! save old CRM cloud fraction - w/o CRM, this is done in cldwat2m.F90
  !-----------------------------------------------------------------------------
  call pbuf_get_field(pbuf,pbuf_get_index('CLDO'),cldo,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))
  call pbuf_get_field(pbuf,pbuf_get_index('CLD') ,cld ,start=(/1,1,itim_old/),kount=(/pcols,pver,1/))

  cldo(1:ncol,1:pver) = cld(1:ncol,1:pver)

  !-----------------------------------------------------------------------------
  ! Aerosol stuff
  !-----------------------------------------------------------------------------
  if (l_tracer_aero) then
    if ( .not. deep_scheme_does_scav_trans() ) then
      !======================================================================
      ! Aerosol wet chemistry determines scavenging and transformations.
      ! This is followed by convective transport of all trace species except
      ! water vapor and condensate. Scavenging needs to occur prior to
      ! transport in order to determine interstitial fraction.
      !======================================================================
      call t_startf('tphysbc_aerosols')

      ! Aerosol wet removal (including aerosol water uptake)
       call aero_model_wetdep( ztodt, dlf, dlf2, cmfmc2, state,  & ! inputs
              sh_e_ed_ratio, mu, md, du, eu, ed, dp, dsubcld,    &
              jt, maxg, ideep, lengath, species_class,           &
              cam_out, pbuf, ptend,                              & ! outputs
              clear_rh=mmf_clear_rh) ! clear air relative humidity for water uptake
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

  !-----------------------------------------------------------------------------
  ! Moist physical parameteriztions complete: 
  ! send dynamical variables, and derived variables to history file
  !-----------------------------------------------------------------------------
  call t_startf('bc_history_write')
  call diag_phys_writeout(state, cam_out%psl)
  call diag_conv(state, ztodt, pbuf)
  call t_stopf('bc_history_write')

  !-----------------------------------------------------------------------------
  ! Write cloud diagnostics on history file
  !-----------------------------------------------------------------------------
  call t_startf('bc_cld_diag_history_write')
  call cloud_diagnostics_calc(state, pbuf)
  call t_stopf('bc_cld_diag_history_write')

  !-----------------------------------------------------------------------------
  ! Radiation computations
  !-----------------------------------------------------------------------------
  if (l_rad) then
    call t_startf('radiation')

    call radiation_tend(state,ptend, pbuf, cam_out, cam_in, &
                        cam_in%landfrac, landm, cam_in%icefrac, cam_in%snowhland, &
                        fsns, fsnt, flns, flnt, fsds, &
                        net_flx, is_cmip6_volc, ztodt, clear_rh=mmf_clear_rh)

    ! Set net flux used by spectral dycores
    do i=1,ncol
      tend%flx_net(i) = net_flx(i)
    end do
    
    ! don't add radiative tendency to GCM temperature in case of MMF
    ! as it was added above as part of crm tendency.
    ptend%s = 0.

    call physics_update(state, ptend, ztodt, tend)
    
    call check_energy_chng(state, tend, "mmfradheat", nstep, ztodt, zero, zero, zero, zero) 

    call t_stopf('radiation')

  end if ! l_rad

  !-----------------------------------------------------------------------------
  ! Diagnostics
  !-----------------------------------------------------------------------------
  if(do_aerocom_ind3) call cloud_top_aerocom(state, pbuf) 


  ! Diagnose the location of the tropopause
  call tropopause_output(state)

  ! Save atmospheric fields to force surface models
  call cam_export(state,cam_out,pbuf)
  
  ! Write export state to history file
  call diag_export(cam_out)

  call check_tracers_fini(tracerint)
  !-----------------------------------------------------------------------------
  !-----------------------------------------------------------------------------
end subroutine tphysbc

!===================================================================================================
!===================================================================================================

subroutine phys_timestep_init(phys_state, cam_out, pbuf2d)
  !-----------------------------------------------------------------------------
  ! Purpose: The place for parameterizations to call per timestep initializations.
  !          This is used to update time interpolated fields from boundary datasets.
  !-----------------------------------------------------------------------------
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
  use mmf_vertical_diffusion,  only: mmf_vertical_diffusion_ts_init
  use radheat,             only: radheat_timestep_init
  use solar_data,          only: solar_data_advance
  use efield,              only: get_efield
  use perf_mod
  use prescribed_ozone,    only: prescribed_ozone_adv
  use prescribed_ghg,      only: prescribed_ghg_adv
  use prescribed_aero,     only: prescribed_aero_adv
  use aerodep_flx,         only: aerodep_flx_adv
  use aircraft_emit,       only: aircraft_emit_adv
  use prescribed_volcaero, only: prescribed_volcaero_adv

  implicit none
  !-----------------------------------------------------------------------------
  type(physics_state), intent(inout), dimension(begchunk:endchunk) :: phys_state
  type(cam_out_t),     intent(inout), dimension(begchunk:endchunk) :: cam_out

  type(physics_buffer_desc), pointer                 :: pbuf2d(:,:)
  !-----------------------------------------------------------------------------

  ! Chemistry surface values
  call chem_surfvals_set()

  ! Solar irradiance
  call solar_data_advance()

  ! Time interpolate for chemistry.
  call chem_timestep_init(phys_state, pbuf2d)

  ! Prescribed tracers
  call prescribed_ozone_adv(phys_state, pbuf2d)
  call prescribed_ghg_adv(phys_state, pbuf2d)
  call prescribed_aero_adv(phys_state, pbuf2d)
  call aircraft_emit_adv(phys_state, pbuf2d)
  call prescribed_volcaero_adv(phys_state, pbuf2d)

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
  call mmf_vertical_diffusion_ts_init(pbuf2d, phys_state)

  ! Time interpolate for tracers, if appropriate
  call tracers_timestep_init(phys_state)

  ! age of air tracers
  call aoa_tracers_timestep_init(phys_state)

end subroutine phys_timestep_init

!===================================================================================================
!===================================================================================================

subroutine add_fld_default_calls()
  use cam_history,        only: addfld, add_default, fieldname_len

  implicit none

  !Add all existing ptend names for the addfld calls
  character(len=20), parameter :: vlist(27) = (/     'topphysbc           '                       ,&
       'chkenergyfix        ','dadadj              ','zm_convr            ','zm_conv_evap        ',&
       'momtran             ','zm_conv_tend        ','UWSHCU              ','convect_shallow     ',&
       'pcwdetrain_mac      ','macro_park          ','macrop              ','micro_mg            ',&
       'cldwat_mic          ','aero_model_wetdep_ma','convtran2           ','cam_radheat         ',&
       'chemistry           ','vdiff               ','rayleigh_friction   ','aero_model_drydep_ma',&
       'Grav_wave_drag      ','convect_shallow_off ','clubb_ice1          ','clubb_det           ',&
       'clubb_ice4          ','clubb_srf           ' /)

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

!===================================================================================================
!===================================================================================================

end module physpkg
