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
  !-----------------------------------------------------------------------

  use shr_kind_mod,     only: r8 => shr_kind_r8
  use spmd_utils,       only: masterproc
  use physconst,        only: latvap, latice, rh2o
  use physics_types,    only: physics_state, physics_tend, physics_state_set_grid, &
       physics_ptend, physics_tend_init, physics_update,    &
       physics_type_alloc, physics_ptend_dealloc
  use phys_grid,        only: get_ncols_p
  use phys_gmean,       only: gmean_mass
  use ppgrid,           only: begchunk, endchunk, pcols, pver, pverp
  use constituents,     only: pcnst, cnst_name, cnst_get_ind
  use camsrfexch,       only: cam_out_t, cam_in_t

  use cam_control_mod,  only: ideal_phys, adiabatic
  use phys_control,     only: phys_do_flux_avg, waccmx_is
  use scamMod,          only: single_column, scm_crm_mode
  use flux_avg,         only: flux_avg_init
  use infnan,           only: posinf, assignment(=)
#ifdef SPMD
  use mpishorthand
#endif
  use perf_mod
  use cam_logfile,     only: iulog
  use camsrfexch,      only: cam_export

  implicit none
  private

  !  Physics buffer index
  integer ::  teout_idx          = 0  

  integer ::  tini_idx           = 0 
  integer ::  qini_idx           = 0 
  integer ::  cldliqini_idx      = 0 
  integer ::  cldiceini_idx      = 0 

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
    use physics_buffer, only: pbuf_init_time
    use shr_kind_mod,       only: r8 => shr_kind_r8
    use spmd_utils,         only: masterproc
    use constituents,       only: pcnst, cnst_add, cnst_chk_dim, cnst_name

    use cam_control_mod,    only: moist_physics
    use phys_control,       only: phys_do_flux_avg, phys_getopts, waccmx_is
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
    use aerosol_intr,       only: aerosol_register_cnst
    use carma_intr,         only: carma_register
    use cam3_aero_data,     only: cam3_aero_data_on, cam3_aero_data_register
    use cam3_ozone_data,    only: cam3_ozone_data_on, cam3_ozone_data_register
    use ghg_data,           only: ghg_data_register
    use vertical_diffusion, only: vd_register
    use convect_deep,       only: convect_deep_register
    use convect_shallow,    only: convect_shallow_register
    use radiation,          only: radiation_register
    use co2_cycle,          only: co2_register
    use flux_avg,           only: flux_avg_register
    use iondrag,            only: iondrag_register
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
    use physics_buffer,     only: pbuf_add_field, dtype_r8

    implicit none
    !---------------------------Local variables-----------------------------
    !
    integer  :: m        ! loop index
    integer  :: mm       ! constituent index 
    !-----------------------------------------------------------------------

    character(len=16) :: microp_scheme
    logical           :: do_clubb_sgs

    call phys_getopts( microp_scheme_out = microp_scheme )
    call phys_getopts( do_clubb_sgs_out  = do_clubb_sgs )

    ! Initialize pbuf_times
    call pbuf_init_time()

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
       elseif( microp_scheme == 'MG' ) then
          if (.not. do_clubb_sgs) call macrop_driver_register()
          call microp_aero_register()
          call microp_driver_register()
       end if
       
       ! Register CLUBB_SGS here
       if (do_clubb_sgs) call clubb_register_cam()
       

       call pbuf_add_field('PREC_STR',  'physpkg',dtype_r8,(/pcols/),prec_str_idx)
       call pbuf_add_field('SNOW_STR',  'physpkg',dtype_r8,(/pcols/),snow_str_idx)
       call pbuf_add_field('PREC_PCW',  'physpkg',dtype_r8,(/pcols/),prec_pcw_idx)
       call pbuf_add_field('SNOW_PCW',  'physpkg',dtype_r8,(/pcols/),snow_pcw_idx)
       call pbuf_add_field('PREC_SED',  'physpkg',dtype_r8,(/pcols/),prec_sed_idx)
       call pbuf_add_field('SNOW_SED',  'physpkg',dtype_r8,(/pcols/),snow_sed_idx)

       call conv_water_register()

       ! chemical constituents
       call chem_register()

       ! co2 constituents
       call co2_register()

       ! register data model ozone with pbuf
       if (cam3_ozone_data_on) then
          call cam3_ozone_data_register()
       end if
       call prescribed_volcaero_register()
       call prescribed_ozone_register()
       call prescribed_aero_register()
       call prescribed_ghg_register()
       call sslt_rebin_register

       ! CAM3 prescribed aerosols
       if (cam3_aero_data_on) then
          call cam3_aero_data_register()
       end if

       ! register various data model gasses with pbuf
       call ghg_data_register()

       ! carma microphysics
       ! 
       ! NOTE: Needs to come before aerosol_register_cnst, so that the CARMA
       ! flags are defined by then.
       call carma_register()

       ! Register iondrag variables with pbuf
       call iondrag_register()

       ! Register ionosphere variables with pbuf if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call ionos_register()
       endif

       ! aerosols
       call aerosol_register_cnst()

       call aircraft_emit_register()

       ! deep convection
       call convect_deep_register

       !  shallow convection
       call convect_shallow_register

       ! radiation
       call radiation_register
       call cloud_diagnostics_register

       ! vertical diffusion
       if (.not. do_clubb_sgs) call vd_register()
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
    use abortutils, only : endrun

    use physics_buffer, only : pbuf_get_index, pbuf_get_field, physics_buffer_desc, pbuf_set_field, pbuf_times


    use cam_initfiles,       only: initial_file_get_id, topo_file_get_id
    use pio,                 only: file_desc_t
    use ncdio_atm,           only: infld
    use dycore,              only: dycore_is
    use polar_avg,           only: polar_average
    use short_lived_species, only: initialize_short_lived_species
    use comsrf,              only: landm, sgh, sgh30
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
    integer :: tpert_idx, qpert_idx, pblh_idx

    logical :: found=.false., found2=.false.
    integer :: ierr
    character(len=4) :: dim1name
    integer :: ixcldice, ixcldliq
    nullify(tptr,tptr3d,tptr3d_2,cldptr,convptr_3d)

    fh_ini=>initial_file_get_id()

    !   dynamics variables are handled in dyn_init - here we read variables needed for physics 
    !   but not dynamics

    if(dycore_is('UNSTRUCTURED')) then  
       dim1name='ncol'
    else
       dim1name='lon'
    end if
    if(aqua_planet) then
       sgh = 0._r8
       sgh30 = 0._r8
       landm = 0._r8
    else
       fh_topo=>topo_file_get_id()
       call infld('SGH', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            sgh, found, grid_map='PHYS')
       if(.not. found) call endrun('ERROR: SGH not found on topo file')

       call infld('SGH30', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            sgh30, found, grid_map='PHYS')
       if(.not. found) then
          if (masterproc) write(iulog,*) 'Warning: Error reading SGH30 from topo file.'
          if (masterproc) write(iulog,*) 'The field SGH30 will be filled using data from SGH.'
          sgh30 = sgh
       end if

       call infld('LANDM_COSLAT', fh_topo, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            landm, found, grid_map='PHYS')

       if(.not.found) call endrun(' ERROR: LANDM_COSLAT not found on topo dataset.')
    end if

    allocate(tptr(1:pcols,begchunk:endchunk))

    call infld('PBLH', fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, grid_map='PHYS')
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'PBLH initialized to 0.'
    end if
    pblh_idx = pbuf_get_index('pblh')

    call pbuf_set_field(pbuf2d, pblh_idx, tptr)

    call infld('TPERT', fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr(:,:), found, grid_map='PHYS')
    if(.not. found) then
       tptr(:,:) = 0._r8
       if (masterproc) write(iulog,*) 'TPERT initialized to 0.'
    end if
    tpert_idx = pbuf_get_index( 'tpert')
    call pbuf_set_field(pbuf2d, tpert_idx, tptr)

    fieldname='QPERT'  
    qpert_idx = pbuf_get_index( 'qpert',ierr)
    if (qpert_idx > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
            tptr, found, grid_map='PHYS')
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
    call infld(fieldname, fh_ini, dim1name, 'lat', 1, pcols, begchunk, endchunk, &
         tptr, found, grid_map='PHYS')
    if(.not.found) then
       if(masterproc) write(iulog,*) trim(fieldname), ' initialized to 1000.'
       tptr=1000._r8
    end if
    do n=1,pbuf_times
       call pbuf_set_field(pbuf2d, m, tptr, start=(/1,n/), kount=(/pcols,1/))
    end do
    deallocate(tptr)

    do lchnk=begchunk,endchunk
       cam_out(lchnk)%tbot(:) = posinf
    end do

    !
    ! 3-D fields
    !

    allocate(tptr3d(pcols,pver,begchunk:endchunk))

    fieldname='CLOUD'
    m = pbuf_get_index('CLD')
    call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, grid_map='PHYS')
    if(found) then
       do n = 1, pbuf_times
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if

    fieldname='QCWAT'
    m = pbuf_get_index(fieldname,ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='PHYS')
       if(.not. found) then
          call infld('Q',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, grid_map='PHYS')
          if (found) then
             if (masterproc) write(iulog,*) trim(fieldname), ' initialized with Q'
             if(dycore_is('LR')) call polar_average(pver, tptr3d) 	
          else
             call endrun('  '//trim(subname)//' Error:  Q must be on Initial File')
          end if
       end if
       do n = 1, pbuf_times
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    fieldname = 'ICCWAT'
    m = pbuf_get_index(fieldname, ierr)
    if (m > 0) then
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
          tptr3d, found, grid_map='phys')
       if(found) then
          do n = 1, pbuf_times
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          call cnst_get_ind('CLDICE', ixcldice)
          call infld('CLDICE',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
             tptr3d, found, grid_map='PHYS')
          if(found) then
             do n = 1, pbuf_times
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
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='phys')
       if(found) then
          do n = 1, pbuf_times
             call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
          end do
       else
          allocate(tptr3d_2(pcols,pver,begchunk:endchunk))     
          call cnst_get_ind('CLDICE', ixcldice)
          call cnst_get_ind('CLDLIQ', ixcldliq)
          call infld('CLDICE',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, grid_map='PHYS')
          call infld('CLDLIQ',fh_ini,dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d_2, found2, grid_map='PHYS')
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
             do n = 1, pbuf_times
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
       call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
            tptr3d, found, grid_map='phys')
       if(.not.found) then
          call infld('T', fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
               tptr3d, found, grid_map='phys')
          if(dycore_is('LR')) call polar_average(pver, tptr3d) 	
          if (masterproc) write(iulog,*) trim(fieldname), ' initialized with T'
       end if
       do n = 1, pbuf_times
          call pbuf_set_field(pbuf2d, m, tptr3d, (/1,1,n/),(/pcols,pver,1/))
       end do
    end if

    deallocate(tptr3d)
    allocate(tptr3d(pcols,pverp,begchunk:endchunk))

    fieldname = 'TKE'
    m = pbuf_get_index( 'tke')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0.01_r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.01'
    end if


    fieldname = 'KVM'
    m = pbuf_get_index('kvm')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    if (found) then
       call pbuf_set_field(pbuf2d, m, tptr3d)
    else
       call pbuf_set_field(pbuf2d, m, 0._r8)
       if (masterproc) write(iulog,*) trim(fieldname), ' initialized to 0.'
    end if


    fieldname = 'KVH'
    m = pbuf_get_index('kvh')
    call infld(fieldname, fh_ini, dim1name, 'ilev', 'lat', 1, pcols, 1, pverp, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
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
    call infld(fieldname, fh_ini, dim1name, 'lev', 'lat', 1, pcols, 1, pver, begchunk, endchunk, &
         tptr3d, found, grid_map='phys')
    if(found) then
       do n = 1, pbuf_times
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
    use physconst,          only: rair, cpair, gravit, stebol, tmelt, &
                                  latvap, latice, rh2o, rhoh2o, pstd, zvir, &
                                  karman, rhodair, physconst_init 
    use ref_pres,           only: pref_edge, pref_mid

    use aerosol_intr,       only: aerosol_init
    use carma_intr,         only: carma_init
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
    use phys_control,       only: phys_getopts, waccmx_is
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
    use qbo,                only: qbo_init
    use iondrag,            only: iondrag_init
#if ( defined OFFLINE_DYN )
    use metdata,            only: metdata_phys_init
#endif
    use ionosphere,	   only: ionos_init  ! Initialization of ionosphere module (WACCM-X)
    use majorsp_diffusion,  only: mspd_init   ! Initialization of major species diffusion module (WACCM-X)
    use clubb_intr,         only: clubb_ini_cam
    use sslt_rebin,         only: sslt_rebin_init
    use tropopause,         only: tropopause_init
    use solar_data,         only: solar_data_init
    use rad_solar_var,      only: rad_solar_var_init

    ! Input/output arguments
    type(physics_state), pointer       :: phys_state(:)
    type(physics_tend ), pointer       :: phys_tend(:)
    type(physics_buffer_desc), pointer :: pbuf2d(:,:)

    type(cam_out_t),intent(inout)      :: cam_out(begchunk:endchunk)

    ! local variables
    integer :: lchnk

    character(len=16) :: microp_scheme 
    logical           :: do_clubb_sgs

    !-----------------------------------------------------------------------

    ! Get microphysics option
    call phys_getopts(microp_scheme_out = microp_scheme)
    call phys_getopts(do_clubb_sgs_out  = do_clubb_sgs )

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

    ! diag_init makes addfld calls for dynamics fields that are output from
    ! the physics decomposition
    call diag_init()

    call check_energy_init()

    call tracers_init()

    ! age of air tracers
    call aoa_tracers_init()

    teout_idx = pbuf_get_index( 'TEOUT')

    ! For adiabatic or ideal physics don't need to initialize any of the
    ! parameterizations below:
    if (adiabatic .or. ideal_phys) return

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

    ! Initialize some aerosol code
    call aerosol_init(pbuf2d)

    ! initialize carma
    call carma_init()

    ! solar irradiance data modules
    call solar_data_init()


    ! Prognostic chemistry.
    call chem_init(phys_state,pbuf2d)

    ! Prescribed tracers
    call prescribed_ozone_init()
    call prescribed_ghg_init()
    call prescribed_aero_init()
    call aerodep_flx_init()
    call aircraft_emit_init()
    call prescribed_volcaero_init()

    ! co2 cycle            
    if (co2_transport()) then
       call co2_init()
    end if

    ! CAM3 prescribed ozone
    if (cam3_ozone_data_on) call cam3_ozone_data_init(phys_state)

    call gw_init()

    call rayleigh_friction_init()

    call pbl_utils_init(gravit, karman, cpair, rair, zvir)
    if (.not. do_clubb_sgs) call vertical_diffusion_init(pbuf2d)

    if ( waccmx_is('ionosphere') .or. waccmx_is('neutral') ) then
       call mspd_init ()
       ! Initialization of ionosphere module if mode set to ionosphere
       if( waccmx_is('ionosphere') ) then
          call ionos_init()
       endif
    endif

    call tsinti(tmelt, latvap, rair, stebol, latice)

    call radiation_init

    call rad_solar_var_init()

    call cloud_diagnostics_init

    call radheat_init(pref_mid)

    call convect_shallow_init(pref_edge)

    call cldfrc_init

    call convect_deep_init(pref_edge)

    if( microp_scheme == 'RK' ) then
       call stratiform_init()
    elseif( microp_scheme == 'MG' ) then 
       if (.not. do_clubb_sgs) call macrop_driver_init()
       call microp_aero_init()
       call microp_driver_init(pbuf2d)
       call conv_water_init
    end if

    ! initiate CLUBB within CAM
    if (do_clubb_sgs) call clubb_ini_cam(pbuf2d)

    call qbo_init

    call iondrag_init(pref_mid)

#if ( defined OFFLINE_DYN )
    call metdata_phys_init()
#endif
    call sslt_rebin_init()
    call tropopause_init()

    prec_dp_idx  = pbuf_get_index('PREC_DP')
    snow_dp_idx  = pbuf_get_index('SNOW_DP')
    prec_sh_idx  = pbuf_get_index('PREC_SH')
    snow_sh_idx  = pbuf_get_index('SNOW_SH')

end subroutine phys_init

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run1(phys_state, ztodt, phys_tend, pbuf2d,  cam_in, cam_out)
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! First part of atmospheric physics package before updating of surface models
    ! 
    !-----------------------------------------------------------------------
    use time_manager,   only: get_nstep
    use cam_diagnostics,only: diag_allocate, diag_physvar_ic
    use check_energy,   only: check_energy_gmean

    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, pbuf_allocate
#if (defined BFB_CAM_SCAM_IOP )
    use cam_history,    only: outfld
#endif
    use comsrf,         only: fsns, fsnt, flns, sgh30, flnt, landm, fsds
    use abortutils,     only: endrun
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
    !-----------------------------------------------------------------------
    !
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! indices
    integer :: ncol                              ! number of columns
    integer :: nstep                             ! current timestep number
#if (! defined SPMD)
    integer  :: mpicom = 0
#endif
    type(physics_buffer_desc), pointer :: phys_buffer_chunk(:)

    call t_startf ('physpkg_st1')
    nstep = get_nstep()

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

    call t_stopf ('physpkg_st1')

    if ( adiabatic .or. ideal_phys )then
       call t_startf ('bc_physics')
       call phys_run1_adiabatic_or_ideal(ztodt, phys_state, phys_tend,  pbuf2d)
       call t_stopf ('bc_physics')
    else
       call t_startf ('physpkg_st1')

       call pbuf_allocate(pbuf2d, 'physpkg')
       call diag_allocate()

       !-----------------------------------------------------------------------
       ! Advance time information
       !-----------------------------------------------------------------------

       call phys_timestep_init( phys_state, cam_out, pbuf2d)

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
       call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C, phys_buffer_chunk)
       do c=begchunk, endchunk
          !
          ! Output physics terms to IC file
          !
          phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)

          call t_startf ('diag_physvar_ic')
          call diag_physvar_ic ( c,  phys_buffer_chunk, cam_out(c), cam_in(c) )
          call t_stopf ('diag_physvar_ic')

          call tphysbc (ztodt, fsns(1,c), fsnt(1,c), flns(1,c), flnt(1,c), phys_state(c),        &
                       phys_tend(c), phys_buffer_chunk,  fsds(1,c), landm(1,c),          &
                       sgh30(1,c), cam_out(c), cam_in(c) )

       end do

       call t_adj_detailf(-1)
       call t_stopf ('bc_physics')

       ! Don't call the rest in CRM mode
       if(single_column.and.scm_crm_mode) return

#ifdef TRACER_CHECK
       call gmean_mass ('between DRY', phys_state)
#endif
    end if

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
    type(physics_ptend) :: ptend           ! indivdual parameterization tendencies
    real(r8)            :: flx_heat(pcols) ! effective sensible heat flux
    real(r8)            :: zero(pcols)     ! array of zeros

    ! physics buffer field for total energy
    integer itim, ncol
    real(r8), pointer, dimension(:) :: teout
    logical, SAVE :: first_exec_of_phys_run1_adiabatic_or_ideal  = .TRUE.
    !-----------------------------------------------------------------------

    nstep = get_nstep()
    zero  = 0._r8

    ! Associate pointers with physics buffer fields
    itim = pbuf_old_tim_idx()
    if (first_exec_of_phys_run1_adiabatic_or_ideal) then
       first_exec_of_phys_run1_adiabatic_or_ideal  = .FALSE.
    endif

!$OMP PARALLEL DO PRIVATE (C, NCOL, PTEND, FLX_HEAT)
    do c=begchunk, endchunk

       ncol = get_ncols_p(c)

       ! Initialize the physics tendencies to zero.
       call physics_tend_init(phys_tend(c))

       ! Dump dynamics variables to history buffers
       call diag_phys_writeout(phys_state(c))

       if (dycore_is('LR') .or. dycore_is('SE') ) then
          call check_energy_fix(phys_state(c), ptend, nstep, flx_heat)
          call physics_update(phys_state(c), ptend, ztodt, phys_tend(c))
          call check_energy_chng(phys_state(c), phys_tend(c), "chkengyfix", nstep, ztodt, &
               zero, zero, zero, flx_heat)
       end if

       if ( ideal_phys )then
          call t_startf('tphysidl')
          call tphysidl(ztodt, phys_state(c), phys_tend(c))
          call t_stopf('tphysidl')
       end if

       ! Save total enery after physics for energy conservation checks
       call pbuf_set_field(pbuf_get_chunk(pbuf2d, c), teout_idx, phys_state(c)%te_cur)


    end do

end subroutine phys_run1_adiabatic_or_ideal

  !
  !-----------------------------------------------------------------------
  !

subroutine phys_run2(phys_state, ztodt, phys_tend, pbuf2d,  cam_out, &
       cam_in )
    !----------------------------------------------------------------------- 
    ! 
    ! Purpose: 
    ! Second part of atmospheric physics package after updating of surface models
    ! 
    !-----------------------------------------------------------------------
    use physics_buffer,         only: physics_buffer_desc, pbuf_get_chunk, pbuf_deallocate, pbuf_update_tim_idx
    use mo_lightning,   only: lightning_no_prod


    use cam_diagnostics,only: diag_deallocate, diag_surf
    use comsrf,         only: trefmxav, trefmnav, sgh, sgh30, fsds 
    use physconst,      only: stebol, latvap
    use carma_intr,     only: carma_accumulate_stats
#if ( defined OFFLINE_DYN )
    use metdata,        only: get_met_srf2
#endif
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
    !
    !-----------------------------------------------------------------------
    !---------------------------Local workspace-----------------------------
    !
    integer :: c                                 ! chunk index
    integer :: ncol                              ! number of columns
#if (! defined SPMD)
    integer  :: mpicom = 0
#endif
    type(physics_buffer_desc),pointer, dimension(:)     :: phys_buffer_chunk
    !
    ! If exit condition just return
    !

    if(single_column.and.scm_crm_mode) return

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
    call t_adj_detailf(+1)

!$OMP PARALLEL DO PRIVATE (C, NCOL, phys_buffer_chunk)

    do c=begchunk,endchunk
       ncol = get_ncols_p(c)
       phys_buffer_chunk => pbuf_get_chunk(pbuf2d, c)
       !
       ! surface diagnostics for history files
       !
       call t_startf('diag_surf')
       call diag_surf(cam_in(c), cam_out(c), phys_state(c)%ps,trefmxav(1,c), trefmnav(1,c))
       call t_stopf('diag_surf')

       call tphysac(ztodt, cam_in(c),  &
            sgh(1,c), sgh30(1,c), cam_out(c),                              &
            phys_state(c), phys_tend(c), phys_buffer_chunk,&
            fsds(1,c))
    end do                    ! Chunk loop

    call t_adj_detailf(-1)
    call t_stopf('ac_physics')

#ifdef TRACER_CHECK
    call gmean_mass ('after tphysac FV:WET)', phys_state)
#endif

    call t_startf ('carma_accumulate_stats')
    call carma_accumulate_stats()
    call t_stopf ('carma_accumulate_stats')

    call t_startf ('physpkg_st2')
    call pbuf_deallocate(pbuf2d, 'physpkg')

    call pbuf_update_tim_idx()
    call diag_deallocate()
    call t_stopf ('physpkg_st2')

end subroutine phys_run2

  !
  !----------------------------------------------------------------------- 
  !

subroutine phys_final( phys_state, phys_tend, pbuf2d )
    use physics_buffer, only : physics_buffer_desc, pbuf_deallocate
    use chemistry, only : chem_final
    use carma_intr, only : carma_final
    use wv_saturation, only : wv_sat_final
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

    if(associated(pbuf2d)) then
       call pbuf_deallocate(pbuf2d,'global')
       deallocate(pbuf2d)
    end if
    deallocate(phys_state)
    deallocate(phys_tend)
    call chem_final
    call carma_final
    call wv_sat_final

end subroutine phys_final


subroutine tphysac (ztodt,   cam_in,  &
       sgh,     sgh30,                                     &
       cam_out,  state,   tend,    pbuf,            &
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
    use chemistry,          only: chem_is_active, chem_timestep_tend
    use cam_diagnostics,    only: diag_phys_tend_writeout
    use gw_drag,            only: gw_tend
    use vertical_diffusion, only: vertical_diffusion_tend
    use rayleigh_friction,  only: rayleigh_friction_tend
    use constituents,       only: cnst_get_ind
    use physics_types,      only: physics_state, physics_tend, physics_ptend, physics_update,    &
         physics_dme_adjust, set_dry_to_wet, physics_state_check
    use majorsp_diffusion,  only: mspd_intr  ! WACCM-X major diffusion
    use ionosphere,         only: ionos_intr ! WACCM-X ionosphere
    use phys_control,       only: phys_getopts
    use tracers,            only: tracers_timestep_tend
    use aoa_tracers,        only: aoa_tracers_timestep_tend
    use physconst,          only: rhoh2o, latvap,latice
    use aerosol_intr,       only: aerosol_emis_intr, aerosol_drydep_intr
    use carma_intr,         only: carma_emission_tend, carma_timestep_tend
    use carma_flags_mod,    only: carma_do_aerosol, carma_do_emission
    use check_energy,       only: check_energy_chng
    use check_energy,       only: check_tracers_data, check_tracers_init, check_tracers_chng
    use time_manager,       only: get_nstep
    use abortutils,         only: endrun
    use dycore,             only: dycore_is
    use cam_control_mod,    only: aqua_planet 
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_set
    use charge_neutrality,  only: charge_fix
    use qbo,                only: qbo_relax
    use iondrag,            only: iondrag_calc, do_waccm_ions
    use clubb_intr,         only: clubb_surface
    use perf_mod
    use phys_control,       only: phys_do_flux_avg, waccmx_is
    use flux_avg,           only: flux_avg_run

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


    type(check_tracers_data):: tracerint             ! tracer mass integrals and cummulative boundary fluxes

    !
    !---------------------------Local workspace-----------------------------
    !
    type(physics_ptend)     :: ptend               ! indivdual parameterization tendencies

    integer  :: nstep                              ! current timestep number
    real(r8) :: zero(pcols)                        ! array of zeros

    integer :: lchnk                                ! chunk identifier
    integer :: ncol                                 ! number of atmospheric columns
    integer i,k,m                 ! Longitude, level indices
    integer :: yr, mon, day, tod       ! components of a date
    integer :: ixcldice, ixcldliq      ! constituent indices for cloud liquid and ice water.

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

    ! physics buffer fields for total energy and mass adjustment
    integer itim, ifld

    real(r8), pointer, dimension(:,:) :: tini
    real(r8), pointer, dimension(:,:) :: cld
    real(r8), pointer, dimension(:,:) :: qini
    real(r8), pointer, dimension(:,:) :: cldliqini
    real(r8), pointer, dimension(:,:) :: cldiceini
    real(r8), pointer, dimension(:,:) :: dtcore
    real(r8), pointer, dimension(:,:) :: ast     ! relative humidity cloud fraction 

    logical :: do_clubb_sgs 

    ! Debug physics_state.
    logical :: state_debug_checks
    !
    !-----------------------------------------------------------------------
    !
    lchnk = state%lchnk
    ncol  = state%ncol

    nstep = get_nstep()
    
    call phys_getopts( do_clubb_sgs_out       = do_clubb_sgs, &
                       state_debug_checks_out = state_debug_checks)

    ! Adjust the surface fluxes to reduce instabilities in near sfc layer
    if (phys_do_flux_avg()) then 
       call flux_avg_run(state, cam_in,  pbuf, nstep, ztodt)
    endif

    ! Validate the physics state.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysac")

    call t_startf('tphysac_init')
    ! Associate pointers with physics buffer fields
    itim = pbuf_old_tim_idx()


    ifld = pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, start=(/1,1,itim/),kount=(/pcols,pver,1/))

    ifld = pbuf_get_index('AST')
    call pbuf_get_field(pbuf, ifld, ast, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

    !
    ! accumulate fluxes into net flux array for spectral dycores
    ! jrm Include latent heat of fusion for snow
    !
    do i=1,ncol
       tend%flx_net(i) = tend%flx_net(i) + cam_in%shf(i) + (cam_out%precc(i) &
            + cam_out%precl(i))*latvap*rhoh2o &
            + (cam_out%precsc(i) + cam_out%precsl(i))*latice*rhoh2o
    end do

    ! emission of aerosols at surface
    call aerosol_emis_intr (state, cam_in)

    if (carma_do_emission) then
       ! carma emissions
       call carma_emission_tend (state, ptend, cam_in, ztodt)
       call physics_update(state, ptend, ztodt, tend)
    end if

    ! get nstep and zero array for energy checker
    zero = 0._r8
    nstep = get_nstep()
    call check_tracers_init(state, tracerint)

    ! Check if latent heat flux exceeds the total moisture content of the
    ! lowest model layer, thereby creating negative moisture.

    call qneg4('TPHYSAC '       ,lchnk               ,ncol  ,ztodt ,               &
         state%q(1,pver,1),state%rpdel(1,pver) ,cam_in%shf ,         &
         cam_in%lhf , cam_in%cflx )

    call t_stopf('tphysac_init')
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

    !===================================================
    ! Vertical diffusion/pbl calculation
    ! Call vertical diffusion code (pbl, free atmosphere and molecular)
    !===================================================

    ! If CLUBB is called, do not call vertical diffusion, but obukov length and
    !   surface friction velocity still need to be computed.  In addition, 
    !   surface fluxes need to be updated here for constituents 
    if (do_clubb_sgs) then

       call clubb_surface ( state, ptend, ztodt, cam_in, surfric, obklen)
       
       ! Update surface flux constituents 
       call physics_update(state, ptend, ztodt, tend)

    else

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
    
    endif


    !===================================================
    ! Rayleigh friction calculation
    !===================================================
    call t_startf('rayleigh_friction')
    call rayleigh_friction_tend( ztodt, state, ptend)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('rayleigh_friction')

    if (do_clubb_sgs) then
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, zero, zero, zero, zero)
    else
      call check_energy_chng(state, tend, "vdiff", nstep, ztodt, cam_in%cflx(:,1), zero, &
           zero, cam_in%shf)
    endif
    
    call check_tracers_chng(state, tracerint, "vdiff", nstep, ztodt, cam_in%cflx)

    !  aerosol dry deposition processes
    call t_startf('aero_drydep')
    call aerosol_drydep_intr (state, ptend, cam_in, cam_out, ztodt,  &
         fsds, obklen, surfric, prect, pbuf)
    call physics_update(state, ptend, ztodt, tend)
    call t_stopf('aero_drydep')

   ! CARMA microphysics
   !
   ! NOTE: This does both the timestep_tend for CARMA aerosols as well as doing the dry
   ! deposition for CARMA aerosols. It needs to follow vertical_diffusion_tend, so that
   ! obklen and surfric have been calculated. It needs to follow aerosol_drydep_intr, so
   ! that cam_out%xxxdryxxx fields have already been set for CAM aerosols and cam_out
   ! can be added to for CARMA aerosols.
   if (carma_do_aerosol) then
     call t_startf('carma_timestep_tend')
     call carma_timestep_tend(state, cam_in, cam_out, ptend, ztodt, pbuf, obklen=obklen, ustar=surfric)
     call physics_update(state, ptend, ztodt, tend)
   
     call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, zero, zero, zero)
     call t_stopf('carma_timestep_tend')
   end if


    !---------------------------------------------------------------------------------
    !	... enforce charge neutrality
    !---------------------------------------------------------------------------------
    call charge_fix( ncol, state%q(:,:,:) )

    !===================================================
    ! Gravity wave drag
    !===================================================
    call t_startf('gw_tend')

    call gw_tend(state, sgh, pbuf, ztodt, ptend, cam_in)

    call physics_update(state, ptend, ztodt, tend)
    ! Check energy integrals
    call check_energy_chng(state, tend, "gwdrag", nstep, ztodt, zero, zero, zero, zero)
    call t_stopf('gw_tend')

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


    !-------------- Energy budget checks vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv

    call pbuf_set_field(pbuf, teout_idx, state%te_cur, (/1,itim/),(/pcols,1/))       


    ! store dse after tphysac in buffer
    do k = 1,pver
       dtcore(:ncol,k) = state%s(:ncol,k)
    end do

    !*** BAB's FV heating kludge *** apply the heating as temperature tendency.
    !*** BAB's FV heating kludge *** modify the temperature in the state structure
    tmp_t(:ncol,:pver) = state%t(:ncol,:pver)
    state%t(:ncol,:pver) = tini(:ncol,:pver) + ztodt*tend%dtdt(:ncol,:pver)

    !
    ! FV: convert dry-type mixing ratios to moist here because physics_dme_adjust
    !     assumes moist. This is done in p_d_coupling for other dynamics. Bundy, Feb 2004.


    if ( dycore_is('LR') .or. dycore_is('SE')) call set_dry_to_wet(state)    ! Physics had dry, dynamics wants moist


    ! Scale dry mass and energy (does nothing if dycore is EUL or SLD)
    call cnst_get_ind('CLDLIQ', ixcldliq)
    call cnst_get_ind('CLDICE', ixcldice)
    tmp_q     (:ncol,:pver) = state%q(:ncol,:pver,1)
    tmp_cldliq(:ncol,:pver) = state%q(:ncol,:pver,ixcldliq)
    tmp_cldice(:ncol,:pver) = state%q(:ncol,:pver,ixcldice)
    call physics_dme_adjust(state, tend, qini, ztodt)
!!!   REMOVE THIS CALL, SINCE ONLY Q IS BEING ADJUSTED. WON'T BALANCE ENERGY. TE IS SAVED BEFORE THIS
!!!   call check_energy_chng(state, tend, "drymass", nstep, ztodt, zero, zero, zero, zero)

    !-------------- Energy budget checks ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
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

    call clybry_fam_set( ncol, lchnk, map2chm, state%q, pbuf )

end subroutine tphysac

subroutine tphysbc (ztodt,               &
       fsns,    fsnt,    flns,    flnt,    state,   &
       tend,    pbuf,     fsds,    landm,            &
       sgh30, cam_out, cam_in )
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
    use physics_buffer,          only : pbuf_get_index, pbuf_old_tim_idx, pbuf_times
    use shr_kind_mod,    only: r8 => shr_kind_r8

    use stratiform,      only: stratiform_tend
    use phys_control,    only: phys_getopts
    use microp_driver,   only: microp_driver_tend
    use microp_aero,     only: microp_aero_run
    use macrop_driver,   only: macrop_driver_tend
    use physics_types,   only: physics_state, physics_tend, physics_ptend, physics_update, &
         physics_ptend_init, physics_ptend_sum, physics_state_check
    use cam_diagnostics, only: diag_conv_tend_ini, diag_phys_writeout, diag_conv, diag_export, diag_state_b4_phys_write
    use cam_history,     only: outfld
    use physconst,       only: cpair, latvap
    use constituents,    only: pcnst, qmin, cnst_get_ind
    use convect_deep,    only: convect_deep_tend, convect_deep_tend_2, deep_scheme_does_scav_trans
    use time_manager,    only: is_first_step, get_nstep
    use convect_shallow, only: convect_shallow_tend
    use check_energy,    only: check_energy_chng, check_energy_fix
    use check_energy,    only: check_tracers_data, check_tracers_init, check_tracers_chng
    use dycore,          only: dycore_is
    use aerosol_intr,    only: aerosol_wet_intr
    use carma_intr,      only: carma_wetdep_tend, carma_timestep_tend
    use carma_flags_mod, only: carma_do_detrain, carma_do_cldice, carma_do_cldliq,  carma_do_wetdep
    use radiation,       only: radiation_tend
    use cloud_diagnostics, only: cloud_diagnostics_calc
    use perf_mod
#ifdef MODAL_AERO
    use modal_aero_data, only: qneg3_worst_thresh_amode
#endif
    use mo_gas_phase_chemdr,only: map2chm
    use clybry_fam,         only: clybry_fam_adj
    use clubb_intr,      only: clubb_tend_cam
    use sslt_rebin,      only: sslt_rebin_adv
    use tropopause,      only: tropopause_output
    use abortutils,      only: endrun

    implicit none

    !
    ! Arguments
    !
    real(r8), intent(in) :: ztodt                          ! 2 delta t (model time increment)
    real(r8), intent(inout) :: fsns(pcols)                   ! Surface solar absorbed flux
    real(r8), intent(inout) :: fsnt(pcols)                   ! Net column abs solar flux at model top
    real(r8), intent(inout) :: flns(pcols)                   ! Srf longwave cooling (up-down) flux
    real(r8), intent(inout) :: flnt(pcols)                   ! Net outgoing lw flux at model top
    real(r8), intent(inout) :: fsds(pcols)                   ! Surface solar down flux
    real(r8), intent(in) :: landm(pcols)                   ! land fraction ramp
    real(r8), intent(in) :: sgh30(pcols)                   ! Std. deviation of 30 s orography for tms

    type(physics_state), intent(inout) :: state
    type(physics_tend ), intent(inout) :: tend
    type(physics_buffer_desc), pointer :: pbuf(:)

    type(cam_out_t),     intent(inout) :: cam_out
    type(cam_in_t),      intent(in)    :: cam_in


    !
    !---------------------------Local workspace-----------------------------
    !

    type(physics_ptend)   :: ptend            ! indivdual parameterization tendencies
    type(physics_state)   :: state_sc         ! state for sub-columns
    type(physics_ptend)   :: ptend_sc         ! ptend for sub-columns
    type(physics_ptend)   :: ptend_aero       ! ptend for microp_aero
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

    integer  i,k,m                             ! Longitude, level, constituent indices
    integer :: ixcldice, ixcldliq              ! constituent indices for cloud liquid and ice water.

    ! physics buffer fields to compute tendencies for stratiform package
    integer itim, ifld
    real(r8), pointer, dimension(:,:) :: cld        ! cloud fraction


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

    ! carma precipitation variables
    real(r8) :: prec_sed_carma(pcols)          ! total precip from cloud sedimentation (CARMA)
    real(r8) :: snow_sed_carma(pcols)          ! snow from cloud ice sedimentation (CARMA)

    ! stratiform precipitation variables
    real(r8),pointer :: prec_str(:)    ! sfc flux of precip from stratiform (m/s)
    real(r8),pointer :: snow_str(:)     ! sfc flux of snow from stratiform   (m/s)
    real(r8),pointer :: prec_pcw(:)     ! total precip from prognostic cloud scheme
    real(r8),pointer :: snow_pcw(:)     ! snow from prognostic cloud scheme
    real(r8),pointer :: prec_sed(:)     ! total precip from cloud sedimentation
    real(r8),pointer :: snow_sed(:)     ! snow from cloud ice sedimentation

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

    real(r8)  :: cmeliq(pcols,pver)                      ! Rate of cond-evap of liq within the cloud

    logical   :: lq(pcnst)

    !  pass macro to micro
    character(len=16) :: microp_scheme 
    character(len=16) :: macrop_scheme

    ! Debug physics_state.
    logical :: state_debug_checks

    call phys_getopts( microp_scheme_out      = microp_scheme, &
                       macrop_scheme_out      = macrop_scheme, &
                       state_debug_checks_out = state_debug_checks)
    
    !-----------------------------------------------------------------------
    call t_startf('bc_init')

    zero = 0._r8
    zero_tracers(:,:) = 0._r8

    lchnk = state%lchnk
    ncol  = state%ncol

    rtdt = 1._r8/ztodt

    nstep = get_nstep()


    ! Associate pointers with physics buffer fields
    itim = pbuf_old_tim_idx()
    ifld = pbuf_get_index('CLD')
    call pbuf_get_field(pbuf, ifld, cld, (/1,1,itim/),(/pcols,pver,1/))

    call pbuf_get_field(pbuf, teout_idx, teout, (/1,itim/), (/pcols,1/))

    call pbuf_get_field(pbuf, tini_idx, tini)
    call pbuf_get_field(pbuf, qini_idx, qini)
    call pbuf_get_field(pbuf, cldliqini_idx, cldliqini)
    call pbuf_get_field(pbuf, cldiceini_idx, cldiceini)

    ifld   =  pbuf_get_index('DTCORE')
    call pbuf_get_field(pbuf, ifld, dtcore, start=(/1,1,itim/), kount=(/pcols,pver,1/) )

    ifld    = pbuf_get_index('FRACIS')
    call pbuf_get_field(pbuf, ifld, fracis, start=(/1,1,1/), kount=(/pcols, pver, pcnst/)  )

    ! Set physics tendencies to 0
    tend %dTdt(:ncol,:pver)  = 0._r8
    tend %dudt(:ncol,:pver)  = 0._r8
    tend %dvdt(:ncol,:pver)  = 0._r8

    !
    ! Make sure that input tracers are all positive (otherwise,
    ! clybry_fam_adj will crash every few years).
    !

#ifdef MODAL_AERO
    call qneg3_modalx1( &
         'TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q, qneg3_worst_thresh_amode )
#else
    call qneg3('TPHYSBCb',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q )
#endif

    ! Validate state coming from the dynamics.
    if (state_debug_checks) &
         call physics_state_check(state, name="before tphysbc (dycore?)")

    call clybry_fam_adj( ncol, lchnk, map2chm, state%q, pbuf )

    ! Since clybry_fam_adj operates directly on the tracers, and has no
    ! physics_update call, re-run qneg3.

#ifdef MODAL_AERO
    call qneg3_modalx1( &
         'TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q, qneg3_worst_thresh_amode )
#else
    call qneg3('TPHYSBCc',lchnk  ,ncol    ,pcols   ,pver    , &
         1, pcnst, qmin  ,state%q )
#endif

    ! Validate output of clybry_fam_adj.
    if (state_debug_checks) &
         call physics_state_check(state, name="clybry_fam_adj")

    fracis (:ncol,:,1:pcnst) = 1._r8
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
    call t_startf('energy_fixer')

    !*** BAB's FV heating kludge *** save the initial temperature
    tini(:ncol,:pver) = state%t(:ncol,:pver)
    if (dycore_is('LR') .or. dycore_is('SE') ) then
       call check_energy_fix(state, ptend, nstep, flx_heat)
       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "chkengyfix", nstep, ztodt, zero, zero, zero, flx_heat)
    end if
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
    if( nstep > pbuf_times-1 ) then
       do k = 1,pver
          dtcore(:ncol,k) = (state%s(:ncol,k) - dtcore(:ncol,k))/(cpair*ztodt)
       end do
       call outfld( 'DTCORE', dtcore, pcols, lchnk )
    end if

    call t_stopf('energy_fixer')
    !
    !===================================================
    ! Dry adjustment
    ! This code block is not a good example of interfacing a parameterization
    !===================================================
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
         rliq,    &
         ztodt,   &
         state,   ptend, cam_in%landfrac, pbuf) 
    call t_stopf('convect_deep_tend')

    call physics_update(state, ptend, ztodt, tend)

    call pbuf_get_field(pbuf, prec_dp_idx, prec_dp )
    call pbuf_get_field(pbuf, snow_dp_idx, snow_dp )
    call pbuf_get_field(pbuf, prec_sh_idx, prec_sh )
    call pbuf_get_field(pbuf, snow_sh_idx, snow_sh )
    call pbuf_get_field(pbuf, prec_str_idx, prec_str )
    call pbuf_get_field(pbuf, snow_str_idx, snow_str )
    call pbuf_get_field(pbuf, prec_sed_idx, prec_sed )
    call pbuf_get_field(pbuf, snow_sed_idx, snow_sed )

    ! Check energy integrals, including "reserved liquid"
    flx_cnd(:ncol) = prec_dp(:ncol) + rliq(:ncol)
    call check_energy_chng(state, tend, "convect_deep", nstep, ztodt, zero, flx_cnd, snow_dp, zero)

    !
    ! Call Hack (1994) convection scheme to deal with shallow/mid-level convection
    !
    call t_startf ('convect_shallow_tend')

    call convect_shallow_tend (ztodt   , cmfmc,  cmfmc2  ,&
         dlf        , dlf2   ,  rliq   , rliq2, & 
         state      , ptend  ,  pbuf)
    call t_stopf ('convect_shallow_tend')

    call physics_update(state, ptend, ztodt, tend)

    flx_cnd(:ncol) = prec_sh(:ncol) + rliq2(:ncol)
    call check_energy_chng(state, tend, "convect_shallow", nstep, ztodt, zero, flx_cnd, snow_sh, zero)

    call check_tracers_chng(state, tracerint, "convect_shallow", nstep, ztodt, zero_tracers)

    call t_stopf('moist_convection')

    ! Rebin the 4-bin version of sea salt into bins for coarse and accumulation
    ! modes that correspond to the available optics data.  This is only necessary
    ! for CAM-RT.  But it's done here so that the microphysics code which is called
    ! from the stratiform interface has access to the same aerosols as the radiation
    ! code.
    call sslt_rebin_adv(pbuf,  state)
    
    !===================================================
    ! Calculate tendencies from CARMA bin microphysics.
    !===================================================
    !
    ! If CARMA is doing detrainment, then on output, rliq no longer represents water reserved
    ! for detrainment, but instead represents potential snow fall. The mass and number of the
    ! snow are stored in the physics buffer and will be incorporated by the MG microphysics.
    !
    ! Currently CARMA cloud microphysics is only supported with the MG microphysics.
    call t_startf('carma_timestep_tend')

    if (carma_do_cldice .or. carma_do_cldliq) then
       call carma_timestep_tend(state, cam_in, cam_out, ptend, ztodt, pbuf, dlf=dlf, rliq=rliq, &
            prec_str=prec_str, snow_str=snow_str, prec_sed=prec_sed_carma, snow_sed=snow_sed_carma)
       call physics_update(state, ptend, ztodt, tend)

       ! Before the detrainment, the reserved condensate is all liquid, but if CARMA is doing
       ! detrainment, then the reserved condensate is snow.
       if (carma_do_detrain) then
          call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, prec_str+rliq, snow_str+rliq, zero)
       else
          call check_energy_chng(state, tend, "carma_tend", nstep, ztodt, zero, prec_str, snow_str, zero)
       end if
    end if

    call t_stopf('carma_timestep_tend')

    if( microp_scheme == 'RK' ) then

       !===================================================
       ! Calculate stratiform tendency (sedimentation, detrain, cloud fraction and microphysics )
       !===================================================
       call t_startf('stratiform_tend')

       call stratiform_tend(state, ptend, pbuf, ztodt, &
            cam_in%icefrac, cam_in%landfrac, cam_in%ocnfrac, &
            landm, cam_in%snowhland, & ! sediment
            dlf, dlf2, & ! detrain
            rliq  , & ! check energy after detrain
            cmfmc,   cmfmc2, &
            cam_in%ts,      cam_in%sst,        zdu)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "cldwat_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

       call t_stopf('stratiform_tend')

    elseif( microp_scheme == 'MG' ) then

       !===================================================
       ! Calculate macrophysical tendency (sedimentation, detrain, cloud fraction)
       !===================================================

       call t_startf('macrop_tend')

       ! don't call Park macrophysics if CLUBB is called
       if (macrop_scheme .ne. 'CLUBB_SGS') then

          call macrop_driver_tend(state, ptend, ztodt, &
               cam_in%landfrac, cam_in%ocnfrac, &
               cam_in%snowhland, & ! sediment
               dlf, dlf2, & ! detrain
               cmfmc,   cmfmc2, &
               cam_in%ts,      cam_in%sst, zdu,  pbuf, &
               cmeliq, det_s, det_ice)

          !  Since we "added" the reserved liquid back in this routine, we need 
	  !    to account for it in the energy checker
          flx_cnd(:ncol) = -1._r8*rliq(:ncol) 
	  flx_heat(:ncol) = det_s(:ncol)
          
          call physics_update(state, ptend, ztodt, tend)
          call check_energy_chng(state, tend, "macrop_tend", nstep, ztodt, zero, flx_cnd, det_ice, flx_heat)
       
       else ! Calculate CLUBB macrophysics

          ! =====================================================
          !    CLUBB call (PBL, shallow convection, macrophysics)
          ! =====================================================  
   
          call clubb_tend_cam(state,ptend,pbuf,1.0_r8*ztodt,&
             cmfmc, cmfmc2, cam_in, sgh30, dlf, det_s, det_ice, cmeliq)

          !  Since we "added" the reserved liquid back in this routine, we need 
	  !    to account for it in the energy checker
          flx_cnd(:ncol) = -1._r8*rliq(:ncol) 
	  flx_heat(:ncol) = cam_in%shf(:ncol) + det_s(:ncol)

          !    Update physics tendencies and copy state to state_eq, because that is 
          !      input for microphysics              
          call physics_update(state, ptend, ztodt, tend)
          call check_energy_chng(state, tend, "clubb_tend", nstep, ztodt, cam_in%lhf/latvap, flx_cnd, det_ice, flx_heat)
 
       endif 

       call t_stopf('macrop_tend') 

       !===================================================
       ! Calculate cloud microphysics 
       !===================================================

       call t_startf('microp_aero_run')
       call microp_aero_run(state, ptend_aero, ztodt, pbuf)
       call t_stopf('microp_aero_run')

       call t_startf('microp_tend')

       call microp_driver_tend( &
            state, ptend, ztodt, pbuf, cmeliq)

       ! combine aero and micro tendencies
       call physics_ptend_sum(ptend_aero, ptend, ncol)

       call physics_update(state, ptend, ztodt, tend)
       call check_energy_chng(state, tend, "microp_tend", nstep, ztodt, zero, prec_str, snow_str, zero)

       call physics_ptend_dealloc(ptend_aero)
       call t_stopf('microp_tend')

    endif

    ! Add the precipitation from CARMA to the precipitation from stratiform.
    if (carma_do_cldice .or. carma_do_cldliq) then
       prec_sed(:ncol) = prec_sed(:ncol) + prec_sed_carma(:ncol)
       snow_sed(:ncol) = snow_sed(:ncol) + snow_sed_carma(:ncol)
    end if

    if ( .not. deep_scheme_does_scav_trans() ) then

       !===================================================
       !  Aerosol wet chemistry determines scavenging fractions, and transformations
       !
       !
       !  Then do convective transport of all trace species except water vapor and
       !     cloud liquid and ice (we needed to do the scavenging first
       !     to determine the interstitial fraction) 
       !===================================================

       call t_startf('bc_aerosols')
       call aerosol_wet_intr (state, ptend, ztodt, pbuf,  cam_out, dlf)
       call physics_update(state, ptend, ztodt, tend)

       if (carma_do_wetdep) then
          ! CARMA wet deposition
          !
          ! NOTE: It needs to follow aerosol_drydep_intr, so that cam_out%xxxwetxxx
          ! fields have already been set for CAM aerosols and cam_out can be added
          ! to for CARMA aerosols.
          call t_startf ('carma_wetdep_tend')
          call carma_wetdep_tend(state, ptend, ztodt, pbuf, dlf, cam_out)
          call physics_update(state, ptend, ztodt, tend)
          call t_stopf ('carma_wetdep_tend')
       end if

       call t_startf ('convect_deep_tend2')
       call convect_deep_tend_2( state,   ptend,  ztodt,  pbuf ) 
       call t_stopf ('convect_deep_tend2')

       call physics_update(state, ptend, ztodt, tend)

       ! check tracer integrals
       call check_tracers_chng(state, tracerint, "cmfmca", nstep, ztodt,  zero_tracers)

       call t_stopf('bc_aerosols')

   endif

    !===================================================
    ! Moist physical parameteriztions complete: 
    ! send dynamical variables, and derived variables to history file
    !===================================================

    call t_startf('bc_history_write')
    call diag_phys_writeout(state, cam_out%psl)
    call diag_conv(state, ztodt, pbuf)

    call t_stopf('bc_history_write')

    !===================================================
    ! Write cloud diagnostics on history file
    !===================================================

    call t_startf('bc_cld_diag_history_write')

    call cloud_diagnostics_calc(state, pbuf)

    call t_stopf('bc_cld_diag_history_write')

    !===================================================
    ! Radiation computations
    !===================================================
    call t_startf('radiation')

    call radiation_tend(state,ptend, pbuf, &
         cam_out, cam_in, &
         cam_in%landfrac,landm,cam_in%icefrac, cam_in%snowhland, &
         fsns,    fsnt, flns,    flnt,  &
         fsds, net_flx)

    ! Set net flux used by spectral dycores
    do i=1,ncol
       tend%flx_net(i) = net_flx(i)
    end do
    call physics_update(state, ptend, ztodt, tend)
    call check_energy_chng(state, tend, "radheat", nstep, ztodt, zero, zero, zero, net_flx)

    call t_stopf('radiation')

    ! Diagnose the location of the tropopause and its location to the history file(s).
    call t_startf('tropopause')
    call tropopause_output(state)
    call t_stopf('tropopause')

    ! Save atmospheric fields to force surface models
    call t_startf('cam_export')
    call cam_export (state,cam_out,pbuf)
    call t_stopf('cam_export')

    ! Write export state to history file
    call t_startf('diag_export')
    call diag_export(cam_out)
    call t_stopf('diag_export')

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
  use carma_intr,          only: carma_timestep_init
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
  use aerodep_flx,         only: aerodep_flx_adv
  use aircraft_emit,       only: aircraft_emit_adv
  use prescribed_volcaero, only: prescribed_volcaero_adv


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

  call carma_timestep_init()

  ! Time interpolate for tracers, if appropriate
  call tracers_timestep_init(phys_state)

  ! age of air tracers
  call aoa_tracers_timestep_init(phys_state)

end subroutine phys_timestep_init

end module physpkg
