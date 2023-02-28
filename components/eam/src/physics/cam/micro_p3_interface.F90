module micro_p3_interface

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!
  !! Interface between E3SM and P3 microphysics
  !!
  !! Author: Peter Caldwell
  !!
  !! Last updated: 2018-09-12
  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  use shr_kind_mod,   only: rtype=>shr_kind_r8
  use ppgrid,         only: pcols,pver,pverp

!comment: I think Kai added handle_errmsg. It would be better to
!use standard E3SM libraries if possible.
  use error_messages, only: handle_errmsg

  use physics_types,  only: physics_state, &
                            physics_ptend, &
                            physics_ptend_init
  use physconst,      only: mwdry, cpair, mwh2o, gravit, rair, cpliq, pi, &
                            rh2o, latvap, latice, tmelt, rhoh2o, rairv
  use constituents,   only: cnst_add, pcnst, sflxnam, apcnst, bpcnst, pcnst,&
                            cnst_name, cnst_get_ind,cnst_longname
  use physics_buffer, only: physics_buffer_desc, dtype_r8, &
                            pbuf_get_field, pbuf_add_field,dyn_time_lvls,dtype_i4, &
                            pbuf_set_field, pbuf_get_index, &
                            pbuf_old_tim_idx
  use ref_pres,       only: top_lev=>trop_cloud_top_lev
  use phys_control,   only: phys_getopts
  use cam_abortutils, only: endrun
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use time_manager,   only: is_first_step, get_curr_date
  use perf_mod,       only: t_startf, t_stopf
  use micro_p3_utils, only: p3_qc_autocon_expon, p3_qc_accret_expon
  use pio,            only: file_desc_t, pio_nowrite
  use cam_pio_utils,    only: cam_pio_openfile,cam_pio_closefile
  use cam_grid_support, only: cam_grid_check, cam_grid_id, cam_grid_get_dim_names
  use ncdio_atm,       only: infld
  use ppgrid,         only: begchunk, endchunk, pcols, pver, pverp,psubcols

  implicit none
  save

  public :: micro_p3_init, micro_p3_register, micro_p3_tend, &
            micro_p3_init_cnst, micro_p3_implements_cnst &
            ,micro_p3_readnl

  character(len=16), parameter :: unset_str = 'UNSET'

  private

  !Define indices for state%q constituents at module level so
  !defining them in micro_p3_register makes them permanently
  !available.
  CHARACTER(len=16) :: precip_frac_method = 'max_overlap'  ! AaronDonahue, Hard-coded for now, should be fixed in the future

  integer, public ::    &
       ixcldliq = -1,   & ! cloud liquid amount index
       ixcldice = -1,      & ! ice index
       ixnumliq = -1,   & ! cloud liquid number index
       ixnumice = -1,   & ! cloud ice number index
       ixrain   = -1,   & ! rain index
       ixnumrain= -1,   & ! rain number index
       ixcldrim = -1,      & ! rime index ??
       ixrimvol  = -1,  & ! rime volume index ??
       ixqm  = -1      ! ?? index ??

!! pbuf
   integer :: &
      cldo_idx,           &
      qme_idx,            &
      precip_total_tend_idx,          &
      nevapr_idx,         &
      dei_idx,            &
      rate1_cw2pr_st_idx, &
      mu_idx,             &
      lambdac_idx,        &
      rei_idx,            &
      rel_idx,            &
      ls_flxprc_idx,      &
      ls_flxsnw_idx,      &
      ls_reffrain_idx,    &
      ls_reffsnow_idx,    &
      cv_reffliq_idx,     &
      cv_reffice_idx,     &
      qr_evap_tend_idx,      &
      cmeliq_idx,         &
      relvar_idx,         &
      qv_prev_idx,        &
      t_prev_idx,         &
      accre_enhan_idx,    &
      ccn3_idx


! Physics buffer indices for fields registered by other modules
   integer :: &
      ast_idx = -1

   integer :: &
      ni_activated_idx = -1,           &
      npccn_idx = -1,          &
      prec_str_idx = -1,       &
      prec_pcw_idx = -1,       &
      prec_sed_idx = -1,       &
      snow_str_idx = -1,       &
      snow_pcw_idx = -1,       &
      snow_sed_idx = -1

   real(rtype) :: &
      micro_mg_accre_enhan_fac = huge(1.0_rtype), & !Accretion enhancement factor from namelist
      prc_coef1_in             = huge(1.0_rtype), &
      prc_exp_in               = huge(1.0_rtype), &
      prc_exp1_in              = huge(1.0_rtype)

   integer :: ncnst

   character(len=8), parameter :: &      ! Constituent names
      cnst_names(8) = (/'CLDLIQ', 'CLDICE','NUMLIQ','NUMICE', &
                      'RAINQM', 'CLDRIM','NUMRAI','BVRIM '/)

   character(len=128) :: micro_p3_lookup_dir     = unset_str ! location of p3 input files
   character(len=16)  :: micro_p3_tableversion   = unset_str ! P3 table version
   logical            :: micro_aerosolactivation = .false.   ! Use aerosol activation
   logical            :: micro_subgrid_cloud     = .false.   ! Use subgrid cloudiness
   logical            :: micro_tend_output       = .false.   ! Default microphysics tendencies to output file
   logical            :: do_prescribed_CCN        = .false.   ! Use prescribed CCN
   contains
!===============================================================================
subroutine micro_p3_readnl(nlfile)

  use namelist_utils,  only: find_group_name
  use units,           only: getunit, freeunit
  use mpishorthand

  character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

  ! Local variables
  integer :: unitn, ierr
  character(len=*), parameter :: subname = 'micro_p3_cam_readnl'

  namelist /micro_nl/ &
       micro_p3_tableversion, micro_p3_lookup_dir, micro_aerosolactivation, micro_subgrid_cloud, &
       micro_tend_output, p3_qc_autocon_expon, p3_qc_accret_expon, do_prescribed_CCN

  !-----------------------------------------------------------------------------

  if (masterproc) then
     unitn = getunit()
     open( unitn, file=trim(nlfile), status='old' )
     call find_group_name(unitn, 'micro_nl', status=ierr)
     if (ierr == 0) then
        read(unitn, micro_nl, iostat=ierr)
        if (ierr /= 0) then
           call endrun(subname // ':: ERROR reading namelist')
        end if
     end if
     close(unitn)
     call freeunit(unitn)

     write(iulog,'(A50)') ' ----- P3 Namelist Values: -----'
     write(iulog,'(A29,1x,A19)')  'micro_p3_tableversion: ',   micro_p3_tableversion
     write(iulog,'(A20,1x,A100)') 'micro_p3_lookup_dir: ',     micro_p3_lookup_dir
     write(iulog,'(A30,1x,L)')    'micro_aerosolactivation: ', micro_aerosolactivation
     write(iulog,'(A30,1x,L)')    'micro_subgrid_cloud: ',     micro_subgrid_cloud
     write(iulog,'(A30,1x,L)')    'micro_tend_output: ',       micro_tend_output
     write(iulog,'(A30,1x,8e12.4)') 'p3_qc_autocon_expon',        p3_qc_autocon_expon
     write(iulog,'(A30,1x,8e12.4)') 'p3_qc_accret_expon',         p3_qc_accret_expon
     write(iulog,'(A30,1x,L)')    'do_prescribed_CCN: ',       do_prescribed_CCN

  end if

#ifdef SPMD
  ! Broadcast namelist variables
  call mpibcast(micro_p3_tableversion,   len(micro_p3_tableversion), mpichar, 0, mpicom)
  call mpibcast(micro_p3_lookup_dir,     len(micro_p3_lookup_dir),   mpichar, 0, mpicom)
  call mpibcast(micro_aerosolactivation, 1,                          mpilog,  0, mpicom)
  call mpibcast(micro_subgrid_cloud,     1,                          mpilog,  0, mpicom)
  call mpibcast(micro_tend_output,       1,                          mpilog,  0, mpicom)
  call mpibcast(p3_qc_autocon_expon,      1,                          mpir8,   0, mpicom)
  call mpibcast(p3_qc_accret_expon,       1,                          mpir8,   0, mpicom)
  call mpibcast(do_prescribed_CCN,     1,                          mpilog,  0,    mpicom)

#endif

  ! Check to make sure p3 table version is valid
  select case (trim(micro_p3_tableversion))
    case ('4')
       ! Version 4 is valid
    case ('4.1')
      ! Version 4.1 is valid
    case ('4.1.1')
      ! Version 4.1.1 is valid
    case default
       print *, micro_p3_tableversion
       call bad_version_endrun()
  end select

  if (masterproc) write(iulog,'(A50)') ' ----- P3 READ NL Finshed: -----'

contains

  subroutine bad_version_endrun
    ! Endrun wrapper with a more useful error message.
    character(len=128) :: errstring
    write(errstring,*) "Invalid version number specified for P3 microphysics: ", &
         micro_p3_tableversion
    call endrun(errstring)
  end subroutine bad_version_endrun

end subroutine micro_p3_readnl
  !================================================================================================

  subroutine micro_p3_register()

  logical :: prog_modal_aero ! prognostic aerosols

  if (masterproc) write(iulog,'(A20)') ' P3 register start ...'

  call phys_getopts( prog_modal_aero_out   = prog_modal_aero )

   ncnst = 0
    ! Register Microphysics Constituents
    ! (i.e. members of state%q) and save indices.
    ! TODO make sure the cnst_names match what we think they are here.
    !================
   call cnst_add(cnst_names(1), mwdry, cpair, 0._rtype, ixcldliq, &
         longname='Grid box averaged cloud liquid amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(2), mwdry, cpair, 0._rtype, ixcldice, &
         longname='Grid box averaged cloud ice amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(3), mwh2o, cpair, 0._rtype, ixnumliq, &
         longname='Grid box averaged cloud liquid number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(4), mwh2o, cpair, 0._rtype, ixnumice, &
         longname='Grid box averaged cloud ice number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(5), mwh2o, cpair, 0._rtype, ixrain, &
         longname='Grid box averaged rain amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(6), mwh2o, cpair, 0._rtype, ixcldrim, &
         longname='Grid box averaged riming amount', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(7), mwh2o, cpair, 0._rtype, ixnumrain, &
         longname='Grid box averaged rain number', &
         is_convtran1=.true.)
   ncnst = ncnst + 1
   call cnst_add(cnst_names(8), mwh2o, cpair, 0._rtype, ixrimvol, &
         longname='Grid box averaged riming volume', &
         is_convtran1=.true.)
   ncnst = ncnst + 1

    ! Add Variables to Pbuf
    !================
   !! module microp_aero
   call pbuf_add_field('CLDO','global', dtype_r8,(/pcols,pver,dyn_time_lvls/),cldo_idx)

   !! module wetdep
   call pbuf_add_field('QME',  'physpkg',dtype_r8,(/pcols,pver/), qme_idx)
   call pbuf_add_field('PRAIN','physpkg',dtype_r8,(/pcols,pver/), precip_total_tend_idx)
   call pbuf_add_field('NEVAPR','physpkg',dtype_r8,(/pcols,pver/), nevapr_idx)

   !! module aero_model
   if (prog_modal_aero) then
      call pbuf_add_field('RATE1_CW2PR_ST','physpkg',dtype_r8,(/pcols,pver/),rate1_cw2pr_st_idx)
   endif

   !! module clubb_intr
   call pbuf_add_field('PRER_EVAP',  'global', dtype_r8,(/pcols,pver/), qr_evap_tend_idx)

   !! module radiation_data & module cloud_rad_props
   call pbuf_add_field('DEI',        'physpkg',dtype_r8,(/pcols,pver/), dei_idx)
   call pbuf_add_field('MU',         'physpkg',dtype_r8,(/pcols,pver/), mu_idx)
   call pbuf_add_field('LAMBDAC',    'physpkg',dtype_r8,(/pcols,pver/), lambdac_idx)

   !! module cospsimulator_intr
   call pbuf_add_field('REL',        'physpkg',dtype_r8,(/pcols,pver/), rel_idx)
   call pbuf_add_field('REI',        'physpkg',dtype_r8,(/pcols,pver/), rei_idx)
   call pbuf_add_field('LS_FLXPRC',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxprc_idx)
   call pbuf_add_field('LS_FLXSNW',  'physpkg',dtype_r8,(/pcols,pverp/), ls_flxsnw_idx)
   call pbuf_add_field('LS_REFFRAIN','physpkg',dtype_r8,(/pcols,pver/), ls_reffrain_idx)
   call pbuf_add_field('LS_REFFSNOW','physpkg',dtype_r8,(/pcols,pver/), ls_reffsnow_idx)
   call pbuf_add_field('CV_REFFLIQ', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffliq_idx)
   call pbuf_add_field('CV_REFFICE', 'physpkg',dtype_r8,(/pcols,pver/), cv_reffice_idx)

   call pbuf_add_field('RELVAR',     'global',dtype_r8,(/pcols,pver/),   relvar_idx)
   call pbuf_add_field('ACCRE_ENHAN','global',dtype_r8,(/pcols,pver/), accre_enhan_idx)

   call pbuf_add_field('QV_PREV_WET',     'global',dtype_r8,(/pcols,pver/), qv_prev_idx)
   call pbuf_add_field('T_PREV',      'global',dtype_r8,(/pcols,pver/), t_prev_idx)

   if (masterproc) write(iulog,'(A20)') '    P3 register finished'
  end subroutine micro_p3_register

  !================================================================================================
  function micro_p3_implements_cnst(name)

    ! Return true if specified constituent is implemented by the
    ! microphysics package

    character(len=*), intent(in) :: name        ! constituent name
    logical :: micro_p3_implements_cnst    ! return value

    micro_p3_implements_cnst = any(name == cnst_names)

  end function micro_p3_implements_cnst


  !================================================================================================

  subroutine micro_p3_init_cnst(name, q)

    ! Initialize the microphysics constituents, if they are
    ! not read from the initial file.

    character(len=*), intent(in) :: name     ! constituent name
    real(rtype), intent(out) :: q(:,:)   ! mass mixing ratio (gcol, plev)

    if (micro_p3_implements_cnst(name)) q = 0.0_rtype

  end subroutine micro_p3_init_cnst

  !================================================================================================

  subroutine micro_p3_init(pbuf2d)
    use micro_p3,       only: p3_init
    use cam_history,    only: addfld, add_default, horiz_only
    use cam_history_support, only: add_hist_coord
    use micro_p3_utils, only: micro_p3_utils_init
    use read_spa_data,  only: ccn_names
    use read_spa_data,  only: is_spa_active
    use shr_log_mod,    only: errMsg => shr_log_errMsg

    type(physics_buffer_desc),  pointer :: pbuf2d(:,:)
    integer        :: m, mm
    integer        :: ierr
    logical :: history_amwg         ! output the variables used by the AMWG diag package
    logical :: history_verbose      ! produce verbose history output
    logical :: history_budget       ! Output tendencies and state variables for CAM4
    integer :: budget_histfile      ! output history file number for budget fields
                                   ! temperature, water vapor, cloud ice and cloud

    !sanity check for spa
    !spa must be active if do_prescribed_CCN is true
    if(do_prescribed_CCN .and. .not. is_spa_active) then
       call endrun('SPA must be active if do_prescribed_CCN is true, '//errmsg(__FILE__,__LINE__))
    endif

    call micro_p3_utils_init(cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice, &
             cpliq,tmelt,pi,iulog,masterproc)

    ! CALL P3 INIT:
    !==============
    !might want to add all E3SM parameter vals to p3_init call...
    ast_idx      = pbuf_get_index('AST') !! from CLUBB
    cmeliq_idx   = pbuf_get_index('CMELIQ') !! from CLUBB Rate of cond-evap of liq within the cloud

    !!
    !! for ice nucleation
    !!

    ni_activated_idx     = pbuf_get_index('NAAI') !! from microp
    npccn_idx    = pbuf_get_index('NPCCN')!! from microp

    prec_str_idx = pbuf_get_index('PREC_STR') !! from physpkg
    snow_str_idx = pbuf_get_index('SNOW_STR') !! from physpkg
    prec_sed_idx = pbuf_get_index('PREC_SED') !! from physpkg
    snow_sed_idx = pbuf_get_index('SNOW_SED') !! from physpkg
    prec_pcw_idx = pbuf_get_index('PREC_PCW') !! from physpkg
    snow_pcw_idx = pbuf_get_index('SNOW_PCW') !! from physpkg

    !Get indices for SPA treatment
    if(do_prescribed_CCN)ccn3_idx = pbuf_get_index(ccn_names(1))

    call p3_init(micro_p3_lookup_dir,micro_p3_tableversion)

    ! Initialize physics buffer grid fields for accumulating precip and
    ! condensation
    if (is_first_step()) then

       call pbuf_set_field(pbuf2d, cldo_idx,   0._rtype)
       call pbuf_set_field(pbuf2d, relvar_idx, 2._rtype)
       call pbuf_set_field(pbuf2d, accre_enhan_idx, micro_mg_accre_enhan_fac)
       call pbuf_set_field(pbuf2d, qr_evap_tend_idx,  0._rtype)
       call pbuf_set_field(pbuf2d, qv_prev_idx,  0._rtype)
       call pbuf_set_field(pbuf2d, t_prev_idx,  0._rtype)

    end if

    ! INITIALIZE OUTPUT
    !==============
    do m = 1, ncnst
       call cnst_get_ind(cnst_names(m), mm)
       if ( any(mm == (/ ixcldliq, ixcldice, ixrain, ixcldrim /)) ) then
          ! mass mixing ratios
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'kg/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', 'kg/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else if ( any(mm == (/ ixnumliq, ixnumice, ixnumrain /)) ) then
          ! number concentrations
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', '1/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', '1/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else if ( mm == ixrimvol ) then
          ! number concentrations
          call addfld(cnst_name(mm), (/ 'lev' /), 'A', 'm3/kg', &
            cnst_longname(mm) )
          call addfld(sflxnam(mm), horiz_only, 'A', 'm3/m2/s', &
            trim(cnst_name(mm))//' surface flux')
       else
          call endrun( "micro_p3_acme_init: &
               &Could not call addfld for constituent with unknown units.")
       endif
    end do
    call addfld(apcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' after physics'  )
    call addfld(apcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' after physics'  )
    call addfld(bpcnst(ixcldliq), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldliq))//' before physics' )
    call addfld(bpcnst(ixcldice), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldice))//' before physics' )
    call addfld(apcnst(ixrain),   (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' after physics'  )
    call addfld(bpcnst(ixrain),   (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixrain))//' before physics' )
    call addfld(apcnst(ixcldrim), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldrim))//' after physics'  )
    call addfld(bpcnst(ixcldrim), (/ 'lev' /), 'A', 'kg/kg', trim(cnst_name(ixcldrim))//' before physics' )

    ! microphysics cloud fraction fields
    call addfld('CLOUDFRAC_LIQ_MICRO', (/ 'lev' /), 'A', 'unitless', 'Grid box liquid cloud fraction in microphysics' )
    call addfld('CLOUDFRAC_ICE_MICRO', (/ 'lev' /), 'A', 'unitless', 'Grid box ice cloud fraction in microphysics' )
    call addfld('CLOUDFRAC_RAIN_MICRO', (/ 'lev' /), 'A', 'unitless', 'Grid box rain cloud fraction in microphysics' )

    call addfld ('CME', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of cond-evap within the cloud'                      )
    call addfld ('FICE', (/ 'lev' /), 'A', 'fraction', 'Fractional ice content within cloud'                     )
    call addfld ('ICWMRST', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-stratus water mixing ratio'                )
    call addfld ('ICIMRST', (/ 'lev' /), 'A', 'kg/kg', 'Prognostic in-stratus ice mixing ratio'                  )

   ! MG microphysics diagnostics
    call addfld ('QV2QI_DEPOS', (/ 'lev' /), 'A', 'kg/kg/s', 'Rate of deposition/sublimation of cloud ice'             )
    call addfld ('QCSEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Cloud water mixing ratio tendency from sedimentation'    )
    call addfld ('QISEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Cloud ice mixing ratio tendency from sedimentation'      )
    call addfld ('QRSEDTEN', (/ 'lev' /), 'A', 'kg/kg/s', 'Rain mixing ratio tendency from sedimentation'           )

   ! History variables for CAM5 microphysics
    call addfld ('ICWNC', (/ 'lev' /), 'A', 'm-3', 'Prognostic in-cloud water number conc'                   )
    call addfld ('ICINC', (/ 'lev' /), 'A', 'm-3', 'Prognostic in-cloud ice number conc'                     )
    call addfld ('CDNUMC', horiz_only,    'A', '1/m2', 'Vertically-integrated droplet concentration'             )
    call addfld ('MPICLWPI', horiz_only,    'A', 'kg/m2', 'Vertically-integrated &
         &in-cloud Initial Liquid WP (Before Micro)' )
    call addfld ('MPICIWPI', horiz_only,    'A', 'kg/m2', 'Vertically-integrated &
         &in-cloud Initial Ice WP (Before Micro)'    )

   ! Averaging for cloud particle number and size
   call addfld ('AWNC', (/ 'lev' /), 'A', 'm-3', 'Average cloud water number conc'                         )
   call addfld ('AWNI', (/ 'lev' /), 'A', 'm-3', 'Average cloud ice number conc'                           )
   call addfld ('AREL', (/ 'lev' /), 'A', 'Micron', 'Average droplet effective radius'                        )
   call addfld ('AREI', (/ 'lev' /), 'A', 'Micron', 'Average ice effective radius'                            )
   ! Frequency arrays for above
   call addfld ('FREQL', (/ 'lev' /), 'A', 'fraction', 'Fractional occurrence of liquid'                          )
   call addfld ('FREQI', (/ 'lev' /), 'A', 'fraction', 'Fractional occurrence of ice'                             )

   ! Average cloud top particle size and number (liq, ice) and frequency
   call addfld ('REL', (/ 'lev' /), 'A', 'micron', 'REL stratiform cloud effective radius liquid')
   call addfld ('REI', (/ 'lev' /), 'A', 'micron', 'REI stratiform cloud effective radius ice')

!!== KZ_DCS
   call addfld ('DCST',(/ 'lev' /), 'A','m','dcs')
!!== KZ_DCS
   ! diagnostic precip
   call addfld ('QRAIN',(/ 'lev' /), 'A','kg/kg','Diagnostic grid-mean rain mixing ratio'         )
   call addfld ('QSNOW',(/ 'lev' /), 'A','kg/kg','Diagnostic grid-mean snow mixing ratio'         )
   call addfld ('NRAIN',(/ 'lev' /), 'A','m-3','Diagnostic grid-mean rain number conc'         )
   call addfld ('NSNOW',(/ 'lev' /), 'A','m-3','Diagnostic grid-mean snow number conc'         )


   ! Aerosol information
   call addfld ('NCAL',(/ 'lev' /), 'A','1/m3','Number Concentation Activated for Liquid')
   call addfld ('NCAI',(/ 'lev' /), 'A','1/m3','Number Concentation Activated for Ice')

   ! Average rain and snow mixing ratio (Q), number (N) and diameter (D), with frequency
   call addfld ('AQRAIN',(/ 'lev' /), 'A','kg/kg','Average rain mixing ratio'         )
   call addfld ('AQSNOW',(/ 'lev' /), 'A','kg/kg','Average snow mixing ratio'         )
   call addfld ('ANRAIN',(/ 'lev' /), 'A','m-3','Average rain number conc'         )
   call addfld ('ANSNOW',(/ 'lev' /), 'A','m-3','Average snow number conc'         )
   call addfld ('ADRAIN',(/ 'lev' /), 'A','Micron','Average rain effective Diameter'         )
   call addfld ('ADSNOW',(/ 'lev' /), 'A','Micron','Average snow effective Diameter'         )
   call addfld ('FREQR',(/ 'lev' /), 'A','fraction','Fractional occurrence of rain'       )
   call addfld ('FREQS',(/ 'lev' /), 'A','fraction','Fractional occurrence of snow'       )

   ! precipitation efficiency & other diagnostic fields
   call addfld('UMR', (/ 'lev' /), 'A',   'm/s', 'Mass-weighted rain  fallspeed'              )

   ! Record of inputs/outputs from p3_main
   call add_hist_coord('P3_input_dim',  16, 'Input field dimension for p3_main subroutine',  'N/A', (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16 /))
   call add_hist_coord('P3_output_dim', 32, 'Output field dimension for p3_main subroutine', 'N/A', (/ 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32 /))
   call addfld('P3_input',  (/ 'ilev         ', 'P3_input_dim ' /),  'I', 'N/A', 'Inputs for p3_main subroutine')
   call addfld('P3_output', (/ 'ilev         ', 'P3_output_dim' /), 'I', 'N/A', 'Outputs for p3_main subroutine')
   ! Record of microphysics tendencies
   ! warm-phase process rates
   call addfld('P3_qrcon',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for rain condensation   (Not in paper?)')
   call addfld('P3_qc2qr_accret_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for cloud droplet accretion by rain')
   call addfld('P3_qc2qr_autoconv_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for cloud droplet autoconversion to rain')
   call addfld('P3_nc_accret_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in cloud droplet number from accretion by rain')
   call addfld('P3_nc2nr_autoconv_tend', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in cloud droplet number from autoconversion')
   call addfld('P3_nc_selfcollect_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in cloud droplet number from self-collection  (Not in paper?)')
   call addfld('P3_nr_selfcollect_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in rain number from self-collection  (Not in paper?)')
   call addfld('P3_nc_nuceat_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in cloud droplet number from activation of CCN')
   call addfld('P3_qccon',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for cloud droplet condensation')
   call addfld('P3_qcnuc',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for activation of cloud droplets from CCN')
   call addfld('P3_qr2qv_evap_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for rain evaporation')
   call addfld('P3_qcevp',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for cloud droplet evaporation')
   call addfld('P3_nr_evap_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in rain number from evaporation')
   call addfld('P3_ncautr', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in rain number from autoconversion of cloud water')
   ! ice-phase process rates
   call addfld('P3_qccol',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for collection of cloud water by ice')
   call addfld('P3_qwgrth', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 wet growth rate')
   call addfld('P3_qidep',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for vapor deposition')
   call addfld('P3_qrcol',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for collection rain mass by ice')
   call addfld('P3_qinuc',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for deposition/condensation freezing nuc')
   call addfld('P3_nc_collect_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in cloud droplet number from collection by ice')
   call addfld('P3_nr_collect_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in rain number from collection by ice')
   call addfld('P3_ni_nucleat_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in ice number from deposition/cond-freezing nucleation')
   call addfld('P3_qi2qv_sublim_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for sublimation of ice')
   call addfld('P3_qi2qr_melt_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for melting of ice')
   call addfld('P3_ni2nr_melt_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for melting of ice')
   call addfld('P3_ni_sublim_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in ice number from sublimation')
   call addfld('P3_ni_selfcollect_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for change in ice number from collection within a category (Not in paper?)')
   call addfld('P3_qc2qi_hetero_frz_tend', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for immersion freezing droplets')
   call addfld('P3_qr2qi_immers_frz_tend', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for immersion freezing rain')
   call addfld('P3_nc2ni_immers_frz_tend', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for immersion freezing droplets')
   call addfld('P3_nr2ni_immers_frz_tend', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for immersion freezing rain')
   call addfld('P3_nr_ice_shed_tend', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for source for rain number from collision of rain/ice above freezing and shedding')
   call addfld('P3_qc2qr_ice_shed_tend',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for source for rain mass due to cloud water/ice collision above freezing and shedding or wet growth and shedding')
   call addfld('P3_ncshdc', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for source for rain number due to cloud water/ice collision above freezing  and shedding (combined with NRSHD in the paper)')
   ! Sedimentation
   call addfld('P3_sed_CLDLIQ',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for liquid cloud content due to sedimentation')
   call addfld('P3_sed_NUMLIQ',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for liquid cloud number due to sedimentation')
   call addfld('P3_sed_CLDRAIN', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for rain cloud content due to sedimentation')
   call addfld('P3_sed_NUMRAIN', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for rain cloud number due to sedimentation')
   call addfld('P3_sed_CLDICE',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for ice cloud content due to sedimentation')
   call addfld('P3_sed_NUMICE',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for ice cloud number due to sedimentation')
   ! Microphysics Processes
   call addfld('P3_mtend_CLDLIQ',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for liquid cloud content due to micro processes')
   call addfld('P3_mtend_NUMLIQ',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for liquid cloud number due to micro processes')
   call addfld('P3_mtend_CLDRAIN', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for rain cloud content due to micro processes')
   call addfld('P3_mtend_NUMRAIN', (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for rain cloud number due to micro processes')
   call addfld('P3_mtend_CLDICE',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for ice cloud content due to micro processes')
   call addfld('P3_mtend_NUMICE',  (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for ice cloud number due to micro processes')
   call addfld('P3_mtend_Q',       (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for water vapor due to micro processes')
   call addfld('P3_mtend_TH',      (/ 'lev' /), 'A', 'kg/kg/s', 'P3 Tendency for potential temp. number due to micro processes')
   ! phase change tendencies
   call addfld('vap_liq_exchange',  (/ 'lev' /), 'A', 'kg/kg/s', 'Tendency for conversion from/to vapor phase to/from liquid phase')
   call addfld('vap_ice_exchange',  (/ 'lev' /), 'A', 'kg/kg/s', 'Tendency for conversion from/to vapor phase to/from frozen phase')
   call addfld('liq_ice_exchange',  (/ 'lev' /), 'A', 'kg/kg/s', 'Tendency for conversion from/to liquid phase to/from frozen phase')

   ! determine the add_default fields
   call phys_getopts(history_amwg_out           = history_amwg         , &
                     history_verbose_out        = history_verbose      , &
                     history_budget_out         = history_budget       , &
                     history_budget_histfile_num_out = budget_histfile)
   if (history_amwg) then
      call add_default ('FICE    ', 1, ' ')
      call add_default ('AQRAIN   ', 1, ' ')
      call add_default ('AQSNOW   ', 1, ' ')
      call add_default ('ANRAIN   ', 1, ' ')
      call add_default ('ANSNOW   ', 1, ' ')
      if (history_verbose) then
         call add_default ('ADRAIN   ', 1, ' ')
         call add_default ('ADSNOW   ', 1, ' ')
      endif
      call add_default ('AREI     ', 1, ' ')
      call add_default ('AREL     ', 1, ' ')
      call add_default ('AWNC     ', 1, ' ')
      call add_default ('AWNI     ', 1, ' ')
      call add_default ('CDNUMC   ', 1, ' ')
      call add_default ('FREQR    ', 1, ' ')
      call add_default ('FREQS    ', 1, ' ')
      call add_default ('FREQL    ', 1, ' ')
      call add_default ('FREQI    ', 1, ' ')
      do m = 1, ncnst
         call cnst_get_ind(cnst_names(m), mm)
         call add_default(cnst_name(mm), 1, ' ')
         ! call add_default(precip_ice_fluxnam(mm),   1, ' ')
      end do
      ! Microphysics cloud fractions
      call add_default ('CLOUDFRAC_LIQ_MICRO', 1, ' ')
      call add_default ('CLOUDFRAC_ICE_MICRO', 1, ' ')
      call add_default ('CLOUDFRAC_RAIN_MICRO', 1, ' ')
      ! Phase change tendencies
      call add_default('vap_liq_exchange',  1, ' ')
      call add_default('vap_ice_exchange',  1, ' ')
      call add_default('liq_ice_exchange',  1, ' ')
      ! Microphysics tendencies
      ! warm-phase process rates
      if (micro_tend_output) then
         call add_default('P3_qrcon',  1, ' ')
         call add_default('P3_qc2qr_accret_tend',  1, ' ')
         call add_default('P3_qc2qr_autoconv_tend',  1, ' ')
         call add_default('P3_nc_accret_tend',  1, ' ')
         call add_default('P3_nc2nr_autoconv_tend', 1, ' ')
         call add_default('P3_nc_selfcollect_tend',  1, ' ')
         call add_default('P3_nr_selfcollect_tend',  1, ' ')
         call add_default('P3_nc_nuceat_tend',  1, ' ')
         call add_default('P3_qccon',  1, ' ')
         call add_default('P3_qcnuc',  1, ' ')
         call add_default('P3_qr2qv_evap_tend',  1, ' ')
         call add_default('P3_qcevp',  1, ' ')
         call add_default('P3_nr_evap_tend',  1, ' ')
         call add_default('P3_ncautr', 1, ' ')
         ! ice-phase process rates
         call add_default('P3_qccol',  1, ' ')
         call add_default('P3_qwgrth', 1, ' ')
         call add_default('P3_qidep',  1, ' ')
         call add_default('P3_qrcol',  1, ' ')
         call add_default('P3_qinuc',  1, ' ')
         call add_default('P3_nc_collect_tend',  1, ' ')
         call add_default('P3_nr_collect_tend',  1, ' ')
         call add_default('P3_ni_nucleat_tend',  1, ' ')
         call add_default('P3_qi2qv_sublim_tend',  1, ' ')
         call add_default('P3_qi2qr_melt_tend',  1, ' ')
         call add_default('P3_ni2nr_melt_tend',  1, ' ')
         call add_default('P3_ni_sublim_tend',  1, ' ')
         call add_default('P3_ni_selfcollect_tend',  1, ' ')
         call add_default('P3_qc2qi_hetero_frz_tend', 1, ' ')
         call add_default('P3_qr2qi_immers_frz_tend', 1, ' ')
         call add_default('P3_nc2ni_immers_frz_tend', 1, ' ')
         call add_default('P3_nr2ni_immers_frz_tend', 1, ' ')
         call add_default('P3_nr_ice_shed_tend', 1, ' ')
         call add_default('P3_qc2qr_ice_shed_tend',  1, ' ')
         call add_default('P3_ncshdc', 1, ' ')
         ! Sedimentation
         call add_default('P3_sed_CLDLIQ',  1, ' ')
         call add_default('P3_sed_NUMLIQ',  1, ' ')
         call add_default('P3_sed_CLDRAIN', 1, ' ')
         call add_default('P3_sed_NUMRAIN', 1, ' ')
         call add_default('P3_sed_CLDICE',  1, ' ')
         call add_default('P3_sed_NUMICE',  1, ' ')
         ! Microphysics Processes
         call add_default('P3_mtend_CLDLIQ',  1, ' ')
         call add_default('P3_mtend_NUMLIQ',  1, ' ')
         call add_default('P3_mtend_CLDRAIN', 1, ' ')
         call add_default('P3_mtend_NUMRAIN', 1, ' ')
         call add_default('P3_mtend_CLDICE',  1, ' ')
         call add_default('P3_mtend_NUMICE',  1, ' ')
         call add_default('P3_mtend_Q',       1, ' ')
         call add_default('P3_mtend_TH',      1, ' ')
      end if
   end if

  end subroutine micro_p3_init

  !================================================================================================
    subroutine get_cloud_fraction(its,ite,kts,kte,ast,qc,qr,qi,method, &
                  cld_frac_i,cld_frac_l,cld_frac_r)

       use micro_p3_utils, only: mincld, qsmall

       integer,intent(in)                                 :: its,ite,kts,kte
       real(rtype),dimension(its:ite,kts:kte),intent(in)  :: ast, qc, qr, qi
       character(len=16),intent(in)                       :: method
       real(rtype),dimension(its:ite,kts:kte),intent(out) :: cld_frac_i, cld_frac_l, cld_frac_r
       real(rtype),dimension(its:ite,kts:kte)             :: cldm

       integer  :: i,k
       integer  :: ktop, kbot, kdir

       call t_startf('micro_p3_get_cloud_fraction')
       ktop = kts        !k of top level
       kbot = kte        !k of bottom level
       kdir = -1         !(k: 1=top, nk=bottom)

       cldm(:,:)  = mincld
       cld_frac_i(:,:) = mincld
       cld_frac_l(:,:) = mincld
       do k = kbot,ktop,kdir
          do i=its,ite
             cldm(i,k)  = max(ast(i,k), mincld)
             cld_frac_i(i,k) = max(ast(i,k), mincld)
             cld_frac_l(i,k) = max(ast(i,k), mincld)
             cld_frac_r(i,k) = cldm(i,k)
          end do
       end do

       !!
       !! precipitation fraction
       !!
       IF (trim(method) == 'in_cloud') THEN
          DO k = ktop-kdir,kbot,-kdir
             DO i=its,ite
                ! in_cloud means that precip_frac (cld_frac_r) = cloud (cldm) frac when cloud mass
                ! is present. Below cloud, precip frac is equal to the cloud
                ! fraction from the last layer that had cloud. Since presence or
                ! absence of cloud is defined as mass > qsmall, sub-cloud precip
                ! frac for the in_cloud method tends to be very small and is
                ! very sensitive to tiny changes in condensate near cloud base.
                IF (qc(i,k) .lt. qsmall .and. qi(i,k) .lt. qsmall) THEN
                   ! max(cld_frac_r above and cld_frac_r for this layer) is taken here
                   ! because code is incapable of handling cld_frac_r<cldm for a
                   ! given grid cell
                   cld_frac_r(i,k) = max(cld_frac_r(i,k+kdir),cld_frac_r(i,k))
                END IF
             END DO !i
          END DO !k
       ELSE IF (trim(method) == 'max_overlap') THEN
       ! max overlap is the max cloud fraction in all layers above which are
       ! connected to this one by a continuous band of precip mass. If
       ! there's no precip mass falling into a cell, it's precip frac is equal
       ! to the cloud frac, which is probably ~zero.

       ! IF rain or ice mix ratios are smaller than threshold,
       ! then leave cld_frac_r as cloud fraction at current level
          DO k = ktop-kdir,kbot,-kdir
             DO i=its,ite
                IF (qr(i,k+kdir) .ge. qsmall .or. qi(i,k+kdir) .ge. qsmall) THEN
                   cld_frac_r(i,k) = max(cld_frac_r(i,k+kdir),cld_frac_r(i,k))
                END IF
             END DO ! i
          END DO ! k
       END IF


       call t_stopf('micro_p3_get_cloud_fraction')
       return
    end subroutine get_cloud_fraction

  !================================================================================================
  subroutine micro_p3_tend(state, ptend, dtime, pbuf)

    use phys_grid,      only: get_rlat_all_p, get_rlon_all_p, get_gcol_all_p
    use time_manager,   only: get_nstep
    use cam_history,    only: outfld
    use time_manager,   only: get_nstep
    use micro_p3,       only: p3_main
    use micro_p3_utils, only: avg_diameter, &
                              rho_h2o, &
                              rho_h2os, &
                              qsmall, &
                              mincld, &
                              inv_cp
    use physics_utils, only: calculate_drymmr_from_wetmmr, calculate_wetmmr_from_drymmr

    !INPUT/OUTPUT VARIABLES
    type(physics_state),         intent(in)    :: state
    type(physics_ptend),         intent(out)   :: ptend
    real(rtype),                 intent(in)    :: dtime
    type(physics_buffer_desc),   pointer       :: pbuf(:)
    logical :: lq(pcnst)   !list of what constituents to update

    !INTERNAL VARIABLES
    real(rtype) :: dz(pcols,pver)        !geometric layer thickness              m
    real(rtype) :: cldliq(pcols,pver)     !cloud liquid water mixing ratio        kg/kg
    real(rtype) :: numliq(pcols,pver)     !cloud liquid water drop concentraiton  #/kg
    real(rtype) :: rain(pcols,pver)       !rain water mixing ratio                kg/kg
    real(rtype) :: numrain(pcols,pver)    !rain water number concentration        #/kg
    real(rtype) :: qv_dry(pcols,pver)     !dry water vapor mixing ratio           kg/kg
    real(rtype) :: qv_prev_dry(pcols,pver)!dry water vapor mixing ratio(previous time step) kg/kg
    real(rtype) :: qv_wet_in(pcols,pver)  !wet water vapor mixing ratio (input for P3) kg/kg
    real(rtype) :: qv_wet_out(pcols,pver) !wet water vapor mixing ratio (output from P3) kg/kg
    real(rtype) :: ice(pcols,pver)        !total ice water mixing ratio           kg/kg
    real(rtype) :: qm(pcols,pver)      !rime ice mixing ratio                  kg/kg
    real(rtype) :: numice(pcols,pver)     !total ice crystal number concentration #/kg
    real(rtype) :: rimvol(pcols,pver)     !rime volume mixing ratio               m3/kg
    real(rtype) :: temp(pcols,pver)       !temperature copy needed for tendency   K
    real(rtype) :: th(pcols,pver)         !potential temperature                  K
    real(rtype) :: precip_liq_surf(pcols)         !precipitation rate, liquid             m s-1
    real(rtype) :: precip_ice_surf(pcols)         !precipitation rate, solid              m s-1

    real(rtype) :: rho_qi(pcols,pver)  !bulk density of ice                    kg m-1
    real(rtype) :: pres(pcols,pver)       !pressure at midlevel                   hPa
    real(rtype) :: qv2qi_depos_tend(pcols,pver)
    real(rtype) :: precip_liq_flux(pcols,pver+1)     !grid-box average rain flux (kg m^-2s^-1) pverp
    real(rtype) :: precip_ice_flux(pcols,pver+1)     !grid-box average ice/snow flux (kg m^-2s^-1) pverp
    real(rtype) :: inv_exner(pcols,pver)       !inverse exner formula for converting between potential and normal temp
    real(rtype) :: cld_frac_r(pcols,pver)      !rain cloud fraction
    real(rtype) :: cld_frac_l(pcols,pver)      !liquid cloud fraction
    real(rtype) :: cld_frac_i(pcols,pver)      !ice cloud fraction
    real(rtype) :: tend_out(pcols,pver,49) !microphysical tendencies
    real(rtype), dimension(pcols,pver) :: liq_ice_exchange ! sum of liq-ice phase change tendenices
    real(rtype), dimension(pcols,pver) :: vap_liq_exchange ! sum of vap-liq phase change tendenices
    real(rtype), dimension(pcols,pver) :: vap_ice_exchange ! sum of vap-ice phase change tendenices
    real(rtype) :: dummy_out(pcols,pver)    ! dummy_output variable for p3_main to replace unused variables.

    !Prescribed CCN concentration
    real(rtype), dimension(pcols,pver) :: nccn_prescribed

    ! PBUF Variables
    real(rtype), pointer :: ast(:,:)      ! Relative humidity cloud fraction
    real(rtype), pointer :: ni_activated(:,:)     ! ice nucleation number
    real(rtype), pointer :: npccn(:,:)    ! liquid activation number tendency
    real(rtype), pointer :: cmeliq(:,:)
    !!
    real(rtype), pointer :: prec_str(:)    ! [Total] Sfc flux of precip from stratiform [ m/s ]
    real(rtype), pointer :: prec_sed(:)    ! Surface flux of total cloud water from sedimentation
    real(rtype), pointer :: prec_pcw(:)    ! Sfc flux of precip from microphysics [ m/s ]
    real(rtype), pointer :: snow_str(:)    ! [Total] Sfc flux of snow from stratiform   [ m/s ]
    real(rtype), pointer :: snow_pcw(:)    ! Sfc flux of snow from microphysics [ m/s ]
    real(rtype), pointer :: snow_sed(:)    ! Surface flux of cloud ice from sedimentation
    real(rtype), pointer :: relvar(:,:)    ! cloud liquid relative variance [-]
    real(rtype), pointer :: cldo(:,:)      ! Old cloud fraction
    real(rtype), pointer :: qr_evap_tend(:,:) ! precipitation evaporation rate
    real(rtype), pointer :: qv_prev_wet(:,:)   ! qv from previous p3_main call
    real(rtype), pointer :: t_prev(:,:)    ! t from previous p3_main call
    !! wetdep
    real(rtype), pointer :: qme(:,:)
    real(rtype), pointer :: precip_total_tend(:,:)        ! Total precipitation (rain + snow)
    real(rtype), pointer :: nevapr(:,:)       ! Evaporation of total precipitation (rain + snow)
    !! COSP simulator
    real(rtype), pointer :: rel(:,:)          ! Liquid effective drop radius (microns)
    real(rtype), pointer :: rei(:,:)          ! Ice effective drop size (microns)
    real(rtype), pointer :: flxprc(:,:)     ! P3 grid-box mean flux_large_scale_cloud_rain+snow at interfaces (kg/m2/s)
    real(rtype), pointer :: flxsnw(:,:)     ! P3 grid-box mean flux_large_scale_cloud_snow at interfaces (kg/m2/s)
    real(rtype), pointer :: reffrain(:,:)   ! P3 diagnostic rain effective radius (um)
    real(rtype), pointer :: reffsnow(:,:)   ! P3 diagnostic snow effective radius (um)
    real(rtype), pointer :: cvreffliq(:,:)    ! convective cloud liquid effective radius (um)
    real(rtype), pointer :: cvreffice(:,:)    ! convective cloud ice effective radius (um)
    !! radiation
    real(rtype), pointer :: dei(:,:)          ! Ice effective diameter (um)
    real(rtype), pointer :: mu(:,:)           ! Size distribution shape parameter for radiation
    real(rtype), pointer :: lambdac(:,:)      ! Size distribution slope parameter for radiation
    ! DONE PBUF
    ! For recording inputs/outputs to p3_main
    real(rtype) :: p3_main_inputs(pcols,pver+1,17) ! Record of inputs for p3_main
    real(rtype) :: p3_main_outputs(pcols,pver+1,31) ! Record of outputs for p3_main

    ! Derived Variables
    real(rtype) :: icimrst(pcols,pver) ! stratus ice mixing ratio - on grid
    real(rtype) :: icwmrst(pcols,pver) ! stratus water mixing ratio - on grid
    real(rtype) :: rho(pcols,pver)
    real(rtype) :: drout2(pcols,pver)
    real(rtype) :: reff_rain(pcols,pver)
    real(rtype) :: col_location(pcols,3),tmp_loc(pcols)  ! Array of column lon (index 1) and lat (index 2)
    integer     :: tmpi_loc(pcols) ! Global column index temp array

    ! Variables used for microphysics output
    real(rtype) :: aqrain(pcols,pver)
    real(rtype) :: anrain(pcols,pver)
    real(rtype) :: nfice(pcols,pver)
    real(rtype) :: efcout(pcols,pver)
    real(rtype) :: efiout(pcols,pver)
    real(rtype) :: ncout(pcols,pver)
    real(rtype) :: niout(pcols,pver)
    real(rtype) :: freqr(pcols,pver)
    real(rtype) :: freql(pcols,pver)
    real(rtype) :: freqi(pcols,pver)
    real(rtype) :: cdnumc(pcols)
    real(rtype) :: icinc(pcols,pver)
    real(rtype) :: icwnc(pcols,pver)


    integer :: it                      !timestep counter                       -
    integer :: its, ite                !horizontal bounds (column start,finish)
    integer :: kts                     !closest level to TOM                   -
    integer :: kte                     !near surface level                     -

    logical :: do_predict_nc           !prognostic droplet concentration or not?
    logical :: do_subgrid_clouds       !use subgrid cloudiness in tendency calculations?
    integer :: icol, ncol, k
    integer :: psetcols, lchnk
    integer :: itim_old
    real(rtype) :: T_virtual

    ! For rrtmg optics. specified distribution.
    real(rtype), parameter :: dcon   = 25.e-6_rtype         ! Convective size distribution effective radius (um)
    real(rtype), parameter :: mucon  = 5.3_rtype            ! Convective size distribution shape parameter
    real(rtype), parameter :: deicon = 50._rtype            ! Convective ice effective diameter (um)
    real(rtype), pointer :: ccn_trcdat(:,:) !BSINGH - receive ccn values in this variable

    call t_startf('micro_p3_tend_init')

    psetcols = state%psetcols
    lchnk = state%lchnk

    !+++ Aaron Donahue
    itim_old = pbuf_old_tim_idx()

    !============================
    ! All external PBUF variables:
    ! INPUTS
    call pbuf_get_field(pbuf, ast_idx,         ast, start=(/1,1,itim_old/), kount=(/psetcols,pver,1/))
    call pbuf_get_field(pbuf, ni_activated_idx,        ni_activated                                                  )
    call pbuf_get_field(pbuf, npccn_idx,       npccn                                                 )
    call pbuf_get_field(pbuf, cmeliq_idx,      cmeliq                                                )
    ! OUTPUTS
    call pbuf_get_field(pbuf,    prec_str_idx,  prec_str)
    call pbuf_get_field(pbuf,    snow_str_idx,  snow_str)
    call pbuf_get_field(pbuf,    prec_sed_idx,  prec_sed)
    call pbuf_get_field(pbuf,    snow_sed_idx,  snow_sed)
    call pbuf_get_field(pbuf,    prec_pcw_idx,  prec_pcw)
    call pbuf_get_field(pbuf,    snow_pcw_idx,  snow_pcw)
    !============================
    ! All internal PBUF variables
    ! INPUTS
    call pbuf_get_field(pbuf,      relvar_idx,    relvar                                                   )
    call pbuf_get_field(pbuf,      t_prev_idx,    t_prev                                                   )
    call pbuf_get_field(pbuf,     qv_prev_idx,    qv_prev_wet                                                  )
    ! OUTPUTS
    call pbuf_get_field(pbuf,        cldo_idx,      cldo, start=(/1,1,itim_old/), kount=(/psetcols,pver,1/))
    call pbuf_get_field(pbuf,         qme_idx,       qme                                                   )
    call pbuf_get_field(pbuf,       precip_total_tend_idx,     precip_total_tend                                                   )
    call pbuf_get_field(pbuf,      nevapr_idx,    nevapr                                                   )
    call pbuf_get_field(pbuf,   qr_evap_tend_idx, qr_evap_tend                                                   )
    call pbuf_get_field(pbuf,         rei_idx,       rei                                                   ) ! ice eff. rad
    call pbuf_get_field(pbuf,         rel_idx,       rel                                                   ) ! liq. eff. rad
    call pbuf_get_field(pbuf,         dei_idx,       dei                                                   )
    call pbuf_get_field(pbuf,          mu_idx,        mu                                                   )
    call pbuf_get_field(pbuf,     lambdac_idx,   lambdac                                                   )
    call pbuf_get_field(pbuf,   ls_flxprc_idx,    flxprc                                                   )
    call pbuf_get_field(pbuf,   ls_flxsnw_idx,    flxsnw                                                   )
    call pbuf_get_field(pbuf, ls_reffrain_idx,  reffrain                                                   )
    call pbuf_get_field(pbuf, ls_reffsnow_idx,  reffsnow                                                   )
    call pbuf_get_field(pbuf,  cv_reffliq_idx, cvreffliq                                                   )
    call pbuf_get_field(pbuf,  cv_reffice_idx, cvreffice                                                   )

    ncol = state%ncol
    !==============
    ! Some pre-microphysics INITIALIZATION
    !==============
    cldo(:ncol,top_lev:pver)=ast(:ncol,top_lev:pver)

    ! INITIALIZE PTEND
    !==============
    !ptend is an output variable. Since not substepping in micro, don't need
    !a local copy.
    lq            = .false. !initialize all constituents to false.
    lq(1)         = .true.
    lq(ixcldliq)  = .true.
    lq(ixcldice)  = .true.
    lq(ixnumliq)  = .true.
    lq(ixnumice)  = .true.
    lq(ixrain)    = .true.
    lq(ixcldrim)  = .true.
    lq(ixnumrain) = .true.
    lq(ixrimvol)  = .true.
    call physics_ptend_init(ptend, psetcols, "micro_p3", ls=.true., lq=lq)

    ! HANDLE AEROSOL ACTIVATION
    !==============
    do_predict_nc = micro_aerosolactivation

    ! ASSIGN TOP AND BOTTOM INDICES FOR GRID
    !==============
    !kts is closest level to top of model. Instead of 1 (top-of-model),
    !we use the previously-defined trop_cloud_top_lev to reduce the number of
    !levels we need to calculate and to avoid upper-atmos regions where this
    !micro-physics is inappropriate. kte is the near-surface level = pver.

    kts=top_lev
    kte=pver

    ! HANDLE TIMESTEP COUNTER, GET LON/LAT values (used by error warning code)
    !==============
    !p3 wants to know the timestep number (which it calls "it") because it
    !handles things differently on the first step, where it doesn't have values
    !yet. E3SM has a handy function for deciding if this is the first step, so
    !we hack "it" with "is_first_step()" for now. Eventually, we should replace
    !"it" with a logical.
    it = get_nstep()
    tmp_loc =-999.0_rtype
    call get_rlon_all_p(lchnk,ncol,tmp_loc)
    col_location(:ncol,2) = tmp_loc(:ncol)*180.0_rtype/pi
    call get_rlat_all_p(lchnk,ncol,tmp_loc)
    col_location(:ncol,3) = tmp_loc(:ncol)*180.0_rtype/pi
    call get_gcol_all_p(lchnk,ncol,tmpi_loc)
    col_location(:ncol,1) = real(tmpi_loc(:ncol))

    ! MAKE LOCAL COPIES OF VARS MODIFIED BY P3
    !==============
    !local copies are needed because state is passed into this routine as intent=in
    !while P3 seeks to modify state variables in-place. Also, we need a copy of
    !old values in order to back out ptend values later. Traditionally, a local copy
    !is created by copying the whole state. It is much cheaper to just copy the
    !variables we need.

    !---------------------------------------------------------------------------------------
    !Wet to dry mixing ratios:
    !-------------------------
    !Since state constituents from the host model are  wet mixing ratios and P3 needs these
    !constituents in dry mixing ratios, we convert the wet mixing ratios to dry mixing ratio
    !while assigning state constituents to the local variables
    !NOTE:Function calculate_drymmr_from_wetmmr takes 3 arguments: (number of columns, wet mmr and
    ! "wet" water vapor mixing ratio)
    !---------------------------------------------------------------------------------------
    qv_wet_in   = state%q(:,:,1) ! Get "wet" water vapor mixing ratio from state
    !Compute dry mixing ratios for all the constituents
    qv_dry(:ncol,:pver)      = calculate_drymmr_from_wetmmr(ncol, pver, qv_wet_in,              qv_wet_in)
    qv_prev_dry(:ncol,:pver) = calculate_drymmr_from_wetmmr(ncol, pver, qv_prev_wet,            qv_wet_in)
    cldliq(:ncol,:pver)      = calculate_drymmr_from_wetmmr(ncol, pver, state%q(:,:,ixcldliq),  qv_wet_in)
    numliq(:ncol,:pver)      = calculate_drymmr_from_wetmmr(ncol, pver, state%q(:,:,ixnumliq),  qv_wet_in)
    rain(:ncol,:pver)        = calculate_drymmr_from_wetmmr(ncol, pver, state%q(:,:,ixrain),    qv_wet_in)
    numrain(:ncol,:pver)     = calculate_drymmr_from_wetmmr(ncol, pver, state%q(:,:,ixnumrain), qv_wet_in)
    ice(:ncol,:pver)         = calculate_drymmr_from_wetmmr(ncol, pver, state%q(:,:,ixcldice),  qv_wet_in)
    qm(:ncol,:pver)          = calculate_drymmr_from_wetmmr(ncol, pver, state%q(:,:,ixcldrim),  qv_wet_in) !Aaron, changed ixqm to ixcldrim to match Kai's code
    numice(:ncol,:pver)      = calculate_drymmr_from_wetmmr(ncol, pver, state%q(:,:,ixnumice),  qv_wet_in)
    rimvol(:ncol,:pver)      = calculate_drymmr_from_wetmmr(ncol, pver, state%q(:,:,ixrimvol),  qv_wet_in)

    ! COMPUTE GEOMETRIC THICKNESS OF GRID & CONVERT T TO POTENTIAL TEMPERATURE
    !==============
    ! TODO: Create a general function to calculate Exner's formula that can be
    ! used by all parameterizations, such as P3 and SHOC.
    ! This would take a bit more work, so we have decided to delay this task
    ! until a later stage of code cleanup.
    inv_exner(:ncol,:pver) = 1._rtype/((state%pmiddry(:ncol,:pver)*1.e-5_rtype)**(rair*inv_cp))
    do icol = 1,ncol
       do k = 1,pver
          ! Note, there is a state%zi variable that could be used to calculate
          ! dz, but that is in a wet coordinate frame rather than dry.  Now that
          ! P3 is using dry MMR we instead calculated dz using virtual
          ! temperature and pressure.
          T_virtual  = state%t(icol,k) * (1.0 + qv_dry(icol,k)*(1.0*mwdry/mwh2o - 1.0))
          dz(icol,k) = (rair/gravit) * state%pdeldry(icol,k) * T_virtual / state%pmiddry(icol,k) 
          th(icol,k) = state%t(icol,k)*inv_exner(icol,k) 
       end do
    end do

    its     = 1
    ite     = state%ncol
    kts     = 1
    kte     = pver
    pres    = state%pmiddry(:,:)
    ! Initialize the raidation dependent variables.
    mu      = 0.0_rtype !mucon
    lambdac = 0.0_rtype !(mucon + 1._rtype)/dcon
    dei     = 50.0_rtype !deicon
    ! Determine the cloud fraction and precip cover
    cld_frac_i(:,:) = 1.0_rtype
    cld_frac_l(:,:) = 1.0_rtype
    cld_frac_r(:,:) = 1.0_rtype
    do_subgrid_clouds = micro_subgrid_cloud
    if (do_subgrid_clouds) &
        call get_cloud_fraction(its,ite,kts,kte,ast(its:ite,kts:kte),cldliq(its:ite,kts:kte), &
                rain(its:ite,kts:kte),ice(its:ite,kts:kte),precip_frac_method, &
                cld_frac_i(its:ite,kts:kte),cld_frac_l(its:ite,kts:kte),cld_frac_r(its:ite,kts:kte))
    call t_stopf('micro_p3_tend_init')

    p3_main_inputs(:,:,:) = -999._rtype
    do k = 1,pver
      p3_main_inputs(1,k,1)  = ast(1,k)
      p3_main_inputs(1,k,2)  = ni_activated(1,k)
      p3_main_inputs(1,k,3)  = npccn(1,k)
      p3_main_inputs(1,k,4)  = pres(1,k)
      p3_main_inputs(1,k,5)  = state%zi(1,k)
      p3_main_inputs(1,k,6)  = state%T(1,k)
      p3_main_inputs(1,k,7)  = qv_dry(1,k)
      p3_main_inputs(1,k,8)  = cldliq(1,k)
      p3_main_inputs(1,k,9)  = ice(1,k)
      p3_main_inputs(1,k,10) = numliq(1,k)
      p3_main_inputs(1,k,11) = numice(1,k)
      p3_main_inputs(1,k,12) = rain(1,k)
      p3_main_inputs(1,k,13) = numrain(1,k)
      p3_main_inputs(1,k,14) = qm(1,k)
      p3_main_inputs(1,k,15) = rimvol(1,k)
      p3_main_inputs(1,k,16) = state%pdeldry(1,k)
      p3_main_inputs(1,k,17) = relvar(1,k)
    end do
    p3_main_inputs(1,pver+1,5) = state%zi(1,pver+1)

    if (do_prescribed_CCN) then
       call pbuf_get_field(pbuf, ccn3_idx, ccn_trcdat) ! now you can use ccn_trcdat anywhere in this code
       ! ccn is uniformly distributed throughout the cell, but P3 computes in-cloud values assuming cell-averages are comprised
       ! of zero values outside cloud. Preemptively multiplying by cldfrac here is needed to get the correct in-cloud ccn value in P3
       nccn_prescribed = ccn_trcdat
    end if

    ! CALL P3
    !==============
    ! TODO: get proper value for 'it' from time module
    dummy_out(:,:) = 0.0_rtype
    precip_liq_surf = 0.0_rtype
    precip_ice_surf = 0.0_rtype
    prec_pcw = 0.0_rtype
    snow_pcw = 0.0_rtype
    vap_liq_exchange = 0.0_rtype

    call t_startf('micro_p3_tend_loop')
    call p3_main( &
         cldliq(its:ite,kts:kte),     & ! INOUT  cloud, mass mixing ratio         kg kg-1
         numliq(its:ite,kts:kte),     & ! INOUT  cloud, number mixing ratio       #  kg-1
         rain(its:ite,kts:kte),       & ! INOUT  rain, mass mixing ratio          kg kg-1
         numrain(its:ite,kts:kte),    & ! INOUT  rain, number mixing ratio        #  kg-1
         th(its:ite,kts:kte),         & ! INOUT  potential temperature            K
         qv_dry(its:ite,kts:kte),     & ! INOUT  water vapor mixing ratio         kg kg-1
         dtime,                       & ! IN     model time step                  s
         ice(its:ite,kts:kte),        & ! INOUT  ice, total mass mixing ratio     kg kg-1
         qm(its:ite,kts:kte),         & ! INOUT  ice, rime mass mixing ratio      kg kg-1
         numice(its:ite,kts:kte),     & ! INOUT  ice, total number mixing ratio   #  kg-1
         rimvol(its:ite,kts:kte),     & ! INOUT  ice, rime volume mixing ratio    m3 kg-1
         pres(its:ite,kts:kte),       & ! IN     pressure at cell midpoints       Pa
         dz(its:ite,kts:kte),         & ! IN     vertical grid spacing            m
         npccn(its:ite,kts:kte),      & ! IN ccn activation number tendency kg-1 s-1
         nccn_prescribed(its:ite,kts:kte), & ! IN ccn prescribed concentration
         ni_activated(its:ite,kts:kte),    & ! IN activated ice nuclei concentration kg-1
         relvar(its:ite,kts:kte),     & ! IN cloud liquid relative variance
         it,                          & ! IN     time step counter NOTE: starts at 1 for first time step
         precip_liq_surf(its:ite),            & ! OUT    surface liquid precip rate       m s-1
         precip_ice_surf(its:ite),            & ! OUT    surface frozen precip rate       m s-1
         its,                         & ! IN     horizontal index lower bound     -
         ite,                         & ! IN     horizontal index upper bound     -
         kts,                         & ! IN     vertical index lower bound       -
         kte,                         & ! IN     vertical index upper bound       -
         rel(its:ite,kts:kte),        & ! OUT    effective radius, cloud          m
         rei(its:ite,kts:kte),        & ! OUT    effective radius, ice            m
         rho_qi(its:ite,kts:kte),  & ! OUT    bulk density of ice              kg m-3
         do_predict_nc,               & ! IN     .true.=prognostic Nc, .false.=specified Nc
         do_prescribed_CCN,           & ! IN
         ! AaronDonahue new stuff
         state%pdeldry(its:ite,kts:kte),  & ! IN pressure level thickness for computing total mass
         inv_exner(its:ite,kts:kte),      & ! IN exner values
         qv2qi_depos_tend(its:ite,kts:kte),    & ! OUT Deposition/sublimation rate of cloud ice
         precip_total_tend(its:ite,kts:kte),      & ! OUT Total precipitation (rain + snow)
         nevapr(its:ite,kts:kte),     & ! OUT evaporation of total precipitation (rain + snow)
         qr_evap_tend(its:ite,kts:kte),  & ! OUT rain evaporation
         precip_liq_flux(its:ite,kts:kte+1),     & ! OUT grid-box average rain flux (kg m^-2s^-1) pverp
         precip_ice_flux(its:ite,kts:kte+1),     & ! OUT grid-box average ice/snow flux (kgm^-2 s^-1) pverp
         cld_frac_r(its:ite,kts:kte),      & ! IN rain cloud fraction
         cld_frac_l(its:ite,kts:kte),      & ! IN liquid cloud fraction
         cld_frac_i(its:ite,kts:kte),      & ! IN ice cloud fraction
         tend_out(its:ite,kts:kte,:), & ! OUT p3 microphysics tendencies
         mu(its:ite,kts:kte),         & ! OUT Size distribution shape parameter for radiation
         lambdac(its:ite,kts:kte),    & ! OUT Size distribution slope parameter for radiation
         liq_ice_exchange(its:ite,kts:kte),& ! OUT sum of liq-ice phase change tendenices
         vap_liq_exchange(its:ite,kts:kte),& ! OUT sun of vap-liq phase change tendencies
         vap_ice_exchange(its:ite,kts:kte),& ! OUT sum of vap-ice phase change tendencies
         qv_prev_dry(its:ite,kts:kte),         & ! IN  qv at end of prev p3_main call   kg kg-1
         t_prev(its:ite,kts:kte),          & ! IN  t at end of prev p3_main call    K
         col_location(its:ite,:3)          & ! IN column locations
         )

    p3_main_outputs(:,:,:) = -999._rtype
    do k = 1,pver
      p3_main_outputs(1,k, 1) = cldliq(1,k)
      p3_main_outputs(1,k, 2) = numliq(1,k)
      p3_main_outputs(1,k, 3) = rain(1,k)
      p3_main_outputs(1,k, 4) = numrain(1,k)
      p3_main_outputs(1,k, 5) = th(1,k)
      p3_main_outputs(1,k, 6) = qv_dry(1,k)
      p3_main_outputs(1,k, 7) = ice(1,k)
      p3_main_outputs(1,k, 8) = qm(1,k)
      p3_main_outputs(1,k, 9) = numice(1,k)
      p3_main_outputs(1,k,10) = rimvol(1,k)
      p3_main_outputs(1,k,14) = rel(1,k)
      p3_main_outputs(1,k,15) = rei(1,k)
      p3_main_outputs(1,k,18) = rho_qi(1,k)
      p3_main_outputs(1,k,19) = qv2qi_depos_tend(1,k)
      p3_main_outputs(1,k,20) = precip_total_tend(1,k)
      p3_main_outputs(1,k,21) = nevapr(1,k)
      p3_main_outputs(1,k,22) = qr_evap_tend(1,k)
      p3_main_outputs(1,k,23) = precip_liq_flux(1,k)
      p3_main_outputs(1,k,24) = precip_ice_flux(1,k)
      p3_main_outputs(1,k,27) = mu(1,k)
      p3_main_outputs(1,k,28) = lambdac(1,k)
      p3_main_outputs(1,k,29) = liq_ice_exchange(1,k)
      p3_main_outputs(1,k,30) = vap_liq_exchange(1,k)
      p3_main_outputs(1,k,31) = vap_ice_exchange(1,k)
    end do
    p3_main_outputs(1,1,11) = precip_liq_surf(1)
    p3_main_outputs(1,1,12) = precip_ice_surf(1)
    p3_main_outputs(1,pver+1,23) = precip_liq_flux(1,pver+1)
    p3_main_outputs(1,pver+1,24) = precip_ice_flux(1,pver+1)
    call outfld('P3_input',  p3_main_inputs,  pcols, lchnk)
    call outfld('P3_output', p3_main_outputs, pcols, lchnk)

    !MASSAGE OUTPUT TO FIT E3SM EXPECTATIONS
    !=============

    !TODO: figure out what else other E3SM parameterizations need from micro and make sure
    !they are assigned here. The comments below are a step in that direction.


    !cloud_rad_props also uses snow radiative properties which aren't available from
    !P3 (perhaps because ice phase in p3 includes *all* ice already?).

    !BACK OUT TENDENCIES FROM STATE CHANGES
    !=============

    !DRY-TO-WET MMRs:
    !================
    !Since the host model needs wet mixing ratio tendencies(state vector has wet mixing ratios),
    !we need to convert dry mixing ratios from P3 to wet mixing ratios before extracting tendencies
    !NOTE: water vapor mixing ratio argument in calculate_wetmmr_from_drymmr function has to be dry water vapor mixing ratio

    qv_wet_out(:ncol,:pver) = calculate_wetmmr_from_drymmr(ncol, pver, qv_dry,  qv_dry)
    cldliq(:ncol,:pver)     = calculate_wetmmr_from_drymmr(ncol, pver, cldliq,  qv_dry)
    numliq(:ncol,:pver)     = calculate_wetmmr_from_drymmr(ncol, pver, numliq,  qv_dry)
    rain(:ncol,:pver)       = calculate_wetmmr_from_drymmr(ncol, pver, rain,    qv_dry)
    numrain(:ncol,:pver)    = calculate_wetmmr_from_drymmr(ncol, pver, numrain, qv_dry)
    ice(:ncol,:pver)        = calculate_wetmmr_from_drymmr(ncol, pver, ice,     qv_dry)
    numice(:ncol,:pver)     = calculate_wetmmr_from_drymmr(ncol, pver, numice,  qv_dry)
    qm(:ncol,:pver)         = calculate_wetmmr_from_drymmr(ncol, pver, qm,      qv_dry)
    rimvol(:ncol,:pver)     = calculate_wetmmr_from_drymmr(ncol, pver, rimvol,  qv_dry)

    temp(:ncol,:pver) = th(:ncol,:pver)/inv_exner(:ncol,:pver)
    ptend%s(:ncol,:pver)           = cpair*( temp(:ncol,:pver) - state%t(:ncol,:pver) )/dtime
    ptend%q(:ncol,:pver,1)         = ( max(0._rtype,qv_wet_out(:ncol,:pver)     ) - state%q(:ncol,:pver,1)         )/dtime
    ptend%q(:ncol,:pver,ixcldliq)  = ( max(0._rtype,cldliq(:ncol,:pver) ) - state%q(:ncol,:pver,ixcldliq)  )/dtime
    ptend%q(:ncol,:pver,ixnumliq)  = ( max(0._rtype,numliq(:ncol,:pver) ) - state%q(:ncol,:pver,ixnumliq)  )/dtime
    ptend%q(:ncol,:pver,ixrain)    = ( max(0._rtype,rain(:ncol,:pver)   ) - state%q(:ncol,:pver,ixrain)    )/dtime
    ptend%q(:ncol,:pver,ixnumrain) = ( max(0._rtype,numrain(:ncol,:pver)) - state%q(:ncol,:pver,ixnumrain) )/dtime
    ptend%q(:ncol,:pver,ixcldice)  = ( max(0._rtype,ice(:ncol,:pver)    ) - state%q(:ncol,:pver,ixcldice)  )/dtime
    ptend%q(:ncol,:pver,ixnumice)  = ( max(0._rtype,numice(:ncol,:pver) ) - state%q(:ncol,:pver,ixnumice)  )/dtime
    ptend%q(:ncol,:pver,ixcldrim)  = ( max(0._rtype,qm(:ncol,:pver)  ) - state%q(:ncol,:pver,ixcldrim)  )/dtime
    ptend%q(:ncol,:pver,ixrimvol)  = ( max(0._rtype,rimvol(:ncol,:pver) ) - state%q(:ncol,:pver,ixrimvol)  )/dtime

    ! Update t_prev and qv_prev to be used by evap_precip
    t_prev(:ncol,:pver) = temp(:ncol,:pver)
    qv_prev_wet(:ncol,:pver) = qv_wet_out(:ncol,:pver)

    call t_stopf('micro_p3_tend_loop')
    call t_startf('micro_p3_tend_finish')
   ! Following MG interface as a template:

    ! Net micro_p3 condensation rate
    ! NOTE: These are probably in dry coordinate, not wet.
    qme(:ncol,top_lev:pver) = cmeliq(:ncol,top_lev:pver) + qv2qi_depos_tend(:ncol,top_lev:pver)  ! qv2qi_depos_tend is output from p3 micro
    ! Add cmeliq to  vap_liq_exchange
    vap_liq_exchange(:ncol,top_lev:pver) = vap_liq_exchange(:ncol,top_lev:pver) + cmeliq(:ncol,top_lev:pver)

!====================== Export variables/Conservation START ======================!
     !For precip, accumulate only total precip in prec_pcw and snow_pcw variables.
    ! Other precip output variables are set to 0
    ! Do not subscript by ncol here, because in physpkg we divide the whole
    ! array and need to avoid an FPE due to uninitialized data.
    prec_pcw = precip_liq_surf + precip_ice_surf
    prec_sed = 0._rtype
    prec_str = prec_pcw + prec_sed

    snow_pcw = precip_ice_surf
    snow_sed = 0._rtype
    snow_str = snow_pcw + snow_sed
!====================== Export variables/Conservation END ======================!

!====================== Radiation Specific Outputs START ======================!

   ! Calculate rho for size distribution
   ! parameter calculations and average it if needed

   rho(:ncol,top_lev:) = &
      state%pmid(:ncol,top_lev:) / (rair*temp(:ncol,top_lev:))
   ! ------------------------------------------------------------ !
   ! Compute in cloud ice and liquid mixing ratios                !
   ! Note that 'iclwp, iciwp' are used for radiation computation. !
   ! ------------------------------------------------------------ !

   icinc = 0._rtype
   icwnc = 0._rtype

   do k = top_lev, pver
      do icol = 1, ncol
         ! Limits for in-cloud mixing ratios consistent with P3 microphysics
         ! in-cloud mixing ratio maximum limit of 0.005 kg/kg
         icimrst(icol,k)   = min( state%q(icol,k,ixcldice) / max(mincld,cld_frac_i(icol,k)),0.005_rtype )
         icwmrst(icol,k)   = min( state%q(icol,k,ixcldliq) / max(mincld,cld_frac_l(icol,k)),0.005_rtype )
         icinc(icol,k)     = state%q(icol,k,ixnumice) / max(mincld,cld_frac_i(icol,k)) * &
              state%pmid(icol,k) / (287.15_rtype*state%t(icol,k))
         icwnc(icol,k)     = state%q(icol,k,ixnumliq) / max(mincld,cld_frac_l(icol,k)) * &
              state%pmid(icol,k) / (287.15_rtype*state%t(icol,k))
      end do
   end do

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !!
   !! derived fields
   !!
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !cloud_rad_props needs ice effective diameter, which Kai calculates as below:
    !   dei = rei*diag_rhopo(i,k,iice)/rho_h2os*2._rtype
    !where rhopo is bulk ice density from table lookup (taken from table_val_ice_bulk_dens, here written as rho_qi) and rho_h2os=917.0 is a constant parameter.
   !! Effective radius for cloud liquid
   rel(:ncol,top_lev:) = rel(:ncol,top_lev:) * 1e6_rtype  ! Rescale rel to be in microns
   !! Effective radius for cloud ice
   rei(:ncol,top_lev:) = rei(:ncol,top_lev:) * 1e6_rtype  ! Rescale rei to be in microns
   !! Effective diameter for cloud ice
   dei(:ncol,top_lev:) = rei(:ncol,top_lev:) * 2._rtype

   !!
   !! Limiters for low cloud fraction
   !!

   do k = top_lev, pver
      do icol = 1, ncol
         if ( ast(icol,k) < 1.e-4_rtype ) then
            mu(icol,k) = mucon
            lambdac(icol,k) = (mucon + 1._rtype)/dcon
            dei(icol,k) = deicon
         end if
      end do
   end do

   !!
   !! New output fields
   !!
   efcout      = 0._rtype
   efiout      = 0._rtype
   ncout       = 0._rtype
   niout       = 0._rtype
   freql       = 0._rtype
   freqi       = 0._rtype
   cdnumc      = 0._rtype
   nfice       = 0._rtype

   ! FICE
   do k = top_lev, pver
      do icol = 1, ncol
      if (ice(icol,k).gt.qsmall .and. (rain(icol,k)+ice(icol,k)+cldliq(icol,k)).gt.qsmall) then
         nfice(icol,k)=min(ice(icol,k)/(rain(icol,k)+ice(icol,k)+cldliq(icol,k)),1._rtype)
      else
         nfice(icol,k)=0._rtype
      end if
      end do
   end do

   ! Column droplet concentration
   cdnumc(:ncol) = sum(numliq(:ncol,top_lev:pver) * &
        state%pdel(:ncol,top_lev:pver)/gravit, dim=2)
   do k = top_lev, pver
      do icol = 1, ncol
         if ( cld_frac_l(icol,k) > 0.01_rtype .and. icwmrst(icol,k) > 5.e-5_rtype ) then
            efcout(icol,k) = rel(icol,k) * cld_frac_l(icol,k)
            ncout(icol,k)  = icwnc(icol,k) * cld_frac_l(icol,k)
            freql(icol,k)  = cld_frac_l(icol,k)
         end if
         if ( cld_frac_i(icol,k) > 0.01_rtype .and. icimrst(icol,k) > 1.e-6_rtype ) then
            efiout(icol,k) = rei(icol,k) * cld_frac_i(icol,k)
            niout(icol,k)  = icinc(icol,k) * cld_frac_i(icol,k)
            freqi(icol,k)  = cld_frac_i(icol,k)
         end if
      end do
   end do

   ! note: 1e-6 kgho2/kgair/s * 1000. pa / (9.81 m/s2) / 1000 kgh2o/m3 = 1e-7 m/s
   ! this is 1ppmv of h2o in 10hpa
   ! alternatively: 0.1 mm/day * 1.e-4 m/mm * 1/86400 day/s = 1.e-9

   !!
   !! Rain/Snow effective diameter
   !!
   drout2    = 0._rtype
   reff_rain = 0._rtype
   aqrain    = 0._rtype
   anrain    = 0._rtype
   freqr     = 0._rtype
   ! Prognostic precipitation
   where (rain(:ncol,top_lev:) >= 1.e-7_rtype)
      drout2(:ncol,top_lev:) = avg_diameter( &
           rain(:ncol,top_lev:), &
           numrain(:ncol,top_lev:) * rho(:ncol,top_lev:), &
           rho(:ncol,top_lev:), rho_h2o)

      aqrain(:ncol,top_lev:) = rain(:ncol,top_lev:) * cld_frac_r(:ncol,top_lev:)
      anrain(:ncol,top_lev:) = numrain(:ncol,top_lev:) * cld_frac_r(:ncol,top_lev:)
      freqr(:ncol,top_lev:) = cld_frac_r(:ncol,top_lev:)
      reff_rain(:ncol,top_lev:) = drout2(:ncol,top_lev:) * &
           1.5_rtype * 1.e6_rtype
   end where

!====================== COSP Specific Outputs START ======================!
! LS_FLXPRC, LS_FLXSNW, LS_REFFRAIN, LS_REFFSNOW, CV_REFFLIQ, CV_REFFICE
   !== Grid-box mean flux_large_scale_cloud at interfaces (kg/m2/s)
! flxprc and flxsnw are used in COSP to compute precipitation fractional
! area and derive precipitation (rain and snow) mixing ratios. Including iflx
! and cflx in precipitation fluxes would result in additional effects of cloud liquid and
! ice on cosp's smiluated lidar and radar reflectivity signal through the rain/snow
! portion of calculations that are handled separately from that of cloud liquid
! and ice. If included, it would not exactly amount to double counting the effect of
! cloud liquid and ice because the mixing ratio derived from iflx and cflx epected to be much smaller
! than the actual grid-mean cldliq and cldice, and rain or snow size distribution
! would be used to compute the lidar/radar signal strength.
!
! Note that it would need to include iflx and cflx to make the values at surface
! interface consistent with large scale precipitation rates.

    ! array must be zeroed beyond trop_cloud_top_pre otherwise undefined values will be used in cosp.
    flxprc(:ncol,1:top_lev) = 0.0_rtype ! Rain+Snow
    flxsnw(:ncol,1:top_lev) = 0.0_rtype ! Snow

    flxprc(:ncol,top_lev:pverp) = precip_liq_flux(:ncol,top_lev:pverp) + precip_ice_flux(:ncol,top_lev:pverp) ! need output from p3
    flxsnw(:ncol,top_lev:pverp) = precip_ice_flux(:ncol,top_lev:pverp) ! need output from p3

    cvreffliq(:ncol,top_lev:pver) = 9.0_rtype
    cvreffice(:ncol,top_lev:pver) = 37.0_rtype

    reffrain(:,:) = 0._rtype
    reffsnow(:,:) = 0._rtype
    reffrain(:ncol,top_lev:pver) = reff_rain(:ncol,top_lev:pver)
    reffsnow(:ncol,top_lev:pver) = 1000._rtype !! dummy value, the choice here impacts the COSP output variable: CFAD_DBZE94_CS.  TODO: Figure out if this is ok, change if needed.

!====================== COSP Specific Outputs  END ======================!

    !WRITE OUTPUT
    !=============
   call outfld('AQRAIN',      aqrain,      pcols,    lchnk)
   call outfld('ANRAIN',      anrain,      pcols,    lchnk)
   call outfld('AREL',        efcout,      pcols,    lchnk)
   call outfld('AREI',        efiout,      pcols,    lchnk)
   call outfld('AWNC' ,       ncout,       pcols,    lchnk)
   call outfld('AWNI' ,       niout,       pcols,    lchnk)
   call outfld('FICE',        nfice,       pcols,    lchnk)
   call outfld('FREQL',       freql,       pcols,    lchnk)
   call outfld('FREQI',       freqi,       pcols,    lchnk)
   call outfld('FREQR',       freqr,       pcols,    lchnk)
   call outfld('CDNUMC',      cdnumc,      pcols,    lchnk)

   call outfld('CLOUDFRAC_LIQ_MICRO',  cld_frac_l,      pcols, lchnk)
   call outfld('CLOUDFRAC_ICE_MICRO',  cld_frac_i,      pcols, lchnk)
   call outfld('CLOUDFRAC_RAIN_MICRO', cld_frac_r,      pcols, lchnk)

   ! Write p3 tendencies as output
   ! warm-phase process rates
   call outfld('P3_qrcon',  tend_out(:,:, 1), pcols, lchnk)
   call outfld('P3_qc2qr_accret_tend',  tend_out(:,:, 2), pcols, lchnk)
   call outfld('P3_qc2qr_autoconv_tend',  tend_out(:,:, 3), pcols, lchnk)
   call outfld('P3_nc_accret_tend',  tend_out(:,:, 4), pcols, lchnk)
   call outfld('P3_nc2nr_autoconv_tend', tend_out(:,:, 5), pcols, lchnk)
   call outfld('P3_nc_selfcollect_tend',  tend_out(:,:, 6), pcols, lchnk)
   call outfld('P3_nr_selfcollect_tend',  tend_out(:,:, 7), pcols, lchnk)
   call outfld('P3_nc_nuceat_tend',  tend_out(:,:, 8), pcols, lchnk)
   call outfld('P3_qccon',  tend_out(:,:, 9), pcols, lchnk)
   call outfld('P3_qcnuc',  tend_out(:,:,10), pcols, lchnk)
   call outfld('P3_qr2qv_evap_tend',  tend_out(:,:,11), pcols, lchnk)
   call outfld('P3_qcevp',  tend_out(:,:,12), pcols, lchnk)
   call outfld('P3_nr_evap_tend',  tend_out(:,:,13), pcols, lchnk)
   call outfld('P3_ncautr', tend_out(:,:,14), pcols, lchnk)
   ! ice-phase process rate
   call outfld('P3_qccol',  tend_out(:,:,15), pcols, lchnk)
   call outfld('P3_qwgrth', tend_out(:,:,16), pcols, lchnk) ! Not a tendency in itself, it is used to build qccol and qrcol
   call outfld('P3_qidep',  tend_out(:,:,17), pcols, lchnk)
   call outfld('P3_qrcol',  tend_out(:,:,18), pcols, lchnk)
   call outfld('P3_qinuc',  tend_out(:,:,19), pcols, lchnk)
   call outfld('P3_nc_collect_tend',  tend_out(:,:,20), pcols, lchnk)
   call outfld('P3_nr_collect_tend',  tend_out(:,:,21), pcols, lchnk)
   call outfld('P3_ni_nucleat_tend',  tend_out(:,:,22), pcols, lchnk)
   call outfld('P3_qi2qv_sublim_tend',  tend_out(:,:,23), pcols, lchnk)
   call outfld('P3_qi2qr_melt_tend',  tend_out(:,:,24), pcols, lchnk)
   call outfld('P3_ni2nr_melt_tend',  tend_out(:,:,25), pcols, lchnk)
   call outfld('P3_ni_sublim_tend',  tend_out(:,:,26), pcols, lchnk)
   call outfld('P3_ni_selfcollect_tend',  tend_out(:,:,27), pcols, lchnk)
   call outfld('P3_qc2qi_hetero_frz_tend', tend_out(:,:,28), pcols, lchnk)
   call outfld('P3_qr2qi_immers_frz_tend', tend_out(:,:,29), pcols, lchnk)
   call outfld('P3_nc2ni_immers_frz_tend', tend_out(:,:,30), pcols, lchnk)
   call outfld('P3_nr2ni_immers_frz_tend', tend_out(:,:,31), pcols, lchnk)
   call outfld('P3_nr_ice_shed_tend', tend_out(:,:,32), pcols, lchnk)
   call outfld('P3_qc2qr_ice_shed_tend',  tend_out(:,:,33), pcols, lchnk)
!   call outfld('P3_qcmul',  tend_out(:,:,34), pcols, lchnk) ! Not actually used, so not actually recorded.  Commented out here for continuity of the array.
   call outfld('P3_ncshdc', tend_out(:,:,35), pcols, lchnk)
   ! sedimentation
   call outfld('P3_sed_CLDLIQ',  tend_out(:,:,36), pcols, lchnk)
   call outfld('P3_sed_NUMLIQ',  tend_out(:,:,37), pcols, lchnk)
   call outfld('P3_sed_CLDRAIN', tend_out(:,:,38), pcols, lchnk)
   call outfld('P3_sed_NUMRAIN', tend_out(:,:,39), pcols, lchnk)
   call outfld('P3_sed_CLDICE',  tend_out(:,:,40), pcols, lchnk)
   call outfld('P3_sed_NUMICE',  tend_out(:,:,41), pcols, lchnk)
   ! microphysics processes
   call outfld('P3_mtend_CLDLIQ',  tend_out(:,:,42), pcols, lchnk)
   call outfld('P3_mtend_NUMLIQ',  tend_out(:,:,43), pcols, lchnk)
   call outfld('P3_mtend_CLDRAIN', tend_out(:,:,44), pcols, lchnk)
   call outfld('P3_mtend_NUMRAIN', tend_out(:,:,45), pcols, lchnk)
   call outfld('P3_mtend_CLDICE',  tend_out(:,:,46), pcols, lchnk)
   call outfld('P3_mtend_NUMICE',  tend_out(:,:,47), pcols, lchnk)
   call outfld('P3_mtend_Q',       tend_out(:,:,48), pcols, lchnk)
   call outfld('P3_mtend_TH',      tend_out(:,:,49), pcols, lchnk)
   ! Phase change tendencies
   call outfld('vap_ice_exchange',      vap_ice_exchange,      pcols, lchnk)
   call outfld('vap_liq_exchange',      vap_liq_exchange,      pcols, lchnk)
   call outfld('liq_ice_exchange',      liq_ice_exchange,      pcols, lchnk)

   call t_stopf('micro_p3_tend_finish')
  end subroutine micro_p3_tend

  !================================================================================================

end module micro_p3_interface
