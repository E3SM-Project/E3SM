module CNRestMod
  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Read/Write to/from CN info to CLM restart file. 
  !
  ! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils  , only : endrun
  use shr_log_mod , only : errMsg => shr_log_errMsg
  use spmdMod     , only : masterproc
  use clm_varctl  , only : iulog, override_bgc_restart_mismatch_dump
  use decompMod   , only : bounds_type
  use restUtilMod
  use ncdio_pio
  use pio
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: CNRest
  !----------------------------------------------------------------------- 

contains

  !-----------------------------------------------------------------------
  subroutine CNRest ( bounds, ncid, flag )
    !
    ! !DESCRIPTION: 
    ! Read/write CN restart data
    !
    ! !USES:
    use shr_infnan_mod,   only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use shr_const_mod,    only : SHR_CONST_PDB
    use clm_varpar,       only : ndecomp_pools, nlevdecomp, crop_prog
    use clm_time_manager, only : is_restart, get_nstep
    use clm_varcon,       only : nlevgrnd, dzsoi_decomp
    use clm_varctl,       only : use_c13, use_c14, use_vertsoilc, use_nitrif_denitrif
    use clm_varctl,       only : use_century_decomp, use_cndv, spinup_state
    use clm_varcon,       only : c13ratio, c14ratio, spval
    use decompMod,        only : bounds_type, get_proc_global
    use clmtype
    !
    ! !ARGUMENTS:
    implicit none
    type(bounds_type), intent(in) :: bounds ! bounds
    type(file_desc_t)             :: ncid   ! netcdf id
    character(len=*), intent(in)  :: flag   !'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    real(r8)           :: c3_del13c            ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(r8)           :: c4_del13c            ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(r8)           :: c3_r1                ! isotope ratio (13c/12c) for C3 photosynthesis
    real(r8)           :: c4_r1                ! isotope ratio (13c/12c) for C4 photosynthesis
    real(r8)           :: c3_r2                ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(r8)           :: c4_r2                ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    integer            :: c,p,j,k,i,l          ! indices 
    integer            :: numg_global          ! total number of grid cells, globally
    integer            :: numl_global          ! total number of landunits, globally
    integer            :: numc_global          ! total number of columns, globally
    integer            :: nump_global          ! total number of pfts, globally
    real(r8)           :: m                    ! multiplier for the exit_spinup code
    logical            :: readvar              ! determine if variable is on initial file
    logical            :: do_io                ! whether to do i/o for the given variable
    integer            :: dimlen               ! dimension length
    integer            :: err_code             ! error code
    character(len=128) :: varname              ! temporary
    integer            :: itemp                ! temporary 
    integer , pointer  :: iptemp(:)            ! pointer to memory to be allocated
    real(r8), pointer  :: ptr1d(:), ptr2d(:,:) ! temporary arrays for slicing larger arrays
    integer            :: ier                  ! error status
    integer            :: nstep                ! time step number
    integer            :: idata
    type(var_desc_t)   :: vardesc              ! local vardesc
    logical            :: exit_spinup = .false.
    logical            :: enter_spinup = .false.
    type(pft_cstate_type), pointer :: pcisos
    type(pft_cstate_type), pointer :: pcbulks
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
    integer            :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer           :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !-----------------------------------------------------------------------

    if ( use_c13 ) then
       pcisos => pc13s
       pcbulks => pcs
    endif

    ! Get expected total number of points, for later error checks
    call get_proc_global(numg_global, numl_global, numc_global, nump_global)

    ! Determine necessary subgrid bounds

    if ( use_c13 ) then
       c3_del13c = -28._r8
       c4_del13c = -13._r8
       c3_r1 = SHR_CONST_PDB + ((c3_del13c*SHR_CONST_PDB)/1000._r8)
       c3_r2 = c3_r1/(1._r8 + c3_r1)
       c4_r1 = SHR_CONST_PDB + ((c4_del13c*SHR_CONST_PDB)/1000._r8)
       c4_r2 = c4_r1/(1._r8 + c4_r1)
    endif

    !--------------------------------
    ! pft ecophysiological variables 
    !--------------------------------

    ! dormant_flag
    call restartvar(ncid=ncid, flag=flag, varname='dormant_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='dormancy flag', units='unitless', &
         interpinic_flag='interp', readvar=readvar, data=pepv%dormant_flag) 

    ! days_active
    call restartvar(ncid=ncid, flag=flag, varname='days_active', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='number of days since last dormancy', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%days_active) 

    ! onset_flag
    call restartvar(ncid=ncid, flag=flag, varname='onset_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='flag if critical growing degree-day sum is exceeded', units='unitless' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%onset_flag) 

    ! onset_counter
    call restartvar(ncid=ncid, flag=flag, varname='onset_counter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset days counter', units='sec' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%onset_counter) 

    ! onset_gddflag
    call restartvar(ncid=ncid, flag=flag, varname='onset_gddflag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset flag for growing degree day sum', units='' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%onset_gddflag) 

    ! onset_fdd
    call restartvar(ncid=ncid, flag=flag, varname='onset_fdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset freezing degree days counter', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%onset_fdd) 

    ! onset_gdd
    call restartvar(ncid=ncid, flag=flag, varname='onset_gdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset growing degree days', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%onset_gdd) 

    ! onset_swi
    call restartvar(ncid=ncid, flag=flag, varname='onset_swi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='onset soil water index', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%onset_swi) 

    ! offset_flag
    call restartvar(ncid=ncid, flag=flag, varname='offset_flag', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset flag', units='unitless' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%offset_flag) 

    ! offset_counter
    call restartvar(ncid=ncid, flag=flag, varname='offset_counter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset days counter', units='sec' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%offset_counter) 

    ! offset_fdd
    call restartvar(ncid=ncid, flag=flag, varname='offset_fdd', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='offset freezing degree days counter', units='days' , &
         interpinic_flag='interp', readvar=readvar, data=pepv%offset_fdd) 

    ! offset_swi
    call restartvar(ncid=ncid, flag=flag, varname='offset_swi', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%offset_swi) 

    if (crop_prog) then
       ! fert_counter
       call restartvar(ncid=ncid, flag=flag, varname='fert_counter', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%fert_counter)

       ! fert
       call restartvar(ncid=ncid, flag=flag, varname='fert', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pnf%fert)
    end if

    ! lgsf
    call restartvar(ncid=ncid, flag=flag, varname='lgsf', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%lgsf) 

    ! bglfr
    call restartvar(ncid=ncid, flag=flag, varname='bglfr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%bglfr) 

    ! bgtr
    call restartvar(ncid=ncid, flag=flag, varname='bgtr', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%bgtr) 

    ! annavg_t2m
    call restartvar(ncid=ncid, flag=flag, varname='annavg_t2m', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%annavg_t2m) 

    ! tempavg_t2m
    call restartvar(ncid=ncid, flag=flag, varname='tempavg_t2m', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%tempavg_t2m) 

    ! gpp
    call restartvar(ncid=ncid, flag=flag, varname='gpp_pepv', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%gpp) 

    ! availc
    call restartvar(ncid=ncid, flag=flag, varname='availc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%availc) 

    ! xsmrpool_recover
    call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_recover', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%xsmrpool_recover) 

    if ( use_c13 ) then
       ! xsmrpool_c13ratio
       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_c13ratio', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%xsmrpool_c13ratio) 
    endif

    ! alloc_pnow
    call restartvar(ncid=ncid, flag=flag, varname='alloc_pnow', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%alloc_pnow) 

    ! c_allometry
    call restartvar(ncid=ncid, flag=flag, varname='c_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%c_allometry) 

    ! n_allometry
    call restartvar(ncid=ncid, flag=flag, varname='n_allometry', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%n_allometry) 

    ! plant_ndemand
    call restartvar(ncid=ncid, flag=flag, varname='plant_ndemand', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%plant_ndemand) 

    ! tempsum_potential_gpp
    call restartvar(ncid=ncid, flag=flag, varname='tempsum_potential_gpp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%tempsum_potential_gpp) 

    !annsum_potential_gpp 
    call restartvar(ncid=ncid, flag=flag, varname='annsum_potential_gpp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%annsum_potential_gpp) 

    ! tempmax_retransn
    call restartvar(ncid=ncid, flag=flag, varname='tempmax_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%tempmax_retransn) 

    ! annmax_retransn
    call restartvar(ncid=ncid, flag=flag, varname='annmax_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%annmax_retransn) 

    ! avail_retransn
    call restartvar(ncid=ncid, flag=flag, varname='avail_retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%avail_retransn) 

    ! plant_nalloc
    call restartvar(ncid=ncid, flag=flag, varname='plant_nalloc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%plant_nalloc) 

    ! plant_calloc
    call restartvar(ncid=ncid, flag=flag, varname='plant_calloc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%plant_calloc) 

    ! excess_cflux
    call restartvar(ncid=ncid, flag=flag, varname='excess_cflux', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%excess_cflux) 

    ! downreg
    call restartvar(ncid=ncid, flag=flag, varname='downreg', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%downreg) 

    ! prev_leafc_to_litter
    call restartvar(ncid=ncid, flag=flag, varname='prev_leafc_to_litter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%prev_leafc_to_litter) 

    ! prev_frootc_to_litter
    call restartvar(ncid=ncid, flag=flag, varname='prev_frootc_to_litter', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%prev_frootc_to_litter) 

    ! tempsum_npp
    call restartvar(ncid=ncid, flag=flag, varname='tempsum_npp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%tempsum_npp) 
 
    ! annsum_npp
    call restartvar(ncid=ncid, flag=flag, varname='annsum_npp', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pepv%annsum_npp) 

    if ( use_c13 ) then
       ! rc13_canair
       call restartvar(ncid=ncid, flag=flag, varname='rc13_canair', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%rc13_canair) 

       ! rc13_psnsun
       call restartvar(ncid=ncid, flag=flag, varname='rc13_psnsun', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%rc13_psnsun) 

       ! rc13_psnsha
       call restartvar(ncid=ncid, flag=flag, varname='rc13_psnsha', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%rc13_psnsha) 
    endif

    if (crop_prog) then
       ! grain_flag
       call restartvar(ncid=ncid, flag=flag, varname='grain_flag', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%grain_flag)
    end if

    !--------------------------------
    ! pft carbon state variables 
    !--------------------------------

    ! leafc
    call restartvar(ncid=ncid, flag=flag, varname='leafc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%leafc) 

    ! leafc_storage
    call restartvar(ncid=ncid, flag=flag, varname='leafc_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%leafc_storage) 

    ! leafc_xfer
    call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%leafc_xfer) 

    ! frootc
    call restartvar(ncid=ncid, flag=flag, varname='frootc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%frootc) 

    ! frootc_storage
    call restartvar(ncid=ncid, flag=flag, varname='frootc_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%frootc_storage) 

    !frootc_xfer 
    call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%frootc_xfer) 

    ! livestemc
    call restartvar(ncid=ncid, flag=flag, varname='livestemc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%livestemc) 

    ! livestemc_storage
    call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%livestemc_storage) 

    ! livestemc_xfer
    call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%livestemc_xfer) 

    ! deadstemc
    call restartvar(ncid=ncid, flag=flag, varname='deadstemc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%deadstemc) 

    ! deadstemc_storage
    call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%deadstemc_storage) 

    ! deadstemc_xfer
    call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%deadstemc_xfer) 

    ! livecrootc
    call restartvar(ncid=ncid, flag=flag, varname='livecrootc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%livecrootc) 

    ! livecrootc_storage
    call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%livecrootc_storage) 

    ! livecrootc_xfer
    call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%livecrootc_xfer) 

    ! deadcrootc
    call restartvar(ncid=ncid, flag=flag, varname='deadcrootc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%deadcrootc) 

    ! deadcrootc_storage
    call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%deadcrootc_storage) 

    ! deadcrootc_xfer
    call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%deadcrootc_xfer) 

    ! gresp_storage
    call restartvar(ncid=ncid, flag=flag, varname='gresp_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%gresp_storage) 

    ! gresp_xfer
    call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%gresp_xfer) 

    ! cpool
    call restartvar(ncid=ncid, flag=flag, varname='cpool', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%cpool) 

    ! xsmrpool
    call restartvar(ncid=ncid, flag=flag, varname='xsmrpool', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%xsmrpool) 

    ! pft_ctrunc
    call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%pft_ctrunc) 

    ! totvegc
    call restartvar(ncid=ncid, flag=flag, varname='totvegc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pcs%totvegc) 

    if ( use_c13 ) then
       !--------------------------------
       ! C13 pft carbon state variables 
       !--------------------------------

       ! leafc
       call restartvar(ncid=ncid, flag=flag, varname='leafc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%leafc)
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%leafc with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%leafc(i) = pcbulks%leafc(i) * c3_r2
          else
             pcisos%leafc(i) = pcbulks%leafc(i) * c4_r2
          endif
       end if

       ! leafc_storage
       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%leafc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%leafc_storage with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%leafc_storage(i) = pcbulks%leafc_storage(i) * c3_r2
          else
             pcisos%leafc_storage(i) = pcbulks%leafc_storage(i) * c4_r2
          endif
       end if

       ! leafc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%leafc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%leafc_xfer with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%leafc_xfer(i) = pcbulks%leafc_xfer(i) * c3_r2
          else
             pcisos%leafc_xfer(i) = pcbulks%leafc_xfer(i) * c4_r2
          endif
       end if

       ! frootc
       call restartvar(ncid=ncid, flag=flag, varname='frootc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%frootc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%frootc with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%frootc(i) = pcbulks%frootc(i) * c3_r2
          else
             pcisos%frootc(i) = pcbulks%frootc(i) * c4_r2
          endif
       end if

       ! frootc_storage
       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%frootc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%frootc_storage with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%frootc_storage(i) = pcbulks%frootc_storage(i) * c3_r2
          else
             pcisos%frootc_storage(i) = pcbulks%frootc_storage(i) * c4_r2
          endif
       end if

       !frootc_xfer 
       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%frootc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%frootc_xfer with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%frootc_xfer(i) = pcbulks%frootc_xfer(i) * c3_r2
          else
             pcisos%frootc_xfer(i) = pcbulks%frootc_xfer(i) * c4_r2
          endif
       end if

       ! livestemc
       call restartvar(ncid=ncid, flag=flag, varname='livestemc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%livestemc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%livestemc with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%livestemc(i) = pcbulks%livestemc(i) * c3_r2
          else
             pcisos%livestemc(i) = pcbulks%livestemc(i) * c4_r2
          endif
       end if

       ! livestemc_storage
       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%livestemc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%livestemc_storage with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%livestemc_storage(i) = pcbulks%livestemc_storage(i) * c3_r2
          else
             pcisos%livestemc_storage(i) = pcbulks%livestemc_storage(i) * c4_r2
          endif
       end if

       ! livestemc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%livestemc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%livestemc_xfer with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%livestemc_xfer(i) = pcbulks%livestemc_xfer(i) * c3_r2
          else
             pcisos%livestemc_xfer(i) = pcbulks%livestemc_xfer(i) * c4_r2
          endif
       end if

       ! deadstemc
       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%deadstemc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%deadstemc with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%deadstemc(i) = pcbulks%deadstemc(i) * c3_r2
          else
             pcisos%deadstemc(i) = pcbulks%deadstemc(i) * c4_r2
          endif
       end if

       ! deadstemc_storage
       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%deadstemc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%deadstemc_storage with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%deadstemc_storage(i) = pcbulks%deadstemc_storage(i) * c3_r2
          else
             pcisos%deadstemc_storage(i) = pcbulks%deadstemc_storage(i) * c4_r2
          endif
       end if

       ! deadstemc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%deadstemc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%deadstemc_xfer with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%deadstemc_xfer(i) = pcbulks%deadstemc_xfer(i) * c3_r2
          else
             pcisos%deadstemc_xfer(i) = pcbulks%deadstemc_xfer(i) * c4_r2
          endif
       end if

       ! livecrootc
       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%livecrootc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%livecrootc with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%livecrootc(i) = pcbulks%livecrootc(i) * c3_r2
          else
             pcisos%livecrootc(i) = pcbulks%livecrootc(i) * c4_r2
          endif
       end if

       ! livecrootc_storage
       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%livecrootc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%livecrootc_storage with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%livecrootc_storage(i) = pcbulks%livecrootc_storage(i) * c3_r2
          else
             pcisos%livecrootc_storage(i) = pcbulks%livecrootc_storage(i) * c4_r2
          endif
       end if

       ! livecrootc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%livecrootc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%livecrootc_xfer with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%livecrootc_xfer(i) = pcbulks%livecrootc_xfer(i) * c3_r2
          else
             pcisos%livecrootc_xfer(i) = pcbulks%livecrootc_xfer(i) * c4_r2
          endif
       end if

       ! deadcrootc
       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%deadcrootc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%deadcrootc with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%deadcrootc(i) = pcbulks%deadcrootc(i) * c3_r2
          else
             pcisos%deadcrootc(i) = pcbulks%deadcrootc(i) * c4_r2
          endif
       end if

       ! deadcrootc_storage
       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%deadcrootc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%deadcrootc_storage with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%deadcrootc_storage(i) = pcbulks%deadcrootc_storage(i) * c3_r2
          else
             pcisos%deadcrootc_storage(i) = pcbulks%deadcrootc_storage(i) * c4_r2
          endif
       end if

       ! deadcrootc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%deadcrootc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%deadcrootc_xfer with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%deadcrootc_xfer(i) = pcbulks%deadcrootc_xfer(i) * c3_r2
          else
             pcisos%deadcrootc_xfer(i) = pcbulks%deadcrootc_xfer(i) * c4_r2
          endif
       end if

       ! gresp_storage
       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%gresp_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%gresp_storage with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%gresp_storage(i) = pcbulks%gresp_storage(i) * c3_r2
          else
             pcisos%gresp_storage(i) = pcbulks%gresp_storage(i) * c4_r2
          endif
       end if

       ! gresp_xfer
       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%gresp_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%gresp_xfer with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%gresp_xfer(i) = pcbulks%gresp_xfer(i) * c3_r2
          else
             pcisos%gresp_xfer(i) = pcbulks%gresp_xfer(i) * c4_r2
          endif
       end if

       ! cpool
       call restartvar(ncid=ncid, flag=flag, varname='cpool_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%cpool) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%cpool with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%cpool(i) = pcbulks%cpool(i) * c3_r2
          else
             pcisos%cpool(i) = pcbulks%cpool(i) * c4_r2
          endif
       end if

       ! xsmrpool
       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%xsmrpool) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%xsmrpool with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%xsmrpool(i) = pcbulks%xsmrpool(i) * c3_r2
          else
             pcisos%xsmrpool(i) = pcbulks%xsmrpool(i) * c4_r2
          endif
       end if

       ! pft_ctrunc
       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%pft_ctrunc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%pft_ctrunc with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%pft_ctrunc(i) = pcbulks%pft_ctrunc(i) * c3_r2
          else
             pcisos%pft_ctrunc(i) = pcbulks%pft_ctrunc(i) * c4_r2
          endif
       end if

       ! totvegc
       call restartvar(ncid=ncid, flag=flag, varname='totvegc_13', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc13s%totvegc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc13s%totvegc with atmospheric c13 value'
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%totvegc(i) = pcbulks%totvegc(i) * c3_r2
          else
             pcisos%totvegc(i) = pcbulks%totvegc(i) * c4_r2
          endif
       end if
    endif

    if ( use_c14 ) then
       !--------------------------------
       ! C14 pft carbon state variables 
       !--------------------------------

       ! leafc
       call restartvar(ncid=ncid, flag=flag, varname='leafc_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%leafc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%leafc with atmospheric c14 value'
          if (pcs%leafc(i) .ne. spval .and. &
               .not. isnan(pcs%leafc(i)) ) then
             pc14s%leafc(i) = pcs%leafc(i) * c14ratio
          endif
       end if

       ! leafc_storage
       call restartvar(ncid=ncid, flag=flag, varname='leafc_storage_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%leafc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%leafc_storage with atmospheric c14 value'
          if (pcs%leafc_storage(i) .ne. spval .and. &
               .not. isnan(pcs%leafc_storage(i)) ) then
             pc14s%leafc_storage(i) = pcs%leafc_storage(i) * c14ratio
          endif
       end if

       ! leafc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='leafc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%leafc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%leafc_xfer with atmospheric c14 value'
          if (pcs%leafc_xfer(i) .ne. spval .and. &
               .not. isnan(pcs%leafc_xfer(i)) ) then
             pc14s%leafc_xfer(i) = pcs%leafc_xfer(i) * c14ratio
          endif
       end if

       ! frootc
       call restartvar(ncid=ncid, flag=flag, varname='frootc_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%frootc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%frootc with atmospheric c14 value'
          if (pcs%frootc(i) .ne. spval .and. &
               .not. isnan(pcs%frootc(i)) ) then
             pc14s%frootc(i) = pcs%frootc(i) * c14ratio
          endif
       end if

       ! frootc_storage
       call restartvar(ncid=ncid, flag=flag, varname='frootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%frootc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%frootc_storage with atmospheric c14 value'
          if (pcs%frootc_storage(i) .ne. spval .and. .not. isnan(pcs%frootc_storage(i)) ) then
             pc14s%frootc_storage(i) = pcs%frootc_storage(i) * c14ratio
          endif
       end if

       !frootc_xfer 
       call restartvar(ncid=ncid, flag=flag, varname='frootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%frootc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%frootc_xfer with atmospheric c14 value'
          if (pcs%frootc_xfer(i) .ne. spval .and. .not. isnan(pcs%frootc_xfer(i)) ) then
             pc14s%frootc_xfer(i) = pcs%frootc_xfer(i) * c14ratio
          endif
       end if

       ! livestemc
       call restartvar(ncid=ncid, flag=flag, varname='livestemc_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%livestemc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%livestemc with atmospheric c14 value'
          if (pcs%livestemc(i) .ne. spval .and. .not. isnan(pcs%livestemc(i)) ) then
             pc14s%livestemc(i) = pcs%livestemc(i) * c14ratio
          endif
       end if

       ! livestemc_storage
       call restartvar(ncid=ncid, flag=flag, varname='livestemc_storage_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%livestemc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%livestemc_storage with atmospheric c14 value'
          if (pcs%livestemc_storage(i) .ne. spval .and. .not. isnan(pcs%livestemc_storage(i)) ) then
             pc14s%livestemc_storage(i) = pcs%livestemc_storage(i) * c14ratio
          endif
       end if

       ! livestemc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='livestemc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%livestemc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%livestemc_xfer with atmospheric c14 value'
          if (pcs%livestemc_xfer(i) .ne. spval .and. .not. isnan(pcs%livestemc_xfer(i)) ) then
             pc14s%livestemc_xfer(i) = pcs%livestemc_xfer(i) * c14ratio
          endif
       end if

       ! deadstemc
       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%deadstemc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%deadstemc with atmospheric c14 value'
          if (pcs%deadstemc(i) .ne. spval .and. .not. isnan(pcs%deadstemc(i)) ) then
             pc14s%deadstemc(i) = pcs%deadstemc(i) * c14ratio
          endif
       end if

       ! deadstemc_storage
       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_storage_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%deadstemc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%deadstemc_storage with atmospheric c14 value'
          if (pcs%deadstemc_storage(i) .ne. spval .and. .not. isnan(pcs%deadstemc_storage(i)) ) then
             pc14s%deadstemc_storage(i) = pcs%deadstemc_storage(i) * c14ratio
          endif
       end if

       ! deadstemc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='deadstemc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%deadstemc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%deadstemc_xfer with atmospheric c14 value'
          if (pcs%deadstemc_xfer(i) .ne. spval .and. .not. isnan(pcs%deadstemc_xfer(i)) ) then
             pc14s%deadstemc_xfer(i) = pcs%deadstemc_xfer(i) * c14ratio
          endif
       end if

       ! livecrootc
       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%livecrootc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%livecrootc with atmospheric c14 value'
          if (pcs%livecrootc(i) .ne. spval .and. .not. isnan(pcs%livecrootc(i)) ) then
             pc14s%livecrootc(i) = pcs%livecrootc(i) * c14ratio
          endif
       end if

       ! livecrootc_storage
       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%livecrootc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%livecrootc_storage with atmospheric c14 value'
          if (pcs%livecrootc_storage(i) .ne. spval .and..not. isnan(pcs%livecrootc_storage(i)) ) then
             pc14s%livecrootc_storage(i) = pcs%livecrootc_storage(i) * c14ratio
          endif
       end if

       ! livecrootc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='livecrootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%livecrootc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%livecrootc_xfer with atmospheric c14 value'
          if (pcs%livecrootc_xfer(i) .ne. spval .and. .not. isnan(pcs%livecrootc_xfer(i)) ) then
             pc14s%livecrootc_xfer(i) = pcs%livecrootc_xfer(i) * c14ratio
          endif
       end if

       ! deadcrootc
       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%deadcrootc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%deadcrootc with atmospheric c14 value'
          if (pcs%deadcrootc(i) .ne. spval .and. .not. isnan(pcs%deadcrootc(i)) ) then
             pc14s%deadcrootc(i) = pcs%deadcrootc(i) * c14ratio
          endif
       end if

       ! deadcrootc_storage
       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_storage_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%deadcrootc_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%deadcrootc_storage with atmospheric c14 value'
          if (pcs%deadcrootc_storage(i) .ne. spval .and. .not. isnan(pcs%deadcrootc_storage(i)) ) then
             pc14s%deadcrootc_storage(i) = pcs%deadcrootc_storage(i) * c14ratio
          endif
       end if

       ! deadcrootc_xfer
       call restartvar(ncid=ncid, flag=flag, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%deadcrootc_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog) 'initializing pc14s%deadcrootc_xfer with atmospheric c14 value'
          if (pcs%deadcrootc_xfer(i) .ne. spval .and. &
               .not. isnan(pcs%deadcrootc_xfer(i)) ) then
             pc14s%deadcrootc_xfer(i) = pcs%deadcrootc_xfer(i) * c14ratio
          endif
       end if

       ! gresp_storage
       call restartvar(ncid=ncid, flag=flag, varname='gresp_storage_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%gresp_storage) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%gresp_storage with atmospheric c14 value'
          if (pcs%gresp_storage(i) .ne. spval .and. .not. isnan(pcs%gresp_storage(i)) ) then
             pc14s%gresp_storage(i) = pcs%gresp_storage(i) * c14ratio
          endif
       end if

       ! gresp_xfer
       call restartvar(ncid=ncid, flag=flag, varname='gresp_xfer_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%gresp_xfer) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%gresp_xfer with atmospheric c14 value'
          if (pcs%gresp_xfer(i) .ne. spval .and. .not. isnan(pcs%gresp_xfer(i)) ) then
             pc14s%gresp_xfer(i) = pcs%gresp_xfer(i) * c14ratio
          endif
       end if

       ! cpool
       call restartvar(ncid=ncid, flag=flag, varname='cpool_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%cpool) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%cpool with atmospheric c14 value'
          if (pcs%cpool(i) .ne. spval .and. .not. isnan(pcs%cpool(i)) ) then
             pc14s%cpool(i) = pcs%cpool(i) * c14ratio
          endif
       end if

       ! xsmrpool
       call restartvar(ncid=ncid, flag=flag, varname='xsmrpool_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%xsmrpool) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%xsmrpool with atmospheric c14 value'
          if (pcs%xsmrpool(i) .ne. spval .and. .not. isnan(pcs%xsmrpool(i)) ) then
             pc14s%xsmrpool(i) = pcs%xsmrpool(i) * c14ratio
          endif
       end if

       ! pft_ctrunc
       call restartvar(ncid=ncid, flag=flag, varname='pft_ctrunc_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%pft_ctrunc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%pft_ctrunc with atmospheric c14 value'
          if (pcs%pft_ctrunc(i) .ne. spval .and. .not. isnan(pcs%pft_ctrunc(i)) ) then
             pc14s%pft_ctrunc(i) = pcs%pft_ctrunc(i) * c14ratio
          endif
       end if

       ! totvegc
       call restartvar(ncid=ncid, flag=flag, varname='totvegc_14', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pc14s%totvegc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing pc14s%totvegc with atmospheric c14 value'
          if (pcs%totvegc(i) .ne. spval .and. .not. isnan(pcs%totvegc(i)) ) then
             pc14s%totvegc(i) = pcs%totvegc(i) * c14ratio
          endif
       end if

       ! rc14_atm
       call restartvar(ncid=ncid, flag=flag, varname='rc14_atm', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%rc14_atm) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing cc14s%rc14_atm with atmospheric c14 value'
          pepv%rc14_atm(i) = c14ratio
       end if

    endif

    !--------------------------------
    ! pft nitrogen state variables
    !--------------------------------

    ! leafn
    call restartvar(ncid=ncid, flag=flag, varname='leafn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%leafn) 

    ! leafn_storage
    call restartvar(ncid=ncid, flag=flag, varname='leafn_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%leafn_storage) 

    ! leafn_xfer
    call restartvar(ncid=ncid, flag=flag, varname='leafn_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%leafn_xfer) 

    ! frootn
    call restartvar(ncid=ncid, flag=flag, varname='frootn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%frootn) 

    ! frootn_storage
    call restartvar(ncid=ncid, flag=flag, varname='frootn_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%frootn_storage) 

    ! frootn_xfer
    call restartvar(ncid=ncid, flag=flag, varname='frootn_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%frootn_xfer) 

    ! livestemn
    call restartvar(ncid=ncid, flag=flag, varname='livestemn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%livestemn) 

    ! livestemn_storage
    call restartvar(ncid=ncid, flag=flag, varname='livestemn_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%livestemn_storage) 

    ! livestemn_xfer
    call restartvar(ncid=ncid, flag=flag, varname='livestemn_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%livestemn_xfer) 

    ! deadstemn
    call restartvar(ncid=ncid, flag=flag, varname='deadstemn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%deadstemn) 

    !deadstemn_storage 
    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%deadstemn_storage) 

    !deadstemn_xfer 
    call restartvar(ncid=ncid, flag=flag, varname='deadstemn_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%deadstemn_xfer) 

    ! livecrootn
    call restartvar(ncid=ncid, flag=flag, varname='livecrootn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%livecrootn) 

    ! livecrootn_storage
    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%livecrootn_storage) 

    !livecrootn_xfer 
    call restartvar(ncid=ncid, flag=flag, varname='livecrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%livecrootn_xfer) 

    ! deadcrootn
    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%deadcrootn) 

    ! deadcrootn_storage
    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_storage', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%deadcrootn_storage) 

    ! deadcrootn_xfer
    call restartvar(ncid=ncid, flag=flag, varname='deadcrootn_xfer', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%deadcrootn_xfer) 

    !retransn 
    call restartvar(ncid=ncid, flag=flag, varname='retransn', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%retransn) 

    ! npool
    call restartvar(ncid=ncid, flag=flag, varname='npool', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%npool) 

    ! pft_ntrunc
    call restartvar(ncid=ncid, flag=flag, varname='pft_ntrunc', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pns%pft_ntrunc) 

    !--------------------------------
    ! column physical state variables
    !--------------------------------

    ! fpi
    if (use_vertsoilc) then
       ptr2d => cps%fpi_vr
       call restartvar(ncid=ncid, flag=flag, varname='fpi_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='fraction of potential immobilization',  units='unitless', &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => cps%fpi_vr(:,1) ! nlevdecomp = 1; so treat as 1D variable
       call restartvar(ncid=ncid, flag=flag, varname='fpi', xtype=ncd_double,  &
            dim1name='column', &
            long_name='fraction of potential immobilization',  units='unitless', &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    ! som_adv_coef
    if (use_vertsoilc) then
       ptr2d => cps%som_adv_coef
       call restartvar(ncid=ncid, flag=flag, varname='som_adv_coef_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='SOM advective flux', units='m/s', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    end if
    
    ! som_diffus_coef
    if (use_vertsoilc) then
       ptr2d => cps%som_diffus_coef
       call restartvar(ncid=ncid, flag=flag, varname='som_diffus_coef_vr', xtype=ncd_double,  &
            dim1name='column',dim2name='levgrnd', switchdim=.true., &
            long_name='SOM diffusivity due to bio/cryo-turbation',  units='m^2/s', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    end if

    ! if (use_nitrif_denitrif) then
    !     ! tmean_monthly_max
    !     call cnrest_addfld_decomp(ncid=ncid, flag=flag, varname='tmean_monthly_max', &
    !         long_name='', units='', &
    !     flag=flag, data_rl=cps%tmean_monthly_max, readvar=readvar, &
    !
    !     ! tmean_monthly
    !     call cnrest_addfld_decomp(ncid=ncid, flag=flag, varname='tmean_monthly', &
    !         long_name='', units='', &
    !     flag=flag, data_rl=cps%tmean_monthly, readvar=readvar, &
    ! endif

    ! fpg
    call restartvar(ncid=ncid, flag=flag, varname='fpg', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%fpg) 

    ! annsum_counter
    call restartvar(ncid=ncid, flag=flag, varname='annsum_counter', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%annsum_counter) 

    ! cannsum_npp
    call restartvar(ncid=ncid, flag=flag, varname='cannsum_npp', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%cannsum_npp) 

    ! col_lag_npp
    call restartvar(ncid=ncid, flag=flag, varname='col_lag_npp', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%col_lag_npp) 

    ! cannavg_t2m
    call restartvar(ncid=ncid, flag=flag, varname='cannavg_t2m', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%cannavg_t2m) 

    ! for fire model changed by F. Li and S. Levis
    ! burndate
    call restartvar(ncid=ncid, flag=flag, varname='burndate', xtype=ncd_int,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%burndate) 

    !lfc
    call restartvar(ncid=ncid, flag=flag, varname='lfc', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%lfc) 

    !wf
    call restartvar(ncid=ncid, flag=flag, varname='wf', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%wf) 

    !btran2
    call restartvar(ncid=ncid, flag=flag, varname='btran2', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=pps%btran2) 

    !--------------------------------
    ! column carbon state variables
    !--------------------------------

    do k = 1, ndecomp_pools
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
       if (use_vertsoilc) then
          ptr2d => ccs%decomp_cpools_vr(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => ccs%decomp_cpools_vr(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
               dim1name='column', &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end do

    ! col_ctrunc
    if (use_vertsoilc) then
       ptr2d => ccs%col_ctrunc_vr
       call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc_vr', xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => ccs%col_ctrunc_vr(:,1) ! nlevdecomp = 1; so treat as 1D variable
       call restartvar(ncid=ncid, flag=flag, varname='col_ctrunc', xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if
    if (flag=='read' .and. .not. readvar) then
       call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! ! nfixation_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='nfixation_prof', &
    !          long_name='', units='', &
    !          flag=flag, data_rl=cps%nfixation_prof, readvar=readvar)
    ! ! ndep_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='ndep_prof', &
    !          long_name='', units='', &
    !          flag=flag, data_rl=cps%ndep_prof, readvar=readvar)

    ! altmax
    call restartvar(ncid=ncid, flag=flag, varname='altmax', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%altmax) 

    ! altmax_lastyear
    call restartvar(ncid=ncid, flag=flag, varname='altmax_lastyear', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%altmax_lastyear) 

    ! altmax_indx
    call restartvar(ncid=ncid, flag=flag, varname='altmax_indx', xtype=ncd_int,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%altmax_indx) 

    ! altmax_lastyear_indx
    call restartvar(ncid=ncid, flag=flag, varname='altmax_lastyear_indx', xtype=ncd_int,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cps%altmax_lastyear_indx) 

    ! seedc
    call restartvar(ncid=ncid, flag=flag, varname='seedc', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=ccs%seedc) 

    ! totlitc
    call restartvar(ncid=ncid, flag=flag, varname='totlitc', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=ccs%totlitc) 

    ! totcolc
    call restartvar(ncid=ncid, flag=flag, varname='totcolc', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=ccs%totcolc) 

    ! prod10c
    call restartvar(ncid=ncid, flag=flag, varname='prod10c', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=ccs%prod10c) 

    ! prod100c
    call restartvar(ncid=ncid, flag=flag, varname='prod100c', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=ccs%prod100c) 

    ! totsomc
    call restartvar(ncid=ncid, flag=flag, varname='totsomc', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=ccs%totsomc) 

    !--------------------------------
    ! C13 column carbon state variables
    !--------------------------------

    if ( use_c13 ) then

       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
          if (use_vertsoilc) then
             ptr2d => cc13s%decomp_cpools_vr(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => cc13s%decomp_cpools_vr(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cc13s%decomp_cpools_vr with atmospheric c13 value for: '//varname
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (ccs%decomp_cpools_vr(i,j,k) .ne. spval .and. .not. isnan(ccs%decomp_cpools_vr(i,j,k)) ) then
                         cc13s%decomp_cpools_vr(i,j,k) = ccs%decomp_cpools_vr(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do

       ! seedc
       call restartvar(ncid=ncid, flag=flag, varname='seedc_13', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc13s%seedc) 
       if (flag=='read' .and. .not. readvar) then
          if (ccs%seedc(i) .ne. spval .and. .not. isnan(ccs%seedc(i)) ) then
             cc13s%seedc(i) = ccs%seedc(i) * c3_r2
          end if
       end if

       ! col_ctrunc_13
       if (use_vertsoilc) then
          ptr2d => cc13s%col_ctrunc_vr
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => cc13s%col_ctrunc_vr(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c13", xtype=ncd_double,  &
               dim1name='column', &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if

       ! totlitc
       call restartvar(ncid=ncid, flag=flag, varname='totlitc_13', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc13s%totlitc) 
       if (flag=='read' .and. .not. readvar) then
          if (ccs%totlitc(i) .ne. spval .and. .not. isnan( ccs%totlitc(i) ) ) then
             cc13s%totlitc(i) = ccs%totlitc(i) * c3_r2
          end if
       end if

       ! totcolc
       call restartvar(ncid=ncid, flag=flag, varname='totcolc_13', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc13s%totcolc) 
       if (flag=='read' .and. .not. readvar) then
          if (ccs%totcolc(i) .ne. spval .and. .not. isnan (ccs%totcolc(i) ) ) then
             cc13s%totcolc(i) = ccs%totcolc(i) * c3_r2
          end if
       end if

       ! prod10c
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_13', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc13s%prod10c) 
       if (flag=='read' .and. .not. readvar) then
          if (ccs%prod10c(i) .ne. spval .and. .not. isnan( ccs%prod10c(i) ) ) then
             cc13s%prod10c(i) = ccs%prod10c(i) * c3_r2
          endif
       end if

       ! prod100c
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_13', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc13s%prod100c) 
       if (flag=='read' .and. .not. readvar) then
          if (ccs%prod100c(i) .ne. spval .and. .not. isnan( ccs%prod100c(i) ) ) then
             cc13s%prod100c(i) = ccs%prod100c(i) * c3_r2
          endif
       end if
    endif

    !--------------------------------
    ! C14 column carbon state variables
    !--------------------------------

    if ( use_c14 ) then

       do k = 1, ndecomp_pools
          varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
          if (use_vertsoilc) then
             ptr2d => cc14s%decomp_cpools_vr(:,:,k)
             call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
                  dim1name='column', dim2name='levgrnd', switchdim=.true., &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp', readvar=readvar, data=ptr2d)
          else
             ptr1d => cc14s%decomp_cpools_vr(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
             call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
                  dim1name='column', &
                  long_name='',  units='', fill_value=spval, &
                  interpinic_flag='interp' , readvar=readvar, data=ptr1d)
          end if
          if (flag=='read' .and. .not. readvar) then
             write(iulog,*) 'initializing cc14s%decomp_cpools_vr with atmospheric c14 value for: '//trim(varname)
             do i = bounds%begc,bounds%endc
                do j = 1, nlevdecomp
                   if (ccs%decomp_cpools_vr(i,j,k) .ne. spval .and. .not. isnan(ccs%decomp_cpools_vr(i,j,k)) ) then
                         cc14s%decomp_cpools_vr(i,j,k) = ccs%decomp_cpools_vr(i,j,k) * c3_r2
                   endif
                end do
             end do
          end if
       end do

       ! seedc
       call restartvar(ncid=ncid, flag=flag, varname='seedc_14', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc14s%seedc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing cc14s%seedc with atmospheric c14 value'
          do i = bounds%begc,bounds%endc
             if (ccs%seedc(i) .ne. spval .and. .not. isnan(ccs%seedc(i)) ) then
                cc14s%seedc(i) = ccs%seedc(i) * c14ratio
             endif
          end do
       end if

       ! col_ctrunc_c14
       if (use_vertsoilc) then
          ptr2d => cc14s%col_ctrunc_vr
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14_vr", xtype=ncd_double,  &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => cc14s%col_ctrunc_vr(:,1)
          call restartvar(ncid=ncid, flag=flag, varname="col_ctrunc_c14", xtype=ncd_double,  &
               dim1name='column', &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if

       ! totlitc
       call restartvar(ncid=ncid, flag=flag, varname='totlitc_14', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc14s%totlitc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing cc14s%totlitc with atmospheric c14 value'
          if (ccs%totlitc(i) .ne. spval .and. .not. isnan(ccs%totlitc(i)) ) then
             cc14s%totlitc(i) = ccs%totlitc(i) * c14ratio
          endif
       end if

       ! totcolc
       call restartvar(ncid=ncid, flag=flag, varname='totcolc_14', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc14s%totcolc) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing cc14s%totcolc with atmospheric c14 value'
          if (ccs%totcolc(i) .ne. spval .and. .not. isnan(ccs%totcolc(i)) ) then
             cc14s%totcolc(i) = ccs%totcolc(i) * c14ratio
          endif
       end if

       ! prod10c
       call restartvar(ncid=ncid, flag=flag, varname='prod10c_14', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc14s%prod10c) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing cc14s%prod10c with atmospheric c14 value'
          if (ccs%prod10c(i) .ne. spval .and. .not. isnan(ccs%prod10c(i)) ) then
             cc14s%prod10c(i) = ccs%prod10c(i) * c14ratio
          endif
       end if

       ! prod100c
       call restartvar(ncid=ncid, flag=flag, varname='prod100c_14', xtype=ncd_double,  &
            dim1name='column', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=cc14s%prod100c) 
       if (flag=='read' .and. .not. readvar) then
          write(iulog,*) 'initializing cc14s%prod100c with atmospheric c14 value'
          if (ccs%prod100c(i) .ne. spval .and. .not. isnan(ccs%prod100c(i)) ) then
             cc14s%prod100c(i) = ccs%prod100c(i) * c14ratio
          endif
       end if
    endif

    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------

    ! sminn
    if (use_vertsoilc) then
       ptr2d => cns%sminn_vr
       call restartvar(ncid=ncid, flag=flag, varname="sminn_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => cns%sminn_vr(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="sminn", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if
    if (flag=='read' .and. .not. readvar) then
       call endrun(msg='ERROR::'//trim(varname)//' is required on an initialization dataset'//&
            errMsg(__FILE__, __LINE__))
    end if

    ! decomposing N pools
    do k = 1, ndecomp_pools
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'n'
       if (use_vertsoilc) then
          ptr2d => cns%decomp_npools_vr(:,:,k)
          call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => cns%decomp_npools_vr(:,1,k)
          call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
               dim1name='column', &
               long_name='',  units='', fill_value=spval, &
               interpinic_flag='interp' , readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: '//trim(varname)//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if
    end do

    ! col_ntrunc
    if (use_vertsoilc) then
       ptr2d => cns%col_ntrunc_vr
       call restartvar(ncid=ncid, flag=flag, varname="col_ntrunc_vr", xtype=ncd_double,  &
            dim1name='column', dim2name='levgrnd', switchdim=.true., &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp', readvar=readvar, data=ptr2d)
    else
       ptr1d => cns%col_ntrunc_vr(:,1)
       call restartvar(ncid=ncid, flag=flag, varname="col_ntrunc", xtype=ncd_double,  &
            dim1name='column', &
            long_name='',  units='', fill_value=spval, &
            interpinic_flag='interp' , readvar=readvar, data=ptr1d)
    end if

    if (use_nitrif_denitrif) then

       ! f_nit_vr
       if (use_vertsoilc) then
          ptr2d => cnf%f_nit_vr(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='f_nit_vr_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => cnf%f_nit_vr(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='f_nit_vr', xtype=ncd_double, &
               dim1name='column', &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg='ERROR:: f_nit_vr'//' is required on an initialization dataset'//&
               errMsg(__FILE__, __LINE__))
       end if

       ! pot_f_nit_vr
       if (use_vertsoilc) then
          ptr2d => cnf%pot_f_nit_vr(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='pot_f_nit_vr_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='potential soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => cnf%pot_f_nit_vr(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='pot_f_nit_vr', xtype=ncd_double, &
               dim1name='column', &
               long_name='soil nitrification flux', units='gN/m3/s', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: pot_f_nit_vr'//' is required on an initialization dataset' )
       end if

       ! smin_no3_vr
       if (use_vertsoilc) then
          ptr2d => cns%smin_no3_vr(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_no3_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d)
       else
          ptr1d => cns%smin_no3_vr(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_no3', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_no3_vr'//' is required on an initialization dataset' )
       end if

       ! smin_nh4
       if (use_vertsoilc) then
          ptr2d => cns%smin_nh4_vr(:,:)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4_vr', xtype=ncd_double, &
               dim1name='column', dim2name='levgrnd', switchdim=.true., &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr2d) 
       else
          ptr1d => cns%smin_nh4_vr(:,1)
          call restartvar(ncid=ncid, flag=flag, varname='smin_nh4', xtype=ncd_double, &
               dim1name='column', &
               long_name='', units='', &
               interpinic_flag='interp', readvar=readvar, data=ptr1d)
       end if
       if (flag=='read' .and. .not. readvar) then
          call endrun(msg= 'ERROR:: smin_nh4_vr'//' is required on an initialization dataset' )
       end if

    end if

    ! Set the integrated sminn based on sminn_vr, as is done in CNSummaryMod (this may
    ! not be the most appropriate method or place to do this)
    cns%sminn(bounds%begc:bounds%endc) = 0._r8
    do j = 1, nlevdecomp
       do c = bounds%begc, bounds%endc
          cns%sminn(c) = cns%sminn(c) + cns%sminn_vr(c,j) * dzsoi_decomp(j)
       end do
    end do

    ! totcoln
    call restartvar(ncid=ncid, flag=flag, varname='totcoln', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cns%totcoln) 

    ! seedn
    call restartvar(ncid=ncid, flag=flag, varname='seedn', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cns%seedn) 

    ! prod10n
    call restartvar(ncid=ncid, flag=flag, varname='prod10n', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cns%prod10n) 

    ! prod100n
    call restartvar(ncid=ncid, flag=flag, varname='prod100n', xtype=ncd_double,  &
         dim1name='column', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=cns%prod100n) 

    ! decomp_cascade_state  
    ! the purpose of this is to check to make sure the bgc used matches what the restart file was generated with.  
    ! add info about the SOM decomposition cascade

    if (use_century_decomp) then
       decomp_cascade_state = 1
    else
       decomp_cascade_state = 0
    end if
    ! add info about the nitrification / denitrification state
    if (use_nitrif_denitrif) then
       decomp_cascade_state = decomp_cascade_state + 10
    end if
    if (flag == 'write') itemp = decomp_cascade_state    
    call restartvar(ncid=ncid, flag=flag, varname='decomp_cascade_state', xtype=ncd_int,  &
         long_name='BGC of the model that wrote this restart file:' &
         // '  1s column: 0 = CLM-CN cascade, 1 = Century cascade;' &
         // ' 10s column: 0 = CLM-CN denitrification, 10 = Century denitrification', units='', &
         interpinic_flag='skip', readvar=readvar, data=itemp)
    if (flag=='read') then
       if (.not. readvar) then
          ! assume, for sake of backwards compatibility, that if decomp_cascade_state 
          ! is not in the restart file, then the current model state is the same as 
          ! the prior model state
          restart_file_decomp_cascade_state = decomp_cascade_state
          if ( masterproc ) write(iulog,*) ' CNRest: WARNING!  Restart file does not ' &
               // ' contain info on decomp_cascade_state used to generate the restart file.  '
          if ( masterproc ) write(iulog,*) '   Assuming the same as current setting: ', decomp_cascade_state
       else
          restart_file_decomp_cascade_state = itemp  
          if (decomp_cascade_state .ne. restart_file_decomp_cascade_state ) then
             if ( masterproc ) then
                write(iulog,*) 'CNRest: ERROR--the decomposition cascade differs between the current ' &
                     // ' model state and the model that wrote the restart file. '
                write(iulog,*) 'The model will be horribly out of equilibrium until after a lengthy spinup. '
                write(iulog,*) 'Stopping here since this is probably an error in configuring the run. '
                write(iulog,*) 'If you really wish to proceed, then override by setting '
                write(iulog,*) 'override_bgc_restart_mismatch_dump to .true. in the namelist'
                if ( .not. override_bgc_restart_mismatch_dump ) then
                   call endrun(msg= ' CNRest: Stopping. Decomposition cascade mismatch error.'//&
                        errMsg(__FILE__, __LINE__))
                endif
             endif
          endif
       end if
    end if

    ! spinup_state
    if (flag == 'write') idata = spinup_state
    call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
         long_name='Spinup state of the model that wrote this restart file: ' &
         // ' 0 = normal model mode, 1 = AD spinup', units='', &
         interpinic_flag='copy', readvar=readvar,  data=idata)
    if (flag == 'read') then
       if (readvar) then
          restart_file_spinup_state = idata
       else
          ! assume, for sake of backwards compatibility, that if spinup_state is not in 
          ! the restart file then current model state is the same as prior model state
          restart_file_spinup_state = spinup_state
          if ( masterproc ) then
             write(iulog,*) ' CNRest: WARNING!  Restart file does not contain info ' &
                  // ' on spinup state used to generate the restart file. '
             write(iulog,*) '   Assuming the same as current setting: ', spinup_state
          end if
       end if
    end if

    ! now compare the model and restart file spinup states, and either take the 
    ! model into spinup mode or out of it if they are not identical
    ! taking model out of spinup mode requires multiplying each decomposing pool 
    ! by the associated AD factor.
    ! putting model into spinup mode requires dividing each decomposing pool 
    ! by the associated AD factor.
    ! only allow this to occur on first timestep of model run.

    nstep = get_nstep()  

    if (flag == 'read' .and. spinup_state .ne. restart_file_spinup_state ) then
       if (spinup_state .eq. 0 .and. restart_file_spinup_state .eq. 1 ) then
          if ( masterproc ) write(iulog,*) ' CNRest: taking SOM pools out of AD spinup mode'
          exit_spinup = .true.
       else if (spinup_state .eq. 1 .and. restart_file_spinup_state .eq. 0 ) then
          if ( masterproc ) write(iulog,*) ' CNRest: taking SOM pools into AD spinup mode'
          enter_spinup = .true.
       else
          call endrun(msg=' CNRest: error in entering/exiting spinup.  spinup_state ' &
               // ' != restart_file_spinup_state, but do not know what to do'//&
               errMsg(__FILE__, __LINE__))
       end if
       if (nstep .ge. 2) then
          call endrun(msg=' CNRest: error in entering/exiting spinup - should occur only when nstep = 1'//&
               errMsg(__FILE__, __LINE__))
       endif
       do k = 1, ndecomp_pools
          if ( exit_spinup ) then
             m = decomp_cascade_con%spinup_factor(k)
          else if ( enter_spinup ) then
             m = 1. / decomp_cascade_con%spinup_factor(k)
          end if
          do c = bounds%begc, bounds%endc
             do j = 1, nlevdecomp
                ccs%decomp_cpools_vr(c,j,k) = ccs%decomp_cpools_vr(c,j,k) * m
                if ( use_c13 ) then
                   cc13s%decomp_cpools_vr(c,j,k) = cc13s%decomp_cpools_vr(c,j,k) * m
                endif
                if ( use_c14 ) then
                   cc14s%decomp_cpools_vr(c,j,k) = cc14s%decomp_cpools_vr(c,j,k) * m
                endif
                cns%decomp_npools_vr(c,j,k) = cns%decomp_npools_vr(c,j,k) * m
             end do
          end do
       end do
    end if

    if ( .not. is_restart() .and. nstep .eq. 1 ) then
       do i = bounds%begp,bounds%endp
          if (pftcon%c3psn(pft%itype(i)) == 1._r8) then
             pcisos%grainc(i) = pcbulks%grainc(i) * c3_r2
             pcisos%grainc_storage(i) = pcbulks%grainc_storage(i) * c3_r2
             pcisos%grainc_xfer(i) = pcbulks%grainc_xfer(i) * c3_r2
             pcisos%dispvegc(i) = pcbulks%dispvegc(i) * c3_r2
             pcisos%storvegc(i) = pcbulks%storvegc(i) * c3_r2
             pcisos%totvegc(i) = pcbulks%totvegc(i) * c3_r2
             pcisos%totpftc(i) = pcbulks%totpftc(i) * c3_r2
             pcisos%woodc(i) = pcbulks%woodc(i) * c3_r2
          else
             pcisos%grainc(i) = pcbulks%grainc(i) * c4_r2
             pcisos%grainc_storage(i) = pcbulks%grainc_storage(i) * c4_r2
             pcisos%grainc_xfer(i) = pcbulks%grainc_xfer(i) * c4_r2
             pcisos%dispvegc(i) = pcbulks%dispvegc(i) * c4_r2
             pcisos%storvegc(i) = pcbulks%storvegc(i) * c4_r2
             pcisos%totvegc(i) = pcbulks%totvegc(i) * c4_r2
             pcisos%totpftc(i) = pcbulks%totpftc(i) * c4_r2
             pcisos%woodc(i) = pcbulks%woodc(i) * c4_r2
          end if
       end do
    end if

    if (use_cndv) then
       ! pft type dgvm physical state - crownarea
       call restartvar(ncid=ncid, flag=flag, varname='CROWNAREA', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pdgvs%crownarea)

       ! tempsum_litfall
       call restartvar(ncid=ncid, flag=flag, varname='tempsum_litfall', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%tempsum_litfall)

       ! annsum_litfall
       call restartvar(ncid=ncid, flag=flag, varname='annsum_litfall', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pepv%annsum_litfall)

       ! nind
       call restartvar(ncid=ncid, flag=flag, varname='nind', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pdgvs%nind)

       ! fpcgrid
       call restartvar(ncid=ncid, flag=flag, varname='fpcgrid', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pdgvs%fpcgrid)

       ! fpcgridold
       call restartvar(ncid=ncid, flag=flag, varname='fpcgridold', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pdgvs%fpcgridold)

       ! pft type dgvm physical state - tmomin20
       do_io = .true.
       if (flag == 'read') then
          ! On a read, confirm that this variable has the expected size; if not, don't
          ! read it (instead leave it at its arbitrary initial value). This is needed to
          ! support older initial conditions for which this variable had a different size.
          call ncd_inqvdlen(ncid, 'TMOMIN20', 1, dimlen, err_code)
          if (dimlen /= nump_global) then
             do_io = .false.
          end if
       end if
       if (do_io) then
          call restartvar(ncid=ncid, flag=flag, varname='TMOMIN20', xtype=ncd_double,  &
               dim1name='pft', &
               long_name='',units='', &
               interpinic_flag='interp', readvar=readvar, data=pdgvs%tmomin20)
       end if

       ! pft type dgvm physical state - agdd20
       do_io = .true.
       if (flag == 'read') then
          ! On a read, confirm that this variable has the expected size; if not, don't
          ! read it (instead leave it at its arbitrary initial value). This is needed to
          ! support older initial conditions for which this variable had a different size.
          call ncd_inqvdlen(ncid, 'AGDD20', 1, dimlen, err_code)
          if (dimlen /= nump_global) then
             do_io = .false.
          end if
       end if
       if (do_io) then
          call restartvar(ncid=ncid, flag=flag, varname='AGDD20', xtype=ncd_double,  &
               dim1name='pft',&
               long_name='',units='', &
               interpinic_flag='interp', readvar=readvar, data=pdgvs%agdd20)
       end if

       ! pft type dgvm physical state - t_mo_min
       call restartvar(ncid=ncid, flag=flag, varname='T_MO_MIN', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pdgvs%t_mo_min)

       ! present  
       if (flag == 'read' .or. flag == 'write') then
          allocate (iptemp(bounds%begp:bounds%endp), stat=ier)
       end if
       if (flag == 'write') then
          do p = bounds%begp,bounds%endp
             iptemp(p) = 0
             if (pdgvs%present(p)) iptemp(p) = 1
          end do
       end if
       call restartvar(ncid=ncid, flag=flag, varname='present', xtype=ncd_int,  &
            dim1name='pft',&
            long_name='',units='', &
            interpinic_flag='interp', readvar=readvar, data=iptemp)
       if (flag=='read' .and. readvar) then
          do p = bounds%begp,bounds%endp
             pdgvs%present(p) = .false.
             if (iptemp(p) == 1) pdgvs%present(p) = .true.
          end do
       end if
       if (flag == 'read' .or. flag == 'write') then
          deallocate (iptemp)
       end if

       ! leafcmax
       call restartvar(ncid=ncid, flag=flag, varname='leafcmax', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pcs%leafcmax)

       ! heatstress
       call restartvar(ncid=ncid, flag=flag, varname='heatstress', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pdgvs%heatstress)

       ! greffic
       call restartvar(ncid=ncid, flag=flag, varname='greffic', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='', units='', &
            interpinic_flag='interp', readvar=readvar, data=pdgvs%greffic)
    end if

  end subroutine CNRest

end module CNRestMod

