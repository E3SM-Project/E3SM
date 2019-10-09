module CNDVType
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing routines to drive the annual dynamic vegetation
  ! that works with CN, reset related variables,
  ! and initialize/reset time invariant variables
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use clm_varctl          , only : use_cndv, iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC DATA TYPES:
  !
  ! DGVM-specific ecophysiological constants structure (patch-level)
  type, public :: dgv_ecophyscon_type

     real(r8), pointer :: crownarea_max(:)   ! patch tree maximum crown area [m2]
     real(r8), pointer :: tcmin(:)           ! patch minimum coldest monthly mean temperature [units?]
     real(r8), pointer :: tcmax(:)           ! patch maximum coldest monthly mean temperature [units?]
     real(r8), pointer :: gddmin(:)          ! patch minimum growing degree days (at or above 5 C)
     real(r8), pointer :: twmax(:)           ! patch upper limit of temperature of the warmest month [units?]
     real(r8), pointer :: reinickerp(:)      ! patch parameter in allometric equation
     real(r8), pointer :: allom1(:)          ! patch parameter in allometric
     real(r8), pointer :: allom2(:)          ! patch parameter in allometric
     real(r8), pointer :: allom3(:)          ! patch parameter in allometric

  end type dgv_ecophyscon_type
  type(dgv_ecophyscon_type), public :: dgv_ecophyscon
  !
  ! DGVM state variables structure
  type, public :: dgvs_type

     real(r8), pointer, public :: agdd_patch        (:) ! patch accumulated growing degree days above 5
     real(r8), pointer, public :: agddtw_patch      (:) ! patch accumulated growing degree days above twmax
     real(r8), pointer, public :: agdd20_patch      (:) ! patch 20-yr running mean of agdd
     real(r8), pointer, public :: tmomin20_patch    (:) ! patch 20-yr running mean of tmomin
     logical , pointer, public :: present_patch     (:) ! patch whether PFT present in patch
     logical , pointer, public :: pftmayexist_patch (:) ! patch if .false. then exclude seasonal decid patches from tropics
     real(r8), pointer, public :: nind_patch        (:) ! patch number of individuals (#/m**2)
     real(r8), pointer, public :: lm_ind_patch      (:) ! patch individual leaf mass
     real(r8), pointer, public :: lai_ind_patch     (:) ! patch LAI per individual
     real(r8), pointer, public :: fpcinc_patch      (:) ! patch foliar projective cover increment (fraction) 
     real(r8), pointer, public :: fpcgrid_patch     (:) ! patch foliar projective cover on gridcell (fraction)
     real(r8), pointer, public :: fpcgridold_patch  (:) ! patch last yr's fpcgrid
     real(r8), pointer, public :: crownarea_patch   (:) ! patch area that each individual tree takes up (m^2)
     real(r8), pointer, public :: greffic_patch     (:)
     real(r8), pointer, public :: heatstress_patch  (:)

   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: InitAccBuffer
     procedure , public  :: InitAccVars
     procedure , public  :: UpdateCNDVAccVars
     procedure , private :: InitAllocate 
     procedure , private :: InitCold     
     procedure , private :: InitHistory
     
  end type dgvs_type
  !-----------------------------------------------------------------------

contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(dgvs_type) :: this
    type(bounds_type), intent(in) :: bounds  

    ! Note - need allocation so that associate statements can be used
    ! at run time for NAG (allocation of variables is needed) - history
    ! should only be initialized if use_cndv is true

    call this%InitAllocate (bounds)

    if (use_cndv) then
       call this%InitCold (bounds)
       call this%InitHistory (bounds)
    end if

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    use clm_varpar     , only : numpft
    use pftvarcon      , only : allom1s, allom2s, allom1, allom2, allom3, reinickerp
    use pftvarcon      , only : ntree, nbrdlf_dcd_brl_shrub
    use pftvarcon      , only : pftpar20, pftpar28, pftpar29, pftpar30, pftpar31
    !
    ! !ARGUMENTS:
    class(dgvs_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: m
    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp
  
    allocate(this%agdd_patch        (begp:endp)) ;     this%agdd_patch        (:) = nan
    allocate(this%agddtw_patch      (begp:endp)) ;     this%agddtw_patch      (:) = nan
    allocate(this%agdd20_patch      (begp:endp)) ;     this%agdd20_patch      (:) = nan
    allocate(this%tmomin20_patch    (begp:endp)) ;     this%tmomin20_patch    (:) = nan
    allocate(this%present_patch     (begp:endp)) ;     this%present_patch     (:) = .false.
    allocate(this%pftmayexist_patch (begp:endp)) ;     this%pftmayexist_patch (:) = .true.
    allocate(this%nind_patch        (begp:endp)) ;     this%nind_patch        (:) = nan
    allocate(this%lm_ind_patch      (begp:endp)) ;     this%lm_ind_patch      (:) = nan
    allocate(this%lai_ind_patch     (begp:endp)) ;     this%lai_ind_patch     (:) = nan
    allocate(this%fpcinc_patch      (begp:endp)) ;     this%fpcinc_patch      (:) = nan
    allocate(this%fpcgrid_patch     (begp:endp)) ;     this%fpcgrid_patch     (:) = nan
    allocate(this%fpcgridold_patch  (begp:endp)) ;     this%fpcgridold_patch  (:) = nan
    allocate(this%crownarea_patch   (begp:endp)) ;     this%crownarea_patch   (:) = nan
    allocate(this%greffic_patch     (begp:endp)) ;     this%greffic_patch     (:) = nan
    allocate(this%heatstress_patch  (begp:endp)) ;     this%heatstress_patch  (:) = nan

    allocate(dgv_ecophyscon%crownarea_max (0:numpft)) 
    allocate(dgv_ecophyscon%tcmin         (0:numpft))         
    allocate(dgv_ecophyscon%tcmax         (0:numpft))         
    allocate(dgv_ecophyscon%gddmin        (0:numpft))        
    allocate(dgv_ecophyscon%twmax         (0:numpft))         
    allocate(dgv_ecophyscon%reinickerp    (0:numpft))    
    allocate(dgv_ecophyscon%allom1        (0:numpft))        
    allocate(dgv_ecophyscon%allom2        (0:numpft))        
    allocate(dgv_ecophyscon%allom3        (0:numpft))        

    do m = 0,numpft
       dgv_ecophyscon%crownarea_max(m) = pftpar20(m)
       dgv_ecophyscon%tcmin(m)         = pftpar28(m)
       dgv_ecophyscon%tcmax(m)         = pftpar29(m)
       dgv_ecophyscon%gddmin(m)        = pftpar30(m)
       dgv_ecophyscon%twmax(m)         = pftpar31(m)
       dgv_ecophyscon%reinickerp(m)    = reinickerp
       dgv_ecophyscon%allom1(m)        = allom1
       dgv_ecophyscon%allom2(m)        = allom2
       dgv_ecophyscon%allom3(m)        = allom3
       ! modification for shrubs by X.D.Z
       if (m > ntree .and. m <= nbrdlf_dcd_brl_shrub ) then 
          dgv_ecophyscon%allom1(m) = allom1s
          dgv_ecophyscon%allom2(m) = allom2s
       end if
    end do

  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !USES:
    use shr_kind_mod  , only : r8 => shr_kind_r8
    use shr_const_mod , only : SHR_CONST_TKFRZ
    use decompMod     , only : bounds_type
    !
    ! !ARGUMENTS:
    class(dgvs_type) :: this 
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: p           ! patch index
    !-----------------------------------------------------------------------

    do p = bounds%begp,bounds%endp
       this%present_patch(p)   = .false.
       this%crownarea_patch(p) = 0._r8
       this%nind_patch(p)      = 0._r8
       this%agdd20_patch(p)    = 0._r8
       this%tmomin20_patch(p)  = SHR_CONST_TKFRZ - 5._r8 !initialize this way for Phenology code
    end do

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize history variables
    !
    ! !USES:
    use histFileMod, only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(dgvs_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------
    
    call hist_addfld1d (fname='AGDD', units='K',  &
         avgflag='A', long_name='growing degree-days base 5C', &
         ptr_patch=this%agdd_patch)

  end subroutine InitHistory


  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    ! 
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use clm_varcon , only : spval  
    use spmdMod    , only : masterproc
    use decompMod  , only : get_proc_global
    use restUtilMod
    use ncdio_pio
    use pio
    !
    ! !ARGUMENTS:
    class(dgvs_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer           :: j,c,p		! indices
    logical           :: readvar	! determine if variable is on initial file
    logical           :: do_io		! whether to do i/o for the given variable
    integer           :: nump_global	! total number of patches, globally
    integer           :: dimlen         ! dimension length
    integer           :: ier            ! error status
    integer           :: itemp		! temporary 
    integer , pointer :: iptemp(:)	! pointer to memory to be allocated
    integer           :: err_code       ! error code
    !-----------------------------------------------------------------------

    ! Get expected total number of points, for later error checks
    call get_proc_global(np=nump_global)

    call restartvar(ncid=ncid, flag=flag, varname='CROWNAREA', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%crownarea_patch)

    call restartvar(ncid=ncid, flag=flag, varname='nind', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%nind_patch)

    call restartvar(ncid=ncid, flag=flag, varname='fpcgrid', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fpcgrid_patch)

    call restartvar(ncid=ncid, flag=flag, varname='fpcgridold', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%fpcgridold_patch)

    ! tmomin20
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
            interpinic_flag='interp', readvar=readvar, data=this%tmomin20_patch)
    end if

    ! agdd20
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
            interpinic_flag='interp', readvar=readvar, data=this%agdd20_patch)
    end if

    ! present  
    if (flag == 'read' .or. flag == 'write') then
       allocate (iptemp(bounds%begp:bounds%endp), stat=ier)
    end if
    if (flag == 'write') then
       do p = bounds%begp,bounds%endp
          iptemp(p) = 0
          if (this%present_patch(p)) iptemp(p) = 1
       end do
    end if
    call restartvar(ncid=ncid, flag=flag, varname='present', xtype=ncd_int,  &
         dim1name='pft',&
         long_name='',units='', &
         interpinic_flag='interp', readvar=readvar, data=iptemp)
    if (flag=='read' .and. readvar) then
       do p = bounds%begp,bounds%endp
          this%present_patch(p) = .false.
          if (iptemp(p) == 1) this%present_patch(p) = .true.
       end do
    end if
    if (flag == 'read' .or. flag == 'write') then
       deallocate (iptemp)
    end if

    call restartvar(ncid=ncid, flag=flag, varname='heatstress', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%heatstress_patch)

    call restartvar(ncid=ncid, flag=flag, varname='greffic', xtype=ncd_double,  &
         dim1name='pft', &
         long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%greffic_patch)

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine InitAccBuffer (this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize accumulation buffer for all required module accumulated fields
    ! This routine set defaults values that are then overwritten by the
    ! restart file for restart or branch runs
    ! Each interval and accumulation type is unique to each field processed.
    ! Routine [initAccBuffer] defines the fields to be processed
    ! and the type of accumulation. 
    ! Routine [updateCNDVAccVars] does the actual accumulation for a given field.
    ! Fields are accumulated by calls to subroutine [update_accum_field]. 
    ! To accumulate a field, it must first be defined in subroutine [initAccVars] 
    ! and then accumulated by calls to [updateCNDVAccVars].
    !
    ! This should only be called if use_cndv is true.
    !
    ! !USES 
    use accumulMod       , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(dgvs_type) :: this
    type(bounds_type), intent(in) :: bounds  

    !
    ! !LOCAL VARIABLES:
    integer, parameter :: not_used = huge(1)

    !---------------------------------------------------------------------

    ! The following are accumulated fields.
    ! These types of fields are accumulated until a trigger value resets
    ! the accumulation to zero (see subroutine update_accum_field).
    ! Hence, [accper] is not valid.

    call init_accum_field (name='AGDDTW', units='K', &
         desc='growing degree-days base twmax', accum_type='runaccum', accum_period=not_used, &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='AGDD', units='K', &
         desc='growing degree-days base 5C', accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

  end subroutine InitAccBuffer

  !-----------------------------------------------------------------------
  subroutine InitAccVars(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module variables that are associated with
    ! time accumulated fields. This routine is called for both an initial run
    ! and a restart run (and must therefore must be called after the restart file 
    ! is read in and the accumulation buffer is obtained)
    !
    ! This should only be called if use_cndv is true.
    !
    ! !USES 
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(dgvs_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: nstep
    integer :: ier            ! error status
    real(r8), pointer :: rbufslp(:)  ! temporary

    !---------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level pft field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg=" allocation error for rbufslp"//&
            errMsg(__FILE__, __LINE__))
    endif

    nstep = get_nstep()

    call extract_accum_field ('AGDDTW', rbufslp, nstep) 
    this%agddtw_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('AGDD', rbufslp, nstep) 
    this%agdd_patch(begp:endp) = rbufslp(begp:endp)

    deallocate(rbufslp)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateCNDVAccVars(this, bounds, temperature_vars)
    !
    ! !DESCRIPTION:
    ! Update accumulated variables. Should be called every time step.
    !
    ! This should only be called if use_cndv is true.
    !
    ! !USES:
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep, get_curr_date
    use pftvarcon        , only : ndllf_dcd_brl_tree
    use TemperatureType  , only : temperature_type
    use accumulMod       , only : update_accum_field, extract_accum_field, accumResetVal
    !
    ! !ARGUMENTS:
    class(dgvs_type)       , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    type(temperature_type) , intent(in)    :: temperature_vars
    !
    ! !LOCAL VARIABLES:
    integer           :: p     ! index
    integer           :: ier   ! error status
    integer           :: dtime ! timestep size [seconds]
    integer           :: nstep ! timestep number
    integer           :: year  ! year (0, ...) for nstep
    integer           :: month ! month (1, ..., 12) for nstep
    integer           :: day   ! day of month (1, ..., 31) for nstep
    integer           :: secs  ! seconds into current date for nstep
    integer           :: begp, endp
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level

    character(len=*), parameter :: subname = 'UpdateCNDVAccVars'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp; endp = bounds%endp

    dtime = get_step_size()
    nstep = get_nstep()
    call get_curr_date (year, month, day, secs)

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate growing degree days based on 10-day running mean temperature.
    ! The trigger to reset the accumulated values to zero is -99999.
    
    ! Accumulate and extract AGDDTW (gdd base twmax, which is 23 deg C
    ! for boreal woody patches)
    
    do p = begp,endp
       rbufslp(p) = max(0._r8, &
            (temperature_vars%t_a10_patch(p) - SHR_CONST_TKFRZ - dgv_ecophyscon%twmax(ndllf_dcd_brl_tree)) &
            * dtime/SHR_CONST_CDAY)
       if (month==1 .and. day==1 .and. secs==int(dtime)) rbufslp(p) = accumResetVal
    end do
    call update_accum_field  ('AGDDTW', rbufslp, nstep)
    call extract_accum_field ('AGDDTW', this%agddtw_patch, nstep)

    ! Accumulate and extract AGDD

    do p = begp,endp
       rbufslp(p) = max(0.0_r8, &
            (temperature_vars%t_ref2m_patch(p) - (SHR_CONST_TKFRZ + 5.0_r8)) * dtime/SHR_CONST_CDAY)
       !
       ! Fix (for bug 1858) from Sam Levis to reset the annual AGDD variable
       ! 
       if (month==1 .and. day==1 .and. secs==int(dtime)) rbufslp(p) = accumResetVal
    end do
    call update_accum_field  ('AGDD', rbufslp, nstep)
    call extract_accum_field ('AGDD', this%agdd_patch, nstep)

    deallocate(rbufslp)

  end subroutine UpdateCNDVAccVars


end module CNDVType
