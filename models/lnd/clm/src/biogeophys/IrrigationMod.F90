module IrrigationMod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculates irrigation flux.
  !
  ! Usage:
  !
  !   - Call CalcIrrigationNeeded in order to compute whether and how much irrigation is
  !     needed for the next call to ApplyIrrigation. This should be called once per
  !     timestep.
  ! 
  !   - Call ApplyIrrigation in order to calculate qflx_irrig. This should be called
  !     exactly once per time step, before the first time qflx_irrig is needed by other
  !     parts of the code. It is acceptable for this to be called earlier in the timestep
  !     than CalcIrrigationNeeded.
  !
  !   - Access the timestep's irrigation flux via qflx_irrig_patch or
  !     qflx_irrig_col. These should be treated as read-only.
  !
  ! Design notes:
  !
  !   In principle, ApplyIrrigation and CalcIrrigationNeeded could be combined into a
  !   single routine. However, right now that is challenging, because qflx_irrig is
  !   needed earlier in the driver loop than when btran becomes available (and
  !   CalcIrrigationNeeded depends on btran). (And qflx_irrig is also used late in the
  !   driver loop - so it wouldn't work, for example, to calculate qflx_irrig after btran
  !   is computed, and then save it on the restart file for the next iteration of the
  !   driver loop: then the uses of qflx_irrig early and late in the driver loop would be
  !   inconsistent.)
  !
  !   If we could have access to btran earlier in the driver loop, so that
  !   CalcIrrigationNeeded could be called before the first time qflx_irrig is needed,
  !   then there might be some advantage to combining ApplyIrrigation and
  !   CalcIrrigationNeeded - or at least calling these two routines from the same place.
  !   In particular: this separation of the irrigation calculation into two routines that
  !   are done at different times in the driver loop makes it harder and less desirable to
  !   nest the irrigation object within some other object: Doing so might make it harder
  !   to do the two separate steps at the right time, and would lead to less clarity about
  !   how these two steps are ordered with respect to the rest of the driver loop. So if
  !   we start trying to create a hierarchy of objects in CLM, we may want to rethink this
  !   design.
  !
  ! !USES:
#include "shr_assert.h"
  use shr_kind_mod     , only : r8 => shr_kind_r8
  use decompMod        , only : bounds_type, get_proc_global
  use shr_log_mod      , only : errMsg => shr_log_errMsg
  use clm_varcon       , only : isecspday, degpsec, denh2o, spval
  use clm_varpar       , only : nlevgrnd
  use clm_time_manager , only : get_step_size
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  use GridcellType     , only : grc                
  use ColumnType       , only : col                
  use PatchType        , only : patch                
  !
  implicit none
  private

  ! !PUBLIC TYPES:
  
  ! This type is public (and its components are public, too) to aid unit testing
  type, public :: irrigation_params_type
     ! Note that we give default initialization values here. Once these parameters are
     ! moved to the params file, we'll probably want to get rid of this default
     ! initialization.
     
     ! Minimum LAI for irrigation
     real(r8) :: irrig_min_lai = 0.0_r8

     ! BTRAN threshold for irrigation
     ! Irrigate when btran falls below 0.999999 rather than 1 to allow for round-off error
     real(r8) :: irrig_btran_thresh = 0.999999_r8

     ! Time of day to check whether we need irrigation, seconds (0 = midnight). 
     ! We start applying the irrigation in the time step FOLLOWING this time, 
     ! since we won't begin irrigating until the next call to ApplyIrrigation
     integer  :: irrig_start_time = isecspday/4

     ! Desired amount of time to irrigate per day (sec). Actual time may 
     ! differ if this is not a multiple of dtime. Irrigation won't work properly 
     ! if dtime > secsperday
     integer  :: irrig_length = isecspday/6       

     ! Determines target soil moisture level for irrigation. If h2osoi_liq_so 
     ! is the soil moisture level at which stomata are fully open and 
     ! h2osoi_liq_sat is the soil moisture level at saturation (eff_porosity), 
     ! then the target soil moisture level is 
     !     (h2osoi_liq_so + irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)). 
     ! A value of 0 means that the target soil moisture level is h2osoi_liq_so. 
     ! A value of 1 means that the target soil moisture level is h2osoi_liq_sat
     real(r8) :: irrig_factor = 0.7_r8            

  end type irrigation_params_type


  type, public :: irrigation_type
     private
     ! Public data members
     ! Note: these should be treated as read-only by other modules
     real(r8), pointer, public :: qflx_irrig_patch(:) ! patch irrigation flux (mm H2O/s)
     real(r8), pointer, public :: qflx_irrig_col  (:) ! col irrigation flux (mm H2O/s)

     ! Private data members; set in initialization:
     type(irrigation_params_type) :: params
     integer :: dtime                ! land model time step (sec)
     integer :: irrig_nsteps_per_day ! number of time steps per day in which we irrigate
     real(r8), pointer :: relsat_so_patch(:,:) ! relative saturation at which smp = smpso (i.e., full stomatal opening) [patch, nlevgrnd]

     ! Private data members; time-varying:
     real(r8), pointer :: irrig_rate_patch         (:)  ! current irrigation rate [mm/s]
     integer , pointer :: n_irrig_steps_left_patch (:)  ! number of time steps for which we still need to irrigate today (if 0, ignore)
     
   contains
     ! Public routines
     ! NOTE(wjs, 2014-10-15) Workaround for pgi bug (pgi 14.7): Add an "Irrigation" prefix to some  generic routines like "Init"
     ! (without this workaround, pgi compilation fails in restFileMod)
     procedure, public :: Init => IrrigationInit
     procedure, public :: InitForTesting ! version of Init meant for unit testing
     procedure, public :: Restart
     procedure, public :: ApplyIrrigation
     procedure, public :: CalcIrrigationNeeded
     procedure, public :: Clean => IrrigationClean ! deallocate memory

     ! Public simply to support unit testing; should not be used from CLM code
     procedure, public, nopass :: IrrigationDeficit  ! compute the irrigation deficit for one layer of one point

     ! Private routines
     procedure, private :: InitAllocate => IrrigationInitAllocate
     procedure, private :: InitHistory => IrrigationInitHistory
     procedure, private :: InitCold => IrrigationInitCold
     procedure, private :: CalcIrrigNstepsPerDay   ! given dtime, calculate irrig_nsteps_per_day
     procedure, private :: PointNeedsCheckForIrrig ! whether a given point needs to be checked for irrigation now
  end type irrigation_type

  interface irrigation_params_type
     module procedure irrigation_params_constructor
  end interface irrigation_params_type
  
contains

  ! ========================================================================
  ! Constructors
  ! ========================================================================

  !-----------------------------------------------------------------------
  function irrigation_params_constructor(irrig_min_lai, irrig_btran_thresh, &
       irrig_start_time, irrig_length, irrig_factor) &
       result(this)
    !
    ! !DESCRIPTION:
    ! Create an irrigation_params instance
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(irrigation_params_type) :: this  ! function result
    real(r8), intent(in) :: irrig_min_lai
    real(r8), intent(in) :: irrig_btran_thresh
    integer , intent(in) :: irrig_start_time
    integer , intent(in) :: irrig_length
    real(r8), intent(in) :: irrig_factor
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'irrigation_params_constructor'
    !-----------------------------------------------------------------------
    
    this%irrig_min_lai = irrig_min_lai
    this%irrig_btran_thresh = irrig_btran_thresh
    this%irrig_start_time = irrig_start_time
    this%irrig_length = irrig_length
    this%irrig_factor = irrig_factor

  end function irrigation_params_constructor


  ! ========================================================================
  ! Infrastructure routines (initialization, restart, etc.)
  ! ========================================================================
  
  !------------------------------------------------------------------------
  subroutine IrrigationInit(this, bounds, soilstate_inst, soil_water_retention_curve)
    use SoilStateType , only : soilstate_type

    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve

    call this%InitAllocate(bounds)
    call this%InitHistory(bounds)
    call this%InitCold(bounds, soilstate_inst, soil_water_retention_curve)
  end subroutine IrrigationInit

  !-----------------------------------------------------------------------
  subroutine InitForTesting(this, bounds, params, dtime, relsat_so)
    !
    ! !DESCRIPTION:
    ! Does initialization needed for unit testing. Allows caller to prescribe values of
    ! some internal variables.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(irrigation_type)       , intent(inout) :: this
    type(bounds_type)            , intent(in)    :: bounds
    type(irrigation_params_type) , intent(in)    :: params
    integer                      , intent(in)    :: dtime ! model time step (sec)
    real(r8)                     , intent(in)    :: relsat_so( bounds%begp: , 1: ) ! relative saturation at which smp = smpso [patch, nlevgrnd]
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'InitForTesting'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL((ubound(relsat_so) == (/bounds%endp, nlevgrnd/)), errMsg(__FILE__, __LINE__))

    call this%InitAllocate(bounds)
    this%params = params
    this%dtime = dtime
    this%irrig_nsteps_per_day = this%CalcIrrigNstepsPerDay(dtime)
    this%relsat_so_patch(:,:) = relsat_so(:,:)

  end subroutine InitForTesting


  !-----------------------------------------------------------------------
  subroutine IrrigationInitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize irrigation data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    integer :: begc, endc

    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp
    begc = bounds%begc; endc= bounds%endc

    allocate(this%qflx_irrig_patch         (begp:endp))          ; this%qflx_irrig_patch         (:)   = nan
    allocate(this%qflx_irrig_col           (begc:endc))          ; this%qflx_irrig_col           (:)   = nan
    allocate(this%relsat_so_patch          (begp:endp,nlevgrnd)) ; this%relsat_so_patch(:,:) = nan
    allocate(this%irrig_rate_patch         (begp:endp))          ; this%irrig_rate_patch         (:)   = nan
    allocate(this%n_irrig_steps_left_patch (begp:endp))          ; this%n_irrig_steps_left_patch (:)   = 0

  end subroutine IrrigationInitAllocate

  !-----------------------------------------------------------------------
  subroutine IrrigationInitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize irrigation history fields
    !
    ! !USES:
    use histFileMod  , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp= bounds%endp

    this%qflx_irrig_patch(begp:endp) = spval
    call hist_addfld1d (fname='QIRRIG', units='mm/s', &
         avgflag='A', long_name='water added through irrigation', &
         ptr_patch=this%qflx_irrig_patch)

  end subroutine IrrigationInitHistory

  !-----------------------------------------------------------------------
  subroutine IrrigationInitCold(this, bounds, soilstate_inst, soil_water_retention_curve)
    !
    ! !DESCRIPTION:
    ! Do cold-start initialization for irrigation data structure
    !
    ! !USES:
    use pftconMod     , only : pftcon
    use SoilStateType , only : soilstate_type
    use pftconMod     , only : noveg
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    type(soilstate_type)   , intent(in)    :: soilstate_inst
    class(soil_water_retention_curve_type), intent(in) :: soil_water_retention_curve
    !
    ! !LOCAL VARIABLES:
    integer :: p ! patch index
    integer :: c ! col index
    integer :: j ! level index

    character(len=*), parameter :: subname = 'InitCold'
    !-----------------------------------------------------------------------

    associate( &
         sucsat => soilstate_inst%sucsat_col , & ! Input:  [real(r8) (:,:) ]  minimum soil suction (mm) (constant)
         bsw    => soilstate_inst%bsw_col    , & ! Input:  [real(r8) (:,:) ]  Clapp and Hornberger "b"  (constant)
         smpso  => pftcon%smpso                & ! Input:  soil water potential at full stomatal opening (mm) (constant)
         )
    
      do j = 1, nlevgrnd
         do p = bounds%begp, bounds%endp
            c = patch%column(p)
            if (patch%itype(p) /= noveg) then
               call soil_water_retention_curve%soil_suction_inverse( &
                    smp_target = smpso(patch%itype(p)), &
                    smpsat = sucsat(c,j), &
                    bsw = bsw(c,j), &
                    s_target = this%relsat_so_patch(p,j))
            end if
         end do
      end do

      this%dtime = get_step_size()
      this%irrig_nsteps_per_day = this%CalcIrrigNstepsPerDay(this%dtime)

    end associate

  end subroutine IrrigationInitCold

  !-----------------------------------------------------------------------
  pure function CalcIrrigNstepsPerDay(this, dtime) result(irrig_nsteps_per_day)
    !
    ! !DESCRIPTION:
    ! Given dtime (sec), determine number of irrigation steps per day
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    integer :: irrig_nsteps_per_day  ! function result
    class(irrigation_type) , intent(in) :: this
    integer                , intent(in) :: dtime ! model time step (sec)
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'CalcIrrigNstepsPerDay'
    !-----------------------------------------------------------------------
    
    irrig_nsteps_per_day = ((this%params%irrig_length + (dtime - 1))/dtime)  ! round up

  end function CalcIrrigNstepsPerDay



  !-----------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Handle restart of irrigation variables
    !
    ! !USES:
    use ncdio_pio        , only : file_desc_t, ncd_inqvdlen, ncd_double, ncd_int
    use restUtilMod
    !
    ! !ARGUMENTS:
    class(irrigation_type) :: this
    type(bounds_type), intent(in)    :: bounds 
    type(file_desc_t), intent(inout) :: ncid   ! netcdf id
    character(len=*) , intent(in)    :: flag   ! 'read', 'write' or 'define'
    !
    ! !LOCAL VARIABLES:
    logical :: do_io
    integer :: dimlen       ! dimension length
    integer :: nump_global  ! total number of patchs, globally
    integer :: err_code     ! error code
    logical :: readvar      ! determine if variable is on initial file

    character(len=*), parameter :: subname = 'Restart'
    !-----------------------------------------------------------------------
    
    ! Get expected total number of points, for later error checks
    call get_proc_global(np=nump_global)

    do_io = .true.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'n_irrig_steps_left', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
          readvar = .false.
       end if
    else if (flag == 'define' .or. do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='n_irrig_steps_left', xtype=ncd_int,  &
            dim1name='pft', &
            long_name='number of irrigation time steps left', units='#', &
            interpinic_flag='interp', readvar=readvar, data=this%n_irrig_steps_left_patch)
       if (flag=='read' .and. .not. readvar) then
          this%n_irrig_steps_left_patch = 0
       end if
    end if

    do_io = .true.
    if (flag == 'read') then
       ! On a read, confirm that this variable has the expected size; if not, don't read
       ! it (instead give it a default value). This is needed to support older initial
       ! conditions for which this variable had a different size.
       call ncd_inqvdlen(ncid, 'irrig_rate', 1, dimlen, err_code)
       if (dimlen /= nump_global) then
          do_io = .false.
          readvar = .false.
       end if
    else if (flag == 'define' .or. do_io) then
       call restartvar(ncid=ncid, flag=flag, varname='irrig_rate', xtype=ncd_double,  &
            dim1name='pft', &
            long_name='irrigation rate', units='mm/s', &
            interpinic_flag='interp', readvar=readvar, data=this%irrig_rate_patch)
       if (flag=='read' .and. .not. readvar) then
          this%irrig_rate_patch = 0.0_r8
       end if
    end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine IrrigationClean(this)
    !
    ! !DESCRIPTION:
    ! Deallocate memory
    !
    ! !ARGUMENTS:
    class(irrigation_type), intent(inout) :: this
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'Clean'
    !-----------------------------------------------------------------------
    
    deallocate(this%qflx_irrig_patch)
    deallocate(this%qflx_irrig_col)
    deallocate(this%relsat_so_patch)
    deallocate(this%irrig_rate_patch)
    deallocate(this%n_irrig_steps_left_patch)

  end subroutine IrrigationClean


  ! ========================================================================
  ! Science routines
  ! ========================================================================
  
  !-----------------------------------------------------------------------
  subroutine ApplyIrrigation(this, bounds)
    !
    ! !DESCRIPTION:
    ! Apply the irrigation computed by CalcIrrigationNeeded to qflx_irrig.
    !
    ! Should be called once, AND ONLY ONCE, per time step. After this is called, you may
    ! access qflx_irrig_patch or qflx_irrig_col.
    !
    ! !USES:
    use subgridAveMod, only : p2c
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: p  ! patch index
    
    character(len=*), parameter :: subname = 'ApplyIrrigation'
    !-----------------------------------------------------------------------
    
    ! This should be called exactly once per time step, so that this counter decrease
    ! works correctly.

    do p = bounds%begp, bounds%endp
       if (this%n_irrig_steps_left_patch(p) > 0) then
          this%qflx_irrig_patch(p)         = this%irrig_rate_patch(p)
          this%n_irrig_steps_left_patch(p) = this%n_irrig_steps_left_patch(p) - 1
       else
          this%qflx_irrig_patch(p) = 0._r8
       end if
    end do

    call p2c (bounds = bounds, &
         parr = this%qflx_irrig_patch(bounds%begp:bounds%endp), &
         carr = this%qflx_irrig_col(bounds%begc:bounds%endc), &
         p2c_scale_type = 'unity')

  end subroutine ApplyIrrigation


  !-----------------------------------------------------------------------
  subroutine CalcIrrigationNeeded(this, bounds, num_exposedvegp, filter_exposedvegp, &
       time_prev, elai, btran, rootfr, t_soisno, eff_porosity, h2osoi_liq)
    !
    ! !DESCRIPTION:
    ! Calculate whether and how much irrigation is needed for each column. However, this
    ! does NOT actually set the irrigation flux.
    !
    ! !USES:
    use shr_const_mod      , only : SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    class(irrigation_type) , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds

    ! time of day (in seconds since 0Z) at start of timestep
    integer, intent(in) :: time_prev

    ! number of points in filter_exposedvegp
    integer, intent(in) :: num_exposedvegp

    ! patch filter for non-snow-covered veg
    integer, intent(in) :: filter_exposedvegp(:)

    ! one-sided leaf area index with burying by snow [patch]
    real(r8), intent(in) :: elai( bounds%begp: )

    ! transpiration wetness factor (0 to 1) [patch]
    real(r8), intent(in) :: btran( bounds%begp: )

    ! fraction of roots in each soil later [patch]
    real(r8), intent(in) :: rootfr( bounds%begp: , 1: )

    ! col soil temperature (K) [col, nlevgrnd] (note that this does NOT contain the snow levels)
    real(r8), intent(in) :: t_soisno( bounds%begc: , 1: )

    ! effective porosity (0 to 1) [col, nlevgrnd]
    real(r8), intent(in) :: eff_porosity( bounds%begc: , 1: )

    ! column liquid water (kg/m2) [col, nlevgrnd] (note that this does NOT contain the snow levels)
    real(r8), intent(in) :: h2osoi_liq( bounds%begc: , 1: )
    
    !
    ! !LOCAL VARIABLES:
    integer :: f    ! filter index
    integer :: p    ! patch index
    integer :: c    ! column index
    integer :: g    ! gridcell index
    integer :: j    ! level

    ! difference between desired soil moisture level for this layer and current soil moisture level [kg/m2]
    real(r8) :: deficit

    ! where do we need to check soil moisture to see if we need to irrigate?
    logical  :: check_for_irrig(bounds%begp:bounds%endp)

    ! set to true if we have encountered a frozen soil layer
    logical  :: frozen_soil(bounds%begp:bounds%endp)
    
    character(len=*), parameter :: subname = 'CalcIrrigationNeeded'
    !-----------------------------------------------------------------------
    
    ! Enforce expected array sizes
    SHR_ASSERT_ALL((ubound(elai) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(btran) == (/bounds%endp/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(rootfr) == (/bounds%endp, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(t_soisno) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(eff_porosity) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))
    SHR_ASSERT_ALL((ubound(h2osoi_liq) == (/bounds%endc, nlevgrnd/)), errMsg(__FILE__, __LINE__))


    ! Determine if irrigation is needed (over irrigated soil columns)
    
    ! First, determine in what grid cells we need to bother 'measuring' soil water, to see if we need irrigation
    ! Also set n_irrig_steps_left for these grid cells
    ! n_irrig_steps_left(p) > 0 is ok even if irrig_rate(p) ends up = 0
    ! in this case, we'll irrigate by 0 for the given number of time steps
    
    do f = 1, num_exposedvegp
       p = filter_exposedvegp(f)
       g = patch%gridcell(p)
       check_for_irrig(p) = this%PointNeedsCheckForIrrig( &
            pft_type=patch%itype(p), elai=elai(p), btran=btran(p), &
            time_prev=time_prev, londeg=grc%londeg(g))

       if (check_for_irrig(p)) then
          this%n_irrig_steps_left_patch(p) = this%irrig_nsteps_per_day
          this%irrig_rate_patch(p)         = 0._r8  ! reset; we'll add to this later
       end if
    end do

    ! Now 'measure' soil water for the grid cells identified above and see if the soil is
    ! dry enough to warrant irrigation
    ! (Note: frozen_soil could probably be a column-level variable, but that would be
    ! slightly less robust to potential future modifications)
    frozen_soil(bounds%begp : bounds%endp) = .false.
    do j = 1,nlevgrnd
       do f = 1, num_exposedvegp
          p = filter_exposedvegp(f)
          c = patch%column(p)
          if (check_for_irrig(p) .and. .not. frozen_soil(p)) then
             ! if level L was frozen, then we don't look at any levels below L
             if (t_soisno(c,j) <= SHR_CONST_TKFRZ) then
                frozen_soil(p) = .true.
             else if (rootfr(p,j) > 0._r8) then
                deficit = this%IrrigationDeficit( &
                     relsat_so  = this%relsat_so_patch(p,j), &
                     h2osoi_liq = h2osoi_liq(c,j), &
                     eff_porosity = eff_porosity(c,j), &
                     dz = col%dz(c,j), &
                     irrig_factor = this%params%irrig_factor)

                ! Add deficit to irrig_rate, converting units from mm to mm/sec
                this%irrig_rate_patch(p)  = this%irrig_rate_patch(p) + &
                     deficit/(this%dtime*this%irrig_nsteps_per_day)

             end if  ! else if (rootfr(p,j) > 0)
          end if     ! if (check_for_irrig(p) .and. .not. frozen_soil(p))
       end do        ! do f
    end do           ! do j


  end subroutine CalcIrrigationNeeded

  !-----------------------------------------------------------------------
  pure function PointNeedsCheckForIrrig(this, pft_type, elai, btran, time_prev, londeg) &
       result(check_for_irrig)
    !
    ! !DESCRIPTION:
    ! Determine whether a given patch needs to be checked for irrigation now.
    !
    ! !USES:
    use pftconMod, only : pftcon
    !
    ! !ARGUMENTS:
    logical :: check_for_irrig  ! function result
    class(irrigation_type), intent(in) :: this
    integer , intent(in) :: pft_type  ! type of pft in this patch
    real(r8), intent(in) :: elai      ! one-sided leaf area index with burying by snow
    real(r8), intent(in) :: btran     ! transpiration wetness factor (0 to 1)
    integer , intent(in) :: time_prev ! time of day (in seconds since 0Z) at start of timestep
    real(r8), intent(in) :: londeg    ! longitude (degrees)
    !
    ! !LOCAL VARIABLES:
    ! local time at start of time step (seconds after solar midnight)
    integer  :: local_time

    ! number of seconds since the prescribed irrigation start time
    integer  :: seconds_since_irrig_start_time

    character(len=*), parameter :: subname = 'PointNeedsCheckForIrrig'
    !-----------------------------------------------------------------------
    
    if (pftcon%irrigated(pft_type) == 1._r8 .and. &
         elai > this%params%irrig_min_lai .and. &
         btran < this%params%irrig_btran_thresh) then
       ! see if it's the right time of day to start irrigating:
       local_time = modulo(time_prev + nint(londeg/degpsec), isecspday)
       seconds_since_irrig_start_time = modulo(local_time - this%params%irrig_start_time, isecspday)
       if (seconds_since_irrig_start_time < this%dtime) then
          check_for_irrig         = .true.
       else
          check_for_irrig    = .false.
       end if
    else
       check_for_irrig       = .false.
    end if

  end function PointNeedsCheckForIrrig



  !-----------------------------------------------------------------------
  pure function IrrigationDeficit(relsat_so, h2osoi_liq, eff_porosity, dz, irrig_factor) &
       result(deficit)
    !
    ! !DESCRIPTION:
    ! Compute irrigation deficit for a given soil layer at a given point. This is the
    ! difference between the desired soil moisture level for this layer and the current
    ! soil moisture level. [kg/m2]
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    real(r8) :: deficit  ! function result
    real(r8), intent(in) :: relsat_so    ! relative saturation at which smp = smpso
    real(r8), intent(in) :: h2osoi_liq   ! current liquid water in layer (kg/m2)
    real(r8), intent(in) :: eff_porosity ! effective porosity (0 to 1)
    real(r8), intent(in) :: dz           ! level thickness (m)
    real(r8), intent(in) :: irrig_factor ! factor determining the target soil moisture level for irrigation (0 to 1)
    !
    ! !LOCAL VARIABLES:

    ! partial volume of liquid water in layer for which smp_node = smpso [0 to 1]
    real(r8) :: vol_liq_so

    ! liquid water corresponding to vol_liq_so for this layer [kg/m2]
    real(r8) :: h2osoi_liq_so

    ! liquid water corresponding to eff_porosity for this layer [kg/m2]
    real(r8) :: h2osoi_liq_sat

    character(len=*), parameter :: subname = 'IrrigationDeficit'
    !-----------------------------------------------------------------------
    
    vol_liq_so     = eff_porosity * relsat_so
    h2osoi_liq_so  = vol_liq_so * denh2o * dz
    h2osoi_liq_sat = eff_porosity * denh2o * dz
    deficit        = max((h2osoi_liq_so + irrig_factor*(h2osoi_liq_sat - h2osoi_liq_so)) - h2osoi_liq, 0._r8)

  end function IrrigationDeficit

end module IrrigationMod
