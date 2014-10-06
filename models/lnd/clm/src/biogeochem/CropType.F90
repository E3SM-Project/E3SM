module CropType
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module containing variables needed for the crop model
  !
  ! TODO(wjs, 2014-08-05) Move more crop-specific variables into here
  !
  ! !USES:
  use shr_kind_mod        , only : r8 => shr_kind_r8
  use shr_log_mod         , only : errMsg => shr_log_errMsg
  use abortutils          , only : endrun
  use decompMod           , only : bounds_type
  use clm_varcon          , only : spval
  use clm_varpar          , only : crop_prog
  use clm_varctl          , only : iulog
  !
  ! !PUBLIC TYPES:
  implicit none
  private
  save
  !
  ! !PUBLIC DATA TYPES:
  !
  ! Crop state variables structure
  type, public :: crop_type
     real(r8), pointer :: gddplant_patch          (:)   ! patch accum gdd past planting date for crop       (ddays)
     real(r8), pointer :: gddtsoi_patch           (:)   ! patch growing degree-days from planting (top two soil layers) (ddays)

   contains
     procedure, public  :: Init
     procedure, private :: InitAllocate 
     procedure, private :: InitHistory
     procedure, public  :: InitAccBuffer
     procedure, public  :: InitAccVars
     procedure, public  :: UpdateAccVars

  end type crop_type

  !------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine Init(this, bounds)
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'Init'
    !-----------------------------------------------------------------------
    
    call this%InitAllocate(bounds)

    if (crop_prog) then
       call this%InitHistory(bounds)
    end if

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitAllocate'
    !-----------------------------------------------------------------------

    begp = bounds%begp; endp = bounds%endp

    allocate(this%gddplant_patch           (begp:endp))                      ; this%gddplant_patch           (:)   = spval
    allocate(this%gddtsoi_patch            (begp:endp))                      ; this%gddtsoi_patch            (:)   = spval
    
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    !
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begp, endp
    
    character(len=*), parameter :: subname = 'InitHistory'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp; endp = bounds%endp

    this%gddplant_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDDPLANT', units='ddays', &
         avgflag='A', long_name='Accumulated growing degree days past planting date for crop', &
         ptr_patch=this%gddplant_patch, default='inactive')

    this%gddtsoi_patch(begp:endp) = spval
    call hist_addfld1d (fname='GDDTSOI', units='ddays', &
         avgflag='A', long_name='Growing degree-days from planting (top two soil layers)', &
         ptr_patch=this%gddtsoi_patch, default='inactive')

  end subroutine InitHistory


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
    ! Routine [updateAccVars] does the actual accumulation for a given field.
    ! Fields are accumulated by calls to subroutine [update_accum_field]. 
    ! To accumulate a field, it must first be defined in subroutine [initAccVars] 
    ! and then accumulated by calls to [updateAccVars].
    !
    ! Should only be called if crop_prog is true
    !
    ! !USES 
    use accumulMod       , only : init_accum_field
    !
    ! !ARGUMENTS:
    class(crop_type) , intent(in) :: this
    type(bounds_type), intent(in) :: bounds  

    !
    ! !LOCAL VARIABLES:
    integer, parameter :: not_used = huge(1)

    !---------------------------------------------------------------------

    call init_accum_field (name='GDDPLANT', units='K', &
         desc='growing degree-days from planting', accum_type='runaccum', accum_period=not_used,  &
         subgrid_type='pft', numlev=1, init_value=0._r8)

    call init_accum_field (name='GDDTSOI', units='K', &
         desc='growing degree-days from planting (top two soil layers)', accum_type='runaccum', accum_period=not_used,  &
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
    ! !USES:
    use accumulMod       , only : extract_accum_field
    use clm_time_manager , only : get_nstep
    !
    ! !ARGUMENTS:
    class(crop_type),  intent(inout) :: this
    type(bounds_type), intent(in)    :: bounds
    !
    ! !LOCAL VARIABLES:
    integer  :: begp, endp
    integer  :: nstep
    integer  :: ier
    real(r8), pointer :: rbufslp(:)  ! temporary
    
    character(len=*), parameter :: subname = 'InitAccVars'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp; endp = bounds%endp

    ! Allocate needed dynamic memory for single level pft field
    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)' in '
       call endrun(msg=" allocation error for rbufslp"//&
            errMsg(__FILE__, __LINE__))
    endif

    nstep = get_nstep()

    call extract_accum_field ('GDDPLANT', rbufslp, nstep) 
    this%gddplant_patch(begp:endp) = rbufslp(begp:endp)

    call extract_accum_field ('GDDTSOI', rbufslp, nstep) 
    this%gddtsoi_patch(begp:endp)  = rbufslp(begp:endp)

    deallocate(rbufslp)

  end subroutine InitAccVars

  !-----------------------------------------------------------------------
  subroutine UpdateAccVars(this, bounds, temperature_vars, cnstate_vars)
    !
    ! !DESCRIPTION:
    ! Update accumulated variables. Should be called every time step.
    !
    ! Should only be called if crop_prog is true.
    !
    ! !USES:
    use accumulMod       , only : update_accum_field, extract_accum_field, accumResetVal
    use shr_const_mod    , only : SHR_CONST_CDAY, SHR_CONST_TKFRZ
    use clm_time_manager , only : get_step_size, get_nstep
    use pftvarcon        , only : nwcereal, nwcerealirrig, mxtmp, baset
    use TemperatureType  , only : temperature_type
    use CNStateType      , only : cnstate_type
    use PatchType        , only : pft                
    use ColumnType       , only : col
    !
    ! !ARGUMENTS:
    class(crop_type)       , intent(inout) :: this
    type(bounds_type)      , intent(in)    :: bounds
    type(temperature_type) , intent(in)    :: temperature_vars
    type(cnstate_type)     , intent(in) :: cnstate_vars
    !
    ! !LOCAL VARIABLES:
    integer :: p,c   ! indices
    integer :: ivt   ! vegetation type
    integer :: dtime ! timestep size [seconds]
    integer :: nstep ! timestep number
    integer :: ier   ! error status
    integer :: begp, endp
    real(r8), pointer :: rbufslp(:)      ! temporary single level - pft level
    
    character(len=*), parameter :: subname = 'UpdateAccVars'
    !-----------------------------------------------------------------------
    
    begp = bounds%begp; endp = bounds%endp

    dtime = get_step_size()
    nstep = get_nstep()

    ! Allocate needed dynamic memory for single level pft field

    allocate(rbufslp(begp:endp), stat=ier)
    if (ier/=0) then
       write(iulog,*)'update_accum_hist allocation error for rbuf1dp'
       call endrun(msg=errMsg(__FILE__, __LINE__))
    endif

    ! Accumulate and extract GDDPLANT
    
    do p = begp,endp
       if (cnstate_vars%croplive_patch(p)) then ! relative to planting date
          ivt = pft%itype(p)
          rbufslp(p) = max(0._r8, min(mxtmp(ivt), &
               temperature_vars%t_ref2m_patch(p)-(SHR_CONST_TKFRZ + baset(ivt)))) &
               * dtime/SHR_CONST_CDAY
          if (ivt == nwcereal .or. ivt == nwcerealirrig) then
             rbufslp(p) = rbufslp(p)*cnstate_vars%vf_patch(p)
          end if
       else
          rbufslp(p) = accumResetVal
       end if
    end do
    call update_accum_field  ('GDDPLANT', rbufslp, nstep)
    call extract_accum_field ('GDDPLANT', this%gddplant_patch, nstep)

    ! Accumulate and extract GDDTSOI
    ! In agroibis this variable is calculated
    ! to 0.05 m, so here we use the top two soil layers

    do p = begp,endp
       if (cnstate_vars%croplive_patch(p)) then ! relative to planting date
          ivt = pft%itype(p)
          c = pft%column(p)
          rbufslp(p) = max(0._r8, min(mxtmp(ivt), &
               ((temperature_vars%t_soisno_col(c,1)*col%dz(c,1) + &
               temperature_vars%t_soisno_col(c,2)*col%dz(c,2))/(col%dz(c,1)+col%dz(c,2))) - &
               (SHR_CONST_TKFRZ + baset(ivt)))) * dtime/SHR_CONST_CDAY
          if (ivt == nwcereal .or. ivt == nwcerealirrig) then
             rbufslp(p) = rbufslp(p)*cnstate_vars%vf_patch(p)
          end if
       else
          rbufslp(p) = accumResetVal
       end if
    end do
    call update_accum_field  ('GDDTSOI', rbufslp, nstep)
    call extract_accum_field ('GDDTSOI', this%gddtsoi_patch, nstep)

    deallocate(rbufslp)

  end subroutine UpdateAccVars

end module CropType

