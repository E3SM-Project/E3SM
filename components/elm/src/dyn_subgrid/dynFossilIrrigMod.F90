module dynFossilIrrigMod

#include "shr_assert.h"

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Time-varying share of irrigation demand that is met from non-renewable (fossil)
  ! groundwater, an external source outside the ELM and MOSART water stores.
  !
  ! The share is read from the file named by namelist fossil_irrig_file: variable
  ! F_NONREN(time, lsmlat, lsmlon), fraction 0-1, with the usual YEAR(time) variable of
  ! the dynamic-landuse files. The value of the current model year is used for the whole
  ! year (no interpolation in time); the first/last year is held outside the series.
  ! Without a file the share is zero everywhere and the model is unchanged.
  !
  ! Used by two-way irrigation (tw_irr): CanopyHydrologyMod applies
  ! f_fossil * demand as always-met irrigation, HydrologyDrainageMod sends only
  ! (1 - f_fossil) * f_surf * demand to MOSART.
  !
  ! !USES:
  use shr_kind_mod          , only : r8 => shr_kind_r8
  use shr_log_mod           , only : errMsg => shr_log_errMsg
  use decompMod             , only : bounds_type, BOUNDS_LEVEL_PROC
  use abortutils            , only : endrun
  use dynFileMod            , only : dyn_file_type
  use dynVarTimeUninterpMod , only : dyn_var_time_uninterp_type
  use elm_varcon            , only : grlnd
  use elm_varctl            , only : iulog, fossil_irrig_file
  use elm_varsur            , only : f_fossil
  use spmdMod               , only : masterproc
  !
  implicit none
  private
  save
  public :: dynFossilIrrig_init    ! open the file and get the first data
  public :: dynFossilIrrig_interp  ! get the share of the current model year
  !
  type(dyn_file_type), target      :: dynFossilIrrig_file
  type(dyn_var_time_uninterp_type) :: fossil_var
  character(len=*), parameter      :: fossil_varname = 'F_NONREN'
  logical                          :: fossil_active = .false.
  !---------------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine dynFossilIrrig_init(bounds)
    !
    ! !DESCRIPTION:
    ! Called once during initialization. f_fossil is allocated (and zero) already.
    !
    ! !USES:
    use dynTimeInfoMod, only : YEAR_POSITION_START_OF_TIMESTEP
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer :: fossil_shape(1)
    character(len=*), parameter :: subname = 'dynFossilIrrig_init'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    f_fossil(bounds%begg:bounds%endg) = 0._r8
    fossil_active = (len_trim(fossil_irrig_file) > 0)
    if (.not. fossil_active) return

    if (masterproc) then
       write(iulog,*) subname, ': non-renewable irrigation share from ', trim(fossil_irrig_file)
    end if

    dynFossilIrrig_file = dyn_file_type(fossil_irrig_file, YEAR_POSITION_START_OF_TIMESTEP)
    fossil_shape(1) = bounds%endg - bounds%begg + 1
    fossil_var = dyn_var_time_uninterp_type( &
         dyn_file=dynFossilIrrig_file, varname=fossil_varname, &
         dim1name=grlnd, conversion_factor=1.0_r8, &
         do_check_sums_equal_1=.false., data_shape=fossil_shape)

    call dynFossilIrrig_interp(bounds)

  end subroutine dynFossilIrrig_init

  !-----------------------------------------------------------------------
  subroutine dynFossilIrrig_interp(bounds)
    !
    ! !DESCRIPTION:
    ! Called every time step, outside threaded regions (may read the file).
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds  ! proc-level bounds
    !
    ! !LOCAL VARIABLES:
    integer               :: g
    real(r8), allocatable :: this_data(:)
    character(len=*), parameter :: subname = 'dynFossilIrrig_interp'
    !-----------------------------------------------------------------------

    SHR_ASSERT_ALL(bounds%level == BOUNDS_LEVEL_PROC, subname // ': argument must be PROC-level bounds')

    if (.not. fossil_active) return

    call dynFossilIrrig_file%time_info%set_current_year()

    allocate(this_data(bounds%begg:bounds%endg))
    call fossil_var%get_current_data(this_data)
    do g = bounds%begg, bounds%endg
       if (this_data(g) < -1.e-6_r8 .or. this_data(g) > 1._r8 + 1.e-6_r8) then
          write(iulog,*) subname, ': F_NONREN outside 0-1 at g = ', g, this_data(g)
          call endrun(msg=subname//' ERROR: F_NONREN outside 0-1'//errMsg(__FILE__, __LINE__))
       end if
       f_fossil(g) = min(1._r8, max(0._r8, this_data(g)))
    end do
    deallocate(this_data)

  end subroutine dynFossilIrrig_interp

end module dynFossilIrrigMod
