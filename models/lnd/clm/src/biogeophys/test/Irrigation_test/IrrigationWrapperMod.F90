module IrrigationWrapperMod

  ! This module provides a wrapper around the irrigation class, for the sake of testing.
  !
  ! It defines a type that holds the parameters needed to set up a test of irrigation.
  !
  ! In addition, it provides a wrapper to the public interface of Irrigation, which calls
  ! the Irrigation routines with the appropriate arguments. This way, if the argument
  ! list changes, we only need to change the call in one place, rather than in every test.
  !
  ! Finally, it provides a routine that does the setup needed for most tests of
  ! irrigation, and a complementary routine to do the teardown.
  
  use shr_kind_mod  , only : r8 => shr_kind_r8
  use IrrigationMod , only : irrigation_type, irrigation_params_type
  use clm_varpar    , only : nlevgrnd
  use PatchType     , only : patch
  use ColumnType    , only : col
  use GridcellType  , only : grc
  use unittestSubgridMod

  implicit none
  save
  private

  type, public :: irrigation_inputs_type
     ! Irrigation parameters
     type(irrigation_params_type) :: irrigation_params
     
     ! State variables
     real(r8), allocatable :: elai(:)
     real(r8), allocatable :: btran(:)
     real(r8), allocatable :: rootfr(:,:)
     real(r8), allocatable :: t_soisno(:,:)
     real(r8), allocatable :: eff_porosity(:,:)
     real(r8), allocatable :: h2osoi_liq(:,:)
     real(r8), allocatable :: relsat_so(:,:)

     ! Previous model time
     integer :: time_prev
   contains
     ! Computes irrigation deficit for every patch and level
     procedure :: computeDeficits

     ! Wrapper that calls both CalcIrrigationNeeded and ApplyIrrigation
     procedure :: calculateAndApplyIrrigation
  end type irrigation_inputs_type

  integer , parameter, public :: dtime = 1800  ! model time step, seconds

  ! Public routines:
  public :: setupIrrigation         ! Do the setup needed for most tests
  public :: teardownIrrigation      ! Teardown stuff set up by setupIrrigation
  
  ! Private routines:

  ! These routines setup and teardown the external environment used by Irrigation - i.e.,
  ! things accessed via 'use' statements
  private :: setupEnvironment
  private :: teardownEnvironment

  interface irrigation_inputs_type
     module procedure constructor
  end interface irrigation_inputs_type
 
contains

  !-----------------------------------------------------------------------
  type(irrigation_inputs_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an irrigation_inputs_type object.
    !
    ! Values are set up such that there is some irrigation deficit everywhere, and
    ! irrigation would start in the following call to CalcIrrigationNeeded (followed by
    ! ApplyIrrigation). Values are set the same for every patch/column, and are the same
    ! at every level EXCEPT for relsat_so, which varies linearly by level and patch number.
    !
    ! Assumes that nlevgrnd has been set, and that all necessary subgrid setup has been
    ! completed.
    !
    ! !ARGUMENTS:
    !
    ! !LOCAL VARIABLES:
    integer :: p,j
    !
    !-----------------------------------------------------------------------

    ! Set parameters
    constructor%irrigation_params = irrigation_params_type( &
         irrig_min_lai = 0.0_r8, &
         irrig_btran_thresh = 0.99_r8, &
         irrig_start_time = 21600, &
         irrig_length = 14400, &
         irrig_factor = 0.7_r8)
         

    ! ------------------------------------------------------------------------
    ! Set state variables
    ! ------------------------------------------------------------------------

    allocate(constructor%elai(bounds%begp:bounds%endp), source=10._r8)
    allocate(constructor%btran(bounds%begp:bounds%endp), source=0._r8)
    allocate(constructor%rootfr(bounds%begp:bounds%endp, nlevgrnd), source=1._r8/nlevgrnd)
    allocate(constructor%t_soisno(bounds%begc:bounds%endc, nlevgrnd), source=1000._r8)
    allocate(constructor%eff_porosity(bounds%begc:bounds%endc, nlevgrnd), source=1._r8)
    allocate(constructor%h2osoi_liq(bounds%begc:bounds%endc, nlevgrnd), source=0._r8)
    allocate(constructor%relsat_so(bounds%begp:bounds%endp, nlevgrnd))

    do j = 1, nlevgrnd
       do p = bounds%begp, bounds%endp
          constructor%relsat_so(p,j) = 0.1_r8 * j * (p - bounds%begp + 1)
       end do
    end do

    ! Set time_prev to the irrig_start_time minus 1 hour (since we're using a longitude
    ! about 1 hour east of 0Z)
    constructor%time_prev = constructor%irrigation_params%irrig_start_time - 3600

  end function constructor

  !-----------------------------------------------------------------------
  subroutine computeDeficits(this, irrigation, deficits)
    !
    ! !DESCRIPTION:
    ! Computes irrigation deficit for each patch and layer.
    !
    ! Allocates the 'deficits' variable, and gives it a lower bound of bounds%begp
    !
    ! The motivation for this function is: For most of the irrigation tests, we assume
    ! that the IrrigationDeficit function is working correctly, and we want to test the
    ! code that builds on top of these computed deficits. By having this function, we can
    ! avoid having to hard-code the deficits in each test.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(irrigation_inputs_type), intent(in) :: this
    class(irrigation_type), intent(in) :: irrigation
    real(r8), allocatable, intent(out) :: deficits(:,:)
    !
    ! !LOCAL VARIABLES:
    integer :: p, c, j

    character(len=*), parameter :: subname = 'computeDeficits'
    !-----------------------------------------------------------------------
    
    allocate(deficits(bounds%begp:bounds%endp, nlevgrnd))
    do j = 1, nlevgrnd
       do p = bounds%begp, bounds%endp
          c = patch%column(p)
          deficits(p,j) = irrigation%IrrigationDeficit(&
               relsat_so = this%relsat_so(p,j), &
               h2osoi_liq = this%h2osoi_liq(c,j), &
               eff_porosity = this%eff_porosity(c,j), &
               dz = col%dz(c,j), &
               irrig_factor = this%irrigation_params%irrig_factor)
       end do
    end do

  end subroutine computeDeficits

  !-----------------------------------------------------------------------
  subroutine calculateAndApplyIrrigation(this, irrigation, numf, filter)
    !
    ! !DESCRIPTION:
    ! Call CalculateIrrigationNeeded with the given irrigation parameters. Then call
    ! ApplyIrrigation.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(irrigation_inputs_type), intent(in) :: this
    class(irrigation_type), intent(inout) :: irrigation
    integer :: numf      ! number of points in filter
    integer :: filter(:) ! filter over which we run irrigation
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'calculateAndApplyIrrigation'
    !-----------------------------------------------------------------------

    call irrigation%CalcIrrigationNeeded(&
         bounds=bounds, &
         num_exposedvegp = numf, &
         filter_exposedvegp = filter, &
         time_prev = this%time_prev, &
         elai = this%elai, &
         btran = this%btran, &
         rootfr = this%rootfr, &
         t_soisno = this%t_soisno, &
         eff_porosity = this%eff_porosity, &
         h2osoi_liq = this%h2osoi_liq)
    
    call irrigation%ApplyIrrigation(bounds)

  end subroutine calculateAndApplyIrrigation


  ! ========================================================================
  ! Procedures not tied to irrigation_inputs_type, but included in this same module for
  ! convenience.
  ! ========================================================================

  !-----------------------------------------------------------------------
  subroutine setupIrrigation(irrigation_inputs, irrigation, maxpft)
    !
    ! !DESCRIPTION:
    ! Do the setup needed for most tests.
    !
    ! Before calling this, you must set up the subgrid structure.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(irrigation_inputs_type), intent(out) :: irrigation_inputs
    type(irrigation_type), intent(out) :: irrigation
    integer, intent(in) :: maxpft ! max pft type
    !
    ! !LOCAL VARIABLES:
    !-----------------------------------------------------------------------
    
    irrigation_inputs = irrigation_inputs_type()
    call setupEnvironment(maxpft=maxpft)
    call irrigation%InitForTesting(bounds, irrigation_inputs%irrigation_params, &
         dtime, irrigation_inputs%relsat_so)

  end subroutine setupIrrigation

  !-----------------------------------------------------------------------
  subroutine teardownIrrigation(irrigation_inputs, irrigation)
    !
    ! !DESCRIPTION:
    ! Teardown stuff set up by setupIrrigation
    !
    ! Note: nothing is done with irrigation_inputs, but it is included in the argument
    ! list for symmetry with the setup routine, in case anything needs to be done with it
    ! in the future.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    type(irrigation_inputs_type), intent(inout) :: irrigation_inputs
    type(irrigation_type), intent(inout) :: irrigation
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'teardownIrrigation'
    !-----------------------------------------------------------------------
    
    call irrigation%Clean()
    call teardownEnvironment()

  end subroutine teardownIrrigation



  !-----------------------------------------------------------------------
  subroutine setupEnvironment(maxpft)
    !
    ! !DESCRIPTION:
    ! Sets up the external environment used by Irrigation - i.e., things accessed via
    ! 'use' statements.
    !
    ! Assumes nlevgrnd has been set, and that all necessary subgrid setup has been
    ! completed.
    !
    ! !USES:
    use pftconMod , only : pftcon
    use clm_varpar, only : mxpft
    !
    ! !ARGUMENTS:
    integer, intent(in) :: maxpft  ! max pft type that needs to be supported
    !
    !-----------------------------------------------------------------------

    allocate(pftcon%irrigated(0:mxpft), source=1.0_r8)

    col%dz(:,1:nlevgrnd) = 1.0_r8

    ! slightly greater than 1 hour offset
    grc%londeg(:) = 15.1_r8
    
  end subroutine setupEnvironment

  !-----------------------------------------------------------------------
  subroutine teardownEnvironment()
    !
    ! !DESCRIPTION:
    ! Tears down the environment set up by setupEnvironment. Should be called after each
    ! test. Note that this does NOT deallocate the subgrid variables - that cleanup
    ! needs to be done separately.
    !
    ! !USES:
    use pftconMod, only : pftcon
    !
    !-----------------------------------------------------------------------
    
    deallocate(pftcon%irrigated)

  end subroutine teardownEnvironment

end module IrrigationWrapperMod
