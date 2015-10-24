module TracerBoundaryCondType
!
! !DESCRIPTION:
! data type to specify boundary conditions for tracer tranpsort
!
! !USES:
   use shr_kind_mod   , only : r8 => shr_kind_r8
   use decompMod      , only : bounds_type

  implicit none
  save
  private
  !
  ! !PUBLIC DATA:
  !

  !--------------------------------------------------------------------------------
  type, public :: tracerboundarycond_type
     real(r8), pointer :: tracer_gwdif_concflux_top_col( : , : , : ) !tracer concentration or incoming flux imposed at top boundary for dual diffusion calculation
     real(r8), pointer :: condc_toplay_col       ( : , : )           !conductance at the column-air interface
     real(r8), pointer :: bot_concflux_col       ( : , : , : )       !bottom boundary condition
     integer,  pointer :: topbc_type             ( : )               !type of top boundary condition, it depends on tracer type
     integer,  pointer :: botbc_type             ( : )               !type of bottom boundary condition, it depends on tracer type
     integer,  pointer :: jtops_col              ( : )               !index of the top numerical node
   contains
     procedure, public  :: Init
     procedure, public  :: Restart
     procedure, public  :: Reset
     procedure, private :: InitAllocate
     procedure, private :: InitHistory
     procedure, private :: InitCold

  end type tracerboundarycond_type
contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds, betrtracer_vars)
    !
    ! !DESCRIPTION:
    ! Initialize the datatype
    !
    ! !USES:
    use BeTRTracerType, only : BeTRTracer_Type
    !
    ! !ARGUMENTS:
    class(tracerboundarycond_type) :: this
    type(bounds_type), intent(in) :: bounds
    type(BeTRTracer_Type), intent(in) :: betrtracer_vars

    call this%InitAllocate(bounds, betrtracer_vars)
    call this%InitHistory(bounds)
    call this%InitCold(bounds)

  end subroutine Init

  !-----------------------------------------------------------------------
  subroutine InitAllocate(this, bounds, betrtracer_vars)
    use BeTRTracerType, only : BeTRTracer_Type
    !
    ! !DESCRIPTION:
    ! allocate memories to relevant variables
    !
    ! !ARGUMENTS:
    class(tracerboundarycond_type)     :: this
    type(bounds_type),      intent(in) :: bounds
    type(BeTRTracer_Type), intent(in)  :: betrtracer_vars
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    !---------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    allocate(this%tracer_gwdif_concflux_top_col  (begc:endc, 1:2, 1:betrtracer_vars%ntracers))    ! 1: values at previous time step, 2: values at current time step
    allocate(this%bot_concflux_col               (begc:endc, 1:2, 1:betrtracer_vars%ntracers))    ! 1: values at previous time step, 2: values at current time step

    allocate(this%condc_toplay_col         (begc:endc, 1:betrtracer_vars%ntracer_groups))
    allocate(this%topbc_type               (1:betrtracer_vars%ntracer_groups))
    allocate(this%botbc_type               (1:betrtracer_vars%ntracer_groups))
    allocate(this%jtops_col                (begc:endc))
  end subroutine InitAllocate

  !-----------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! History fields initialization
    !
    ! !USES:
    use clm_varcon    , only: spval
    use clm_varpar    , only: nlevsno

    !
    ! !ARGUMENTS:
    class(tracerboundarycond_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc
    real(r8), pointer :: data2dptr_col(:,:) ! temp. pointers for slicing larger arrays


  end subroutine InitHistory
  
  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! do cold initialization
    ! !USES:
    !
    use shr_infnan_mod, only : nan => shr_infnan_nan, assignment(=)
    implicit none
    ! !ARGUMENTS:
    class(tracerboundarycond_type) :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:


    !-----------------------------------------------------------------------
    this%topbc_type(:)                        = -1
    this%botbc_type(:)                        = -1
    this%tracer_gwdif_concflux_top_col(:,:,:) = nan
    this%condc_toplay_col(:,:)                = nan
    this%bot_concflux_col(:,:,:)              = 0._r8
    this%jtops_col(:)                         = 1
  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
    !
    ! !DESCRIPTION:
    ! Read/Write module information to/from restart file.
    !
    ! !USES:
    use clm_varpar , only : nlevsno, nlevsoi
    use clm_varcon , only : spval
    use clm_varctl , only : iulog
    use ncdio_pio
    !
    ! !ARGUMENTS:
    class(tracerboundarycond_type)      :: this
    type(bounds_type)   , intent(in)    :: bounds
    type(file_desc_t)   , intent(inout) :: ncid                                         ! netcdf id
    character(len=*)    , intent(in)    :: flag                                         ! 'read' or 'write'
    !
    ! !LOCAL VARIABLES:
    integer :: j,c     ! indices
    logical :: readvar ! determine if variable is on initial file


  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine Reset(this, column)
    !
    ! !DESCRIPTION:
    ! Intitialize SNICAR variables for fresh snow column
    !
    ! !ARGUMENTS:
    class(tracerboundarycond_type) :: this
    integer           , intent(in) :: column     ! column index



  end subroutine Reset
 end module TracerBoundaryCondType
