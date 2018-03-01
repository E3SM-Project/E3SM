module centuryecaOutType
  use bshr_kind_mod       , only : r8 => shr_kind_r8
  use bshr_log_mod  , only : errMsg => shr_log_errMsg
implicit none
  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: centuryeca_out_type
    real(r8), pointer :: ystates0(:)  => null()
    real(r8), pointer :: ystatesf(:)  => null()
  contains
    procedure, public :: Init
    procedure, private :: InitAllocate

  end type centuryeca_out_type

  public :: create_century_out_type
contains
  !-------------------------------------------------------------------------------
  function create_century_out_type()
  ! DESCRIPTION
  ! constructor
    implicit none
    class(centuryeca_out_type), pointer :: create_century_out_type
    class(centuryeca_out_type), pointer :: centout

    allocate(centout)
    create_century_out_type => centout

  end function create_century_out_type
  !-------------------------------------------------------------------------------

  subroutine Init(this, nstates)
  implicit none
  class(centuryeca_out_type), intent(inout) :: this
  integer, intent(in) :: nstates


  call this%InitAllocate(nstates)
  end subroutine Init
  !-------------------------------------------------------------------------------

  subroutine InitAllocate(this, nstates)
  implicit none
  class(centuryeca_out_type), intent(inout) :: this
  integer, intent(in) :: nstates


  allocate(this%ystates0(nstates))
  allocate(this%ystatesf(nstates))
  end subroutine InitAllocate


end module centuryecaOutType
