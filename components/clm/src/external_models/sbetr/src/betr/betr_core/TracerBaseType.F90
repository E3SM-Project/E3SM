module TracerBaseType
use betr_ctrl, only : max_betr_hist_type
use BeTRHistVarType, only : betr_hist_var_type
implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
  integer, parameter :: loc_str_len=255


  type, public :: tracerbase_type
     integer :: num_hist1d
     integer :: num_hist2d
     type(betr_hist_var_type), allocatable :: hist1d_var(:)
     type(betr_hist_var_type), allocatable :: hist2d_var(:)
  contains
     procedure, public :: tracer_base_init
     procedure, public :: add_hist_var2d
     procedure, public :: add_hist_var1d
     procedure, public :: alloc_hist_list
  end type tracerbase_type

contains

  !-----------------------------------------------------------------------
  subroutine tracer_base_init(this)
  implicit none
  class(tracerbase_type), intent(inout) :: this

  !set number of history files to zero
  this%num_hist1d = 0
  this%num_hist2d = 0
  end subroutine tracer_base_init

  !-----------------------------------------------------------------------
  subroutine add_hist_var2d(this, it, num2d, fname, units, type2d,  &
     avgflag, long_name, default)
  !
  !DESCRIPTION
  !build the namelist for an output 2d variable
  implicit none
  class(tracerbase_type), intent(inout) :: this
  integer         , intent(in):: it
  integer         , intent(inout):: num2d
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: type2d
  character(len=*), intent(in) :: avgflag
  character(len=*), intent(in) :: long_name
  character(len=*),optional, intent(in) :: default

  character(len=12) :: default_loc = "active"

  if(present(default))then
    default_loc=trim(default)
  endif

  num2d = num2d + 1
  if(it/=1)then
    this%hist2d_var(num2d)%varname=trim(fname)
    this%hist2d_var(num2d)%units=trim(units)
    this%hist2d_var(num2d)%avg_flag=trim(avgflag)
    this%hist2d_var(num2d)%long_name=trim(long_name)
    this%hist2d_var(num2d)%use_default=trim(default_loc)
  endif

  end subroutine add_hist_var2d

  !-----------------------------------------------------------------------
  subroutine add_hist_var1d(this, it, num1d, fname, units,  &
     avgflag, long_name, default)
  !
  !DESCRIPTION
  !build the namelist for an output 2d variable
  implicit none
  class(tracerbase_type), intent(inout) :: this
  integer         , intent(in):: it
  integer         , intent(inout):: num1d
  character(len=*), intent(in) :: fname
  character(len=*), intent(in) :: units
  character(len=*), intent(in) :: avgflag
  character(len=*), intent(in) :: long_name
  character(len=*),optional, intent(in) :: default

  character(len=12) :: default_loc = "active"

  if(present(default))then
    default_loc=trim(default)
  endif

  num1d = num1d + 1
  if(it/=1)then
    this%hist1d_var(num1d)%varname=trim(fname)
    this%hist1d_var(num1d)%units=trim(units)
    this%hist1d_var(num1d)%avg_flag=trim(avgflag)
    this%hist1d_var(num1d)%long_name=trim(long_name)
    this%hist1d_var(num1d)%use_default=trim(default_loc)
  endif



  end subroutine add_hist_var1d
  !-----------------------------------------------------------------------
  subroutine alloc_hist_list(this, num1d, num2d)

  implicit none
  class(tracerbase_type), intent(inout) :: this
  integer, intent(in) :: num1d, num2d
  integer :: jj

  !allocate memories

  this%num_hist1d = num1d
  this%num_hist2d = num2d
  allocate(this%hist1d_var(this%num_hist1d))
  allocate(this%hist2d_var(this%num_hist2d))

  end subroutine alloc_hist_list

end module TracerBaseType
