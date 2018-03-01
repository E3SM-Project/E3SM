module betr_regression_module

  use betr_constants, only : betr_filename_length

  implicit none

  private
  character(len=*), private, parameter :: mod_filename = &
       __FILE__

  type, public :: betr_regression_type
     character(len=betr_filename_length), private :: filename
     integer :: output

     logical, public  :: write_regression_output
     integer, private :: num_cells
     ! FIXME(bja, 201603) specifying cell ids requires more careful
     !thought. Maybe just hard code a max....

     !X!integer, allocatable, private :: cell_ids(:)
   contains
     procedure, public  :: Init
     procedure, public  :: OpenOutput
     procedure, public  :: CloseOutput
     procedure, public  :: WriteData
     procedure, private :: ReadNamelist
     procedure, private :: CheckInput
  end type betr_regression_type

contains

  subroutine Init(this, base_filename, namelist_buffer, bstatus)

    use betr_constants, only : betr_namelist_buffer_size
    use BetrStatusType , only : betr_status_type
    implicit none

    class(betr_regression_type)              , intent(inout) :: this
    character(len=betr_filename_length)      , intent(in)    :: base_filename
    character(len=betr_namelist_buffer_size) , intent(in)    :: namelist_buffer
    class(betr_status_type)        , intent(out)   :: bstatus

    this%output = 16
    this%filename = trim(base_filename) // '.regression'

    call this%ReadNamelist(namelist_buffer, bstatus)
    if(bstatus%check_status())return
    call this%CheckInput(bstatus)
  end subroutine Init

  !--------------------------------------------------------------------

  subroutine ReadNamelist(this, namelist_buffer, bstatus)
    use bshr_log_mod   , only : errMsg => shr_log_errMsg
    use betr_constants , only : stdout, betr_string_length_long, betr_namelist_buffer_size
    use BetrStatusType , only : betr_status_type
    implicit none

    class(betr_regression_type), intent(inout) :: this
    character(len=betr_namelist_buffer_size), intent(in) :: namelist_buffer
    class(betr_status_type)        , intent(out)   :: bstatus

    character(len=*), parameter :: subname = 'betr_regression:ReadNamelist'
    ! !LOCAL VARIABLES:
    integer :: nml_error
    integer :: cells
    character(len=betr_string_length_long) :: ioerror_msg
    logical :: write_regression_output
    !-----------------------------------------------------------------------

    namelist / regression_test / cells, write_regression_output
    call bstatus%reset()
    cells = 0
    write_regression_output = .false.

    ! ----------------------------------------------------------------------
    ! Read namelist from standard input.
    ! ----------------------------------------------------------------------

    if ( .true. )then
       ioerror_msg=''
       read(namelist_buffer, nml=regression_test, iostat=nml_error, iomsg=ioerror_msg)
       if (nml_error /= 0) then
          call bstatus%set_msg(msg="ERROR reading betr_regression_test namelist "//errmsg(mod_filename, __LINE__),err=-1)
          return
       end if
    end if

    if (.true.) then
       write(stdout, *)
       write(stdout, *) '--------------------'
       write(stdout, *)
       write(stdout, *) ' betr regression test type :'
       write(stdout, *)
       write(stdout, *) ' regression_test namelist settings :'
       write(stdout, *)
       write(stdout, regression_test)
       write(stdout, *)
       write(stdout, *) '--------------------'
    endif

    this%write_regression_output = write_regression_output
    this%num_cells = cells
  end subroutine ReadNamelist

  !---------------------------------------------------------------------------------

  subroutine CheckInput(this, bstatus)

    use betr_constants , only : betr_string_length_long
    use bshr_log_mod   , only : errMsg => shr_log_errMsg
    use BetrStatusType , only : betr_status_type
    implicit none

    class(betr_regression_type), intent(in) :: this
    class(betr_status_type)        , intent(out)   :: bstatus
    character(len=betr_string_length_long) :: msg

    if (this%write_regression_output) then
       ! sanity check on input
       if (this%num_cells < 0) then
          msg = 'ERROR num cells must be >= 0. '
          call bstatus%set_msg(msg=msg//errmsg(mod_filename, __LINE__),err=-1)
          return
       end if
    end if
  end subroutine CheckInput

  !---------------------------------------------------------------------------------

  subroutine OpenOutput(this)

    use betr_constants, only : stdout

    implicit none

    class(betr_regression_type), intent(inout) :: this

    integer :: output = 16

    write(stdout, '(a, a)') 'Writing regression output to ', this%filename

    open(this%output, file=this%filename, status='REPLACE')

  end subroutine OpenOutput


  !---------------------------------------------------------------------------------

  subroutine CloseOutput(this)

    use betr_constants, only : stdout

    implicit none

    class(betr_regression_type), intent(inout) :: this

    close(this%output)

  end subroutine CloseOutput


  !---------------------------------------------------------------------------------

  subroutine WriteData(this, category, name, data)
    !
    ! Write regression data in a cfg/ini format that can easily be
    ! parsed by external processors.
    !
    ! sections names are the species name (and possible
    ! phase/location?)
    !
    ! section data is just keyword = value for things like min, max,
    ! mean, vector norms, cell point data.

    use bshr_kind_mod  , only : r8 => shr_kind_r8
    use betr_constants , only : stdout
    use betr_constants , only : betr_var_name_length, betr_string_length

    implicit none

    class(betr_regression_type)         , intent(inout) :: this
    character(len=betr_string_length)   , intent(in)    :: category
    character(len=betr_var_name_length) , intent(in)    :: name
    real(r8)                            , intent(in)    :: data(:)

    integer :: cell_increment, num_cells, cell
    real(r8) :: val, local_val

    ! FIXME(bja, 201603) cfg/ini format limits the characters that can
    ! be use in section names. We need to sanitize the names!
    write(this%output, '("[",a,"]")') trim(name)

    write(this%output, '("category = ",a)') trim(category)

    val = minval(data(:))
    if(abs(val)<1.e-50_r8)val=0._r8
    write(this%output, '("min = ",e21.13)') val

    val = maxval(data(:))
    if(abs(val)<1.e-50_r8)val=0._r8
    write(this%output, '("max = ",e21.13)') val

    val = sum(data(:)) / size(data)
    if(abs(val)<1.e-50_r8)val=0._r8
    write(this%output, '("mean = ",e21.13)') val

    if (this%num_cells > 0) then
       num_cells = this%num_cells
       if (num_cells > size(data)) then
          ! need to truncate
          num_cells = size(data)
       end if

       cell_increment = int(size(data) / num_cells)

       do cell = 1, size(data), cell_increment
          local_val = data(cell)
          if(abs(local_val)<1.e-50_r8)local_val=0._r8
          write(this%output, '("cell ", i4, " = ", e21.13)') cell, local_val
       end do
    write(this%output, '(a)')
    end if
  end subroutine WriteData


end module betr_regression_module
