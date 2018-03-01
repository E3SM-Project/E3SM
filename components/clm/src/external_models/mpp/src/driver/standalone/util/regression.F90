module regression_mod

  implicit none

#include <petsc/finclude/petsc.h>

  type, public :: regression_type
     character(len=300), private :: filename
     integer                     :: output
     integer, private            :: num_cells
   contains
     procedure, public :: Init
     procedure, public :: OpenOutput
     procedure, public :: CloseOutput
     procedure, public :: WriteData
  end type regression_type

contains

  !--------------------------------------------------------------------

  subroutine Init(this, base_filename, num_cells)
    !
    implicit none
    !
    class(regression_type), intent(inout) :: this
    character(len=256)    , intent(in)    :: base_filename
    PetscInt                              :: num_cells

    this%output    = 16
    this%filename  = trim(base_filename) // '.regression'
    this%num_cells = num_cells

  end subroutine Init

  !--------------------------------------------------------------------

  subroutine OpenOutput(this)
    !
    implicit none
    !
    class(regression_type) :: this
    !
    PetscInt :: status

    open(this%output, file=this%filename, status='REPLACE', iostat=status)

    if (status /= 0 ) then
       write(*,*)'Unable to open regression file for output'
       stop
    endif

  end subroutine OpenOutput

  !--------------------------------------------------------------------

  subroutine CloseOutput(this)
    !
    implicit none
    !
    class(regression_type) :: this
    !
    PetscInt :: status

    close(this%output,  iostat=status)

    if (status /= 0 ) then
       write(*,*)'Unable to close regression file'
       stop
    endif

  end subroutine CloseOutput

  !--------------------------------------------------------------------

  subroutine WriteData(this, name, category, data)
    !
    implicit none
    !
    class(regression_type) , intent(in) :: this
    character(len=64)      , intent(in)  :: category
    character(len=64)      , intent(in)  :: name
    PetscReal              , intent(in), pointer :: data(:)
    !
    !
    PetscReal                      :: val
    PetscInt                       :: cell_increment
    PetscInt                       :: cell
    PetscInt                       :: num_cells_local

    write(this%output, '("[",a,"]")') trim(name)

    write(this%output, '("category = ",a)') trim(category)

    val = minval(data(:))
    if(abs(val)<1.d-50) val=0.d0
    write(this%output, '("min = ",e21.13)') val

    val = maxval(data(:))
    if(abs(val)<1.d-50)val=0.d0
    write(this%output, '("max = ",e21.13)') val

    val = sum(data(:)) / size(data)
    if(abs(val)<1.d-50)val=0.d0
    write(this%output, '("mean = ",e21.13)') val

    if (this%num_cells > 0) then
       num_cells_local = this%num_cells
       if (num_cells_local > size(data)) then
          ! need to truncate
          num_cells_local = size(data)
       end if

       cell_increment = int(size(data) / num_cells_local)

       do cell = 1, size(data), cell_increment
          val = data(cell)
          if(abs(val)<1.d-50)val=0.d0
          write(this%output, '("cell ", i4, " = ", e21.13)') cell, val
       end do
       write(this%output, '(a)')
    end if

  end subroutine WriteData

end module regression_mod
