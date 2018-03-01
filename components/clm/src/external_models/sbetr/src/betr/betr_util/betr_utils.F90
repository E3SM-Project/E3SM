module betr_utils

  implicit none
  character(len=*), private, parameter :: mod_filename = &
       __FILE__
contains
  ! ----------------------------------------------------------------------
  function remove_filename_extension(filename) result(basename)
    !
    ! Remove the extension from a file name to get a base filename.
    !
    ! We start at the end of the filename and assume that the extension
    ! is marked by a period.
    !

    use betr_constants, only : betr_filename_length

    implicit none

    character(len=betr_filename_length), intent(in) :: filename
    character(len=betr_filename_length) :: basename
    integer :: ext_index

    ext_index = scan(filename, '.', .true.)
    if (ext_index == 0) then
       ! no period marking an extension...
       ext_index = len(trim(filename)) + 1

    end if
    basename = filename(1:ext_index-1)
  end function remove_filename_extension

  ! ----------------------------------------------------------------------
  function num2str(num, fmt)result(ans)
  use betr_constants, only : betr_string_length_long

  implicit none
  integer :: num
  character(len=*), intent(in) :: fmt
  character(len=betr_string_length_long) :: ans

  write(ans,fmt)num
  return
  end function num2str
  ! ----------------------------------------------------------------------
  function log2str(logval)result(str)
  implicit none
  logical, intent(in) :: logval
  character(len=8) :: str

  if(logval)then
    str='.true.'
  else
    str='.false.'
  endif
  end function log2str
end module betr_utils
