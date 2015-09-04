module control

   implicit none
   private

   include 'netcdf.inc'
   

   logical, public :: verbose
   logical, public :: silent
   logical, public :: compute_gauss
   logical, public :: prec_override = .false.
   integer, public :: prec_out

   public :: &
      set_user_skip,    &
      set_user_include, &
      is_user_skip,    &
      is_user_include

   integer, parameter :: fexcl_max = 100
   character(len=nf_max_name) :: fexcl(fexcl_max) = ' '
   integer :: num_fexcl = 0

   integer, parameter :: fincl_max = 100
   character(len=nf_max_name) :: fincl(fincl_max) = ' '
   integer :: num_fincl = 0

!=======================================================================================
contains
!=======================================================================================

subroutine set_user_skip(var_list)

   character(len=*), intent(in) :: var_list

   integer :: i
   integer :: start, str_len, next_comma
   !-------------------------------------------------------------------

   call delimited_string_to_array(var_list, ',', fexcl, num_fexcl)

   if (verbose) write(6,*)'set_user_skip: fexcl:', trim(var_list)

end subroutine set_user_skip

!---------------------------------------------------------------------------------------

subroutine set_user_include(var_list)

   character(len=*), intent(in) :: var_list

   integer :: i
   integer :: start, str_len, next_comma
   !-------------------------------------------------------------------

   call delimited_string_to_array(var_list, ',', fincl, num_fincl)

   if (verbose) write(6,*)'set_user_include: fincl:', trim(var_list)

end subroutine set_user_include

!---------------------------------------------------------------------------------------

logical function is_user_skip(name)

   character(len=*), intent(in) :: name

   integer :: i

   ! look through list of user specified names to skip (i.e.,
   ! don't copy/interpolate these variables on input file to
   ! the output file

   is_user_skip = .false.
   do i = 1, num_fexcl
      if (trim(name) == trim(fexcl(i))) then
         is_user_skip = .true.
         return
      end if
   end do

end function is_user_skip

!---------------------------------------------------------------------------------------

logical function is_user_include(name)

   character(len=*), intent(in) :: name

   integer :: i

   ! look through list of user specified names to include, i.e.,
   ! copy these variables from the template file to the output file.

   is_user_include = .false.
   do i = 1, num_fincl
      if (trim(name) == trim(fincl(i))) then
         is_user_include = .true.
         return
      end if
   end do

end function is_user_include

!---------------------------------------------------------------------------------------

subroutine delimited_string_to_array(delim_str, delim, array, nelem)

   character(len=*), intent(in) :: delim_str  ! string containing delimiter separated tokens
   character(len=1), intent(in) :: delim      ! value of delimiter

   character(len=nf_max_name), intent(out) :: array(:) ! names
   integer,                    intent(out) :: nelem    ! number of tokens

   integer :: i, nelem_max
   integer :: start, str_len, next_delim
   !-------------------------------------------------------------------

   nelem_max = size(array)

   start = 1
   str_len = len_trim(delim_str)
   
   ! check that there's at least one variable specified
   if (str_len == 0) then
      write(6,*) 'delimited_string_to_array: ERROR - input string argument is empty'
      stop
   end if

   do i = 1, nelem_max

      ! scan for delimiter
      next_delim = scan(delim_str(start:str_len), delim)

      ! copy next variable into array
      if (next_delim == 0) then
         ! if scan returns zero there are no delimiters.  So use entire remaining
         ! string as a variable name and we're done
         array(i) = delim_str(start:str_len)
         exit
      else
         array(i) = delim_str(start:start+next_delim-2)
         start = start+next_delim
         if (start > str_len) exit
      end if

   end do

   nelem = i

end subroutine delimited_string_to_array

!---------------------------------------------------------------------------------------

end module control
