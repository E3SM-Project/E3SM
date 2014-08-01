! Utility functions in support of grid types (e.g., grid box average, sub-column)
module grid_flag_utils

use abortutils,    only: endrun

  implicit none
  private

  public clear_grid_types, is_col_type_set, set_col_type, add_col_type, clear_col_type, bit_field_kind

  ! --col_type parameters -- 
  ! A number from 0-31 which identifies which field type is being referenced
  ! Internally used to identify which bit of grid_type is used for the field

  integer, parameter, public :: flag_gridcol = 0, flag_subcol = 1

  ! --grid_type parameters -- 
  ! Each bit of grid_types variables turns on/off that field type
  ! predefined grid-type parameters

  integer,parameter :: bit_field_kind = kind(1)  ! integer field to store bit field info

contains

subroutine clear_grid_types(grid_types_flag)
!---------------------------
! Clear all of the flags
!---------------------------

  integer(bit_field_kind), intent(out) :: grid_types_flag

  grid_types_flag = 0

end subroutine clear_grid_types

logical function is_col_type_set(grid_types_flag, col_type)
!---------------------------
! Returns true if the requested col_type flag is set
!---------------------------
  integer(bit_field_kind), intent(in) :: grid_types_flag
  integer, intent(in) :: col_type

  if (col_type > bit_size(grid_types_flag)) then
     call endrun('is_col_type error: col_type is greater than the number of bits in grid_types_flag') 
  end if

  if (col_type < 0)  then
     call endrun('is_col_type error: col_type must be a positive number') 
  end if

  is_col_type_set = btest(grid_types_flag, col_type)

end function is_col_type_set

integer function set_col_type(col_type)
!---------------------------
! Sets the flag for the requested col_type
!---------------------------
  integer, intent(in) :: col_type

  if (col_type > bit_size(0)) then
     call endrun('set_col_type error: col_type is greater than the number of bits in grid_types_flag') 
  end if

  if (col_type < 0)  then
     call endrun('set_col_type error: col_type must be a positive number') 
  end if

  set_col_type = ibset(0, col_type)

end function set_col_type

integer function add_col_type(grid_types_flag, col_type)
!---------------------------
! Given an existing flag, adds the requested col_type to it
!---------------------------
  integer(bit_field_kind), intent(in) :: grid_types_flag
  integer, intent(in) :: col_type

  if (col_type > bit_size(grid_types_flag)) then
     call endrun('add_col_type error: col_type is greater than the number of bits in grid_types_flag') 
  end if

  if (col_type < 0)  then
     call endrun('add_col_type error: col_type must be a positive number') 
  end if

  add_col_type = ibset(grid_types_flag, col_type)

end function add_col_type

integer function clear_col_type(grid_types_flag, col_type)
!---------------------------
! Given an existing flag, clears the requested col_type from it
!---------------------------
  integer(bit_field_kind), intent(in) :: grid_types_flag
  integer, intent(in) :: col_type

  if (col_type > bit_size(grid_types_flag)) then
     call endrun('clear_col_type error: col_type is greater than the number of bits in grid_types_flag') 
  end if

  if (col_type < 0)  then
     call endrun('clear_col_type error: col_type must be a positive number') 
  end if

  clear_col_type = ibclr(grid_types_flag, col_type)

end function clear_col_type

end module grid_flag_utils
