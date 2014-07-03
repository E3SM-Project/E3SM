module test_dynVarShared
  ! Shared code to set up tests of dyn_var_type and its extensions

  use shr_kind_mod, only : r8 => shr_kind_r8
  use dynFileMod, only : dyn_file_type
  use ncdio_pio, only : ncd_set_var

  implicit none
  private
  save

  public :: create_dyn_file

contains
  
  function create_dyn_file(cur_year) result(dyn_file)
    ! Set up a dyn_file variable for tests. Assumes we're using the mock version of
    ! dynFileMod.
    !
    ! The years in the mock "file" go from 11 - 14. 

    ! The "file" contains two variables: foo_1d, which is a 1-d variable (i.e., just space
    ! & time, no level dimension); and foo_2d, which is a 2-d variable (i.e., includes a
    ! level dimension)

    type(dyn_file_type) :: dyn_file
    integer, intent(in) :: cur_year  ! current model year

    real(r8) :: data1d(3,4)   ! space & time only
    real(r8) :: data2d(6,4)   ! space & level & time; first two dimensions are [2,3]
    integer :: i, lev, time

    dyn_file = dyn_file_type([11,12,13,14], cur_year)
    
    data1d = reshape([1._r8, 2._r8, 3._r8, &   ! year 11
                    4._r8, 5._r8, 6._r8, &     ! year 12
                    7._r8, 8._r8, 9._r8, &     ! year 13
                    10._r8,11._r8,12._r8], &   ! year 14
                    [3, 4])
    call ncd_set_var(dyn_file, 'foo_1d', data1d, [3])

    data2d = reshape([ 1._r8,  2._r8,  3._r8,  4._r8,  5._r8,  6._r8, &  ! year 11
                       7._r8,  8._r8,  9._r8, 10._r8, 11._r8, 12._r8, &  ! year 12
                      13._r8, 14._r8, 15._r8, 16._r8, 17._r8, 18._r8, &  ! year 13
                      19._r8, 20._r8, 21._r8, 22._r8, 23._r8, 24._r8],&  ! year 14
                      [6, 4])
    call ncd_set_var(dyn_file, 'foo_2d', data2d, [2, 3])

  end function create_dyn_file
    
end module test_dynVarShared
