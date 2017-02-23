module kissvec_mod
! Fortran binding for the KISS vectorizable random number generator.

implicit none

integer, parameter :: r8 = selected_real_kind(12)

private
public :: kissvec

contains

  subroutine kissvec( seed1, seed2, seed3, seed4, ran_arr, length)

    ! We can assume that an r8 is a double and an i4 is an int32_t, but we can't
    ! make any guarantees about the relative sizes of a Fortran default integer
    ! and C's size_t.
    use iso_c_binding, only: c_size_t

    integer,                     intent(in)    :: length
    real(r8), dimension(length), intent(inout) :: ran_arr
    integer,  dimension(length), intent(inout) :: seed1, seed2, seed3, seed4

    ! C implementation
    interface
       subroutine kiss_rng(seed1, seed2, seed3, seed4, ran_arr, length) bind(c)
         ! Note that the definition of kiss_rng uses unsigned int, but the
         ! Fortran standard largely ignores signed/unsigned because there
         ! are no unsigned integers in Fortran.
         use iso_c_binding, only: c_int32_t, c_double, c_size_t
         integer(c_int32_t), intent(inout) :: seed1(*), seed2(*), seed3(*), seed4(*)
         real(c_double), intent(inout) :: ran_arr(*)
         integer(c_size_t), value :: length
       end subroutine kiss_rng
    end interface

    call kiss_rng(seed1, seed2, seed3, seed4, ran_arr, int(length, c_size_t))

  end subroutine kissvec

end module kissvec_mod
