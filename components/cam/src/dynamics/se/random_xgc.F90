! random number generator from XGC1 code
! imported November 19, 2014 by P. Worley (worleyph@ornl.gov)

MODULE Ecuyer_random
! L'Ecuyer's 1996 random number generator.
! Fortran version by Alan.Miller @ vic.cmis.csiro.au
! N.B. This version is compatible with Lahey's ELF90
! http://www.ozemail.com.au/~milleraj
! Latest revision - 30 March 1999

  IMPLICIT NONE
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(14, 60)

  ! These are unsigned integers in the C version
  INTEGER, SAVE :: s1 = 1234, s2 = -4567, s3 = 7890

CONTAINS

  SUBROUTINE init_seeds(i1, i2, i3)
    IMPLICIT NONE

    INTEGER, INTENT(IN) :: i1, i2, i3

    s1 = i1
    s2 = i2
    s3 = i3
    IF (IAND(s1,-2) == 0) s1 = i1 - 1023
    IF (IAND(s2,-8) == 0) s2 = i2 - 1023
    IF (IAND(s3,-16) == 0) s3 = i3 - 1023

    RETURN
  END SUBROUTINE init_seeds

  FUNCTION taus88() RESULT(random_numb)
    ! Generates a random number between 0 and 1.  Translated from C function in:
    ! Reference:
    ! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
    ! generators', Math. of Comput., 65, 203-213.

    ! The cycle length is claimed to be about 2^(88) or about 3E+26.
    ! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

    IMPLICIT NONE
    REAL (dp) :: random_numb

    INTEGER   :: b

    ! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
    !      to the left if j > 0, otherwise to the right.

    b  = ISHFT( IEOR( ISHFT(s1,13), s1), -19)
    s1 = IEOR( ISHFT( IAND(s1,-2), 12), b)
    b  = ISHFT( IEOR( ISHFT(s2,2), s2), -25)
    s2 = IEOR( ISHFT( IAND(s2,-8), 4), b)
    b  = ISHFT( IEOR( ISHFT(s3,3), s3), -11)
    s3 = IEOR( ISHFT( IAND(s3,-16), 17), b)
    random_numb = IEOR( IEOR(s1,s2), s3) * 2.3283064365E-10_dp + 0.5_dp

    RETURN
  END FUNCTION taus88

  SUBROUTINE init_seeds_ext(sv)
    IMPLICIT NONE

    INTEGER, INTENT(INOUT) :: sv(3)

    IF (IAND(sv(1),-2) == 0)  sv(1) = sv(1) - 1023
    IF (IAND(sv(2),-8) == 0)  sv(2) = sv(2) - 1023
    IF (IAND(sv(3),-16) == 0) sv(3) = sv(3) - 1023

    RETURN
  END SUBROUTINE init_seeds_ext

  FUNCTION taus88_ext(sv) RESULT(random_numb)
    ! Generates a random number between 0 and 1.  Translated from C function in:
    ! Reference:
    ! L'Ecuyer, P. (1996) `Maximally equidistributed combined Tausworthe
    ! generators', Math. of Comput., 65, 203-213.

    ! The cycle length is claimed to be about 2^(88) or about 3E+26.
    ! Actually - (2^31 - 1).(2^29 - 1).(2^28 - 1).

    IMPLICIT NONE
    INTEGER, INTENT(INOUT) :: sv(3)
    REAL (dp) :: random_numb

    INTEGER   :: b, i1, i2, i3

    ! N.B. ISHFT(i,j) is a bitwise (non-circular) shift operation;
    !      to the left if j > 0, otherwise to the right.

    i1 = sv(1)
    i2 = sv(2)
    i3 = sv(3)
    b  = ISHFT( IEOR( ISHFT(i1,13), i1), -19)
    i1 = IEOR( ISHFT( IAND(i1,-2), 12), b)
    b  = ISHFT( IEOR( ISHFT(i2,2), i2), -25)
    i2 = IEOR( ISHFT( IAND(i2,-8), 4), b)
    b  = ISHFT( IEOR( ISHFT(i3,3), i3), -11)
    i3 = IEOR( ISHFT( IAND(i3,-16), 17), b)
    random_numb = IEOR( IEOR(i1,i2), i3) * 2.3283064365E-10_dp + 0.5_dp
    sv(1) = i1
    sv(2) = i2
    sv(3) = i3

    RETURN
  END FUNCTION taus88_ext

END MODULE Ecuyer_random

!! random number generator module
module random_xgc
  use Ecuyer_random, only: init_seeds, init_seeds_ext, taus88, taus88_ext
  use shr_kind_mod, only: r8 => shr_kind_r8
  implicit none
  type seeds_type
     integer, pointer :: s(:)
  end type seeds_type
  type(seeds_type), allocatable :: sv(:)

contains
  !random number - seed initialize
  ! if ranx is to be executed in a threaded region,
  ! then each thread must set its own seed, otherwise
  ! each thread generates the identical sequence
  subroutine init_ranx(seed)
    implicit none

    integer, intent(in) :: seed
#ifdef _OPENMP
    integer :: omp_get_num_threads
    integer :: omp_get_thread_num
    integer :: num_threads, thread_id
    integer :: i

    thread_id   = omp_get_thread_num()

!$OMP CRITICAL
    num_threads = omp_get_num_threads()
    if (.not. allocated(sv)) then
       allocate( sv(0:num_threads-1) )
       do i=0,num_threads-1
          allocate(sv(i)%s(3))
          sv(i)%s(1) = 1234*seed
          sv(i)%s(2) = 2345*seed + 6789
          sv(i)%s(3) = 4321*seed + 10
          call init_seeds_ext(sv(i)%s)
       end do
    endif
!$OMP END CRITICAL

    sv(thread_id)%s(1) = 1234*seed
    sv(thread_id)%s(2) = 2345*seed + 6789
    sv(thread_id)%s(3) = 4321*seed + 10
    call init_seeds_ext(sv(thread_id)%s)

#else
    call init_seeds(1234*seed,2345*seed+6789,4321*seed+10)
#endif
  end subroutine init_ranx

  ! random number generator
  function ranx()
    implicit none
    real(r8) :: ranx
#ifdef _OPENMP
    integer :: thread_id
    integer :: omp_get_thread_num
    thread_id = omp_get_thread_num()
    ranx=taus88_ext( sv(thread_id)%s )
#else
    ranx=taus88()
#endif
  end function ranx

end module random_xgc
