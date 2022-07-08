program test_vmath

!
! This is a test for the shr_vmath_mod module.
!

use shr_kind_mod, only: r8 => shr_kind_r8
use shr_kind_mod, only: r4 => shr_kind_r4
use shr_kind_mod, only: i8 => shr_kind_i8
use shr_kind_mod, only: i4 => shr_kind_i4
use shr_const_mod, only: pi => shr_const_pi
use shr_vmath_mod

implicit none
integer, parameter :: vlen = 128
real(r8) :: ivec(vlen), rvec(vlen), ovec(vlen), nvec(vlen)
real(r8), parameter :: bigval = 1.0E300_r8
real(r8), parameter :: smallval = 1.0E-300_r8
real(r8), parameter :: tolerance = 1.0E-15_r8
integer :: i
call random_number(ivec)  ! numbers between 0 and 1

ivec = ivec * bigval       ! numbers between 0 and 1e308

call shr_vmath_sqrt(ivec, rvec, vlen)

ovec = dsqrt(ivec)
do i=1,vlen
   if(abs(rvec(i)-ovec(i)) > tolerance) then
      print *,__LINE__,i, ivec(i),rvec(i),ovec(i)
   endif
enddo

rvec = (rvec - ovec)/ovec

call assert(all(abs(rvec) < tolerance),"shr_vmath_sqrt test failed")

call shr_vmath_rsqrt(ivec, rvec, vlen)

ovec = 1.0_r8/ovec

do i=1,vlen
   if(abs((rvec(i)-ovec(i))/ovec(i)) > tolerance) then
      print *,__LINE__,i, ivec(i),rvec(i),ovec(i)
   endif
enddo

rvec = (rvec - ovec)/ovec

call assert(all(abs(rvec) < tolerance),"shr_vmath_rsqrt test failed")

call random_number(nvec)
nvec = (nvec - 0.5_r8)*bigval

call shr_vmath_div(ivec, nvec, rvec, vlen)

ovec = ivec/nvec

rvec = (rvec - ovec)/ovec

call assert(all(abs(rvec) < tolerance),"shr_vmath_div test failed")

call random_number(ivec)
ivec = ivec*1400_r8 - 700_r8

call shr_vmath_exp(ivec, rvec, vlen)

ovec = exp(ivec)
!print *,minval(abs(rvec)),maxval(rvec)

rvec = (rvec - ovec)/ovec

call assert(all(abs(rvec) < tolerance),"shr_vmath_exp test failed")

ivec = ovec
call shr_vmath_log(ivec, rvec, vlen)
ovec = log(ivec)

rvec = (rvec - ovec)/ovec

call assert(all(abs(rvec) < tolerance),"shr_vmath_log test failed")

call random_number(ivec)
ivec = (ivec-0.5_r8)*2.0_r8*pi
call shr_vmath_sin(ivec, rvec, vlen)
ovec = sin(ivec)
rvec = (rvec - ovec)/ovec

call assert(all(abs(rvec) < tolerance),"shr_vmath_sin test failed")

call shr_vmath_cos(ivec, rvec, vlen)
ovec = cos(ivec)
rvec = (rvec - ovec)/ovec

call assert(all(abs(rvec) < tolerance),"shr_vmath_cos test failed")

contains

  subroutine assert(val, msg)
    logical, intent(in) :: val
    character(len=*), intent(in) :: msg

    if (.not. val) then
       print *, msg
       stop 1
    end if

  end subroutine assert

end program test_vmath
