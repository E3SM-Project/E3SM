program test_shr_spfn

use shr_kind_mod, only: r8 => shr_kind_r8, &
                        r4 => shr_kind_r4

use shr_spfn_mod
use test_mod

implicit none

real(r4) :: tmpr4
real(r8) :: tmpr8

character(len=2) :: num_str = ""

integer :: n, fac

call test_init(29)

!--------------------------------------------------
! Test erf at several values
!--------------------------------------------------
call test_is( 0._r4 == shr_spfn_erf(0.0_r4) , &
     "erf(0) is 0 (single precision)." )

call test_is( 1._r4 == shr_spfn_erf(15.0_r4) , &
     "erf(15) is 1 (single precision)." )

call test_is( -1._r4 == shr_spfn_erf(-15.0_r4) , &
     "erf(-15) is -1 (single precision)." )

call test_is( 0._r8 == shr_spfn_erf(0.0_r8) , &
     "erf(0) is 0 (double precision)." )

call test_is( 1._r8 == shr_spfn_erf(30.0_r8) , &
     "erf(30) is 1 (double precision)." )

call test_is( -1._r8 == shr_spfn_erf(-30.0_r8) , &
     "erf(-30) is -1 (double precision)." )

call test_is( shr_spfn_erf(1._r4) == shr_spfn_erf(1.0_r4) , &
     "erf is reproducible (single precision)." )

call test_is( shr_spfn_erf(1._r8) == shr_spfn_erf(1.0_r8) , &
     "erf is reproducible (double precision)." )

!--------------------------------------------------
! Test erfc at several values
!--------------------------------------------------
call test_is( 1._r4 == shr_spfn_erfc(0.0_r4) , &
     "erfc(0) is 1 (single precision)." )

call test_is( 0._r4 == shr_spfn_erfc(15.0_r4) , &
     "erfc(15) is 0 (single precision)." )

call test_is( 2._r4 == shr_spfn_erfc(-15.0_r4) , &
     "erfc(-15) is 2 (single precision)." )

call test_is( 1._r8 == shr_spfn_erfc(0.0_r8) , &
     "erfc(0) is 1 (double precision)." )

call test_is( 0._r8 == shr_spfn_erfc(30.0_r8) , &
     "erfc(30) is 0 (double precision)." )

call test_is( 2._r8 == shr_spfn_erfc(-30.0_r8) , &
     "erfc(-30) is 2 (double precision)." )

call test_is( shr_spfn_erfc(1._r4) == shr_spfn_erfc(1.0_r4) , &
     "erfc is reproducible (single precision)." )

call test_is( shr_spfn_erfc(1._r8) == shr_spfn_erfc(1.0_r8) , &
     "erfc is reproducible (double precision)." )

!--------------------------------------------------
! Test erfc_scaled at several values
!--------------------------------------------------

call test_is( 1._r4 == shr_spfn_erfc_scaled(0._r4) , &
     "erfc_scaled(0) = 1 (single precision)." )

call test_is( 1._r8 == shr_spfn_erfc_scaled(0._r8) , &
     "erfc_scaled(0) = 1 (double precision)." )

! Check relative error at 1.
tmpr4 = exp(1._r4**2)*shr_spfn_erfc(1._r4)
tmpr4 = (shr_spfn_erfc_scaled(1._r4)-tmpr4) / tmpr4

! Tolerance here is kind of an arbitrary criterion.
call test_is( abs(tmpr4) < 1.e-4_r4 , &
     "erfc_scaled(1) ~= exp(1)*erfc(1) (single precision)." )

! Check relative error at 1.
tmpr8 = exp(1._r8**2)*shr_spfn_erfc(1._r8)
tmpr8 = (shr_spfn_erfc_scaled(1._r8)-tmpr8) / tmpr8

! Precision here is kind of an arbitrary criterion.
call test_is( abs(tmpr8) < 1.e-8_r8 , &
     "erfc_scaled(1) ~= exp(1)*erfc(1) (double precision)." )

!--------------------------------------------------
! Test gamma equal to factorials for small ints.
!--------------------------------------------------

! Define 0! = 1.
fac = 1
do n = 1, 4
   ! Get string for this argument.
   write(num_str,"(I2)") n

   ! Only have double-precision gamma.
   tmpr8 = shr_spfn_gamma(real(n,r8))

   call test_is( tmpr8 == real(fac,r8) , &
        "gamma(n) = n+1! for n="//trim(num_str)//"." )

   call test_is( tmpr8 == shr_spfn_igamma(real(n,r8),0._r8) , &
        "gamma(n) = igamma(n,0) for n="//trim(num_str)//"." )

   ! Get next factorial.
   fac = fac*(n)
end do

!--------------------------------------------------
! Test one value for igamma with two nonzero arguments.
!--------------------------------------------------

! igamma(1,x) = exp(-x)
tmpr8 = exp(-1._r8)
tmpr8 = (shr_spfn_igamma(1._r8, 1._r8)-tmpr8) / tmpr8
call test_is( tmpr8 < 1.e-8_r8 , &
     "igamma(1,1) = exp(-1)." )

call test_final()

end program test_shr_spfn
