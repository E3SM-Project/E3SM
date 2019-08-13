module nanMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: nanMod
!
! !DESCRIPTION:
! Set parameters for the floating point flags "inf" Infinity
! and "nan" not-a-number. As well as "bigint" the point
! at which integers start to overflow. These values are used
! to initialize arrays with as a way to detect if arrays
! are being used before being set.
! Note that bigint is the largest possible 32-bit integer.
!
! !USES:
  use shr_kind_mod, only: r8 => shr_kind_r8
!
! !PUBLIC TYPES:
  implicit none
  save
#ifdef __PGI
! quiet nan for portland group compilers
  real(r8), parameter :: inf = O'0777600000000000000000'
  real(r8), parameter :: nan = O'0777700000000000000000'
  integer,  parameter :: bigint = O'17777777777'
#else
! signaling nan otherwise
  real(r8), parameter :: inf = O'0777600000000000000000'
  real(r8), parameter :: nan = O'0777610000000000000000'
  integer,  parameter :: bigint = O'17777777777'
#endif
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein based on cam module created by
! CCM core group
!
!EOP
!-----------------------------------------------------------------------

end module nanMod
