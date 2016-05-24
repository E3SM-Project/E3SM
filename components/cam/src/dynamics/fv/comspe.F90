
module comspe
!BOP
!
! !MODULE: comspe
!
! !USES:
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   implicit none

!
! !PUBLIC DATA MEMBERS:
!
!    Spectral space arrays
!
   real(r8) vz(psp,plev)      ! Vorticity spectral coefficients
   real(r8) d(psp,plev)       ! Divergence spectral coefficients
   real(r8) t(psp,plev)       ! Temperature spectral coefficients
   real(r8) q(psp,plev)       ! Moisture     spectral coefficients
   real(r8) alps(psp)         ! Log-pressure spectral coefficients
   real(r8) hs(psp,plev)      ! hydrostatic matrix for "real" atmosphere
   real(r8) hsnm(psp,plev)    ! vertical normal modes of "hs"
   real(r8) dsnm(psp,plev)    ! vertical normal modes of "ds"
   real(r8) dnm(psp,plev)     ! vertical normal modes of "d"
   real(r8) vznm(psp,plev)    ! vertical normal modes of "vz"
   real(r8) a0nm(psp)         ! wave # coefs (use in vert normal mode space)
   real(r8) bmnm(psp)         ! wave # coefs (use in vert normal mode space)
   real(r8) bpnm(psp)         ! wave # coefs (use in vert normal mode space)
   real(r8) atri(psp)         ! wave # coefs (use in vert normal mode space)
   real(r8) btri(psp)         ! wave # coefs (use in vert normal mode space)
   real(r8) ctri(psp)         ! wave # coefs (use in vert normal mode space)

   integer ncutoff   ! Break-even point for vector lengths in GRCALC
   integer nalp(pmax)        ! Pointer into polynomial arrays
#if ( defined SPMD )
   integer begm(0:plat-1)
   integer endm(0:plat-1)
#else
   integer begm(0:0)
   integer endm(0:0)
   data begm,endm/1,pmmax/
#endif

   integer ncoefi(pmaxp)     ! Pointer to start of coefficient diagonals
   integer nm(pmax)          ! Number of coeffs stored on a given diagonal
   integer nco2(pmax)        ! Complex form of ncoefi
   integer nstart(pmmax)     ! Starting indices for spectral arrays (real)
   integer nlen(pmmax)

   real(r8) :: alp(pspt,plat/2)
   real(r8) :: dalp(pspt,plat/2)

!
! !DESCRIPTION:
!
!   This module defines variables for the spectral algorithm 
!   (probably not used in LRDC).
!
! !REVISION HISTORY:
!   00.10.19   Rosinski   Creation
!   01.03.22   Sawyer     ProTeX Documentation
!
!EOP

end module comspe




