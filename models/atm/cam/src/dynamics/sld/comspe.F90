
module comspe

!----------------------------------------------------------------------- 
! 
! Purpose: Spectral space arrays
! 
! Method: 
! 
! Author: CCM Core Group
! $Author$
! $Id$
! 
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid, only: plev, plat
  use pspect, only: psp, pmax, pmmax, pspt, pmaxp

  implicit none

  private

  public vz, d, t, q, alps, hs, hsnm, dsnm, dnm, vznm, lnpstar
  public a0nm, bmnm, bpnm, atri
  public btri, ctri, ncutoff, nalp, numm, maxm, lpspt, locm, locrm, lnstart
  public ncoefi, nm, nco2, nstart, nlen, alp, dalp

  real(r8) :: vz(psp,plev)        ! Vorticity spectral coefficients
  real(r8) :: d(psp,plev)         ! Divergence spectral coefficients
  real(r8) :: t(psp,plev)         ! Temperature spectral coefficients
  real(r8) :: q(psp,plev)         ! Moisture     spectral coefficients
  real(r8) :: alps(psp)           ! Log-pressure spectral coefficients
  real(r8) :: hs(psp,plev)        ! hydrostatic matrix for "real" atmosphere
  real(r8) :: hsnm(psp,plev)    ! vertical normal modes of "hs"
  real(r8) :: dsnm(psp,plev)    ! vertical normal modes of "ds"
  real(r8) :: dnm(psp,plev)     ! vertical normal modes of "d"
  real(r8) :: vznm(psp,plev)    ! vertical normal modes of "vz"
  real(r8) :: lnpstar(psp)      ! ln (Ps*) (SLD term; Ritchie & Tanguay, 1995)
  real(r8) :: a0nm(psp)         ! wave # coefs (use in vert normal mode space)
  real(r8) :: bmnm(psp)         ! wave # coefs (use in vert normal mode space)
  real(r8) :: bpnm(psp)         ! wave # coefs (use in vert normal mode space)
  real(r8) :: atri(psp)         ! wave # coefs (use in vert normal mode space)
  real(r8) :: btri(psp)         ! wave # coefs (use in vert normal mode space)
  real(r8) :: ctri(psp)         ! wave # coefs (use in vert normal mode space)

  integer :: ncutoff         = huge(1)   ! Break-even point for vector lengths in GRCALC
  integer :: nalp(pmax)      = huge(1)   ! Pointer into polynomial arrays
#if ( defined SPMD )
  integer :: maxm            = huge(1)  ! max number of Fourier wavenumbers per MPI task
  integer :: lpspt           = huge(1)  ! number of local spectral coefficients (NOT USED YET)
  integer, dimension(:), allocatable   :: numm
                                       ! number of Fourier wavenumbers owned per task
  integer, dimension(:,:), allocatable :: locm, locrm
                                       ! assignment of wavenumbers to MPI tasks
  integer, dimension(:), allocatable   :: lnstart 
                                       ! Starting indices for local spectral arrays (real) (NOT USED YET)
#else
  integer :: numm(0:0)       = pmmax
  integer :: maxm            = pmmax
  integer :: lpspt           = pspt
  integer :: locm(1:pmmax, 0:0) = huge(1)
  integer :: locrm(1:2*pmmax, 0:0) = huge(1)
  integer :: lnstart(1:pmmax) = huge(1)
#endif

  integer :: ncoefi(pmaxp)   = huge(1)   ! Pointer to start of coefficient diagonals
  integer :: nm(pmax)        = huge(1)   ! Number of coeffs stored on a given diagonal
  integer :: nco2(pmax)      = huge(1)   ! Complex form of ncoefi
  integer :: nstart(pmmax)   = huge(1)   ! Starting indices for spectral arrays (real)
  integer :: nlen(pmmax)     = huge(1)   ! Length vectors for spectral arrays

  real(r8) :: alp(pspt,plat/2)     ! Legendre polynomials
  real(r8) :: dalp(pspt,plat/2)    ! Legendre polynomial derivatives

end module comspe
