module gw_diffusion

!
! This module contains code computing the effective diffusion of
! constituents and dry static energy due to gravity wave breaking.
!

use gw_utils, only: r8
use vdiff_lu_solver, only: lu_decomp

implicit none
private
save

public :: gw_ediff
public :: gw_diff_tend

contains

!==========================================================================

subroutine gw_ediff(ncol, pver, ngwv, kbot, ktop, tend_level, &
     gwut, ubm, nm, rho, dt, gravit, pmid, rdpm, c, &
     egwdffi, decomp)
!
! Calculate effective diffusivity associated with GW forcing.
!
! Author: F. Sassi, Jan 31, 2001
!
  use gw_utils, only: midpoint_interp
  use vdiff_lu_solver, only: vd_lu_decomp

!-------------------------------Input Arguments----------------------------

  ! Column, level, and gravity wave spectrum dimensions.
  integer, intent(in) :: ncol, pver, ngwv
  ! Bottom and top levels to operate on.
  integer, intent(in) :: kbot, ktop
  ! Per-column bottom index where tendencies are applied.
  integer, intent(in) :: tend_level(ncol)
  ! GW zonal wind tendencies at midpoint.
  real(r8), intent(in) :: gwut(ncol,pver,-ngwv:ngwv)
  ! Projection of wind at midpoints.
  real(r8), intent(in) :: ubm(ncol,pver)
  ! Brunt-Vaisalla frequency.
  real(r8), intent(in) :: nm(ncol,pver)

  ! Density at interfaces.
  real(r8), intent(in) :: rho(ncol,pver+1)
  ! Time step.
  real(r8), intent(in) :: dt
  ! Acceleration due to gravity.
  real(r8), intent(in) :: gravit
  ! Midpoint pressure.
  real(r8), intent(in) :: pmid(ncol,pver)
  ! 1./pdel (pdel is thickness between interfaces).
  real(r8), intent(in) :: rdpm(ncol,pver)
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-ngwv:ngwv)

!-----------------------------Output Arguments-----------------------------
  ! Effective gw diffusivity at interfaces.
  real(r8), intent(out) :: egwdffi(ncol,0:pver)
  ! LU decomposition.
  type(lu_decomp), intent(out) :: decomp

!-----------------------------Local Workspace------------------------------

  ! 0._r8 (array to pass to vd_lu_decomp).
  real(r8) :: zero(ncol)
  ! Effective gw diffusivity at midpoints.
  real(r8) :: egwdffm(ncol,pver)
  ! dt * (gravit*rho)^2/dp
  real(r8) :: tmpi2(ncol,pver+1)
  ! Level and wave indices.
  integer :: k, l
  ! Inverse Prandtl number.
  real(r8), parameter :: prndl=0.25_r8
  ! Density scale height.
  real(r8), parameter :: dscale=7000._r8

!--------------------------------------------------------------------------

  zero = 0.0_r8

  egwdffi = 0._r8
  egwdffm = 0._r8

  ! Calculate effective diffusivity at midpoints.
  do l = -ngwv, ngwv
     do k = ktop+1, kbot

        egwdffm(:,k) = egwdffm(:,k) + &
             prndl * 0.5_r8 * gwut(:,k,l) * (c(:,l)-ubm(:,k)) / nm(:,k)**2

     end do
  end do


  ! Interpolate effective diffusivity to interfaces.
  ! Assume zero at top and bottom interfaces.
  egwdffi(:,ktop+1:kbot-1) = midpoint_interp(egwdffm(:,ktop+1:kbot))

  ! Limit diffusivity to some reasonable value.
  egwdffi = min(150._r8, egwdffi)

  ! Do not calculate diffusivities below level where tendencies are
  ! actually allowed.
  do k= ktop, kbot
     where (k >= tend_level) egwdffi(:,k) = 0.0_r8
  enddo

  ! Calculate dt * (gravit*rho)^2/dp at interior interfaces.
  tmpi2 = 0.0_r8
  do k = ktop+2 , kbot+1
     tmpi2(:,k) = dt * (gravit * rho(:,k))**2 / (pmid(:,k) - pmid(:,k-1))
  end do

  ! Decompose the diffusion matrix.
  ! Note that [ktop,kbot] are model interfaces (beginning at zero), whereas
  ! in vd_lu_decomp they are expected as midpoints.
  call vd_lu_decomp(ncol, pver, ncol, &
       zero, egwdffi, tmpi2, rdpm, dt, gravit, zero, ktop+1, kbot+1, decomp)

end subroutine gw_ediff

!==========================================================================

subroutine gw_diff_tend(ncol, pver, kbot, ktop, q, dt, decomp, dq)

!
! Calculates tendencies from effective diffusion due to gravity wave
! breaking.
!
! Method:
! A constituent flux on interfaces is given by:
!
!              rho * (w'q') = rho * Deff qz
!
! where (all evaluated on interfaces):
!
!        rho   = density
!        qz    = constituent vertical gradient
!        Deff  = effective diffusivity
!
! An effective diffusivity is calculated by adding up the diffusivities
! from all waves (see gw_ediff). The tendency is calculated by invoking LU
! decomposition and solving as for a regular diffusion equation.
!
! Author: Sassi - Jan 2001
!--------------------------------------------------------------------------

  use vdiff_lu_solver, only: vd_lu_solve

!---------------------------Input Arguments--------------------------------

  ! Column and level dimensions.
  integer, intent(in) :: ncol, pver
  ! Bottom and top levels to operate on.
  integer, intent(in) :: kbot, ktop

  ! Constituent to diffuse.
  real(r8), intent(in) :: q(ncol,pver)
  ! Time step.
  real(r8), intent(in) :: dt

  ! LU decomposition.
  type(lu_decomp), intent(in) :: decomp

!--------------------------Output Arguments--------------------------------

  ! Constituent tendencies.
  real(r8), intent(out) :: dq(ncol,pver)

!--------------------------Local Workspace---------------------------------

  ! Temporary storage for constituent.
  real(r8) :: qnew(ncol,pver)
  ! 0._r8 (array to pass to vd_lu_decomp).
  real(r8) :: zero(ncol)

!--------------------------------------------------------------------------
  zero = 0.0_r8
  dq   = 0.0_r8
  qnew = q

  ! Solve the diffusion matrix.
  call vd_lu_solve(ncol, pver, ncol, qnew, decomp, ktop+1, kbot+1, zero)

  ! Evaluate tendency to be reported back.
  dq = (qnew-q) / dt

end subroutine gw_diff_tend

end module gw_diffusion
