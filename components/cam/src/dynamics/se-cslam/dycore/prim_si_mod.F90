module prim_si_mod
  use shr_kind_mod,   only: r8=>shr_kind_r8

  implicit none
  private

  public :: preq_hydrostatic, geopotential_t
  public :: preq_pressure
contains


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!
!  CCM3 hydrostatic integral
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
  subroutine preq_hydrostatic(phi,phis,T_v,p,dp)
    use dimensions_mod, only : np, nlev
    use physconst, only: rair


    !------------------------------Arguments---------------------------------------------------------------
    real(kind=r8), intent(out) :: phi(np,np,nlev)
    real(kind=r8), intent(in) :: phis(np,np)
    real(kind=r8), intent(in) :: T_v(np,np,nlev)
    real(kind=r8), intent(in) :: p(np,np,nlev)
    real(kind=r8), intent(in) :: dp(np,np,nlev)
    !------------------------------------------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k                         ! longitude, level indices
    real(kind=r8) Hkk,Hkl          ! diagonal term of energy conversion matrix
    real(kind=r8), dimension(np,np,nlev) :: phii       ! Geopotential at interfaces
    !-----------------------------------------------------------------------

#if (defined COLUMN_OPENMP)
!$omp parallel do private(k,j,i,hkk,hkl)
#endif
       do j=1,np   !   Loop inversion (AAM)

          do i=1,np
             hkk = dp(i,j,nlev)*0.5_r8/p(i,j,nlev)
             hkl = 2*hkk
             phii(i,j,nlev)  = Rair*T_v(i,j,nlev)*hkl
             phi(i,j,nlev) = phis(i,j) + Rair*T_v(i,j,nlev)*hkk
          end do

          do k=nlev-1,2,-1
             do i=1,np
                ! hkk = dp*ckk
                hkk = dp(i,j,k)*0.5_r8/p(i,j,k)
                hkl = 2*hkk
                phii(i,j,k) = phii(i,j,k+1) + Rair*T_v(i,j,k)*hkl
                phi(i,j,k) = phis(i,j) + phii(i,j,k+1) + Rair*T_v(i,j,k)*hkk
             end do
          end do

          do i=1,np
             ! hkk = dp*ckk
             hkk = 0.5_r8*dp(i,j,1)/p(i,j,1)
             phi(i,j,1) = phis(i,j) + phii(i,j,2) + Rair*T_v(i,j,1)*hkk
          end do

       end do


end subroutine preq_hydrostatic



!
!  The hydrostatic routine from CAM physics.
!  (FV stuff removed)
!  t,q input changed to take t_v
!  removed gravit, so this routine returns PHI, not zm
subroutine geopotential_t(                                 &
       pmid   , pdel   ,  tv      , rair   ,  zm)

!-----------------------------------------------------------------------
!
! Purpose:
! Compute the geopotential height (above the surface) at the midpoints and
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------
    use dimensions_mod,     only : nlev, nlevp, np
    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments



    real(r8), intent(in) :: pmid (np*np,nlev)    ! Midpoint pressures
    real(r8), intent(in) :: pdel (np*np,nlev)    ! layer thickness
    real(r8), intent(in) :: tv    (np*np,nlev)    ! temperature
    real(r8), intent(in) :: rair                 ! Gas constant for dry air
    ! real(r8), intent(in) :: gravit               ! Acceleration of gravity
    ! real(r8), intent(in) :: zvir                 ! rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: zm(np*np,nlev)      ! Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
    integer :: ncol=np*np             ! Number of longitudes

    integer  :: i,k                ! Lon, level indices
    real(r8) :: hkk(np*np)         ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(np*np)         ! off-diagonal element
    real(r8) :: rog                ! Rair / gravit
    real(r8) :: zi(np*np,nlevp)     ! Height above surface at interfaces
!
!-----------------------------------------------------------------------
!
!    rog = rair/gravit
    rog = rair

! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,nlevp) = 0.0_r8
    end do

! Compute zi, zm from bottom up.
! Note, zi(i,k) is the interface above zm(i,k)
    do k = nlev, 1, -1
! First set hydrostatic elements consistent with dynamics
       do i = 1,ncol
          hkl(i) = pdel(i,k) / pmid(i,k)
          hkk(i) = 0.5_r8 * hkl(i)
       end do

! Now compute tv, zm, zi
       do i = 1,ncol
          ! tvfac   = 1._r8 + zvir * q(i,k)
          ! tv      = t(i,k) * tvfac
          zm(i,k) = zi(i,k+1) + rog * tv(i,k) * hkk(i)
          zi(i,k) = zi(i,k+1) + rog * tv(i,k) * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_t





!-----------------------------------------------------------------------
! preq_pressure:
!
! Purpose:
! Define the pressures of the interfaces and midpoints from the
! coordinate definitions and the surface pressure. Originally plevs0!
!
! Method:
!
! Author: B. Boville/ Adapted for HOMME by Rich Loft
!
!-----------------------------------------------------------------------
!
! $Id: prim_si_mod.F90,v 2.10 2005/10/14 20:17:22 jedwards Exp $
! $Author: jedwards $
!
!-----------------------------------------------------------------------

  subroutine preq_pressure (ps0,  ps,               &
       hyai, hybi, hyam, hybm, &
       pint, pmid, pdel)
    use dimensions_mod, only : np, nlev, nlevp
    implicit none

    !-----------------------------------------------------------------------

    real(kind=r8), intent(in)  :: ps0                ! Hybrid coordinate reference pressure (pascals)
    real(kind=r8), intent(in)  :: ps(np,np)          ! Surface pressure (pascals)
    real(kind=r8), intent(in)  :: hyai(nlevp)        ! Hybrid interface A coefficients
    real(kind=r8), intent(in)  :: hybi(nlevp)        ! Hybrid interface B coefficients
    real(kind=r8), intent(in)  :: hyam(nlev)         ! Hybrid midpoint  A coefficients
    real(kind=r8), intent(in)  :: hybm(nlev)         ! Hybrid midpoint  B coefficients
    real(kind=r8), intent(out) :: pint(np,np,nlevp)  ! Pressure at model interfaces
    real(kind=r8), intent(out) :: pmid(np,np,nlev)   ! Pressure at model levels
    real(kind=r8), intent(out) :: pdel(np,np,nlev)   ! Layer thickness (pint(k+1) - pint(k))
    !-----------------------------------------------------------------------

    !---------------------------Local workspace-----------------------------
    integer i,j,k             ! Horizontal, level indices
    !-----------------------------------------------------------------------
    !
    ! Set interface pressures
    !
    do k=1,nlevp
       do j=1,np
          do i=1,np
             pint(i,j,k) = hyai(k)*ps0 + hybi(k)*ps(i,j)
          end do
       end do
    end do
    !
    ! Set midpoint pressures and layer thicknesses
    !
    do k=1,nlev
       do j=1,np
          do i=1,np
             pmid(i,j,k) = hyam(k)*ps0 + hybm(k)*ps(i,j)
             pdel(i,j,k) = pint(i,j,k+1) - pint(i,j,k)
          end do
       end do
    end do

  end subroutine preq_pressure




end module prim_si_mod
