
module geopotential

!---------------------------------------------------------------------------------
! Compute geopotential from temperature or
! compute geopotential and temperature from dry static energy.
!
! The hydrostatic matrix elements must be consistent with the dynamics algorithm.
! The diagonal element is the itegration weight from interface k+1 to midpoint k.
! The offdiagonal element is the weight between interfaces.
! 
! Author: B.Boville, Feb 2001 from earlier code by Boville and S.J. Lin
!---------------------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use ppgrid,       only: pver, pverp
  use dycore,       only: dycore_is

  implicit none
  private
  save

  public geopotential_dse
  public geopotential_t

contains
!===============================================================================
  subroutine geopotential_dse(                                &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  ,  &
       dse    , q      , phis   , rair   , gravit , cpair  ,  &
       zvir   , t      , zi     , zm     , ncol             )
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the temperature  and geopotential height (above the surface) at the
! midpoints and interfaces from the input dry static energy and pressures.
!
!-----------------------------------------------------------------------
!------------------------------Arguments--------------------------------
!
! Input arguments
    integer, intent(in) :: ncol                  ! Number of longitudes

    ! rair, and cpair are passed in as slices of rank 3 arrays allocated
    ! at runtime. Don't specify size to avoid temporary copy.
    real(r8), intent(in) :: piln (:,:)    ! (pcols,pverp) - Log interface pressures
    real(r8), intent(in) :: pmln (:,:)    ! (pcols,pver)  - Log midpoint pressures
    real(r8), intent(in) :: pint (:,:)    ! (pcols,pverp) - Interface pressures
    real(r8), intent(in) :: pmid (:,:)    ! (pcols,pver)  - Midpoint pressures
    real(r8), intent(in) :: pdel (:,:)    ! (pcols,pver)  - layer thickness
    real(r8), intent(in) :: rpdel(:,:)    ! (pcols,pver)  - inverse of layer thickness
    real(r8), intent(in) :: dse  (:,:)    ! (pcols,pver)  - dry static energy
    real(r8), intent(in) :: q    (:,:)    ! (pcols,pver)  - specific humidity
    real(r8), intent(in) :: phis (:)      ! (pcols)       - surface geopotential
    real(r8), intent(in) :: rair (:,:)    !               - Gas constant for dry air
    real(r8), intent(in) :: gravit        !               - Acceleration of gravity
    real(r8), intent(in) :: cpair(:,:)    !               - specific heat at constant p for dry air
    real(r8), intent(in) :: zvir (:,:)    ! (pcols,pver)  - rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: t(:,:)       ! (pcols,pver)  - temperature
    real(r8), intent(out) :: zi(:,:)      ! (pcols,pverp) - Height above surface at interfaces
    real(r8), intent(out) :: zm(:,:)      ! (pcols,pver)  - Geopotential height at mid level
!
!---------------------------Local variables-----------------------------------------
!
    logical  :: fvdyn                   ! finite volume dynamics
    integer  :: i,k                     ! Lon, level, level indices
    real(r8) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(ncol)               ! off-diagonal element
    real(r8) :: rog(ncol,pver)          ! Rair / gravit
    real(r8) :: tv                      ! virtual temperature
    real(r8) :: tvfac                   ! Tv/T
!
!----------------------------------------------------------------------------------
    rog(:ncol,:) = rair(:ncol,:) / gravit

! Set dynamics flag
    fvdyn = dycore_is ('LR')

! The surface height is zero by definition.
    do i = 1,ncol
       zi(i,pverp) = 0.0_r8
    end do

! Compute the virtual temperature, zi, zm from bottom up
! Note, zi(i,k) is the interface above zm(i,k)
    do k = pver, 1, -1

! First set hydrostatic elements consistent with dynamics
       if (fvdyn) then
          do i = 1,ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1._r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          end do
       else
          do i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          end do
       end if

! Now compute tv, t, zm, zi
       do i = 1,ncol
          tvfac   = 1._r8 + zvir(i,k) * q(i,k)
          tv      = (dse(i,k) - phis(i) - gravit*zi(i,k+1)) / ((cpair(i,k) / tvfac) + &
	                                                               rair(i,k)*hkk(i))

          t (i,k) = tv / tvfac

          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_dse

!===============================================================================
  subroutine geopotential_t(                                 &
       piln   , pmln   , pint   , pmid   , pdel   , rpdel  , &
       t      , q      , rair   , gravit , zvir   ,          &
       zi     , zm     , ncol   )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute the geopotential height (above the surface) at the midpoints and 
! interfaces using the input temperatures and pressures.
!
!-----------------------------------------------------------------------

use ppgrid, only : pcols

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer, intent(in) :: ncol                  ! Number of longitudes

    real(r8), intent(in) :: piln (:,:)    ! (pcols,pverp) - Log interface pressures
    real(r8), intent(in) :: pmln (:,:)    ! (pcols,pver)  - Log midpoint pressures
    real(r8), intent(in) :: pint (:,:)    ! (pcols,pverp) - Interface pressures
    real(r8), intent(in) :: pmid (:,:)    ! (pcols,pver)  - Midpoint pressures
    real(r8), intent(in) :: pdel (:,:)    ! (pcols,pver)  - layer thickness
    real(r8), intent(in) :: rpdel(:,:)    ! (pcols,pver)  - inverse of layer thickness
    real(r8), intent(in) :: t    (:,:)    ! (pcols,pver)  - temperature
    real(r8), intent(in) :: q    (:,:)    ! (pcols,pver)  - specific humidity
    real(r8), intent(in) :: rair (:,:)    ! (pcols,pver)  - Gas constant for dry air
    real(r8), intent(in) :: gravit        !               - Acceleration of gravity
    real(r8), intent(in) :: zvir (:,:)    ! (pcols,pver)  - rh2o/rair - 1

! Output arguments

    real(r8), intent(out) :: zi(:,:)      ! (pcols,pverp) - Height above surface at interfaces
    real(r8), intent(out) :: zm(:,:)      ! (pcols,pver)  - Geopotential height at mid level
!
!---------------------------Local variables-----------------------------
!
    logical  :: fvdyn                   ! finite volume dynamics
    integer  :: i,k                     ! Lon, level indices
    real(r8) :: hkk(ncol)               ! diagonal element of hydrostatic matrix
    real(r8) :: hkl(ncol)               ! off-diagonal element
    real(r8) :: rog(ncol,pver)          ! Rair / gravit
    real(r8) :: tv                      ! virtual temperature
    real(r8) :: tvfac                   ! Tv/T
!
!-----------------------------------------------------------------------
!
    rog(:ncol,:) = rair(:ncol,:) / gravit

! Set dynamics flag

    fvdyn = dycore_is ('LR')

! The surface height is zero by definition.

    do i = 1,ncol
       zi(i,pverp) = 0.0_r8
    end do

! Compute zi, zm from bottom up. 
! Note, zi(i,k) is the interface above zm(i,k)

    do k = pver, 1, -1

! First set hydrostatic elements consistent with dynamics

       if (fvdyn) then
          do i = 1,ncol
             hkl(i) = piln(i,k+1) - piln(i,k)
             hkk(i) = 1._r8 - pint(i,k) * hkl(i) * rpdel(i,k)
          end do
       else
          do i = 1,ncol
             hkl(i) = pdel(i,k) / pmid(i,k)
             hkk(i) = 0.5_r8 * hkl(i)
          end do
       end if

! Now compute tv, zm, zi

       do i = 1,ncol
          tvfac   = 1._r8 + zvir(i,k) * q(i,k)
          tv      = t(i,k) * tvfac

          zm(i,k) = zi(i,k+1) + rog(i,k) * tv * hkk(i)
          zi(i,k) = zi(i,k+1) + rog(i,k) * tv * hkl(i)
       end do
    end do

    return
  end subroutine geopotential_t
end module geopotential
