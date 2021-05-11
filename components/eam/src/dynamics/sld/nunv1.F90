subroutine nunv1(lam     ,phi     ,lamdp   ,phidp   ,ud      , &
                 vd      ,coslat  ,grfu    ,grfv    ,grfulat , &
                 grfvlat ,nlon    )
!
!-----------------------------------------------------------------------
!
!  Purpose:
!  
!  Advection of vector terms (u/v prognostic eqn)
!  NOTE:  "gamma" factor is used to renormalize the length of the
!         arrival point vector
!
! Author:  J. Olson
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use scanslt,      only: plond
  implicit none

!------------------------------Arguments--------------------------------
!
  real(r8), intent(in)   :: lam    (plon)        ! vector of arrival point longitudes
  real(r8), intent(in)   :: phi                  ! arrival point latitude
  real(r8), intent(in)   :: lamdp  (plon,plev)   ! endpoint longitude coordinates
  real(r8), intent(in)   :: phidp  (plon,plev)   ! endpoint latitude  coordinates
  real(r8), intent(in)   :: ud     (plon,plev)   ! U(-) evaluated at departure point
  real(r8), intent(in)   :: vd     (plon,plev)   ! V(-) evaluated at departure point
  real(r8), intent(in)   :: coslat               ! cosine of latitude slice
  real(r8), intent(in)   :: grfu   (plon,plev)   ! Nu
  real(r8), intent(in)   :: grfv   (plon,plev)   ! Nv
  real(r8), intent(out)  :: grfulat(plon,plev)   ! Nu
  real(r8), intent(out)  :: grfvlat(plon,plev)   ! Nv
  integer , intent(in)   :: nlon                 ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  integer i,k                   ! indices
  real(r8) dl                   ! lam+ - lam-
  real(r8) sindl                ! sin (lam+ - lam-)
  real(r8) cosdl                ! cos (lam+ - lam-)
  real(r8) sinpa                ! sin (phi+)
  real(r8) sinpd                ! sin (phi-)
  real(r8) cospa                ! cos (phi+)
  real(r8) cospd                ! cos (phi-)
  real(r8) alphud               ! tmp space
  real(r8) alphvd               ! tmp space
  real(r8) betaud               ! tmp space
  real(r8) betavd               ! tmp space
  real(r8) termu                ! tmp space
  real(r8) termv                ! tmp space
  real(r8) gamma                ! tmp space
!
!-----------------------------------------------------------------------
!
  sinpa  = sin(phi)
  cospa  = cos(phi)
!
  do k = 1,plev
     do i = 1,nlon
        dl     = lam(i) - lamdp(i,k)
        sindl  = sin( dl )
        cosdl  = cos( dl )
        sinpd  = sin( phidp(i,k) )
        cospd  = cos( phidp(i,k) )
        alphud = cosdl
        alphvd = sinpd*sindl
        betaud = -sinpa*sindl
        betavd = cospd*cospa + sinpd*sinpa*cosdl
        termu  = alphud*ud(i,k) + alphvd*vd(i,k)
        termv  = betaud*ud(i,k) + betavd*vd(i,k)
        gamma  = (ud(i,k)*ud(i,k) + vd(i,k)*vd(i,k))/(termu*termu + termv*termv)
        gamma  = gamma**0.5_r8
        grfulat(i,k) = grfu(i,k) + gamma*termu*coslat
        grfvlat(i,k) = grfv(i,k) + gamma*termv*coslat
     end do
  end do
!
  return
end subroutine nunv1

