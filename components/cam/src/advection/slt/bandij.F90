
subroutine bandij(dlam    ,phib    ,lamp    ,phip    ,iband   , &
                  jband   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate longitude and latitude indices that identify the
! intervals on the extended grid that contain the departure points.
! Upon entry, all dep. points should be within jintmx intervals of the
! Northern- and Southern-most model latitudes. Note: the algorithm
! relies on certain relationships of the intervals in the Gaussian grid.
! 
! Method: 
!  dlam    Length of increment in equally spaced longitude grid (rad.)
!  phib    Latitude values for the extended grid.
!  lamp    Longitude coordinates of the points.  It is assumed that
!                     0.0 .le. lamp(i) .lt. 2*pi .
!  phip    Latitude coordinates of the points.
!  iband   Longitude index of the points.  This index points into
!          the extended arrays, e.g.,
!                   lam(iband(i)) .le. lamp(i) .lt. lam(iband(i)+1) .
!  jband   Latitude index of the points.  This index points into
!          the extended arrays, e.g.,
!                   phib(jband(i)) .le. phip(i) .lt. phib(jband(i)+1) .
! 
! Author: J. Olson
! 
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid,       only: plon, plev
  use scanslt,      only: platd, i1
  use rgrid,        only: fullgrid
  implicit none

!------------------------------Arguments--------------------------------
  real(r8), intent(in)  :: dlam(platd)        ! longitude increment
  real(r8), intent(in)  :: phib(platd)        ! latitude  coordinates of model grid
  real(r8), intent(in)  :: lamp(plon,plev)    ! longitude coordinates of dep. points
  real(r8), intent(in)  :: phip(plon,plev)    ! latitude  coordinates of dep. points
  integer , intent(in)  :: nlon               ! number of longitudes
  integer , intent(out) :: iband(plon,plev,4) ! longitude index of dep. points
  integer , intent(out) :: jband(plon,plev)   ! latitude  index of dep. points
!-----------------------------------------------------------------------
!
!---------------------------Local workspace-----------------------------
!
  integer i,j,k             ! indices
  real(r8) dphibr           ! reciprocal of an approximate del phi
  real(r8) phibs            ! latitude of southern-most latitude
  real(r8) rdlam(platd)     ! reciprocal of longitude increment
!
!-----------------------------------------------------------------------
!
  dphibr = 1._r8/( phib(platd/2+1) - phib(platd/2) )
  phibs  = phib(1)
  do j = 1,platd
     rdlam(j) = 1._r8/dlam(j)
  end do
!
! Loop over level and longitude

!$OMP PARALLEL DO PRIVATE (K, I)
  do k=1,plev
     do i = 1,nlon
!
! Latitude indices.
!
        jband(i,k) = int ( (phip(i,k) - phibs)*dphibr + 1._r8 )
        if( phip(i,k) >= phib(jband(i,k)+1) ) then
           jband(i,k) = jband(i,k) + 1
        end if
!
! Longitude indices.
!
        iband(i,k,1) = i1 + int( lamp(i,k)*rdlam(jband(i,k)-1))
        if (fullgrid) then
           iband(i,k,2) = iband(i,k,1)
           iband(i,k,3) = iband(i,k,1)
           iband(i,k,4) = iband(i,k,1)
        else
           iband(i,k,2) = i1 + int( lamp(i,k)*rdlam(jband(i,k)  ))
           iband(i,k,3) = i1 + int( lamp(i,k)*rdlam(jband(i,k)+1))
           iband(i,k,4) = i1 + int( lamp(i,k)*rdlam(jband(i,k)+2))
        end if
     end do
  end do

  return
end subroutine bandij
