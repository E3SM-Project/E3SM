
subroutine sphdep(jcen    ,jgc     ,dt      ,ra      ,iterdp  , &
                  locgeo  ,ub      ,uxl     ,uxr     ,lam     , &
                  phib    ,lbasiy  ,lammp   ,phimp   ,lamdp   , &
                  phidp   ,idp     ,jdp     ,vb      ,vxl     , &
                  vxr     ,nlon    ,nlonex  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Compute departure points for semi-Lagrangian transport on surface of
! sphere using midpoint quadrature.  Computations are done in:
!
!   1) "local geodesic"   coordinates for "locgeo" = .true.
!   2) "global spherical" coordinates for "locgeo" = .false.
! 
! Method: 
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
  use pmgrid,       only: plev, plon, plat
  use scanslt,      only: platd, plond, beglatex, endlatex, i1, nxpt, j1
  use abortutils,   only: endrun
  use cam_logfile,  only: iulog
  implicit none
#include <parslt.h>

!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon              ! longitude dimension
  integer , intent(in) :: nlonex(platd)     ! extended longitude dimension
  integer , intent(in) :: jcen              ! index of lat slice (extnd)
  integer , intent(in) :: jgc               ! index of lat slice (model)
  real(r8), intent(in) :: dt                ! time step (seconds)
  real(r8), intent(in) :: ra                ! 1./(radius of earth)
  integer , intent(in) :: iterdp            ! number of iterations
  logical , intent(in) :: locgeo            ! computation type flag
  real(r8), intent(in) :: ub (plond,plev,beglatex:endlatex) ! x-deriv
  real(r8), intent(in) :: vb (plond,plev,beglatex:endlatex) ! x-deriv
  real(r8), intent(in) :: uxl(plond,plev,beglatex:endlatex) ! left  x-deriv (u)
  real(r8), intent(in) :: uxr(plond,plev,beglatex:endlatex) ! right x-deriv 
  real(r8), intent(in) :: vxl(plond,plev,beglatex:endlatex) ! left  x-deriv (v)
  real(r8), intent(in) :: vxr(plond,plev,beglatex:endlatex) ! right x-deriv
  real(r8), intent(in) :: lam(plond,platd)  ! long. coord. of model grid
  real(r8), intent(in) :: phib(platd)       ! lat.  coord. of model grid
  real(r8), intent(in) :: lbasiy(4,2,platd) ! lat interpolation weights
  real(r8), intent(inout) :: lammp(plon,plev)  ! long coord of midpoint
  real(r8), intent(inout) :: phimp(plon,plev)  ! lat  coord of midpoint
  real(r8), intent(out)   :: lamdp(plon,plev)  ! long coord of dep. point
  real(r8), intent(out)   :: phidp(plon,plev)  ! lat  coord of dep. point
  integer , intent(out)   :: idp(plon,plev,4)  ! long index of dep. point
  integer , intent(out)   :: jdp(plon,plev)    ! lat  index of dep. point
!
!  jcen    Index in extended grid corresponding to latitude being
!          forecast.
!  jgc     Index in model    grid corresponding to latitude being
!          forecast.
!  dt      Time interval that parameterizes the parcel trajectory.
!  ra      Reciprocal of radius of earth.
!  iterdp  Number of iterations used for departure point calculation.
!  locgeo  Logical flag to indicate computation in "local geodesic" or
!          "global spherical" space.
!  ub      Longitudinal velocity components in spherical coordinates.
!  uxl     x-derivatives of u at the left  (west) edge of given interval
!  vxl     x-derivatives of v at the left  (west) edge of given interval
!  uxr     x-derivatives of u at the right (east) edge of given interval
!  vxr     x-derivatives of v at the right (east) edge of given interval
!  lam     Longitude values for the extended grid.
!  phib    Latitude  values for the extended grid.
!  lbasiy  Weights for Lagrange cubic interpolation on the unequally
!          spaced latitude grid.
!  lammp   Longitude coordinates of the trajectory mid-points of the
!          parcels that correspond to the global grid points contained
!          in the latitude slice being forecast.  On entry lammp
!          is an initial guess.
!  phimp   Latitude coordinates of the trajectory mid-points of the
!          parcels that correspond to the global grid points contained
!          in the latitude slice being forecast.  On entry phimp
!          is an initial guess.
!  lamdp   Longitude coordinates of the departure points that correspond
!          to the global grid points contained in the latitude slice
!          being forecast.  lamdp is constrained so that
!                     0.0 .le. lamdp(i) .lt. 2*pi .
!  phidp   Latitude coordinates of the departure points that correspond
!          to the global grid points contained in the latitude slice
!          being forecast.  If phidp is computed outside the latitudinal
!          domain of the extended grid, then an abort will be called by
!          subroutine "trjgl".
!  idp     Longitude index of departure points.  This index points into
!          the extended arrays, e.g.,
!              lam (idp(i,k)) .le. lamdp(i,k) .lt. lam (idp(i,k)+1).
!  jdp     Latitude  index of departure points.  This index points into
!          the extended arrays, e.g.,
!              phib(jdp(i,k)) .le. phidp(i,k) .lt. phib(jdp(i,k)+1).
!-----------------------------------------------------------------------

 !------------------------ local variables ------------------------------
  integer iter                     ! index
  integer i, j, k                  ! indices
  integer imax, imin, kmin, kmax   ! indices
  real(r8) finc                    ! time step factor
  real(r8) dttmp                   ! time step (seconds)
  real(r8) dlam(platd)             ! increment of grid in x-direction
  real(r8) phicen                  ! latitude coord of current lat slice
  real(r8) cphic                   ! cos(phicen)
  real(r8) sphic                   ! sin(phicen)
  real(r8) upr  (plon,plev)        ! u in local geodesic coords
  real(r8) vpr  (plon,plev)        ! v in local geodesic coords
  real(r8) lampr(plon,plev)        ! relative long coord of dep pt
  real(r8) phipr(plon,plev)        ! relative lat  coord of dep pt
  real(r8) uvmp (plon,plev,2)      ! u/v (spherical) interpltd to dep pt
  real(r8) fint (plon,plev,ppdy,2) ! u/v x-interpolants
  real(r8) phidpmax
  real(r8) phidpmin
  real(r8) phimpmax
  real(r8) phimpmin
!-----------------------------------------------------------------------
!
  do j=1,platd
     dlam(j) = lam(nxpt+2,j) - lam(nxpt+1,j)
  end do
  phicen = phib(jcen)
  cphic  = cos( phicen )
  sphic  = sin( phicen )
!
! Convert latitude coordinates of trajectory midpoints from spherical
! to local geodesic basis.
!
  if( locgeo ) call s2gphi(lam(i1,jcen) ,cphic   ,sphic   ,lammp   ,phimp, &
                           phipr        ,nlon    )
!
! Loop over departure point iterates.
!
  do 30 iter = 1,iterdp
!
! Compute midpoint indicies.
!
     call bandij(dlam    ,phib    ,lammp   ,phimp   ,idp     , &
                 jdp     ,nlon    )
!
! Hermite cubic interpolation to the x-coordinate of each
! departure point at each y-coordinate required to compute the
! y-interpolants.
!
     call herxin(1      ,1       ,ub      ,uxl   ,uxr          , &
                 lam    ,lammp   ,idp     ,jdp   ,fint(1,1,1,1), &
                 nlon   ,nlonex  )

     call herxin(1      ,1       ,vb      ,vxl   ,vxr          , &
                 lam    ,lammp   ,idp     ,jdp   ,fint(1,1,1,2), &
                 nlon   ,nlonex  )

     call lagyin(2       ,fint     ,lbasiy  ,phimp   ,jdp     , &
                 jcen    ,uvmp     ,nlon    )
!
! Put u/v on unit sphere
!
!$OMP PARALLEL DO PRIVATE (K, I)
     do k = 1,plev
        do i = 1,nlon
           uvmp(i,k,1) = uvmp(i,k,1)*ra
           uvmp(i,k,2) = uvmp(i,k,2)*ra
        end do
     end do
!
! For local geodesic:
!
!   a) Convert velocity coordinates at trajectory midpoints from
!      spherical coordinates to local geodesic coordinates,
!   b) Estimate midpoint parcel trajectory,
!   c) Convert back to spherical coordinates
!
! Else, for global spherical
!
!   Estimate midpoint trajectory with no conversions
!
     if ( locgeo ) then
        call s2gvel(uvmp(1,1,1),uvmp(1,1,2)   ,lam(i1,jcen) ,cphic ,sphic   , &
                    lammp      ,phimp         ,upr          ,vpr   ,nlon    )
        call trajmp(dt      ,upr     ,vpr     ,phipr   ,lampr   , &
                    nlon    )
        dttmp = 0.5_r8*dt
        call g2spos(dttmp   ,lam(i1,jcen) ,phib    ,phicen  ,cphic , &
                    sphic   ,upr          ,vpr     ,lampr   ,phipr , &
                    lammp   ,phimp        ,nlon    )
     else
        call trjmps(dt      ,uvmp(1,1,1) ,uvmp(1,1,2),  phimp   ,lampr   , &
                    phipr   ,nlon    )
        finc = 1._r8
        call trjgl (finc    ,phicen  ,lam(i1,jcen) ,lampr   ,phipr , &
                    lammp   ,phimp   ,nlon    )
     end if
!
! Test that the latitudinal extent of trajectory is NOT over the poles
! Distributed memory case: check that the latitudinal extent of the 
! trajectory is not more than "jintmx" gridpoints away.
!
     phimpmax = -1.e36_r8
     phimpmin =  1.e36_r8
     do k=1,plev
        do i=1,nlon
           if (phimp(i,k)>phimpmax) then
              phimpmax = phimp(i,k)
              imax = i
              kmax = k
           end if
           if (phimp(i,k)<phimpmin) then
              phimpmin = phimp(i,k)
              imin = i
              kmin = k
           end if
        end do
     end do
#if ( defined SPMD )
     if ( phimp(imax,kmax) >= phib(endlatex-nxpt) ) then
#else
     if ( phimp(imax,kmax) >= phib(j1+plat) ) then
#endif
        write(iulog,*)'SPHDEP: ****** MODEL IS BLOWING UP:  CFL condition likely violated *********'
        write(iulog,9000) imax,kmax,jgc
        write(iulog,*)' Possible solutions:  a)  reduce time step'
        write(iulog,*)'                      b)  if initial run, set "DIVDAMPN = 1." in namelist and rerun'
        write(iulog,*)'                      c)  modified code may be in error'
        call endrun
#if ( defined SPMD )
     else if( phimp(imin,kmin) < phib(beglatex+nxpt) ) then
#else
     else if( phimp(imin,kmin) < phib(j1-1) ) then
#endif
        write(iulog,*)'SPHDEP: ****** MODEL IS BLOWING UP:  CFL condition likely violated *********'
        write(iulog,9000) imin,kmin,jgc
        write(iulog,*)' Possible solutions:  a)  reduce time step'
        write(iulog,*)'                      b)  if initial run, set "DIVDAMPN = 1." in namelist and rerun'
        write(iulog,*)'                      c)  modified code may be in error'
        call endrun
     end if

30 continue          ! End of iter=1,iterdp loop
!
! Compute departure points in geodesic coordinates, and convert back
! to spherical coordinates.
!
! Else, compute departure points directly in spherical coordinates
!
     if (locgeo) then
!$OMP PARALLEL DO PRIVATE (K, I)
        do k = 1,plev
           do i = 1,nlon
              lampr(i,k) = 2._r8*lampr(i,k)
              phipr(i,k) = 2._r8*phipr(i,k)
           end do
        end do
        dttmp = dt
        call g2spos(dttmp   ,lam(i1,jcen) ,phib    ,phicen  ,cphic   , &
                    sphic   ,upr          ,vpr     ,lampr   ,phipr   , &
                    lamdp   ,phidp        ,nlon    )
     else
        finc = 2._r8
        call trjgl (finc    ,phicen  ,lam(i1,jcen) ,lampr   ,phipr   , &
                    lamdp   ,phidp   ,nlon    )
     end if
!
! Test that the latitudinal extent of trajectory is NOT over the poles
! Distributed memory case: check that the latitudinal extent of the 
! trajectory is not more than "jintmx" gridpoints away.
!
     phidpmax = -1.e36_r8
     phidpmin =  1.e36_r8
     do k=1,plev
        do i=1,nlon
           if (phidp(i,k)>phidpmax) then
              phidpmax = phidp(i,k)
              imax = i
              kmax = k
           end if
           if (phidp(i,k)<phidpmin) then
              phidpmin = phidp(i,k)
              imin = i
              kmin = k
           end if
        end do
     end do
#if ( defined SPMD )
     if ( phidp(imax,kmax) >= phib(endlatex-nxpt) ) then
#else
     if ( phidp(imax,kmax) >= phib(j1+plat) ) then
#endif
        write(iulog,*)'SPHDEP: ****** MODEL IS BLOWING UP:  CFL condition likely violated *********'
        write(iulog,9000) imax,kmax,jgc
        write(iulog,*)' Possible solutions:  a)  reduce time step'
        write(iulog,*)'                      b)  if initial run, set "DIVDAMPN = 1." in namelist and rerun'
        write(iulog,*)'                      c)  modified code may be in error'
        call endrun
#if ( defined SPMD )
     else if( phidp(imin,kmin) < phib(beglatex+nxpt) ) then
#else
     else if( phidp(imin,kmin) < phib(j1-1) ) then
#endif
        write(iulog,*)'SPHDEP: ****** MODEL IS BLOWING UP:  CFL condition likely violated *********'
        write(iulog,9000) imin,kmin,jgc
        write(iulog,*)' Possible solutions:  a)  reduce time step'
        write(iulog,*)'                      b)  if initial run, set "DIVDAMPN = 1." in namelist and rerun'
        write(iulog,*)'                      c)  modified code may be in error'
        call endrun
     end if
!
! Compute departure point indicies.
!
     call bandij(dlam    ,phib    ,lamdp   ,phidp   ,idp     , &
                 jdp     ,nlon    )

9000 format(//'Parcel associated with longitude ',i5,', level ',i5, &
          ' and latitude ',i5,' is outside the model domain.')

  return
end subroutine sphdep

!============================================================================================

subroutine g2spos(dttmp   ,lam     ,phib    ,phi     ,cosphi  , &
                  sinphi  ,upr     ,vpr     ,lamgc   ,phigc   , &
                  lamsc   ,phisc   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Transform position coordinates for a set of points, each of which is
! associated with a grid point in a global latitude slice, from local
! geodesic to spherical coordinates.
! 
! Method: 
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
  use pmgrid,       only: plon, plev, plat
  use scanslt,      only: plond1, platd, j1
  implicit none

!------------------------------Arguments--------------------------------
  real(r8), intent(in) :: dttmp                ! time step
  real(r8), intent(in) :: lam(plond1)          ! model longitude coordinates
  real(r8), intent(in) :: phib(platd)          ! extended grid latitude coordinates
  real(r8), intent(in) :: phi                  ! current latitude coordinate (radians)
  real(r8), intent(in) :: cosphi               ! cos of current latitude
  real(r8), intent(in) :: sinphi               ! sin of current latitude
  real(r8), intent(in) :: upr  (plon,plev)     ! u-wind in geodesic coord
  real(r8), intent(in) :: vpr  (plon,plev)     ! v-wind in geodesic coord
  real(r8), intent(in) :: lamgc(plon,plev)     ! geodesic long coord. of dep. point
  real(r8), intent(in) :: phigc(plon,plev)     ! geodesic lat  coord. of dep. point
  integer , intent(in) :: nlon                 ! longitude dimension 
  real(r8), intent(out):: lamsc(plon,plev)     ! spherical long coord. of dep. point
  real(r8), intent(out):: phisc(plon,plev)     ! spherical lat  coord. of dep. point
!
!
!  dttmp  Time step over which midpoint/endpoint trajectory is
!         calculated (seconds).
!  lam    Longitude coordinates of the global grid points in spherical
!         system.  The grid points in the global array are the reference
!         points for the local geodesic systems.
!  phib   Latitude values for the extended grid.
!  phi    Latitude coordinate (in the global grid) of the current
!         latitude slice.
!  cosphi cos( phi )
!  sinphi sin( phi )
!  upr    zonal      velocity at departure point in local geodesic coord
!  vpr    Meridional velocity at departure point in local geodesic coord
!  lamgc  Longitude coordinate of points in geodesic coordinates.
!  phigc  Latitude coordinate of points in geodesic coordinates.
!  lamsc  Longitude coordinate of points in spherical coordinates.
!  phisc  Latitude coordinate of points in spherical coordinates.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,ii,k                ! indices
  integer nval(plev)            ! number of values returned from whenfgt
  integer indx(plon,plev)       ! index holder
  real(r8) pi                   ! 4.*atan(1.)
  real(r8) twopi                ! 2.*pi
  real(r8) pi2                  ! pi/2
  real(r8) sgnphi               ! holds sign of phi
  real(r8) sphigc               ! sin(phigc)
  real(r8) cphigc               ! cos(phigc)
  real(r8) clamgc               ! cos(lamgc)
  real(r8) slam2                ! sin(lamgc)**2
  real(r8) phipi2               ! tmp variable
  real(r8) slamgc(plon,plev)    ! sin(lamgc)
  real(r8) dlam(plon,plev)      ! zonal extent of trajectory
  real(r8) coeff                ! tmp variable
  real(r8) distmx               ! max distance
  real(r8) dist(plon,plev)      ! approx. distance traveled along traj.
  real(r8) fac                  ! 1. - 10*eps, eps from mach. precision
  integer s_nval
!-----------------------------------------------------------------------
!
  fac    = 1._r8 - 10._r8*epsilon (fac)
  pi     = 4._r8*atan(1._r8)
  twopi  = pi*2._r8
  pi2    = pi/2._r8
  coeff  = (1.1_r8*dttmp)**2
  distmx = (sign(pi2,phi) - phi)**2/coeff
  sgnphi = sign( 1._r8, phi )

!$OMP PARALLEL DO PRIVATE (K, I, SPHIGC, CPHIGC, CLAMGC, S_NVAL)
  do k=1,plev
     do i=1,nlon
        sphigc      = sin( phigc(i,k) )
        cphigc      = cos( phigc(i,k) )
        slamgc(i,k) = sin( lamgc(i,k) )
        clamgc      = cos( lamgc(i,k) )
        phisc(i,k)  = asin((sphigc*cosphi + cphigc*sinphi*clamgc)*fac)
        if ( abs(phisc(i,k)) .ge. phib(j1+plat)*fac ) then
           phisc(i,k) = sign( phib(j1+plat),phisc(i,k) )*fac
        end if
        dlam(i,k) = asin((slamgc(i,k)*cphigc/cos(phisc(i,k)))*fac)
!
! Compute estimated trajectory distance based upon winds alone
!
        dist(i,k) = upr(i,k)**2 + vpr(i,k)**2
     end do
!
! Determine which trajectories may have crossed over pole
!
     s_nval = 0
     do i=1,nlon
        if (dist(i,k) > distmx) then
           s_nval = s_nval + 1
           indx(s_nval,k) = i
        end if
     end do
     nval(k) = s_nval
  end do
!
! Check that proper branch of arcsine is used for calculation of
! dlam for those trajectories which may have crossed over pole.
!
!$OMP PARALLEL DO PRIVATE (K, II, I, SLAM2, PHIPI2)
  do k=1,plev
     do ii=1,nval(k)
        i = indx(ii,k)
        slam2 = slamgc(i,k)**2
        phipi2 = asin((sqrt((slam2 - 1._r8)/(slam2 - 1._r8/cosphi**2)))*fac)
        if (sgnphi*phigc(i,k) > phipi2) then
           dlam(i,k) = sign(pi,lamgc(i,k)) - dlam(i,k)
        end if
     end do

     do i=1,nlon
        lamsc(i,k) = lam(i) + dlam(i,k)
!
! Restrict longitude to be in the range [0, twopi).
!  
        if( lamsc(i,k) >= twopi ) lamsc(i,k) = lamsc(i,k) - twopi
        if( lamsc(i,k) < 0.0_r8    ) lamsc(i,k) = lamsc(i,k) + twopi
     end do
  end do

  return
end subroutine g2spos

!============================================================================================

subroutine s2gphi(lam     ,cosphi  ,sinphi  ,lamsc   ,phisc   , &
                  phigc   ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Calculate transformed local geodesic latitude coordinates for a set
! of points, each of which is associated with a grid point in a global
! latitude slice. Transformation is spherical to local geodesic.
! (Williamson and Rasch, 1991)
! 
! Method: 
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
  use scanslt,      only: plond1
  implicit none

!------------------------------Arguments--------------------------------
  real(r8), intent(in)  :: lam(plond1)          ! long coordinates of model grid
  real(r8), intent(in)  :: cosphi               ! cos(latitude)
  real(r8), intent(in)  :: sinphi               ! sin(latitude)
  real(r8), intent(in)  :: lamsc(plon,plev)     ! spher. long coords of dep points
  real(r8), intent(in)  :: phisc(plon,plev)     ! spher. lat  coords of dep points
  integer , intent(in)  :: nlon                 ! longitude dimension 
  real(r8), intent(out) :: phigc(plon,plev)     ! loc geod. lat coords of dep points
!
! lam    longitude coordinates of the global grid points in spherical
!        system. The grid points in the global array are the reference
!        points for the local geodesic systems.
! cosphi cosine of the latitude of the global latitude slice.
! sinphi sine of the latitude of the global latitude slice.
! lamsc  longitude coordinate of dep. points in spherical coordinates.
! phisc  latitude  coordinate of dep. points in spherical coordinates.
! phigc  latitude  coordinate of dep. points in local geodesic coords.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,k                   ! longitude, level indices
  real(r8) sphisc               ! |
  real(r8) cphisc               ! | -- temporary variables
  real(r8) clamsc               ! |
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE (K, I, SPHISC, CPHISC, CLAMSC)
  do k = 1,plev
     do i = 1,nlon
        sphisc = sin( phisc(i,k) )
        cphisc = cos( phisc(i,k) )
        clamsc = cos( lam(i) - lamsc(i,k) )
        phigc(i,k) = asin( sphisc*cosphi - cphisc*sinphi*clamsc )
     end do
  end do

  return
end subroutine s2gphi

!============================================================================================

subroutine s2gvel(udp     ,vdp     ,lam     ,cosphi  ,sinphi  , &
                  lamdp   ,phidp   ,upr     ,vpr     ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Transform velocity components at departure points associated with a
! single latitude slice from spherical coordinates to local geodesic
! coordinates. (Williamson and Rasch, 1991)
! 
! Method: 
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
  use scanslt,      only: plond1
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon                 ! longitude dimension 
  real(r8), intent(in) :: udp(plon,plev)       ! u in spherical coords.
  real(r8), intent(in) :: vdp(plon,plev)       ! v in spherical coords.
  real(r8), intent(in) :: lam(plond1)          ! x-coordinates of model grid
  real(r8), intent(in) :: cosphi               ! cos(latitude)
  real(r8), intent(in) :: sinphi               ! sin(latitude)
  real(r8), intent(in) :: lamdp(plon,plev)     ! spherical longitude coord of dep pt.
  real(r8), intent(in) :: phidp(plon,plev)     ! spherical latitude  coord of dep pt.
  real(r8), intent(out) :: upr(plon,plev)      ! u in local geodesic coords.
  real(r8), intent(out) :: vpr(plon,plev)      ! v in local geodesic coords.
!
! udp    u-component of departure point velocity in spherical coords.
! vdp    v-component of departure point velocity in spherical coords.
! lam    Longitude of arrival point position (model grid point) in spherical coordinates.
! cosphi Cos of latitude of arrival point positions (model grid pt).
! sinphi Sin of latitude of arrival point positions (model grid pt).
! lamdp  Longitude of departure point position in spherical coordinates.
! phidp  Latitude  of departure point position in spherical coordinates.
! upr    u-component of departure point velocity in geodesic coords.
! vpr    v-component of departure point velocity in geodesic coords.
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,k                   ! longitude, level indices
  real(r8) cdlam                ! |
  real(r8) clamp                ! |
  real(r8) cphid                ! |
  real(r8) cphip                ! |
  real(r8) dlam                 ! | -- temporary variables
  real(r8) sdlam                ! |
  real(r8) slamp                ! |
  real(r8) sphid                ! |
  real(r8) sphip                ! |
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE (K, I, DLAM, SDLAM, CDLAM, SPHID, CPHID, SPHIP, &
!$OMP                      CPHIP, SLAMP, CLAMP)
  do k = 1,plev
     do i = 1,nlon
        dlam  = lam(i) - lamdp(i,k)
        sdlam = sin( dlam )
        cdlam = cos( dlam )
        sphid = sin( phidp(i,k) )
        cphid = cos( phidp(i,k) )
        sphip = sphid*cosphi - cphid*sinphi*cdlam
        cphip = cos( asin( sphip ) )
        slamp = -sdlam*cphid/cphip
        clamp = cos( asin( slamp ) )
        vpr(i,k) = (vdp(i,k)*(cphid*cosphi + sphid*sinphi*cdlam) - &
                    udp(i,k)*sinphi*sdlam)/cphip
        upr(i,k) = (udp(i,k)*cdlam + vdp(i,k)*sphid*sdlam + &
                    vpr(i,k)*slamp*sphip)/clamp
     end do
  end do

  return
end subroutine s2gvel

!============================================================================================

subroutine trajmp(dt      ,upr     ,vpr     ,phipr   ,lampr   , &
                  nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Estimate mid-point of parcel trajectory (geodesic coordinates) based
! upon horizontal wind field.
! 
! Method: 
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
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in)    :: nlon             ! longitude dimension
  real(r8), intent(in)    :: dt               ! time step (seconds)
  real(r8), intent(in)    :: upr(plon,plev)   ! u-component of wind in local geodesic
  real(r8), intent(in)    :: vpr(plon,plev)   ! v-component of wind in local geodesic
  real(r8), intent(inout) :: phipr(plon,plev) ! latitude coord of trajectory mid-point
  real(r8), intent(out)   :: lampr(plon,plev) ! longitude coord of traj. mid-point
!
!  dt      Time interval that corresponds to the parcel trajectory.
!  upr     u-coordinate of velocity corresponding to the most recent
!          estimate of the trajectory mid-point (in geodesic system).
!  vpr     v-coordinate of velocity corresponding to the most recent
!          estimate of the trajectory mid-point (in geodesic system).
!  phipr   Phi value at trajectory mid-point (geodesic coordinates).
!          On entry this is the most recent estimate.
!  lampr   Lambda value at trajectory mid-point (geodesic coordinates).
!-----------------------------------------------------------------------

!---------------------------Local variables-----------------------------
  integer i,k                 ! index
!-----------------------------------------------------------------------
!
!$OMP PARALLEL DO PRIVATE (K, I)
  do k=1,plev
     do i = 1,nlon
        lampr(i,k) = -.5_r8*dt* upr(i,k) / cos( phipr(i,k) )
        phipr(i,k) = -.5_r8*dt* vpr(i,k)
     end do
  end do

  return
end subroutine trajmp

!============================================================================================

subroutine trjgl(finc    ,phicen  ,lam     ,lampr   ,phipr   , &
                 lamp    ,phip    ,nlon    )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Map relative trajectory mid/departure point coordinates to global
! latitude/longitude coordinates and test limits
! 
! Method: 
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
  use scanslt,      only: plond1
  implicit none

!------------------------------Arguments--------------------------------
  integer , intent(in) :: nlon                 ! longitude dimension
  real(r8), intent(in) :: finc                 ! number of time increments
  real(r8), intent(in) :: phicen               ! current latitude value in extnded grid
  real(r8), intent(in) :: lam(plond1)          ! longitude values for the extended grid
  real(r8), intent(in) :: lampr(plon,plev)     ! relative x-coordinate of departure pt.
  real(r8), intent(in) :: phipr(plon,plev)     ! relative y-coordinate of departure pt.
  real(r8), intent(out) :: lamp (plon,plev)    ! long coords of traj midpoints
  real(r8), intent(out) :: phip (plon,plev)    ! lat  coords of traj midpoints
!
! finc    Time step factor (1. for midpoint, 2. for dep. point)
! phicen  Latitude value for current latitude being forecast.
! lam     Longitude values for the extended grid.
! lampr   Longitude coordinates (relative to the arrival point) of the
!         trajectory mid-points of the parcels that correspond to the
!         global grid points contained in the latitude slice being forecast.
! phipr   Latitude coordinates (relative to the arrival point) of the
!         trajectory mid-points of the parcels that correspond to the
!         global grid points contained in the latitude slice being forecast.
! lamp    Longitude coordinates of the trajectory mid-points of the
!         parcels that correspond to the global grid points contained
!         in the latitude slice being forecast.
! phip    Latitude  coordinates of the trajectory mid-points of the
!         parcels that correspond to the global grid points contained
!         in the latitude slice being forecast.
!-----------------------------------------------------------------------

!--------------------------Local variables------------------------------
  integer i                 ! longitude index
  integer k                 ! level index
  real(r8) pi               ! 3.14.......
  real(r8) twopi            ! 2*pi
!-----------------------------------------------------------------------
!
  pi = 4._r8*atan(1._r8)
  twopi = pi*2._r8
!$OMP PARALLEL DO PRIVATE (K, I)
  do k = 1,plev
     do i = 1,nlon
        lamp(i,k) = lam(i) + finc*lampr(i,k)
        phip(i,k) = phicen + finc*phipr(i,k)
        if(lamp(i,k) >= twopi) lamp(i,k) = lamp(i,k) - twopi
        if(lamp(i,k) <    0.0_r8) lamp(i,k) = lamp(i,k) + twopi
     end do
  end do

  return
end subroutine trjgl

