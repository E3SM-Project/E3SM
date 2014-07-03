subroutine sltint(jcen    ,fb      ,lam     ,rdphi   , &
                  rdz     ,lbasdy  ,wdz     ,xl      ,xr      , &
                  wgt1x   ,wgt2x   ,wgt3x   ,wgt4x   ,hl      , &
                  hr      ,dhl     ,dhr     ,ys      ,yn      , &
                  wgt1y   ,wgt2y   ,wgt3y   ,wgt4y   ,hs      , &
                  hn      ,dhs     ,dhn     ,wgt1z   ,wgt2z   , &
                  wgt3z   ,wgt4z   ,hb      ,ht      ,dhb     , &
                  dht     ,idp     ,jdp     ,kdp     ,kkdp    , &
                  lhrzint ,lvrtint ,limdrh  ,limdrv  ,fdp     , &
                  nlon    )
!
!-----------------------------------------------------------------------
!
! Purpose:
! Interpolate field to departure points using Hermite or Lagrange
! Cubic interpolation
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
  use pmgrid,       only: plon, plev
  use abortutils, only: endrun
  use scanslt,      only: beglatex, endlatex, i1, plond, platd, nxpt
#if (!defined UNICOSMP)
  use srchutil, only: whenieq
#endif

  implicit none

!----------------------------- Parameters ------------------------------
!
  integer, parameter :: ppdy = 4 ! length of interpolation grid stencil
  integer, parameter :: ppdz = 4 ! length of interp. grid stencil in z
!
!------------------------------Arguments--------------------------------
!
  integer , intent(in)   :: jcen                ! index of latitude in extended grid

  real(r8), intent(in)   :: fb(plond,plev,beglatex:endlatex) ! input field
  real(r8), intent(in)   :: lam   (plond,platd) ! longitude coordinates of model grid
  real(r8), intent(in)   :: rdphi (plon,plev)   ! reciprocal of y-interval
  real(r8), intent(in)   :: rdz   (plon,plev)   ! reciprocal of z-interval
  real(r8), intent(in)   :: lbasdy(4,2,platd)   ! basis functions for lat deriv est.
  real(r8), intent(in)   :: wdz   (4,2,plev)    ! basis functions for vert deriv est.
  real(r8), intent(in)   :: xl    (plon,plev,4) ! weight for x-interpolants (left)
  real(r8), intent(in)   :: xr    (plon,plev,4) ! weight for x-interpolants (right)
  real(r8), intent(in)   :: wgt1x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4x (plon,plev,4) ! weight for x-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hl    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: hr    (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: dhl   (plon,plev,4) ! weight for x-interpolants (Hermite)
  real(r8), intent(in)   :: dhr   (plon,plev,4) ! weight for x-interpolants (Hermite)

  real(r8), intent(in)   :: ys    (plon,plev)   ! weight for y-interpolants (south)
  real(r8), intent(in)   :: yn    (plon,plev)   ! weight for y-interpolants (north)
  real(r8), intent(in)   :: wgt1y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4y (plon,plev)   ! weight for y-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hs    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: hn    (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: dhs   (plon,plev)   ! weight for y-interpolants (Hermite)
  real(r8), intent(in)   :: dhn   (plon,plev)   ! weight for y-interpolants (Hermite)

  real(r8), intent(in)   :: wgt1z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt2z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt3z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: wgt4z (plon,plev)   ! weight for z-interpolants (Lag Cubic)
  real(r8), intent(in)   :: hb    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: ht    (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: dhb   (plon,plev)   ! weight for z-interpolants (Hermite)
  real(r8), intent(in)   :: dht   (plon,plev)   ! weight for z-interpolants (Hermite)

  integer , intent(in)   :: idp   (plon,plev,4) ! index of x-coordinate of dep pt
  integer , intent(in)   :: jdp   (plon,plev)   ! index of y-coordinate of dep pt
  integer , intent(in)   :: kdp   (plon,plev)   ! index of z-coordinate of dep pt
  integer , intent(in)   :: kkdp  (plon,plev)   ! index of z-coordinate of dep pt (alt)

  logical , intent(in)   :: lhrzint             ! flag to do horizontal interpolation
  logical , intent(in)   :: lvrtint             ! flag to do vertical   interpolation
  logical , intent(in)   :: limdrh              ! horizontal derivative limiter flag
  logical , intent(in)   :: limdrv              ! vertical   derivative limiter flag
  real(r8), intent(out)  :: fdp(plon,plev)      ! interpolant
  integer , intent(in)   :: nlon                ! number of longitudes for this latitude
!
!---------------------------Local workspace-----------------------------
!
  integer i                           ! |
  integer ii1,ii2,ii3,ii4             ! |
  integer ii,jj,j                     ! |
  integer k                           ! |
  integer kk                          ! |
  integer jmin                        ! |
  integer jmax                        ! | -- indices
  integer jdpval                      ! |
  integer kdpval                      ! |
  integer icount                      ! |
  integer indx(plon)                  ! |
  integer nval                        ! |
  integer kdimm1                      ! |
  integer kdimm2                      ! |
  integer kdimm3                      ! |
!
  real(r8) fac                        ! factor applied in limiter
  real(r8) tmp1                       ! derivative factor
  real(r8) tmp2                       ! abs(tmp1)
  real(r8) deli                       ! linear derivative
  real(r8) dx                         ! delta x
  real(r8) rdx (platd)                ! 1./dx
  real(r8) rdx6(platd)                ! 1./(6*dx)
  real(r8) fxl                        ! left  derivative estimate
  real(r8) fxr                        ! right derivative estimate
!
  real(r8) f1                         ! |
  real(r8) f2                         ! |
  real(r8) f3                         ! |
  real(r8) f4                         ! |
  real(r8) tmptop                     ! | -- work arrays
  real(r8) tmpbot                     ! |
  real(r8) fintx(plon,plev,ppdy,ppdz) ! |
  real(r8) finty(plon,plev     ,ppdz) ! |
  real(r8) fbot (plon,plev,ppdz)      ! |
  real(r8) ftop (plon,plev,ppdz)      ! |
  integer :: kdim                     ! Vertical levels operating over
!
!-----------------------------------------------------------------------
!
  kdim = plev
  fac  = 3._r8*(1._r8 - 10._r8*epsilon(fac))
  kdimm1 = kdim - 1
  kdimm2 = kdim - 2
  kdimm3 = kdim - 3
!
  do j = 1,platd
     dx      = lam(nxpt+2,j) - lam(nxpt+1,j)
     rdx (j) = 1._r8/dx
     rdx6(j) = 1._r8/(6._r8*dx)
  end do
!
!-----------------------------------------------------------------------
!------------------------- Code Description ----------------------------
!
! Each block of code performs a specific interpolation as follows:
!
!  For 3-D (horizontal AND vertical) interpolation:
!
!    10XX loops: Hermite  cubic/linear interpolation in the horizontal
!    20XX loops: Lagrange cubic/linear interpolation in the horizontal
!    30XX loops: Hermite  cubic/linear interpolation in the vertical
!    40XX loops: Lagrange cubic/linear interpolation in the vertical
!
!  For horizontal interpolation only:
!
!    50XX loops: an optimized Lagrange cubic/linear algorithm
!    60XX loops: Hermite cubic/linear interpolation in the horizontal
!
!  For vertical interpolation only:
!
!    70XX loops: an optimized Lagrange cubic/linear algorithm (no
!                Hermite interpolator available)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  if( lhrzint .and. lvrtint ) then
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!    10XX loops: Hermite  cubic/linear interpolation in the horizontal
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
     if( limdrh ) then
!
! PART 1:  x-interpolation
!
! Loop over fields.
! ..x interpolation at each height needed for z interpolation.
! ...x interpolation at each latitude needed for y interpolation.
!
        do k=1,plev
           do i=1,nlon
              ii1 = idp(i,k,1)
              ii2 = idp(i,k,2)
              ii3 = idp(i,k,3)
              ii4 = idp(i,k,4)
              jj = jdp(i,k)
              kk = kkdp(i,k)
!
! Height level 1:  Linear interpolation on inner two latitudes only
!
!!!           fintx(i,k,1,1) = not used
              fintx(i,k,2,1) = fb (ii2  ,kk-1,jj  )*xl (i,k,2) &
                             + fb (ii2+1,kk-1,jj  )*xr (i,k,2)
              fintx(i,k,3,1) = fb (ii3  ,kk-1,jj+1)*xl (i,k,3) &
                             + fb (ii3+1,kk-1,jj+1)*xr (i,k,3)
!!!           fintx(i,k,4,1) = not used
!
! Height level 2
!
!   Latitude 1:  Linear interpolation
!
              fintx(i,k,1,2) = fb (ii1  ,kk,jj-1)*xl (i,k,1) &
                             + fb (ii1+1,kk,jj-1)*xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
              fxl = (   - 2._r8*fb (ii2-1,kk,jj) &
                        - 3._r8*fb (ii2  ,kk,jj) &
                        + 6._r8*fb (ii2+1,kk,jj) &
                        -    fb (ii2+2,kk,jj) )*rdx6(jj)
              fxr = (        fb (ii2-1,kk,jj) &
                        - 6._r8*fb (ii2  ,kk,jj) &
                        + 3._r8*fb (ii2+1,kk,jj) &
                        + 2._r8*fb (ii2+2,kk,jj) )*rdx6(jj)
!
              deli = (       fb (ii2+1,kk,jj) - &
                             fb (ii2  ,kk,jj) )*rdx(jj)
              tmp1 = fac*deli
              tmp2 = abs( tmp1 )
              if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
              if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
              if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
              if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
              fintx(i,k,2,2) = fb (ii2  ,kk,jj)*hl (i,k,2) &
                             + fb (ii2+1,kk,jj)*hr (i,k,2) &
                             + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
              fxl = (   - 2._r8*fb (ii3-1,kk  ,jj+1) &
                        - 3._r8*fb (ii3  ,kk  ,jj+1) &
                        + 6._r8*fb (ii3+1,kk  ,jj+1) &
                        -    fb (ii3+2,kk  ,jj+1) )*rdx6(jj+1)
              fxr = (        fb (ii3-1,kk  ,jj+1) &
                        - 6._r8*fb (ii3  ,kk  ,jj+1) &
                        + 3._r8*fb (ii3+1,kk  ,jj+1) &
                        + 2._r8*fb (ii3+2,kk  ,jj+1) )*rdx6(jj+1)
!
              deli = (       fb (ii3+1,kk  ,jj+1) - &
                             fb (ii3  ,kk  ,jj+1) )*rdx(jj+1)
              tmp1 = fac*deli
              tmp2 = abs( tmp1 )
              if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
              if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
              if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
              if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
              fintx(i,k,3,2) = fb (ii3  ,kk  ,jj+1)*hl (i,k,3) &
                             + fb (ii3+1,kk  ,jj+1)*hr (i,k,3) &
                             + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
              fintx(i,k,4,2) = fb (ii4  ,kk,jj+2)*xl (i,k,4) &
                             + fb (ii4+1,kk,jj+2)*xr (i,k,4)
!
! Height level 3
!
!   Latitude 1:  Linear interpolation
!
              fintx(i,k,1,3) = fb (ii1  ,kk+1,jj-1)*xl (i,k,1) &
                             + fb (ii1+1,kk+1,jj-1)*xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
              fxl = (   - 2._r8*fb (ii2-1,kk+1,jj  ) &
                        - 3._r8*fb (ii2  ,kk+1,jj  ) &
                        + 6._r8*fb (ii2+1,kk+1,jj  ) &
                        -    fb (ii2+2,kk+1,jj  ) )*rdx6(jj)
              fxr = (        fb (ii2-1,kk+1,jj  ) &
                        - 6._r8*fb (ii2  ,kk+1,jj  ) &
                        + 3._r8*fb (ii2+1,kk+1,jj  ) &
                        + 2._r8*fb (ii2+2,kk+1,jj  ) )*rdx6(jj)
!
              deli = (       fb (ii2+1,kk+1,jj  ) - &
                             fb (ii2  ,kk+1,jj  ) )*rdx(jj)
              tmp1 = fac*deli
              tmp2 = abs( tmp1 )
              if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
              if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
              if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
              if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
              fintx(i,k,2,3) = fb (ii2  ,kk+1,jj  )*hl (i,k,2) &
                             + fb (ii2+1,kk+1,jj  )*hr (i,k,2) &
                             + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
              fxl = (   - 2._r8*fb (ii3-1,kk+1,jj+1) &
                        - 3._r8*fb (ii3  ,kk+1,jj+1) &
                        + 6._r8*fb (ii3+1,kk+1,jj+1) &
                        -    fb (ii3+2,kk+1,jj+1) )*rdx6(jj+1)
              fxr = (        fb (ii3-1,kk+1,jj+1) &
                        - 6._r8*fb (ii3  ,kk+1,jj+1) &
                        + 3._r8*fb (ii3+1,kk+1,jj+1) &
                        + 2._r8*fb (ii3+2,kk+1,jj+1) )*rdx6(jj+1)
!
              deli = (       fb (ii3+1,kk+1,jj+1) - &
                             fb (ii3  ,kk+1,jj+1) )*rdx(jj+1)
              tmp1 = fac*deli
              tmp2 = abs( tmp1 )
              if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
              if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
              if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
              if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
              fintx(i,k,3,3) = fb (ii3  ,kk+1,jj+1)*hl (i,k,3) &
                             + fb (ii3+1,kk+1,jj+1)*hr (i,k,3) &
                             + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
              fintx(i,k,4,3) = fb (ii4  ,kk+1,jj+2)*xl (i,k,4) &
                             + fb (ii4+1,kk+1,jj+2)*xr (i,k,4)
!
! Height level 4:  Linear interpolation on inner two latitudes only
!
!!!           fintx(i,k,1,4) = not used
              fintx(i,k,2,4) = fb (ii2  ,kk+2,jj  )*xl (i,k,2) &
                             + fb (ii2+1,kk+2,jj  )*xr (i,k,2)
              fintx(i,k,3,4) = fb (ii3  ,kk+2,jj+1)*xl (i,k,3) &
                             + fb (ii3+1,kk+2,jj+1)*xr (i,k,3)
!!!           fintx(i,k,4,4) = not used
           end do
        end do
!
! The following loop computes x-derivatives for those cases when the
! departure point lies in either the top or bottom interval of the 
! model grid.  In this special case, data are shifted up or down to
! keep the departure point in the middle interval of the 4-point
! stencil.  Therefore, some derivatives that were computed above will 
! be over-written.
!
        do k=1,plev
           do i=1,nlon
              ii1 = idp(i,k,1)
              ii2 = idp(i,k,2)
              ii3 = idp(i,k,3)
              ii4 = idp(i,k,4)
              jj = jdp(i,k)
              kk = kkdp(i,k)
!
! TOP interval
!
              if(kdp (i,k) .eq. 1) then
!
! shift levels 4 and 2 data to levels 1 and 3, respectively
!
                 fintx(i,k,2,1) = fintx(i,k,2,4)
                 fintx(i,k,3,1) = fintx(i,k,3,4)
!
                 fintx(i,k,1,3) = fintx(i,k,1,2)
                 fintx(i,k,2,3) = fintx(i,k,2,2)
                 fintx(i,k,3,3) = fintx(i,k,3,2)
                 fintx(i,k,4,3) = fintx(i,k,4,2)
!
! Height level 1 (placed in level 2 of stencil):
!
!   Latitude 1:  Linear interpolation
!
                 fintx(i,k,1,2) = fb (ii1  ,1,jj-1)*xl (i,k,1) &
                                + fb (ii1+1,1,jj-1)*xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
                 fxl = (   - 2._r8*fb (ii2-1,1,jj  ) &
                           - 3._r8*fb (ii2  ,1,jj  ) &
                           + 6._r8*fb (ii2+1,1,jj  ) &
                           -    fb (ii2+2,1,jj  ) )*rdx6(jj)
                 fxr = (        fb (ii2-1,1,jj  ) &
                           - 6._r8*fb (ii2  ,1,jj  ) &
                           + 3._r8*fb (ii2+1,1,jj  ) &
                           + 2._r8*fb (ii2+2,1,jj  ) )*rdx6(jj)
!                                      
                 deli = (       fb (ii2+1,1,jj  ) - &
                                fb (ii2  ,1,jj  ) )*rdx(jj)
                 tmp1 = fac*deli
                 tmp2 = abs( tmp1 )
                 if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
                 if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
                 if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
                 if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
                 fintx(i,k,2,2) = fb (ii2  ,1,jj  )*hl (i,k,2) &
                                + fb (ii2+1,1,jj  )*hr (i,k,2) &
                                + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
                 fxl = (   - 2._r8*fb (ii3-1,1,jj+1) &
                           - 3._r8*fb (ii3  ,1,jj+1) &
                           + 6._r8*fb (ii3+1,1,jj+1) &
                           -    fb (ii3+2,1,jj+1) )*rdx6(jj+1)
                 fxr = (        fb (ii3-1,1,jj+1) &
                           - 6._r8*fb (ii3  ,1,jj+1) &
                           + 3._r8*fb (ii3+1,1,jj+1) &
                           + 2._r8*fb (ii3+2,1,jj+1) )*rdx6(jj+1)
!                                      
                 deli = (       fb (ii3+1,1,jj+1) - &
                                fb (ii3  ,1,jj+1) )*rdx(jj+1)
                 tmp1 = fac*deli
                 tmp2 = abs( tmp1 )
                 if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
                 if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
                 if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
                 if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
                 fintx(i,k,3,2) = fb (ii3  ,1,jj+1)*hl (i,k,3) &
                                + fb (ii3+1,1,jj+1)*hr (i,k,3) &
                                + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
                 fintx(i,k,4,2) = fb (ii4  ,1,jj+2)*xl (i,k,4) &
                                + fb (ii4+1,1,jj+2)*xr (i,k,4)
!
! Height level 3 (placed in level 4 of stencil):
!  Linear interpolation on inner two latitudes only
!
!!!              fintx(i,k,1,4) = not used
                 fintx(i,k,2,4) = fb (ii2  ,3,jj  )*xl (i,k,2) &
                                + fb (ii2+1,3,jj  )*xr (i,k,2)
                 fintx(i,k,3,4) = fb (ii3  ,3,jj+1)*xl (i,k,3) &
                                + fb (ii3+1,3,jj+1)*xr (i,k,3)
!!!              fintx(i,k,4,4) = not used
!
! BOT interval
!
              else if(kdp (i,k) .eq. kdimm1) then
!
! shift levels 1 and 3 data to levels 4 and 2, respectively
!
                 fintx(i,k,2,4) = fintx(i,k,2,1)
                 fintx(i,k,3,4) = fintx(i,k,3,1)
!
                 fintx(i,k,1,2) = fintx(i,k,1,3)
                 fintx(i,k,2,2) = fintx(i,k,2,3)
                 fintx(i,k,3,2) = fintx(i,k,3,3)
                 fintx(i,k,4,2) = fintx(i,k,4,3)
!
! Height level 2 (placed in level 1 of stencil):
!  Linear interpolation on inner two latitudes only
!
!!!              fintx(i,k,1,1) =  not used
                 fintx(i,k,2,1) = fb (ii2  ,kdimm2,jj  )*xl (i,k,2) &
                                + fb (ii2+1,kdimm2,jj  )*xr (i,k,2)
                 fintx(i,k,3,1) = fb (ii3  ,kdimm2,jj+1)*xl (i,k,3) &
                                + fb (ii3+1,kdimm2,jj+1)*xr (i,k,3)
!!!              fintx(i,k,4,1) =  not used
!
! Height level 4 (placed in level 3 of stencil):
!
!   Latitude 1:  Linear interpolation
!
                 fintx(i,k,1,3) = fb (ii1  ,kdim,jj-1)*xl (i,k,1) &
                                + fb (ii1+1,kdim,jj-1)*xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
                 fxl = (   - 2._r8*fb (ii2-1,kdim,jj  ) &
                           - 3._r8*fb (ii2  ,kdim,jj  ) &
                           + 6._r8*fb (ii2+1,kdim,jj  ) &
                           -    fb (ii2+2,kdim,jj  ) )*rdx6(jj)
                 fxr = (        fb (ii2-1,kdim,jj  ) &
                           - 6._r8*fb (ii2  ,kdim,jj  ) &
                           + 3._r8*fb (ii2+1,kdim,jj  ) &
                           + 2._r8*fb (ii2+2,kdim,jj  ) )*rdx6(jj)
!                                          
                 deli = (       fb (ii2+1,kdim,jj  ) - &
                                fb (ii2  ,kdim,jj  ) )*rdx(jj)
                 tmp1 = fac*deli
                 tmp2 = abs( tmp1 )
                 if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
                 if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
                 if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
                 if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
                 fintx(i,k,2,3) = fb (ii2  ,kdim,jj  )*hl (i,k,2) &
                                + fb (ii2+1,kdim,jj  )*hr (i,k,2) &
                                + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
                 fxl = (   - 2._r8*fb (ii3-1,kdim,jj+1) &
                           - 3._r8*fb (ii3  ,kdim,jj+1) &
                           + 6._r8*fb (ii3+1,kdim,jj+1) &
                           -    fb (ii3+2,kdim,jj+1) )*rdx6(jj+1)
                 fxr = (        fb (ii3-1,kdim,jj+1) &
                           - 6._r8*fb (ii3  ,kdim,jj+1) &
                           + 3._r8*fb (ii3+1,kdim,jj+1) &
                           + 2._r8*fb (ii3+2,kdim,jj+1) )*rdx6(jj+1)
!                                          
                 deli = (       fb (ii3+1,kdim,jj+1) - &
                                fb (ii3  ,kdim,jj+1) )*rdx(jj+1)
                 tmp1 = fac*deli
                 tmp2 = abs( tmp1 )
                 if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
                 if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
                 if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
                 if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
                 fintx(i,k,3,3) = fb (ii3  ,kdim,jj+1)*hl (i,k,3) &
                                + fb (ii3+1,kdim,jj+1)*hr (i,k,3) &
                                + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
                 fintx(i,k,4,3) = fb (ii4  ,kdim,jj+2)*xl (i,k,4) &
                                + fb (ii4+1,kdim,jj+2)*xr (i,k,4)
              end if
           end do
        end do
!
! PART 2:  y-derivatives
!
        jmin =  1000000
        jmax = -1000000
        do k=1,plev
           do i=1,nlon
              if(jdp(i,k) .lt. jmin) jmin = jdp(i,k)
              if(jdp(i,k) .gt. jmax) jmax = jdp(i,k)
           end do
        end do
!
! Loop over departure latitudes
!
        icount = 0
        do jdpval = jmin,jmax
           do k=1,plev
              call whenieq(nlon    ,jdp(1,k),1       ,jdpval  , &
                           indx    ,nval    )
              icount = icount + nval
!
! y derivatives at the inner height levels (kk = 2,3) needed for
! z-interpolation
!
              do kk  = 2,3
                 do ii = 1,nval
                    i = indx(ii)
                    fbot(i,k,kk) = lbasdy(1,1,jdpval)*fintx(i,k,1,kk) &
                                 + lbasdy(2,1,jdpval)*fintx(i,k,2,kk) &
                                 + lbasdy(3,1,jdpval)*fintx(i,k,3,kk) &
                                 + lbasdy(4,1,jdpval)*fintx(i,k,4,kk)
                    ftop(i,k,kk) = lbasdy(1,2,jdpval)*fintx(i,k,1,kk) &
                                 + lbasdy(2,2,jdpval)*fintx(i,k,2,kk) &
                                 + lbasdy(3,2,jdpval)*fintx(i,k,3,kk) &
                                 + lbasdy(4,2,jdpval)*fintx(i,k,4,kk)
                 end do
              end do
           end do
        end do
        if (icount.ne.nlon*plev) then
           call endrun ('SLTINT:  Did not complete computations for all departure points')
        end if
!
! Apply SCM0 limiter to derivative estimates.
!
        do kk  = 2,3
           do k=1,plev
              do i=1,nlon
                 deli = ( fintx(i,k,3,kk) - fintx(i,k,2,kk) )*rdphi(i,k)
                 tmp1 = fac*deli
                 tmp2 = abs( tmp1 )
                 if( deli*fbot(i,k,kk)   .le. 0.0_r8  ) fbot(i,k,kk) = 0._r8
                 if( deli*ftop(i,k,kk)   .le. 0.0_r8  ) ftop(i,k,kk) = 0._r8
                 if( abs( fbot(i,k,kk) ) .gt. tmp2 ) fbot(i,k,kk) = tmp1
                 if( abs( ftop(i,k,kk) ) .gt. tmp2 ) ftop(i,k,kk) = tmp1
              end do
           end do
        end do
!
! PART 3:  y-interpolants
!
        do k=1,plev
           do i=1,nlon
              finty(i,k,1) = fintx(i,k,2,1)*ys (i,k) &
                           + fintx(i,k,3,1)*yn (i,k)
              finty(i,k,2) = fintx(i,k,2,2)*hs (i,k) + fbot (i,k  ,2)*dhs(i,k) &
                           + fintx(i,k,3,2)*hn (i,k) + ftop (i,k  ,2)*dhn(i,k)
              finty(i,k,3) = fintx(i,k,2,3)*hs (i,k) + fbot (i,k  ,3)*dhs(i,k) &
                           + fintx(i,k,3,3)*hn (i,k) + ftop (i,k  ,3)*dhn(i,k)
              finty(i,k,4) = fintx(i,k,2,4)*ys (i,k) &
                           + fintx(i,k,3,4)*yn (i,k)
           end do
        end do
     endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!    20XX loops: Lagrange cubic/linear interpolation in the horizontal
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
     if( .not. limdrh ) then
!
! PART 1:  X-INTERPOLATION
!
! Loop over fields.
! ..x interpolation at each height needed for z interpolation.
! ...x interpolation at each latitude needed for y interpolation.
!
        do k=1,plev
           do i=1,nlon
              ii1 = idp(i,k,1)
              ii2 = idp(i,k,2)
              ii3 = idp(i,k,3)
              ii4 = idp(i,k,4)
              jj = jdp(i,k)
              kk = kkdp(i,k)
!
! Height level 1:  Linear interpolation on inner two latitudes only
!
!!!           fintx(i,k,1,1) = not used
              fintx(i,k,2,1) = fb (ii2+1,kk-1,jj  )*xr   (i,k,2) &
                             + fb (ii2  ,kk-1,jj  )*xl   (i,k,2)
              fintx(i,k,3,1) = fb (ii3+1,kk-1,jj+1)*xr   (i,k,3) &
                             + fb (ii3  ,kk-1,jj+1)*xl   (i,k,3)
!!!           fintx(i,k,4,1) = not used
!
! Height level 2:  Linear interpolation on outer two latitudes;
!                  Cubic  interpolation on inner two latitudes.
!
              fintx(i,k,1,2) = fb (ii1+1,kk  ,jj-1)*xr   (i,k,1) &
                             + fb (ii1  ,kk  ,jj-1)*xl   (i,k,1)
              fintx(i,k,2,2) = fb (ii2-1,kk  ,jj  )*wgt1x(i,k,2) &
                             + fb (ii2  ,kk  ,jj  )*wgt2x(i,k,2) &
                             + fb (ii2+1,kk  ,jj  )*wgt3x(i,k,2) &
                             + fb (ii2+2,kk  ,jj  )*wgt4x(i,k,2)
              fintx(i,k,3,2) = fb (ii3-1,kk  ,jj+1)*wgt1x(i,k,3) &
                             + fb (ii3  ,kk  ,jj+1)*wgt2x(i,k,3) &
                             + fb (ii3+1,kk  ,jj+1)*wgt3x(i,k,3) &
                             + fb (ii3+2,kk  ,jj+1)*wgt4x(i,k,3)
              fintx(i,k,4,2) = fb (ii4+1,kk  ,jj+2)*xr   (i,k,4) &
                             + fb (ii4  ,kk  ,jj+2)*xl   (i,k,4)
!
! Height level 3:  Linear interpolation on outer two latitudes;
!                  Cubic  interpolation on inner two latitudes.
!
              fintx(i,k,1,3) = fb (ii1+1,kk+1,jj-1)*xr   (i,k,1) &
                             + fb (ii1  ,kk+1,jj-1)*xl   (i,k,1)
              fintx(i,k,2,3) = fb (ii2-1,kk+1,jj  )*wgt1x(i,k,2) &
                             + fb (ii2  ,kk+1,jj  )*wgt2x(i,k,2) &
                             + fb (ii2+1,kk+1,jj  )*wgt3x(i,k,2) &
                             + fb (ii2+2,kk+1,jj  )*wgt4x(i,k,2)
              fintx(i,k,3,3) = fb (ii3-1,kk+1,jj+1)*wgt1x(i,k,3) &
                             + fb (ii3  ,kk+1,jj+1)*wgt2x(i,k,3) &
                             + fb (ii3+1,kk+1,jj+1)*wgt3x(i,k,3) &
                             + fb (ii3+2,kk+1,jj+1)*wgt4x(i,k,3)
              fintx(i,k,4,3) = fb (ii4+1,kk+1,jj+2)*xr   (i,k,4) &
                             + fb (ii4  ,kk+1,jj+2)*xl   (i,k,4)
!
! Height level 4:  Linear interpolation on inner two latitudes only
!
!!!           fintx(i,k,1,4) = not used
              fintx(i,k,2,4) = fb (ii2+1,kk+2,jj  )*xr   (i,k,2) &
                             + fb (ii2  ,kk+2,jj  )*xl   (i,k,2)
              fintx(i,k,3,4) = fb (ii3+1,kk+2,jj+1)*xr   (i,k,3) &
                             + fb (ii3  ,kk+2,jj+1)*xl   (i,k,3)
!!!           fintx(i,k,4,4) = not used
           end do
        end do
!
! The following loop computes x-derivatives for those cases when the
! departure point lies in either the top or bottom interval of the 
! model grid.  In this special case, data are shifted up or down to
! keep the departure point in the middle interval of the 4-point
! stencil.  Therefore, some derivatives that were computed above will 
! be over-written.
!
        do k=1,plev
           do i=1,nlon
              ii1 = idp(i,k,1)
              ii2 = idp(i,k,2)
              ii3 = idp(i,k,3)
              ii4 = idp(i,k,4)
              jj = jdp(i,k)
              kk = kkdp(i,k)
!
! TOP interval
!
              if(kdp (i,k) .eq. 1) then
!
! shift levels 4 and 2 data to levels 1 and 3, respectively
!
                 fintx(i,k,2,1) = fintx(i,k,2,4)
                 fintx(i,k,3,1) = fintx(i,k,3,4)
!
                 fintx(i,k,1,3) = fintx(i,k,1,2)
                 fintx(i,k,2,3) = fintx(i,k,2,2)
                 fintx(i,k,3,3) = fintx(i,k,3,2)
                 fintx(i,k,4,3) = fintx(i,k,4,2)
!
! Height level 1 (placed in level 2 of stencil):
!  Linear interpolation on outer two latitudes;
!  Cubic  interpolation on inner two latitudes.
!
                 fintx(i,k,1,2) = fb (ii1+1,1,jj-1)*xr   (i,k,1) &
                                + fb (ii1  ,1,jj-1)*xl   (i,k,1)
                 fintx(i,k,2,2) = fb (ii2-1,1,jj  )*wgt1x(i,k,2) &
                                + fb (ii2  ,1,jj  )*wgt2x(i,k,2) &
                                + fb (ii2+1,1,jj  )*wgt3x(i,k,2) &
                                + fb (ii2+2,1,jj  )*wgt4x(i,k,2)
                 fintx(i,k,3,2) = fb (ii3-1,1,jj+1)*wgt1x(i,k,3) &
                                + fb (ii3  ,1,jj+1)*wgt2x(i,k,3) &
                                + fb (ii3+1,1,jj+1)*wgt3x(i,k,3) &
                                + fb (ii3+2,1,jj+1)*wgt4x(i,k,3)
                 fintx(i,k,4,2) = fb (ii4+1,1,jj+2)*xr   (i,k,4) &
                                + fb (ii4  ,1,jj+2)*xl   (i,k,4)
!
! Height level 3 (placed in level 4 of stencil):
!  Linear interpolation on inner two latitudes only
!
!!!              fintx(i,k,1,4) = not used
                 fintx(i,k,2,4) = fb (ii2+1,3,jj  )*xr   (i,k,2) &
                                + fb (ii2  ,3,jj  )*xl   (i,k,2)
                 fintx(i,k,3,4) = fb (ii3+1,3,jj+1)*xr   (i,k,3) &
                                + fb (ii3  ,3,jj+1)*xl   (i,k,3)
!!!              fintx(i,k,4,4) = not used
!
! BOT interval
!
              else if(kdp (i,k) .eq. kdimm1) then
!
! shift levels 1 and 3 data to levels 4 and 2, respectively
!
                 fintx(i,k,2,4) = fintx(i,k,2,1)
                 fintx(i,k,3,4) = fintx(i,k,3,1)
!
                 fintx(i,k,1,2) = fintx(i,k,1,3)
                 fintx(i,k,2,2) = fintx(i,k,2,3)
                 fintx(i,k,3,2) = fintx(i,k,3,3)
                 fintx(i,k,4,2) = fintx(i,k,4,3)
!
! Height level 2 (placed in level 1 of stencil):
!  Linear interpolation on inner two latitudes only
!
!!!              fintx(i,k,1,1) = not used
                 fintx(i,k,2,1) = fb (ii2+1,kdimm2,jj  )*xr   (i,k,2) &
                                + fb (ii2  ,kdimm2,jj  )*xl   (i,k,2)
                 fintx(i,k,3,1) = fb (ii3+1,kdimm2,jj+1)*xr   (i,k,3) &
                                + fb (ii3  ,kdimm2,jj+1)*xl   (i,k,3)
!!!              fintx(i,k,4,1) = not used
!
! Height level 4 (placed in level 3 of stencil):
!  Linear interpolation on outer two latitudes;
!  Cubic  interpolation on inner two latitudes.
!
                 fintx(i,k,1,3) = fb (ii1+1,kdim,jj-1)*xr   (i,k,1) &
                                + fb (ii1  ,kdim,jj-1)*xl   (i,k,1)
                 fintx(i,k,2,3) = fb (ii2-1,kdim,jj  )*wgt1x(i,k,2) &
                                + fb (ii2  ,kdim,jj  )*wgt2x(i,k,2) &
                                + fb (ii2+1,kdim,jj  )*wgt3x(i,k,2) &
                                + fb (ii2+2,kdim,jj  )*wgt4x(i,k,2)
                 fintx(i,k,3,3) = fb (ii3-1,kdim,jj+1)*wgt1x(i,k,3) &
                                + fb (ii3  ,kdim,jj+1)*wgt2x(i,k,3) &
                                + fb (ii3+1,kdim,jj+1)*wgt3x(i,k,3) &
                                + fb (ii3+2,kdim,jj+1)*wgt4x(i,k,3)
                 fintx(i,k,4,3) = fb (ii4+1,kdim,jj+2)*xr   (i,k,4) &
                                + fb (ii4  ,kdim,jj+2)*xl   (i,k,4)
              end if
           end do
        end do
!
! PART 2:  Y-INTERPOLATION
!
! Linear on outside of stencil; Lagrange cubic on inside.
!
        do k=1,plev
           do i=1,nlon
              finty(i,k,1) = fintx(i,k,2,1)*ys   (i,k) &
                           + fintx(i,k,3,1)*yn   (i,k)
              finty(i,k,2) = fintx(i,k,1,2)*wgt1y(i,k) &
                           + fintx(i,k,2,2)*wgt2y(i,k) &
                           + fintx(i,k,3,2)*wgt3y(i,k) &
                           + fintx(i,k,4,2)*wgt4y(i,k)
              finty(i,k,3) = fintx(i,k,1,3)*wgt1y(i,k) &
                           + fintx(i,k,2,3)*wgt2y(i,k) &
                           + fintx(i,k,3,3)*wgt3y(i,k) &
                           + fintx(i,k,4,3)*wgt4y(i,k)
              finty(i,k,4) = fintx(i,k,2,4)*ys   (i,k) &
                           + fintx(i,k,3,4)*yn   (i,k)
           end do
        end do
     endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!    30XX loops: Hermite  cubic/linear interpolation in the vertical
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
     if( limdrv ) then
        icount = 0
        do kdpval = 1,kdimm1
           do k=1,plev
              call whenieq(nlon    ,kdp(1,k),1       ,kdpval  , &
                           indx    ,nval    )
              icount = icount + nval
              do ii = 1,nval
                 i = indx(ii)
                 if(kdpval .eq. 1) then
                    ftop(i,k,1) = ( finty(i,k,3) - finty(i,k,2) )*rdz(i,k)
                    fbot(i,k,1) = wdz(1,1,      2)*finty(i,k,2) + &
                                  wdz(2,1,      2)*finty(i,k,3) + &
                                  wdz(3,1,      2)*finty(i,k,4) + &
                                  wdz(4,1,      2)*finty(i,k,1)
                 else if(kdpval .eq. kdimm1) then
                    ftop(i,k,1) = wdz(1,2,kdimm2 )*finty(i,k,4) + &
                                  wdz(2,2,kdimm2 )*finty(i,k,1) + &
                                  wdz(3,2,kdimm2 )*finty(i,k,2) + &
                                  wdz(4,2,kdimm2 )*finty(i,k,3)
                    fbot(i,k,1) = 0._r8
                 else
                    ftop(i,k,1) = wdz(1,1,kdpval )*finty(i,k,1) + &
                                  wdz(2,1,kdpval )*finty(i,k,2) + &
                                  wdz(3,1,kdpval )*finty(i,k,3) + &
                                  wdz(4,1,kdpval )*finty(i,k,4)
                    fbot(i,k,1) = wdz(1,2,kdpval )*finty(i,k,1) + &
                                  wdz(2,2,kdpval )*finty(i,k,2) + &
                                  wdz(3,2,kdpval )*finty(i,k,3) + &
                                  wdz(4,2,kdpval )*finty(i,k,4)
                 endif
              end do
           end do
        end do
!
! Apply SCM0 limiter to derivative estimates.
!
        do k=1,plev
           do i=1,nlon
              deli = ( finty(i,k,3) - finty(i,k,2) )*rdz(i,k)
              tmp1 = fac*deli
              tmp2 = abs( tmp1 )
              if( deli*fbot(i,k,1)   .le. 0.0_r8  ) fbot(i,k,1) = 0._r8
              if( deli*ftop(i,k,1)   .le. 0.0_r8  ) ftop(i,k,1) = 0._r8
              if( abs( fbot(i,k,1) ) .gt. tmp2 ) fbot(i,k,1) = tmp1
              if( abs( ftop(i,k,1) ) .gt. tmp2 ) ftop(i,k,1) = tmp1
              fdp(i,k) = finty(i,k,2)*ht(i,k) + ftop(i,k,1)*dht(i,k) + &
                         finty(i,k,3)*hb(i,k) + fbot(i,k,1)*dhb(i,k)
           end do
        end do
        if (icount.ne.nlon*plev) then
           call endrun ('SLTINT:  Did not complete computations for all departure points')
        endif
     endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!    40XX loops: Lagrange cubic/linear interpolation in the vertical
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
     if( .not. limdrv ) then
        do k=1,plev
           do i=1,nlon
              fdp(i,k) = finty(i,k,1)*wgt1z(i,k) &
                       + finty(i,k,2)*wgt2z(i,k) &
                       + finty(i,k,3)*wgt3z(i,k) &
                       + finty(i,k,4)*wgt4z(i,k)
           end do
        end do
!
! IF the departure point is in either the top or bottom interval of the
! model grid:  THEN
! Perform Hermite cubic interpolation.  The data are shifted up or down
! such that the departure point sits in the middle interval of the
! 4 point stencil (the shift originally took place in routine "LAGXIN").
! Therefore the derivative weights must be applied appropriately to
! account for this shift.  The following overwrites some results from
! the previous loop.
!
        do k=1,plev
           do i=1,nlon
              if(kdp (i,k) .eq. 1) then
                 tmptop = ( finty(i,k,3) - finty(i,k,2) )*rdz(i,k)
                 tmpbot = wdz(1,1,     2)*finty(i,k,2) &
                        + wdz(2,1,     2)*finty(i,k,3) &
                        + wdz(3,1,     2)*finty(i,k,4) &
                        + wdz(4,1,     2)*finty(i,k,1)
                 fdp(i,k) = finty(i,k,2)*ht (i,k) + tmptop      *dht(i,k) &
                          + finty(i,k,3)*hb (i,k) + tmpbot      *dhb(i,k)
              else if(kdp (i,k) .eq. kdimm1) then
                 tmptop = wdz(1,2,kdimm2)*finty(i,k,4) &
                        + wdz(2,2,kdimm2)*finty(i,k,1) &
                        + wdz(3,2,kdimm2)*finty(i,k,2) &
                        + wdz(4,2,kdimm2)*finty(i,k,3)
!!!!!            tmpbot = 0.
                 fdp(i,k) = finty(i,k,2)*ht (i,k) + tmptop       *dht(i,k) &
                          + finty(i,k,3)*hb (i,k)
!!!!!                     + tmpbot      *dhb(i,k)
              end if
           end do
        end do
     end if
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  Horizontal interpolation only
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  elseif(lhrzint) then
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!    50XX loops: an optimized Lagrange cubic/linear algorithm (no
!                Hermite interpolator available)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
     if(limdrh) then
!
! PART 1:  x-interpolation
!
        do k=1,plev
           do i = 1,nlon
              ii1 = idp(i,k,1)
              ii2 = idp(i,k,2)
              ii3 = idp(i,k,3)
              ii4 = idp(i,k,4)
              jj = jdp(i,k)
              kk = kdp(i,k)
!
! Height level 2
!
!   Latitude 1:  Linear interpolation
!
              fintx(i,k,1,2) = fb (ii1  ,kk  ,jj-1)*xl (i,k,1) &
                             + fb (ii1+1,kk  ,jj-1)*xr (i,k,1)
!
!   Latitude 2:  Cubic interpolation
!
              fxl = (   - 2._r8*fb (ii2-1,kk  ,jj  ) &
                        - 3._r8*fb (ii2  ,kk  ,jj  ) &
                        + 6._r8*fb (ii2+1,kk  ,jj  ) &
                        -    fb (ii2+2,kk  ,jj  ) )*rdx6(jj)
              fxr = (        fb (ii2-1,kk  ,jj  ) &
                        - 6._r8*fb (ii2  ,kk  ,jj  ) &
                        + 3._r8*fb (ii2+1,kk  ,jj  ) &
                        + 2._r8*fb (ii2+2,kk  ,jj  ) )*rdx6(jj)
!
              deli = (       fb (ii2+1,kk  ,jj  ) - &
                             fb (ii2  ,kk  ,jj  ) )*rdx(jj)
              tmp1 = fac*deli
              tmp2 = abs( tmp1 )
              if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
              if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
              if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
              if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
              fintx(i,k,2,2) = fb (ii2  ,kk  ,jj  )*hl (i,k,2) &
                             + fb (ii2+1,kk  ,jj  )*hr (i,k,2) &
                             + fxl*dhl(i,k,2) + fxr*dhr(i,k,2)
!
!   Latitude 3:  Cubic interpolation
!
              fxl = (   - 2._r8*fb (ii3-1,kk  ,jj+1) &
                        - 3._r8*fb (ii3  ,kk  ,jj+1) &
                        + 6._r8*fb (ii3+1,kk  ,jj+1) &
                        -    fb (ii3+2,kk  ,jj+1) )*rdx6(jj+1)
              fxr = (        fb (ii3-1,kk  ,jj+1) &
                        - 6._r8*fb (ii3  ,kk  ,jj+1) &
                        + 3._r8*fb (ii3+1,kk  ,jj+1) &
                        + 2._r8*fb (ii3+2,kk  ,jj+1) )*rdx6(jj+1)
!
              deli = (       fb (ii3+1,kk  ,jj+1) - &
                             fb (ii3  ,kk  ,jj+1) )*rdx(jj+1)
              tmp1 = fac*deli
              tmp2 = abs( tmp1 )
              if( deli*fxl   .le. 0.0_r8  ) fxl = 0._r8
              if( deli*fxr   .le. 0.0_r8  ) fxr = 0._r8
              if( abs( fxl ) .gt. tmp2 ) fxl = tmp1
              if( abs( fxr ) .gt. tmp2 ) fxr = tmp1
!
              fintx(i,k,3,2) = fb (ii3  ,kk  ,jj+1)*hl (i,k,3) &
                             + fb (ii3+1,kk  ,jj+1)*hr (i,k,3) &
                             + fxl*dhl(i,k,3) + fxr*dhr(i,k,3)
!
!   Latitude 4:  Linear interpolation
!
              fintx(i,k,4,2) = fb (ii4  ,kk  ,jj+2)*xl (i,k,4) &
                             + fb (ii4+1,kk  ,jj+2)*xr (i,k,4)
           end do
        end do
!
! PART 2:  y-derivatives
!
        jmin =  1000000
        jmax = -1000000
        do k=1,plev
           do i = 1,nlon
              if(jdp(i,k) .lt. jmin) jmin = jdp(i,k)
              if(jdp(i,k) .gt. jmax) jmax = jdp(i,k)
           end do
        end do
!
! Loop over departure latitudes
!
        icount = 0
        do jdpval = jmin,jmax
           do k=1,plev
              call whenieq(nlon    ,jdp(1,k),1       ,jdpval  , &
                           indx    ,nval    )
              icount = icount + nval
!
! y derivatives at the inner height levels (kk = 2,3) needed for
! z-interpolation
!
              do kk  = 2,2
                 do ii = 1,nval
                    i = indx(ii)
                    fbot(i,k,kk) = lbasdy(1,1,jdpval)*fintx(i,k,1,kk) &
                                 + lbasdy(2,1,jdpval)*fintx(i,k,2,kk) &
                                 + lbasdy(3,1,jdpval)*fintx(i,k,3,kk) &
                                 + lbasdy(4,1,jdpval)*fintx(i,k,4,kk)
                    ftop(i,k,kk) = lbasdy(1,2,jdpval)*fintx(i,k,1,kk) &
                                 + lbasdy(2,2,jdpval)*fintx(i,k,2,kk) &
                                 + lbasdy(3,2,jdpval)*fintx(i,k,3,kk) &
                                 + lbasdy(4,2,jdpval)*fintx(i,k,4,kk)
                 end do
              end do
           end do
        end do
        if (icount.ne.nlon*plev) then
           call endrun ('SLTINT:  Did not complete computations for all departure points')
        end if
!
! Apply SCM0 limiter to derivative estimates.
!
        do kk  = 2,2
           do k=1,plev
              do i = 1,nlon
                 deli = ( fintx(i,k,3,kk) - fintx(i,k,2,kk) )*rdphi(i,k)
                 tmp1 = fac*deli
                 tmp2 = abs( tmp1 )
                 if( deli*fbot(i,k,kk)   .le. 0.0_r8  ) fbot(i,k,kk) = 0._r8
                 if( deli*ftop(i,k,kk)   .le. 0.0_r8  ) ftop(i,k,kk) = 0._r8
                 if( abs( fbot(i,k,kk) ) .gt. tmp2 ) fbot(i,k,kk) = tmp1
                 if( abs( ftop(i,k,kk) ) .gt. tmp2 ) ftop(i,k,kk) = tmp1
              end do
           end do
        end do
!
! PART 3:  y-interpolants
!
        do k=1,plev
           do i = 1,nlon
              fdp(i,k) = fintx(i,k,2,2)*hs (i,k) + fbot (i,k,2)*dhs(i,k) &
                       + fintx(i,k,3,2)*hn (i,k) + ftop (i,k,2)*dhn(i,k)
           end do
        end do
     endif
!
     if( .not. limdrh ) then
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!    60XX loops: Hermite cubic/linear interpolation in the horizontal
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
        do k=1,plev
           do i=1,nlon
              ii1 = idp(i,k,1)
              ii2 = idp(i,k,2)
              ii3 = idp(i,k,3)
              ii4 = idp(i,k,4)
              jj = jdp(i,k)
              kk = kdp(i,k)
!
! x-interpolants for the 4 latitudes
!
              f1 = fb(ii1+1,kk,jj-1)*xr   (i,k,1) &
                 + fb(ii1  ,kk,jj-1)*xl   (i,k,1)
              f2 = fb(ii2-1,kk,jj  )*wgt1x(i,k,2) &
                 + fb(ii2  ,kk,jj  )*wgt2x(i,k,2) &
                 + fb(ii2+1,kk,jj  )*wgt3x(i,k,2) &
                 + fb(ii2+2,kk,jj  )*wgt4x(i,k,2)
              f3 = fb(ii3-1,kk,jj+1)*wgt1x(i,k,3) &
                 + fb(ii3  ,kk,jj+1)*wgt2x(i,k,3) &
                 + fb(ii3+1,kk,jj+1)*wgt3x(i,k,3) &
                 + fb(ii3+2,kk,jj+1)*wgt4x(i,k,3)
              f4 = fb(ii4+1,kk,jj+2)*xr   (i,k,4) &
                 + fb(ii4  ,kk,jj+2)*xl   (i,k,4)
!
! y-interpolant
!
              fdp(i,k) = f1*wgt1y(i,k) + f2*wgt2y(i,k) + &
                         f3*wgt3y(i,k) + f4*wgt4y(i,k)
           end do
        end do
     end if
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
!  Vertical interpolation only
!
!    70XX loops: an optimized Lagrange cubic/linear algorithm (no
!                Hermite interpolator available)
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  else if(lvrtint) then
     if(limdrh .or. limdrv) then
        call endrun ('SLTINT: Shape preserving capability not provided for vertical-only interpolation')
     end if

     do k=1,plev
        do i=1,nlon
           kk = kkdp(i,k)
           ii = i1+i-1
           fdp(i,k) = fb(ii,kk-1,jcen)*wgt1z(i,k) &
                    + fb(ii,kk  ,jcen)*wgt2z(i,k) &
                    + fb(ii,kk+1,jcen)*wgt3z(i,k) &
                    + fb(ii,kk+2,jcen)*wgt4z(i,k)
        end do
     end do
!
! IF the departure point is in either the top or bottom interval of the
! model grid:  THEN perform Hermite cubic interpolation.  The following
! overwrites some results from the previous loop.
!
     do k=1,plev
        do i=1,nlon
           ii = i1+i-1
           if(kdp (i,k) .eq. 1) then
!!!!!         tmptop   = 0.
              tmpbot   = wdz(1,1,     2)*fb(ii,     1,jcen) &
                       + wdz(2,1,     2)*fb(ii,     2,jcen) &
                       + wdz(3,1,     2)*fb(ii,     3,jcen) &
                       + wdz(4,1,     2)*fb(ii,     4,jcen)
              fdp(i,k) = fb(ii,1     ,jcen)*ht (i,k) &
                    &  + fb(ii,2     ,jcen)*hb (i,k) &
                    &  + tmpbot            *dhb(i,k)
!!!!!               &  + tmptop            *dht(i,k)
           else if(kdp (i,k) .eq. kdimm1) then
              tmptop = wdz(1,2,kdimm2)*fb(ii,kdimm3,jcen) &
                     + wdz(2,2,kdimm2)*fb(ii,kdimm2,jcen) &
                     + wdz(3,2,kdimm2)*fb(ii,kdimm1,jcen) &
                     + wdz(4,2,kdimm2)*fb(ii,kdim  ,jcen)
!!!!!         tmpbot = 0.
              fdp(i,k) = fb(ii,kdimm1,jcen)*ht (i,k) &
                       + tmptop            *dht(i,k) &
                       + fb(ii,kdim  ,jcen)*hb (i,k)
!!!!!                  + tmpbot            *dhb(i,k)
           end if
        end do
     end do
  else
     call endrun ('SLTINT:  Error: must specify at least one of lhrzint or lvrtint to be .true.')
  end if
!
  return
end subroutine sltint

