
subroutine hrintp(pf      ,pkcnst  ,fb      ,fxl     ,fxr     , &
                  x       ,y       ,dy      ,wdy     ,xdp     , &
                  ydp     ,idp     ,jdp     ,jcen    ,limitd  , &
                  fint    ,fyb     ,fyt     ,fdp     ,nlon    , &
                  nlonex  )

!----------------------------------------------------------------------- 
! 
! Purpose: 
! Interpolate 2-d field to departure point using tensor product
! Hermite cubic interpolation.
! 
! Method: 
! 
! Author: 
! Original version:  J. Olson
! Standardized:      J. Rosinski, June 1992
! Reviewed:          D. Williamson, P. Rasch, August 1992
! Reviewed:          D. Williamson, P. Rasch, March 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plev, plon
   use scanslt,      only: plond, platd, beglatex, endlatex
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
#include <parslt.h>
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: pf                   ! dimension (number of fields)
   integer, intent(in) :: pkcnst               ! dimension (see ext. document)
!
   real(r8), intent(in) :: fb (plond,plev,pkcnst,beglatex:endlatex) ! input fields
   real(r8), intent(in) :: fxl(plond,plev,pf    ,beglatex:endlatex) ! left  x-derivs
   real(r8), intent(in) :: fxr(plond,plev,pf    ,beglatex:endlatex) ! right x-derivs
   real(r8), intent(in) :: x  (plond,platd)        ! long. grid coordinates
   real(r8), intent(in) :: y  (platd)              ! lat.  grid coordinates
   real(r8), intent(in) :: dy (platd)              ! intervals betwn lat grid pts.
   real(r8), intent(in) :: wdy(4,2,platd)          ! lat. derivative weights
   real(r8), intent(in) :: xdp(plon,plev)          ! x-coord of dep. pt.
   real(r8), intent(in) :: ydp(plon,plev)          ! y-coord of dep. pt.
!
   integer, intent(in) :: idp(plon,plev,4)     ! i index of dep. pt.
   integer, intent(in) :: jdp(plon,plev)       ! j index of dep. pt.
   integer, intent(in) :: jcen
!
   logical, intent(in) :: limitd               ! flag for shape-preservation
!
! Output arguments
!
   real(r8), intent(out) :: fint(plon,plev,ppdy,pf) ! x interpolants
   real(r8), intent(out) :: fyb (plon,plev,pf)      ! y-derivatives at bot of int.
   real(r8), intent(out) :: fyt (plon,plev,pf)      ! y-derivatives at top of int.
   real(r8), intent(out) :: fdp (plon,plev,pf)      ! horizontal interpolants

   integer, intent(in) :: nlon
   integer, intent(in) :: nlonex(platd)
!
!-----------------------------------------------------------------------
!
!  pf      Number of fields being interpolated.
!  pkcnst  dimensioning construct for 3-D arrays. (see ext. document)
!  fb      Extended array of data to be interpolated.
!  fxl     x-derivatives at the left  edge of each interval containing
!          the departure point.
!  fxr     x-derivatives at the right edge of each interval containing
!          the departure point.
!  x       Equally spaced x grid values in extended arrays.
!  y       y-coordinate (latitude) values in the extended array.
!  dy      Increment in the y-coordinate value for each interval in the
!          extended array.
!  wdy     Weights for Lagrange cubic derivative estimates on the
!          unequally spaced y-grid.  If grid interval j (in extended
!          array is surrounded by a 4 point stencil, then the
!          derivative at the "bottom" of the interval uses the weights
!          wdy(1,1,j),wdy(2,1,j), wdy(3,1,j), and wdy(4,1,j).  The
!          derivative at the "top" of the interval uses wdy(1,2,j),
!          wdy(2,2,j), wdy(3,2,j) and wdy(4,2,j).
!  xdp     xdp(i,k) is the x-coordinate of the departure point that
!          corresponds to global grid point (i,k) in the latitude slice
!          being forecasted.
!  ydp     ydp(i,k) is the y-coordinate of the departure point that
!          corresponds to global grid point (i,k) in the latitude slice
!          being forecasted.
!  idp     idp(i,k) is the index of the x-interval that contains the
!          departure point corresponding to global grid point (i,k) in
!          the latitude slice being forecasted.
!          Note that
!                x(idp(i,k)) .le. xdp(i,k) .lt. x(idp(i,k)+1) .
!  jdp     jdp(i,k) is the index of the y-interval that contains the
!          departure point corresponding to global grid point (i,k) in
!          the latitude slice being forecasted.
!          Suppose yb contains the y-coordinates of the extended array
!          and ydp(i,k) is the y-coordinate of the departure point
!          corresponding to grid point (i,k).  Then,
!                yb(jdp(i,k)) .le. ydp(i,k) .lt. yb(jdp(i,k)+1) .
!  limitd  Logical flag to specify whether or not the y-derivatives will
!          be limited.
!  fint    WORK ARRAY, results not used on return
!  fyb     WORK ARRAY, results not used on return
!  fyt     WORK ARRAY, results not used on return
!  fdp     Value of field at the horizontal departure points.
!
!-----------------------------------------------------------------------
!
! Hermite cubic interpolation to the x-coordinate of each
! departure point at each y-coordinate required to compute the
! y-derivatives.
!
   call herxin(pf      ,pkcnst  ,fb      ,fxl     ,fxr     , &
               x       ,xdp     ,idp     ,jdp     ,fint    , &
               nlon    ,nlonex  )
!
! Compute y-derivatives.
!
   call cubydr(pf      ,fint    ,wdy     ,jdp     ,jcen    , &
               fyb     ,fyt     ,nlon    )
   if( limitd )then
      call limdy(pf    ,fint    ,dy      ,jdp     ,fyb     , &
                 fyt   ,nlon    )
   end if
!
! Hermite cubic interpolation in the y-coordinate.
!
   call heryin(pf      ,fint    ,fyb     ,fyt     ,y       , &
               dy      ,ydp     ,jdp     ,fdp     ,nlon    )
!
   return
end subroutine hrintp
