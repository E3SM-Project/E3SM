module gauaw_mod
!-----------------------------------------------------------------------
!
! Purpose:
!
!    Module to calculate the Gaussian Weights. Public interface is
!    the subroutine "gauaw( a, w, k )".
!
! Method: 
!
!	The algorithm is described in Davis and Rabinowitz,
!      Journal of Research of the NBS, V 56, Jan 1956.
!
! Author: David Williamson, Jim Hack
!
!-----------------------------------------------------------------------
   use abortutils, only: endrun
   use shr_kind_mod, only: r8 => shr_kind_r8

#ifdef NO_R16
   integer,parameter :: r16= selected_real_kind(12) ! 8 byte real
#else
   integer,parameter :: r16= selected_real_kind(24) ! 16 byte real
#endif

   save

!
! Public subroutines
!

   public gauaw

!
! Variables private to routines inside this module
!

   real(r16), private ::    pi           ! value of pi
   real(r16), private, parameter :: one = 1.0_r16     ! 1. in real(r16).  Needed by atan

!
! Functions private to routines inside this module
!

   private bsslzr

contains

   subroutine gauaw(a, w, k)
!-----------------------------------------------------------------------
!
! Calculate sine of latitudes a(k) and weights w(k) for the gaussian
! quadrature. The algorithm is described in Davis and Rabinowitz,
! Journal of Research of the NBS, V 56, Jan 1956.
! The zeros of the bessel function j0, which are obtained from bsslzr,
! are used as a first guess for the abscissa.
!
! Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic in order to
! achieve (nearly) identical weights and latitudes on all machines.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      L. Bath, Jun 1992
!                    L. Buja, Feb 1996
! Reviewed:          D. Williamson, J. Hack, Aug 1992
!                    D. Williamson, J. Hack, Feb 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   implicit none
!------------------------------Arguments--------------------------------
!
   integer , intent(in)  :: k      ! number of latitudes pole to pole
   real(r8), intent(out) :: a(k)   ! sine of latitudes
   real(r8), intent(out) :: w(k)   ! gaussian weights
!
!---------------------------Local workspace-----------------------------
!
   real(r16) sinlat(k)    ! sine of latitudes
   real(r16) wgt(k)       ! gaussian weights
 
   real(r16) eps          ! convergence criterion
   real(r16) c            ! constant combination
   real(r16) fk           ! real k
   real(r16) xz           ! abscissa estimate
   real(r16) pkm1         ! |
   real(r16) pkm2         ! |-polynomials
   real(r16) pkmrk        ! |
   real(r16) pk           ! |
   real(r16) sp           ! current iteration latitude increment
   real(r16) avsp         ! |sp|
   real(r16) fn           ! real n
#ifdef NO_R16
   parameter (eps = 1.e-15_r16)
#else
   parameter (eps = 1.e-27_r16)
#endif
 
   integer kk           ! k/2 (number of latitudes in hemisphere)
   integer is           ! latitude index
   integer iter         ! iteration counter
   integer n,l          ! indices
!
!-----------------------------------------------------------------------
!
   pi  = 4._r16*atan(one)
!
! The value eps, used for convergence tests in the iterations,
! can be changed.  Newton iteration is used to find the abscissas.
!
   c = (1._r16-(2._r16/pi)**2)*0.25_r16
   fk = k
   kk = k/2
   call bsslzr(sinlat,kk)
   do is=1,kk
     xz = cos(sinlat(is)/(((fk+0.5_r16)**2+c)**0.5_r16))
!
! This is the first approximation to xz
!
     iter = 0
 10  continue
     pkm2 = 1._r16
     pkm1 = xz
     iter = iter + 1
     if (iter.gt.10) then
!
! Error exit
!
       call endrun ('GAUAW: no convergence in 10 iterations')
     end if
!
! Computation of the legendre polynomial
!
     do n=2,k
       fn = n
       pk = ((2._r16*fn-1._r16)*xz*pkm1-(fn-1._r16)*pkm2)/fn
       pkm2 = pkm1
       pkm1 = pk
     enddo
     pkm1 = pkm2
     pkmrk = (fk*(pkm1-xz*pk))/(1._r16-xz**2)
     sp = pk/pkmrk
     xz = xz - sp
     avsp = abs(sp)
     if (avsp.gt.eps) go to 10
     sinlat(is) = xz
     wgt(is) = (2._r16*(1._r16-xz**2))/(fk*pkm1)**2
   end do
!
   if (k.ne.kk*2) then
!
! For odd k computation of weight at the equator
!
     sinlat(kk+1) = 0._r16
     pk = 2._r16/fk**2
     do n=2,k,2
       fn = n
       pk = pk*fn**2/(fn-1._r16)**2
     end do
     wgt(kk+1) = pk
   end if
!
! Complete the sets of abscissas and weights, using the symmetry.
! Also note truncation from real(r16) to real*8
!
   do n=1,kk
     l = k + 1 - n
     a(n) = sinlat(n)
     a(l) = -sinlat(n)
 
     w(n) = wgt(n)
     w(l) = wgt(n)
   end do
 
   return
   end subroutine gauaw


 
!===========================================================================
 
   subroutine bsslzr(bes, n) 
!-----------------------------------------------
!   M o d u l e s 
!-----------------------------------------------
!
! Return n zeros (or if n>50, approximate zeros), of the Bessel function
! j0,in the array bes. The first 50 zeros will be given exactly, and the
! remaining zeros are computed by extrapolation,and therefore not exact.
!
! Modified 1/23/97 by Jim Rosinski to use real*16 arithmetic
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, D. Williamson, August 1992
! Reviewed:          J. Hack, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
   integer , intent(in) :: n 
   real(r16) , intent(inout) :: bes(n) 
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
   integer :: j, nn 
   real(r16), dimension(50) :: bz 

   save bz 
!-----------------------------------------------
!------------------------------Arguments--------------------------------
!
!
!---------------------------Local workspace-----------------------------
!
!
   data bz/ 2.4048255577_r16, 5.5200781103_r16, 8.6537279129_r16, 11.7915344391_r16, &
      14.9309177086_r16, 18.0710639679_r16, 21.2116366299_r16, 24.3524715308_r16, &
      27.4934791320_r16, 30.6346064684_r16, 33.7758202136_r16, 36.9170983537_r16, &
      40.0584257646_r16, 43.1997917132_r16, 46.3411883717_r16, 49.4826098974_r16, &
      52.6240518411_r16, 55.7655107550_r16, 58.9069839261_r16, 62.0484691902_r16, &
      65.1899648002_r16, 68.3314693299_r16, 71.4729816036_r16, 74.6145006437_r16, &
      77.7560256304_r16, 80.8975558711_r16, 84.0390907769_r16, 87.1806298436_r16, &
      90.3221726372_r16, 93.4637187819_r16, 96.6052679510_r16, 99.7468198587_r16, &
      102.8883742542_r16, 106.0299309165_r16, 109.1714896498_r16, 112.3130502805_r16, &
      115.4546126537_r16, 118.5961766309_r16, 121.7377420880_r16, 124.8793089132_r16, &
      128.0208770059_r16, 131.1624462752_r16, 134.3040166383_r16, 137.4455880203_r16, &
      140.5871603528_r16, 143.7287335737_r16, 146.8703076258_r16, 150.0118824570_r16, &
      153.1534580192_r16, 156.2950342685_r16/  
!
   nn = n 
   if (n > 50) then 
      bes(50) = bz(50) 
      do j = 51, n 
         bes(j) = bes(j-1) + pi 
      end do 
      nn = 49 
   endif 
   bes(:nn) = bz(:nn) 
   return  
   end subroutine bsslzr 

end module gauaw_mod
