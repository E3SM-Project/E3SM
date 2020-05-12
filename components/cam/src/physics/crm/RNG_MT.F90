! Mersenne Twister Uniform Random Number Generator ((0,1) Real Number)
!
! call RNG_MT_set_seed(seed) to set the seed (integer)
! If you do not call this, the seed is 4357.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This library is free software; you can redistribute it and/or
! modify it under the terms of the GNU Library General Public
! License as published by the Free Software Foundation; either
! version 2 of the License, or (at your option) any later
! version.
! This library is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
! See the GNU Library General Public License for more details.
! You should have received a copy of the GNU Library General
! Public License along with this library; if not, write to the
! Free Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
! 02111-1307  USA
!
! Copyright (C) 1997 Makoto Matsumoto and Takuji Nishimura.
! When you use this, send an email to: matumoto@math.keio.ac.jp
! with an appropriate reference to your work.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fortran translation by Hiroshi Takano.  Jan. 13, 1999.
!
!   genrand()      -> double precision function grnd()
!   sgenrand(seed) -> subroutine sgrnd(integer seed)
!
! This program uses the following non-standard intrinsics.
!   ishft(i,n): If n>0, shifts bits in i by n positions to left.
!               If n<0, shifts bits in i by n positions to right.
!   iand (i,j): Performs logical AND on corresponding bits of i and j.
!   ior  (i,j): Performs inclusive OR on corresponding bits of i and j.
!   ieor (i,j): Performs exclusive OR on corresponding bits of i and j.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Fortran90 translation by Yusuke Konishi. May. 13, 2013.
!
! The interval is changed from [0, 1] to [0, 1).
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Update by Walter Hannah - March, 2018
! - restructured code to be more consistent with E3SM 
! - changed the interval to (0,1)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
module RNG_MT
   use shr_kind_mod,  only: r8 => shr_kind_r8

   implicit none

   integer, private, parameter :: N      = 624      ! why are these parameters so big?
   integer, private, parameter :: N1     = 625
   integer, private, parameter :: M      = 397
   integer, private, parameter :: MATA   = -1727483681
   integer, private, parameter :: UMASK  = -1147483648 ! originally = -2147483648, but this is too big
   integer, private, parameter :: LMASK  = 2147483647
   integer, private, parameter :: TMASKB = -1658038656
   integer, private, parameter :: TMASKC = -272236544
   
   integer, private :: mti
   integer, private, dimension(2)     :: mag01(0:1) = (/0, MATA/)
   integer, private, dimension(0:N-1) :: mt 

   public RNG_MT_set_seed
   public RNG_MT_gen_rand
   
   contains

!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
! set initial seeds to mt[N] using the generator - Line 25 of Table 1 in
! Knuth 1981, The Art of Computer Programming Vol. 2 (2nd Ed.), pp102
!
subroutine RNG_MT_set_seed(seed)

   integer,intent(in) :: seed

   !!! initialize mt
   mt = N1

   !!! set seed
   mt(0) = iand(seed, -1)
   do mti = 1, N - 1
      mt(mti) = iand(69069 * mt(mti - 1), -1)
   end do

end subroutine RNG_MT_set_seed
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
subroutine RNG_MT_gen_rand(grnd)
   real(r8), intent(out) :: grnd
   integer :: y, kk
   
   if(mti >= N) then
      !!! generate N words at one time
      if(mti == N + 1) then
         !!! if RNG_MT_set_seed() has not been called, a default initial seed is used
         call RNG_MT_set_seed(4357)
         !!! that also means we need to initialize mt
         mt = N1
      endif

      do kk = 0, N - M - 1
         y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
         mt(kk) = ieor(ieor(mt(kk + M), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      do kk = N - M, N - 2
         y = ior(iand(mt(kk), UMASK), iand(mt(kk + 1), LMASK))
         mt(kk) = ieor(ieor(mt(kk + (M - N)), ishft(y, -1)), mag01(iand(y, 1)))
      end do

      y = ior(iand(mt(N - 1), UMASK), iand(mt(0), LMASK))
      mt(N - 1) = ieor(ieor(mt(M - 1), ishft(y, -1)), mag01(iand(y, 1)))
      mti = 0
   endif

   y = mt(mti)
   mti = mti + 1
   y = ieor(y,      ishft(y, -11)         )
   y = ieor(y, iand(ishft(y,   7), TMASKB))
   y = ieor(y, iand(ishft(y,  15), TMASKC))
   y = ieor(y,      ishft(y, -18)         )

   if(y < 0) then
      grnd = (dble(y) + 2.0d0 ** 32) / (2.0d0 ** 32)
   else
      ! grnd = dble(y) / (2.0d0 ** 32)      ! for interval [0,1)
      grnd = (dble(y) + 0.5d0)/(2.0d0**32)  ! for interval (0,1)
   endif

end subroutine RNG_MT_gen_rand
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
end module RNG_MT
