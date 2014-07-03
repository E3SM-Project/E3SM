subroutine trunc
!----------------------------------------------------------------------- 
! 
! Purpose: 
! Check consistency of truncation parameters and evaluate pointers
! and displacements for spectral arrays
! 
! Method: 
! 
! Author: 
! Original version:  CCM1
! Standardized:      L. Bath, June 1992
!                    T. Acker, March 1996
! Reviewed:          J. Hack, D. Williamson, August 1992
! Reviewed:          J. Hack, D. Williamson, April 1996
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use abortutils, only: endrun
!-----------------------------------------------------------------------
   implicit none
!-----------------------------------------------------------------------
!
!---------------------------Local variables-----------------------------
!
   integer m              ! loop index
!
!-----------------------------------------------------------------------
!
! trunc first evaluates truncation parameters for a general pentagonal 
! truncation for which the following parameter relationships are true
!
! 0 .le. |m| .le. ptrm
!
! |m| .le. n .le. |m|+ptrn for |m| .le. ptrk-ptrn
!
! |m| .le. n .le. ptrk     for (ptrk-ptrn) .le. |m| .le. ptrm
!
! Most commonly utilized truncations include:
!  1: triangular  truncation for which ptrk=ptrm=ptrn
!  2: rhomboidal  truncation for which ptrk=ptrm+ptrn
!  3: trapezoidal truncation for which ptrn=ptrk .gt. ptrm
!
! Simple sanity check
! It is necessary that ptrm .ge. ptrk-ptrn .ge. 0
!
   if (ptrm.lt.(ptrk-ptrn)) then
      call endrun ('TRUNC: Error in truncation parameters.  ntrm < (ptrk-ptrn)')
   end if
   if (ptrk.lt.ptrn) then
      call endrun ('TRUNC: Error in truncation parameters.  ptrk < ptrn')
   end if
!
! Evaluate pointers and displacement info based on truncation params
!
! The following ifdef logic seems to have something do with SPMD 
! implementation, although it's not clear how this info is used
! Dave, can you check this with JR?
!
   nstart(1) = 0
   nlen(1) = ptrn + 1
   do m=2,pmmax
      nstart(m) = nstart(m-1) + nlen(m-1)
      nlen(m) = min0(ptrn+1,ptrk+2-m)
   end do
!      write(iulog,*)'Starting index  length'
!      do m=1,ptrm+1
!         write(iulog,'(1x,i14,i8)')nstart(m),nlen(m)
!      end do
!
! Assign wavenumbers  and spectral offsets if not SPMD
!
#if ( ! defined SPMD )
   do m=1,pmmax
      locm(m,0) = m
      lnstart(m) = nstart(m)
   enddo
#endif
!
   return
end subroutine trunc
