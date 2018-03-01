subroutine hordif(k,ztdt)

!----------------------------------------------------------------------- 
! 
! Purpose: 
! 
! Method: 
! Horizontal diffusion of z,d,t,q
! 1. implicit del**2 form above level kmnhd4
! 2. implicit del**4 form at level kmnhd4 and below
! 3. courant number based truncation at level kmxhdc and above
! 4. increased del**2 coefficient at level kmxhd2 and above
!
! Computational note: this routine is multitasked by level, hence it 
! is called once for each k
!
! Author: 
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, J. Hack, August 1992
! Reviewed:          B. Boville, April 1996
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
   use time_manager, only: get_step_size, is_first_step, get_nstep
   use eul_control_mod
   use spmd_utils, only : iam
!-----------------------------------------------------------------------
   implicit none
!------------------------------Arguments--------------------------------
!
! Input arguments
!
   integer, intent(in) :: k            ! level index

   real(r8), intent(in) :: ztdt        ! 2 times time step unless nstep=0
!
!---------------------------Local workspace-----------------------------
!
   integer ir,ii        ! spectral indices       
   integer lmr,lmc      ! spectral indices
   real(r8) dfac            ! large coefficient on del^n multipliers to
!                            strongly damp waves req'd by Courant limiter
   integer  lm,m,n          ! spectral indices
   real(r8) ztodt           ! 2 delta t
   real(r8) zdt             ! model time step
   real(r8) dmpini          ! used to compute divergence damp rate
   real(r8) dmptim          ! used to compute divergence damp rate
   real(r8) dmprat          ! divergence damping rate
   real(r8) coef            ! coeff. used to apply damping rate to divergence
   real(r8) two
!
!-----------------------------------------------------------------------
!DIR$ NOSTREAM
   two=2._r8
!
! Set the horizontal diffusion factors for each wavenumer at this level
! depending on: whether del^2 or del^4 diffusion is to be applied; and 
! whether the courant number limit is to be applied.
!
   if (k .ge. kmnhd4) then        ! Del^4 diffusion factors
      do n=1,pnmax
         hdiftq(n,k) = hdfst4(n)
         hdifzd(n,k) = hdfsd4(n)
      end do
!
! Spectrally truncate selected levels (if courant number too large)
!
      if (k.le. kmxhdc .and. nindex(k).le.pnmax) then
         dfac = 1000._r8
         do n=nindex(k),pnmax
            hdiftq(n,k) = dfac*hdfst4(n)
            hdifzd(n,k) = dfac*hdfsd4(n)
         end do
      end if
   else                      ! Del^2 diffusion factors
      if (k.le.kmxhd2) then
!
! Buggy sun compiler gives wrong answer for following line when
! using  -Qoption f90comp -r8const flags
!          dfac = 2.**(real(kmxhd2-k+1,r8))
         dfac = two**(real(kmxhd2-k+1,r8))
      else
         dfac = 1.0_r8
      end if
      do n=1,pnmax
         hdiftq(n,k) = dfac*hdfst2(n)
         hdifzd(n,k) = dfac*hdfsd2(n)
      end do
!
! Spectrally truncate selected levels (if courant number too large)
!
      if ((k.le.kmxhdc).and.(nindex(k).le.pnmax)) then
         dfac = 1000._r8
         do n=nindex(k),pnmax
            hdiftq(n,k) = dfac*hdfst2(n)
            hdifzd(n,k) = dfac*hdfsd2(n)
         end do
      end if
   end if
!
! Define damping rate for divergence damper
!
   zdt = get_step_size()

!   ztodt = 2._r8*zdt
!   if (is_first_step()) ztodt = .5_r8*ztodt
    ztodt = ztdt	
!
! Initial damping rate (e-folding time = zdt) and then linearly decrease
! to 0. over number of days specified by "divdampn".
!
   coef = 1._r8
   if (divdampn .gt. 0.0_r8) then
      dmpini = 1._r8/(zdt)
      dmptim = divdampn*86400._r8
      dmprat = dmpini * (dmptim - real(get_nstep(),r8)*zdt) / dmptim
      if (dmprat .gt. 0.0_r8) coef = 1.0_r8 / (1.0_r8+ztodt*dmprat)
   endif
!
! Compute time-split implicit factors for this level
!
   do lm=1,numm(iam)
      m=locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
      do n=1,nlen(m)
         ir = lmc + 2*n - 1
         ii = ir + 1
!
! time-split implicit factors
!
         t(ir,k) = t(ir,k)/(1._r8 + ztdt*hdiftq(n+m-1,k))
         t(ii,k) = t(ii,k)/(1._r8 + ztdt*hdiftq(n+m-1,k))
!
         d(ir,k) = d(ir,k)*coef/(1._r8 + ztdt*hdifzd(n+m-1,k))
         d(ii,k) = d(ii,k)*coef/(1._r8 + ztdt*hdifzd(n+m-1,k))
!
         vz(ir,k) = vz(ir,k)/(1._r8 + ztdt*hdifzd(n+m-1,k))
         vz(ii,k) = vz(ii,k)/(1._r8 + ztdt*hdifzd(n+m-1,k))
      end do
   end do
!
   return
end subroutine hordif

