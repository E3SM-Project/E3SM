
subroutine grcalcs (irow    ,ztodt   ,grts    ,grths   ,grds    ,&
                    grzs    ,grus    ,gruhs   ,grvs    ,grvhs   ,&
                    grpss   ,grdpss  ,grpms   ,grpls   ,tmpSPEcoef)
!-----------------------------------------------------------------------
!
! Complete inverse Legendre transforms from spectral to Fourier space at 
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric Fourier coeffs of temperature and
! "grtha" contains the antisymmetric Fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Modified:          P. Worley, October 2002
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use physconst, only: ez, ra
   use eul_control_mod
   use spmd_utils, only : iam
   implicit none

!
! Input arguments
!
   integer, intent(in) :: irow         ! latitude pair index
   real(r8), intent(in) :: ztodt       ! twice the timestep unless nstep = 0
   real(r8), intent(in) :: tmpSPEcoef(plev*24,pnmax,maxm) ! rearranged variables array
!
! Output arguments: symmetric fourier coefficients
!
   real(r8), intent(out) :: grts(2*maxm,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(out) :: grths(2*maxm,plev)   ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(out) :: grds(2*maxm,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(out) :: grzs(2*maxm,plev)    ! sum(n) of z(n,m)*P(n,m)
   real(r8), intent(out) :: grus(2*maxm,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: gruhs(2*maxm,plev)   ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1)) 
   real(r8), intent(out) :: grvs(2*maxm,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvhs(2*maxm,plev)   ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grpss(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grdpss(2*maxm)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grpms(2*maxm)        ! sum(n) of lnps(n,m)*H(n,m)
   real(r8), intent(out) :: grpls(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)*m/a
!
!---------------------------Local workspace-----------------------------
!
   real(r8) dalpn(pspt)         ! (a/(n(n+1)))*derivative of Legendre functions (complex)
   real(r8) zurcor              ! conversion term relating abs. & rel. vort.
   real(r8) tmpGRcoef(plev*24,maxm) ! temporal storage for Fourier coeffs

   integer k                ! level index
   integer lm, m            ! local and global Fourier wavenumber indices of spectral array
   integer mlength          ! number of local wavenumbers
   integer n                ! meridional wavenumber index
   integer ir,ii            ! spectral indices
   integer lmr,lmc          ! spectral indices
   integer lmwave0          ! local index for wavenumber 0
   integer lmrwave0         ! local offset for wavenumber 0
   integer kv               ! level x variable index
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
!
!DIR$ NOSTREAM
   lmwave0 = -1
   lmrwave0 = 0
   dalpn(2) = 0.0_r8
   mlength = numm(iam)
!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      if (m .eq. 1) then
         lmwave0 = lm
         lmrwave0 = lmr
      endif
      do n=1,nlen(m)
         dalpn(lmr+n) = ldalp(lmr+n,irow)*rsq(m+n-1)*ra
      end do
   end do
   zurcor = ez*dalpn(lmrwave0 + 2)
!
! Initialize sums
!
   grpss (:)   = 0._r8
   grpls (:)   = 0._r8
   grpms (:)   = 0._r8
   grdpss(:)   = 0._r8
!cdir collapse
   tmpGRcoef (:,:) = 0._r8
!
! Loop over n for t,q,d,and end of u and v
!
!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      do n=2,nlen(m),2
         do kv=1,plev*8
            tmpGRcoef(kv,lm) = tmpGRcoef(kv,lm) + tmpSPEcoef(kv,n,lm)*dalpn(lmr+n)
         end do
      end do
   end do
!
!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      do n=1,nlen(m),2
         do kv=plev*8+1,plev*24
            tmpGRcoef(kv,lm) = tmpGRcoef(kv,lm) + tmpSPEcoef(kv,n,lm)*lalp(lmr+n,irow)
         end do
      end do
   end do
!
! Combine the two parts of u(m) and v(m)
!
   do lm=1,mlength
      do kv=1,plev*8
         tmpGRcoef(kv,lm) = tmpGRcoef(kv,lm) + tmpGRcoef(kv+plev*16,lm)
      end do
   end do
!
! Save accumulated results to gr* arrays
!
   do lm=1,mlength
!DIR$ PREFERVECTOR
      do k=1,plev
         grus (2*lm-1,k) = tmpGRcoef(k        ,lm)
         grus (2*lm  ,k) = tmpGRcoef(k+plev   ,lm)
         grvs (2*lm-1,k) = tmpGRcoef(k+plev*2 ,lm)
         grvs (2*lm  ,k) = tmpGRcoef(k+plev*3 ,lm)
         gruhs(2*lm-1,k) = tmpGRcoef(k+plev*4 ,lm)
         gruhs(2*lm  ,k) = tmpGRcoef(k+plev*5 ,lm)
         grvhs(2*lm-1,k) = tmpGRcoef(k+plev*6 ,lm)
         grvhs(2*lm  ,k) = tmpGRcoef(k+plev*7 ,lm)

         grts (2*lm-1,k) = tmpGRcoef(k+plev*8 ,lm)
         grts (2*lm  ,k) = tmpGRcoef(k+plev*9 ,lm)
         grths(2*lm-1,k) = tmpGRcoef(k+plev*10,lm)
         grths(2*lm  ,k) = tmpGRcoef(k+plev*11,lm)
         grds (2*lm-1,k) = tmpGRcoef(k+plev*12,lm)
         grds (2*lm  ,k) = tmpGRcoef(k+plev*13,lm)
         grzs (2*lm-1,k) = tmpGRcoef(k+plev*14,lm)
         grzs (2*lm  ,k) = tmpGRcoef(k+plev*15,lm)
      end do
   end do
!
! Remove Coriolis contribution to absolute vorticity from u(m)
! Correction for u:zeta=vz-ez=(zeta+f)-f
!
   if (lmwave0 .ne. -1) then
      do k=1,plev
!        grus(1,k) = grus(1,k) - zurcor
         grus(2*lmwave0-1,k) = grus(2*lmwave0-1,k) - zurcor
      end do
   endif
!
!-----------------------------------------------------------------------
!
! Computation for 1-level variables (ln(p*) and derivatives).
!
!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
!cdir shortloop
      do n=1,nlen(m),2
         ir = lmc + 2*n - 1
         ii = ir + 1
!         
         grpss (2*lm-1) = grpss (2*lm-1) + alps(ir)*lalp(lmr+n,irow)
         grpss (2*lm  ) = grpss (2*lm  ) + alps(ii)*lalp(lmr+n,irow)
!
         grdpss(2*lm-1) = grdpss(2*lm-1) + alps(ir)*lalp(lmr+n,irow)*hdfst4(m+n-1)*ztodt
         grdpss(2*lm  ) = grdpss(2*lm  ) + alps(ii)*lalp(lmr+n,irow)*hdfst4(m+n-1)*ztodt
      end do
   end do

!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
!cdir shortloop
      do n=2,nlen(m),2
         ir = lmc + 2*n - 1
         ii = ir + 1
!
         grpms(2*lm-1) = grpms(2*lm-1) + alps(ir)*ldalp(lmr+n,irow)*ra
         grpms(2*lm  ) = grpms(2*lm  ) + alps(ii)*ldalp(lmr+n,irow)*ra
      end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
      grpls(2*lm-1) = -grpss(2*lm  )*ra*xm(m)
      grpls(2*lm  ) =  grpss(2*lm-1)*ra*xm(m)
   end do
!
   return
end subroutine grcalcs

subroutine grcalca (irow    ,ztodt   ,grta    ,grtha   ,grda    ,&
                    grza    ,grua    ,gruha   ,grva    ,grvha   ,&
                    grpsa   ,grdpsa  ,grpma   ,grpla   ,tmpSPEcoef)

!-----------------------------------------------------------------------
!
! Complete inverse Legendre transforms from spectral to Fourier space at 
! the the given latitude. Only positive latitudes are considered and 
! symmetric and antisymmetric (about equator) components are computed. 
! The sum and difference of these components give the actual fourier 
! coefficients for the latitude circle in the northern and southern 
! hemispheres respectively.
!
! The naming convention is as follows:
!  - The fourier coefficient arrays all begin with "gr";
!  - "t, q, d, z, ps" refer to temperature, specific humidity, 
!     divergence, vorticity, and surface pressure;
!  - "h" refers to the horizontal diffusive tendency for the field.
!  - "s" suffix to an array => symmetric component;
!  - "a" suffix to an array => antisymmetric component.
! Thus "grts" contains the symmetric Fourier coeffs of temperature and
! "grtha" contains the antisymmetric Fourier coeffs of the temperature
! tendency due to horizontal diffusion.
! Three additional surface pressure related quantities are returned:
!  1. "grdpss" and "grdpsa" contain the surface pressure factor
!      (proportional to del^4 ps) used for the partial correction of 
!      the horizontal diffusion to pressure surfaces.
!  2. "grpms" and "grpma" contain the longitudinal component of the 
!      surface pressure gradient.
!  3. "grpls" and "grpla" contain the latitudinal component of the 
!      surface pressure gradient.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Modified:          P. Worley, October 2002
!
!-----------------------------------------------------------------------
!
! $Id$
! $Author$
!
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use physconst, only: ra
   use eul_control_mod
   use spmd_utils, only : iam
   implicit none

!
! Input arguments
!
   integer,  intent(in) :: irow        ! latitude pair index
   real(r8), intent(in) :: ztodt       ! twice the timestep unless nstep = 0
   real(r8), intent(in) :: tmpSPEcoef(plev*24,pnmax,maxm) ! array for rearranged variables
!
!
! Output arguments: antisymmetric fourier coefficients
!
   real(r8), intent(out) :: grta(2*maxm,plev)    ! sum(n) of t(n,m)*P(n,m)
   real(r8), intent(out) :: grtha(2*maxm,plev)   ! sum(n) of K(2i)*t(n,m)*P(n,m)
   real(r8), intent(out) :: grda(2*maxm,plev)    ! sum(n) of d(n,m)*P(n,m)
   real(r8), intent(out) :: grza(2*maxm,plev)    ! sum(n) of z(n,m)*P(n,m)
   real(r8), intent(out) :: grua(2*maxm,plev)    ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: gruha(2*maxm,plev)   ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grva(2*maxm,plev)    ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grvha(2*maxm,plev)   ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
   real(r8), intent(out) :: grpsa(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grdpsa(2*maxm)       ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
   real(r8), intent(out) :: grpma(2*maxm)        ! sum(n) of lnps(n,m)*H(n,m)
   real(r8), intent(out) :: grpla(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)*m/a
!
!---------------------------Local workspace-----------------------------
!
   real(r8) dalpn(pspt)          ! (a/(n(n+1)))*derivative of Legendre functions (complex)
   real(r8) tmpGRcoef(plev*24,maxm) ! temporal storage for Fourier coefficients

   integer k                ! level index
   integer lm, m            ! local and global Fourier wavenumber indices of spectral array
   integer mlength          ! number of local wavenumbers
   integer n                ! meridional wavenumber index
   integer ir,ii            ! spectral indices
   integer lmr,lmc          ! spectral indices
   integer kv               ! level x variable index
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
!
!DIR$ NOSTREAM
   mlength = numm(iam)
!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      do n=1,nlen(m)
         dalpn(lmr+n) = ldalp(lmr+n,irow)*rsq(m+n-1)*ra
      end do
   end do
!
! Initialize sums
!
   grpsa (:) = 0._r8
   grpla (:) = 0._r8
   grpma (:) = 0._r8
   grdpsa(:) = 0._r8
!cdir collapse
   tmpGRcoef(:,:) = 0._r8
!
! Loop over n for t,q,d,and end of u and v
!
!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      do n=1,nlen(m),2
         do kv=1,plev*8
            tmpGRcoef(kv,lm) = tmpGRcoef(kv,lm) + tmpSPEcoef(kv,n,lm)*dalpn(lmr+n)
         end do
      end do
   end do

!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      do n=2,nlen(m),2
         do kv=plev*8+1,plev*24
            tmpGRcoef(kv,lm) = tmpGRcoef(kv,lm) + tmpSPEcoef(kv,n,lm)*lalp(lmr+n,irow)
         end do
      end do
   end do
!
! Combine the two parts of u(m) and v(m)
!
   do lm=1,mlength
      do kv=1,plev*8
         tmpGRcoef(kv,lm) = tmpGRcoef(kv,lm) + tmpGRcoef(kv+plev*16,lm)
      end do
   end do
!
! Save accumulated results to gr* arrays
!
   do lm=1,mlength
!DIR$ PREFERVECTOR
      do k=1,plev
         grua (2*lm-1,k) = tmpGRcoef(k        ,lm)
         grua (2*lm  ,k) = tmpGRcoef(k+plev   ,lm)
         grva (2*lm-1,k) = tmpGRcoef(k+plev*2 ,lm)
         grva (2*lm  ,k) = tmpGRcoef(k+plev*3 ,lm)
         gruha(2*lm-1,k) = tmpGRcoef(k+plev*4 ,lm)
         gruha(2*lm  ,k) = tmpGRcoef(k+plev*5 ,lm)
         grvha(2*lm-1,k) = tmpGRcoef(k+plev*6 ,lm)
         grvha(2*lm  ,k) = tmpGRcoef(k+plev*7 ,lm)

         grta (2*lm-1,k) = tmpGRcoef(k+plev*8 ,lm)
         grta (2*lm  ,k) = tmpGRcoef(k+plev*9 ,lm)
         grtha(2*lm-1,k) = tmpGRcoef(k+plev*10,lm)
         grtha(2*lm  ,k) = tmpGRcoef(k+plev*11,lm)
         grda (2*lm-1,k) = tmpGRcoef(k+plev*12,lm)
         grda (2*lm  ,k) = tmpGRcoef(k+plev*13,lm)
         grza (2*lm-1,k) = tmpGRcoef(k+plev*14,lm)
         grza (2*lm  ,k) = tmpGRcoef(k+plev*15,lm)
      end do
   end do
!
!-----------------------------------------------------------------------
!
! Computation for 1-level variables (ln(p*) and derivatives).
!
!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
!cdir shortloop
      do n=1,nlen(m),2
         ir = lmc + 2*n - 1
         ii = ir + 1

         grpma(2*lm-1) = grpma(2*lm-1) + alps(ir)*ldalp(lmr+n,irow)*ra
         grpma(2*lm  ) = grpma(2*lm  ) + alps(ii)*ldalp(lmr+n,irow)*ra
      end do
   end do

!cdir novector
   do lm=1,mlength
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
!cdir shortloop
      do n=2,nlen(m),2
         ir = lmc + 2*n - 1
         ii = ir + 1
!
         grpsa (2*lm-1) = grpsa (2*lm-1) + alps(ir)*lalp(lmr+n,irow)
         grpsa (2*lm  ) = grpsa (2*lm  ) + alps(ii)*lalp(lmr+n,irow)
!
         grdpsa(2*lm-1) = grdpsa(2*lm-1) + alps(ir)*lalp(lmr+n,irow)*hdfst4(m+n-1)*ztodt
         grdpsa(2*lm  ) = grdpsa(2*lm  ) + alps(ii)*lalp(lmr+n,irow)*hdfst4(m+n-1)*ztodt
      end do
!
! Multiply by m/a to get d(ln(p*))/dlamda
! and by 1/a to get (1-mu**2)d(ln(p*))/dmu
!
      grpla(2*lm-1) = -grpsa(2*lm  )*ra*xm(m)
      grpla(2*lm  ) =  grpsa(2*lm-1)*ra*xm(m)
   end do
!
   return
end subroutine grcalca

subroutine prepGRcalc(tmpSPEcoef)

!-----------------------------------------------------------------------
!
! Rearrange multi-level spectral coefficients for vectorization.
! The results are saved to "tmpSPEcoef" and will be used in
! "grcalcs" and "grcalca". 
!
!-----------------------------------------------------------------------
!
  use shr_kind_mod, only: r8 => shr_kind_r8
  use pmgrid
  use pspect
  use comspe
  use rgrid
  use commap
  use physconst, only: ra
  use eul_control_mod,    only: hdiftq, hdifzd
  use spmd_utils, only : iam
!
  implicit none
!
!
!---------------------------Output argument-----------------------------
!
   real(r8), intent(out) :: tmpSPEcoef(plev*24,pnmax,maxm) ! array for rearranged variables
!
!---------------------------Local workspace-----------------------------
!
   real(r8) raxm
!
   integer lm, m, n, k
   integer lmr, lmc
   integer ir ,ii
!
!-----------------------------------------------------------------------
!
   do lm=1,numm(iam)
      m = locm(lm,iam)
      lmr = lnstart(lm)
      lmc = 2*lmr
      raxm = ra*xm(m)
      do n=1,nlen(m)
         ir = lmc + 2*n - 1
         ii = ir + 1
         do k=1,plev
            tmpSPEcoef(k        ,n,lm) = vz(ir,k)
            tmpSPEcoef(k+plev   ,n,lm) = vz(ii,k)
            tmpSPEcoef(k+plev*2 ,n,lm) = -d(ir,k)
            tmpSPEcoef(k+plev*3 ,n,lm) = -d(ii,k)
            tmpSPEcoef(k+plev*4 ,n,lm) = -vz(ir,k)*hdifzd(n+m-1,k)
            tmpSPEcoef(k+plev*5 ,n,lm) = -vz(ii,k)*hdifzd(n+m-1,k)
            tmpSPEcoef(k+plev*6 ,n,lm) = d(ir,k)*hdifzd(n+m-1,k)
            tmpSPEcoef(k+plev*7 ,n,lm) = d(ii,k)*hdifzd(n+m-1,k)

            tmpSPEcoef(k+plev*8 ,n,lm) = t(ir,k)
            tmpSPEcoef(k+plev*9 ,n,lm) = t(ii,k)
            tmpSPEcoef(k+plev*10,n,lm) = -t(ir,k)*hdiftq(n+m-1,k)
            tmpSPEcoef(k+plev*11,n,lm) = -t(ii,k)*hdiftq(n+m-1,k)
            tmpSPEcoef(k+plev*12,n,lm) = d(ir,k)
            tmpSPEcoef(k+plev*13,n,lm) = d(ii,k)
            tmpSPEcoef(k+plev*14,n,lm) = vz(ir,k)
            tmpSPEcoef(k+plev*15,n,lm) = vz(ii,k)

            tmpSPEcoef(k+plev*16,n,lm) =  d (ii,k)*rsq(m+n-1)*raxm
            tmpSPEcoef(k+plev*17,n,lm) = -d (ir,k)*rsq(m+n-1)*raxm
            tmpSPEcoef(k+plev*18,n,lm) =  vz(ii,k)*rsq(m+n-1)*raxm
            tmpSPEcoef(k+plev*19,n,lm) = -vz(ir,k)*rsq(m+n-1)*raxm
            tmpSPEcoef(k+plev*20,n,lm) = -d (ii,k)*hdifzd(n+m-1,k)*rsq(m+n-1)*raxm
            tmpSPEcoef(k+plev*21,n,lm) =  d (ir,k)*hdifzd(n+m-1,k)*rsq(m+n-1)*raxm
            tmpSPEcoef(k+plev*22,n,lm) = -vz(ii,k)*hdifzd(n+m-1,k)*rsq(m+n-1)*raxm
            tmpSPEcoef(k+plev*23,n,lm) =  vz(ir,k)*hdifzd(n+m-1,k)*rsq(m+n-1)*raxm
         end do
      end do
   end do
!
   return
end subroutine prepGRcalc
