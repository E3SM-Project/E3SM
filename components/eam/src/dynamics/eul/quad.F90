
subroutine quad(lm      ,zdt     ,ztdtsq  ,grlps1  ,grlps2  ,&
                grt1    ,grz1    ,grd1    ,grfu1   ,grfv1   ,&
                grvt1   ,grrh1   ,grt2    ,grz2    ,grd2    ,&
                grfu2   ,grfv2   ,grvt2   ,grrh2   )
!-----------------------------------------------------------------------
!
! Perform gaussian quadrature for 1 Fourier wavenumber (m) to obtain the 
! spectral coefficients of ln(ps), temperature, vorticity, and divergence.
! Add the tendency terms requiring meridional derivatives during the
! transform.
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Rosinski
! Standardized:      J. Rosinski, June 1992
! Reviewed:          B. Boville, D. Williamson, J. Hack, August 1992
! Reviewed:          B. Boville, D. Williamson, April 1996
! Modified:          P. Worley, September 2002
! Modified:          NEC, April 2004
!
!-----------------------------------------------------------------------
   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid
   use pspect
   use comspe
   use rgrid
   use commap
   use physconst, only: rearth
   use spmd_utils, only : iam
   implicit none
!
! Input arguments
!
   integer, intent(in) :: lm                         ! local Fourier wavenumber index
   real(r8), intent(in) :: zdt                       ! timestep(dt) unless nstep = 0
   real(r8), intent(in) :: ztdtsq(pnmax)             ! 2*zdt*n(n+1)/(a^2)
!                                            where n IS the 2-d wavenumber
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMS and and in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
! Suffixes 1 and 2 refer to symmetric and antisymmetric components
! respectively.
!
   real(r8), intent(in) :: grlps1(2*maxm,(plat+1)/2)        ! ln(ps) - symmetric
   real(r8), intent(in) :: grlps2(2*maxm,(plat+1)/2)        ! ln(ps) - antisymmetric
!
! symmetric components
!
   real(r8), intent(in) :: grt1(2*maxm,plev,(plat+1)/2)     ! temperature
   real(r8), intent(in) :: grz1(2*maxm,plev,(plat+1)/2)     ! vorticity
   real(r8), intent(in) :: grd1(2*maxm,plev,(plat+1)/2)     ! divergence
   real(r8), intent(in) :: grfu1(2*maxm,plev,(plat+1)/2)    ! partial u momentum tendency (fu)
   real(r8), intent(in) :: grfv1(2*maxm,plev,(plat+1)/2)    ! partial v momentum tendency (fv)
   real(r8), intent(in) :: grvt1(2*maxm,plev,(plat+1)/2)    ! heat flux
   real(r8), intent(in) :: grrh1(2*maxm,plev,(plat+1)/2)    ! rhs of div eqn (del^2 term)
!
! antisymmetric components
!
   real(r8), intent(in) :: grt2(2*maxm,plev,(plat+1)/2)     ! temperature
   real(r8), intent(in) :: grz2(2*maxm,plev,(plat+1)/2)     ! vorticity
   real(r8), intent(in) :: grd2(2*maxm,plev,(plat+1)/2)     ! divergence
   real(r8), intent(in) :: grfu2(2*maxm,plev,(plat+1)/2)    ! partial u momentum tend (fu)
   real(r8), intent(in) :: grfv2(2*maxm,plev,(plat+1)/2)    ! partial v momentum tend (fv)
   real(r8), intent(in) :: grvt2(2*maxm,plev,(plat+1)/2)    ! heat flux
   real(r8), intent(in) :: grrh2(2*maxm,plev,(plat+1)/2)    ! rhs of div eqn (del^2 term)
!
!---------------------------Local workspace-----------------------------
!
   integer j                          ! latitude pair index
   integer m                          ! global wavenumber index
   integer n                          ! total wavenumber index
   integer ir,ii                      ! spectral indices
   integer lmr,lmc                    ! spectral indices
   integer k                          ! level index
   integer kv                         ! index for vectorization

   real(r8) zcsj                          ! cos**2(lat)*radius of earth
   real(r8) zrcsj                         ! 1./(a*cos^2(lat))
   real(r8) zdtrc                         ! dt/(a*cos^2(lat))
   real(r8) ztdtrc                        ! 2dt/(a*cos^2(lat))
   real(r8) zw((plat+1)/2)                    ! 2*w
   real(r8) ztdtrw((plat+1)/2)                ! 2w*2dt/(a*cos^2(lat))
   real(r8) zwalp                         ! zw*alp
   real(r8) zwdalp                        ! zw*dalp
   real(r8) sqzwalp                       ! ztdtsq*zw*alp

   real(r8) tmpGR1odd(plev*6,(plat+1)/2)      ! temporary space for Fourier coeffs
   real(r8) tmpGR2odd(plev*6,(plat+1)/2)      !
   real(r8) tmpGR3odd(plev*6,(plat+1)/2)      !
   real(r8) tmpGR1evn(plev*6,(plat+1)/2)      !
   real(r8) tmpGR2evn(plev*6,(plat+1)/2)      !
   real(r8) tmpGR3evn(plev*6,(plat+1)/2)      !

   real(r8) tmpSPEodd(plev*6,2*ptrn)      ! temporary space for spectral coeffs
   real(r8) tmpSPEevn(plev*6,2*ptrn)      !
!
!-----------------------------------------------------------------------
!
! Compute constants
!
!$OMP PARALLEL DO PRIVATE(J, ZCSJ, ZRCSJ, ZDTRC, ZTDTRC)
   do j=1,plat/2
      zcsj = cs(j)*rearth
      zrcsj = 1._r8/zcsj
      zdtrc = zdt*zrcsj
      ztdtrc = 2._r8*zdtrc
      zw(j) = w(j)*2._r8
      ztdtrw(j) = ztdtrc*zw(j)
   end do
!
! Accumulate contributions to spectral coefficients of ln(p*), the only
! single level field. Use symmetric or antisymmetric fourier cofficients
! depending on whether the total wavenumber is even or odd.
!
   m  = locm(lm,iam)
   lmr = lnstart(lm)
   lmc = 2*lmr
   do n=1,2*nlen(m)
      alps(lmc+n) = 0._r8
   end do
!$OMP PARALLEL DO PRIVATE(N, J, IR, II, ZWALP)
   do n=1,nlen(m),2
      ir = lmc + 2*n - 1
      ii = ir + 1
      do j=beglatpair(m),plat/2
         zwalp = zw(j)*lalp(lmr+n,j)
         alps(ir) = alps(ir) + grlps1(2*lm-1,j)*zwalp
         alps(ii) = alps(ii) + grlps1(2*lm  ,j)*zwalp
      end do
   end do
!$OMP PARALLEL DO PRIVATE(N, J, IR, II, ZWALP)
   do n=2,nlen(m),2
      ir = lmc + 2*n - 1
      ii = ir + 1
      do j=beglatpair(m),plat/2
         zwalp = zw(j)*lalp(lmr+n,j)
         alps(ir) = alps(ir) + grlps2(2*lm-1,j)*zwalp
         alps(ii) = alps(ii) + grlps2(2*lm  ,j)*zwalp
      end do
   end do
!
! Accumulate contributions to spectral coefficients of the multilevel fields.
! Use symmetric or antisymmetric fourier coefficients depending on whether
! the total wavenumber is even or odd.
!
!
! Initialize temporary storage for spectral coefficients
!
   do n=1,nlen(m)
      do kv=1,plev*6
         tmpSPEodd(kv,n) = 0._r8
         tmpSPEevn(kv,n) = 0._r8
      end do
   end do
!
! Rearrange Fourier coefficients to temporal storage
!
!$OMP PARALLEL DO PRIVATE(J, K)
   do j = beglatpair(m),plat/2
      do k=1,plev

         tmpGR1odd(k       ,j) =  grt1 (2*lm-1,k,j) ! first term for odd n
         tmpGR1odd(k+plev  ,j) =  grt1 (2*lm  ,k,j)
         tmpGR1odd(k+plev*2,j) =  grd1 (2*lm-1,k,j)
         tmpGR1odd(k+plev*3,j) =  grd1 (2*lm  ,k,j)
         tmpGR1odd(k+plev*4,j) =  grz1 (2*lm-1,k,j)
         tmpGR1odd(k+plev*5,j) =  grz1 (2*lm  ,k,j)

         tmpGR2odd(k       ,j) =  grvt2(2*lm-1,k,j) ! second term for odd n
         tmpGR2odd(k+plev  ,j) =  grvt2(2*lm  ,k,j)
         tmpGR2odd(k+plev*2,j) = -grfv2(2*lm-1,k,j)
         tmpGR2odd(k+plev*3,j) = -grfv2(2*lm  ,k,j)
         tmpGR2odd(k+plev*4,j) =  grfu2(2*lm-1,k,j)
         tmpGR2odd(k+plev*5,j) =  grfu2(2*lm  ,k,j)

         tmpGR3odd(k+plev*2,j) =  grrh1(2*lm-1,k,j) ! additional term for odd n
         tmpGR3odd(k+plev*3,j) =  grrh1(2*lm  ,k,j)

         tmpGR1evn(k       ,j) =  grt2 (2*lm-1,k,j) ! first term for even n
         tmpGR1evn(k+plev  ,j) =  grt2 (2*lm  ,k,j)
         tmpGR1evn(k+plev*2,j) =  grd2 (2*lm-1,k,j)
         tmpGR1evn(k+plev*3,j) =  grd2 (2*lm  ,k,j)
         tmpGR1evn(k+plev*4,j) =  grz2 (2*lm-1,k,j)
         tmpGR1evn(k+plev*5,j) =  grz2 (2*lm  ,k,j)

         tmpGR2evn(k       ,j) =  grvt1(2*lm-1,k,j) ! first term for even n
         tmpGR2evn(k+plev  ,j) =  grvt1(2*lm  ,k,j)
         tmpGR2evn(k+plev*2,j) = -grfv1(2*lm-1,k,j)
         tmpGR2evn(k+plev*3,j) = -grfv1(2*lm  ,k,j)
         tmpGR2evn(k+plev*4,j) =  grfu1(2*lm-1,k,j)
         tmpGR2evn(k+plev*5,j) =  grfu1(2*lm  ,k,j)

         tmpGR3evn(k+plev*2,j) =  grrh2(2*lm-1,k,j) ! additional term for even n
         tmpGR3evn(k+plev*3,j) =  grrh2(2*lm  ,k,j)

      end do
   end do
!
! Accumulate first and second terms
!
!$OMP PARALLEL DO PRIVATE(N, J, ZWDALP, ZWALP, KV)
   do n=1,nlen(m),2
      do j=beglatpair(m),plat/2
         zwdalp = ztdtrw(j)*ldalp(lmr+n,j)
         zwalp  = zw(j)    *lalp (lmr+n,j)
         do kv=1,plev*6
            tmpSPEodd(kv,n) = tmpSPEodd(kv,n) + &
                 zwalp*tmpGR1odd(kv,j) + zwdalp*tmpGR2odd(kv,j)
         end do
      end do
   end do
!$OMP PARALLEL DO PRIVATE(N, J, ZWDALP, ZWALP, KV)
   do n=2,nlen(m),2
      do j=beglatpair(m),plat/2
         zwdalp = ztdtrw(j)*ldalp(lmr+n,j)
         zwalp  = zw(j)    *lalp (lmr+n,j)
         do kv=1,plev*6
            tmpSPEevn(kv,n) = tmpSPEevn(kv,n) + &
                 zwalp*tmpGR1evn(kv,j) + zwdalp*tmpGR2evn(kv,j)
         end do
      end do
   end do
!
! Add additional term for divergence
!
!$OMP PARALLEL DO PRIVATE(N, J, SQZWALP, KV)
   do n=1,nlen(m),2
      do j=beglatpair(m),plat/2
         sqzwalp = ztdtsq(n+m-1)*zw(j)*lalp (lmr+n,j)
         do kv=plev*2+1,plev*4
            tmpSPEodd(kv,n) = tmpSPEodd(kv,n) + sqzwalp*tmpGR3odd(kv,j)
         end do
      end do
   end do
!$OMP PARALLEL DO PRIVATE(N, J, SQZWALP, KV)
   do n=2,nlen(m),2
      do j=beglatpair(m),plat/2
         sqzwalp  = ztdtsq(n+m-1)*zw(j)*lalp (lmr+n,j)
         do kv=plev*2+1,plev*4
            tmpSPEevn(kv,n) = tmpSPEevn(kv,n) + sqzwalp*tmpGR3evn(kv,j)
         end do
      end do
   end do
!
! Save accumulated results 
!
!$OMP PARALLEL DO PRIVATE(N, IR, II, K)
   do n=1,nlen(m),2
      ir = lmc+2*n-1
      ii = ir+1
      do k=1,plev
         t (ir,k) = tmpSPEodd(k       ,n)
         t (ii,k) = tmpSPEodd(k+plev  ,n)
         d (ir,k) = tmpSPEodd(k+plev*2,n)
         d (ii,k) = tmpSPEodd(k+plev*3,n)
         vz(ir,k) = tmpSPEodd(k+plev*4,n)
         vz(ii,k) = tmpSPEodd(k+plev*5,n)
      end do
   end do
!$OMP PARALLEL DO PRIVATE(N, IR, II, K)
   do n=2,nlen(m),2
      ir = lmc+2*n-1
      ii = ir+1
      do k=1,plev
         t (ir,k) = tmpSPEevn(k       ,n)
         t (ii,k) = tmpSPEevn(k+plev  ,n)
         d (ir,k) = tmpSPEevn(k+plev*2,n)
         d (ii,k) = tmpSPEevn(k+plev*3,n)
         vz(ir,k) = tmpSPEevn(k+plev*4,n)
         vz(ii,k) = tmpSPEevn(k+plev*5,n)
      end do
   end do
!
   return
end subroutine quad
