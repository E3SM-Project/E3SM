subroutine quad(lm      ,grlps1  ,grlps2  ,grt1    ,grq1    , &
                grz1    ,grd1    ,grfu1   ,grfv1   ,grt2    , &
                grq2    ,grz2    ,grd2    ,grfu2   ,grfv2   )
!-----------------------------------------------------------------------
!
! Purpose:
! Perform gaussian quadrature for 1 Fourier wavenumber (m) to obtain the
! spectral coefficients of ln(ps), temperature, vorticity, divergence
! and specific humidity.  Add the tendency terms requiring meridional
! derivatives during the transform.
!
! Author:   J. Rosinski
! Modified: P. Worley, October 2002
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

!------------------------------Arguments--------------------------------
!
! Fourier coefficient arrays which have a latitude index on them for
! multitasking. These arrays are defined in LINEMS and and used in QUAD
! to compute spectral coefficients. They contain a latitude index so
! that the sums over latitude can be performed in a specified order.
!
! Suffixes 1 and 2 refer to symmetric and antisymmetric components
! respectively.
!
  integer , intent(in)   :: lm                          ! local Fourier wavenumber index
  real(r8), intent(in)   :: grlps1(2*maxm,plat/2)      ! ln(ps) - symmetric
  real(r8), intent(in)   :: grlps2(2*maxm,plat/2)      ! ln(ps) - antisymmetric
!
! symmetric components
!
  real(r8), intent(in)   :: grt1  (2*maxm,plev,plat/2) ! temperature
  real(r8), intent(in)   :: grq1  (2*maxm,plev,plat/2) ! moisture
  real(r8), intent(in)   :: grz1  (2*maxm,plev,plat/2) ! vorticity
  real(r8), intent(in)   :: grd1  (2*maxm,plev,plat/2) ! divergence
  real(r8), intent(in)   :: grfu1 (2*maxm,plev,plat/2) ! partial u momentum tendency (fu)
  real(r8), intent(in)   :: grfv1 (2*maxm,plev,plat/2) ! partial v momentum tendency (fv)
!
! antisymmetric components
!
  real(r8), intent(in)   :: grt2  (2*maxm,plev,plat/2) ! temperature
  real(r8), intent(in)   :: grq2  (2*maxm,plev,plat/2) ! moisture
  real(r8), intent(in)   :: grz2  (2*maxm,plev,plat/2) ! vorticity
  real(r8), intent(in)   :: grd2  (2*maxm,plev,plat/2) ! divergence
  real(r8), intent(in)   :: grfu2 (2*maxm,plev,plat/2) ! partial u momentum tend (fu)
  real(r8), intent(in)   :: grfv2 (2*maxm,plev,plat/2) ! partial v momentum tend (fv)
!
!---------------------------Local workspace-----------------------------
!
  integer j                 ! latitude pair index
  integer m                 ! Fourier wavenumber index
  integer n                 ! index
  integer ir,ii             ! indices
  integer mr,mc             ! indices
  integer k                 ! level index
  real(r8) zcsj             ! cos**2(lat)*radius of earth
  real(r8) zrcsj            ! 1./(a*cos^2(lat))
  real(r8) zw    (plat/2)   ! 2*w
  real(r8) ztdtrw(plat/2)   ! 2w*2dt/(a*cos^2(lat))
  real(r8) zwalp
  real(r8) zwdalp
!
!-----------------------------------------------------------------------
!
! Compute constants
!
  do j = 1,plat/2
     zcsj      = cs(j)*rearth
     zrcsj     = 1._r8/zcsj
     zw(j)     = w(j)*2._r8
     ztdtrw(j) = zrcsj*zw(j)
  end do
!
! Accumulate contributions to spectral coefficients of ln(p*), the only
! single level field. Use symmetric or antisymmetric fourier cofficients
! depending on whether the total wavenumber is even or odd.
!
  m  = locm(lm,iam)
  mr = nstart(m)
  mc = 2*mr
  do n = 1,2*nlen(m)
     alps(mc+n) = 0._r8
  end do
  do j=beglatpair(m),plat/2
     do n=1,nlen(m),2
        ir = mc + 2*n - 1
        ii = ir + 1
        zwalp    = zw(j)*alp(mr+n,j)
        alps(ir) = alps(ir) + grlps1(2*lm-1,j)*zwalp
        alps(ii) = alps(ii) + grlps1(2*lm  ,j)*zwalp
     end do
     do n = 2,nlen(m),2
        ir = mc + 2*n - 1
        ii = ir + 1
        zwalp    = zw(j)*alp(mr+n,j)
        alps(ir) = alps(ir) + grlps2(2*lm-1,j)*zwalp
        alps(ii) = alps(ii) + grlps2(2*lm  ,j)*zwalp
     end do
  end do
!
! Accumulate contributions to spectral coefficients of the multilevel
! fields.  Use symmetric or antisymmetric fourier coefficients depending
! on whether the total wavenumber is even or odd.
!
  do k = 1,plev
     do n = 1,2*nlen(m)
        t (mc+n,k) = 0._r8
        q (mc+n,k) = 0._r8
        d (mc+n,k) = 0._r8
        vz(mc+n,k) = 0._r8
     end do
     do j=beglatpair(m),plat/2
        do n=1,nlen(m),2
           zwdalp = ztdtrw(j)*dalp(mr+n,j)
           zwalp  = zw(j)    *alp (mr+n,j)
           ir = mc + 2*n - 1
           ii = ir + 1
           t (ir,k) = t (ir,k) + zwalp*grt1 (2*lm-1,k,j)
           t (ii,k) = t (ii,k) + zwalp*grt1 (2*lm  ,k,j)
           q (ir,k) = q (ir,k) + zwalp*grq1 (2*lm-1,k,j)
           q (ii,k) = q (ii,k) + zwalp*grq1 (2*lm  ,k,j)
           d (ir,k) = d (ir,k) + grd1 (2*lm-1,k,j)*zwalp - grfv2(2*lm-1,k,j)*zwdalp
           d (ii,k) = d (ii,k) + grd1 (2*lm  ,k,j)*zwalp - grfv2(2*lm  ,k,j)*zwdalp
           vz(ir,k) = vz(ir,k) + grz1 (2*lm-1,k,j)*zwalp + grfu2(2*lm-1,k,j)*zwdalp
           vz(ii,k) = vz(ii,k) + grz1 (2*lm  ,k,j)*zwalp + grfu2(2*lm  ,k,j)*zwdalp
        end do
     end do
     do j=beglatpair(m),plat/2
        do n = 2,nlen(m),2
           zwdalp = ztdtrw(j)*dalp(mr+n,j)
           zwalp  = zw(j)    *alp (mr+n,j)
           ir = mc + 2*n - 1
           ii = ir + 1
           t (ir,k) = t (ir,k) + zwalp*grt2(2*lm-1,k,j)
           t (ii,k) = t (ii,k) + zwalp*grt2(2*lm  ,k,j)
           q (ir,k) = q (ir,k) + zwalp*grq2(2*lm-1,k,j)
           q (ii,k) = q (ii,k) + zwalp*grq2(2*lm  ,k,j)
           d (ir,k) = d (ir,k) + grd2 (2*lm-1,k,j)*zwalp - grfv1(2*lm-1,k,j)*zwdalp
           d (ii,k) = d (ii,k) + grd2 (2*lm  ,k,j)*zwalp - grfv1(2*lm  ,k,j)*zwdalp
           vz(ir,k) = vz(ir,k) + grz2 (2*lm-1,k,j)*zwalp + grfu1(2*lm-1,k,j)*zwdalp
           vz(ii,k) = vz(ii,k) + grz2 (2*lm  ,k,j)*zwalp + grfu1(2*lm  ,k,j)*zwdalp
        end do
     end do
  end do
!
  return
end subroutine quad

