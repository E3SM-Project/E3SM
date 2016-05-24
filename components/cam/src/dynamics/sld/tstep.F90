  subroutine tstep(lm      ,ztdtsq  )
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the vertically coupled system of equations arising
! from the semi-implicit equations for each spectral element along
! two dimensional wavenumber n.  The inverse matrix depends
! only on two dimensional wavenumber and the reference atmosphere.
! It is precomputed and stored for use during the forecast. The routine
! overwrites the d,T and lnps coefficients with the new values.
!
! Original version:  CCM1
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
    use scanslt,      only: epssld
    use commap
    use spmd_utils, only : iam
    implicit none

!------------------------------Arguments--------------------------------
!
! Input arguments
!
    integer , intent(in)   :: lm            ! local Fourier wavenumber index
    real(r8), intent(in)   :: ztdtsq(pnmax) ! 2*dt*(n(n+1)/a^2 where n is 2-d wavenumber
!
!---------------------------Local workspace-----------------------------
!
    real(r8) hhref           ! href/2 (reference hydrostatic matrix / 2)
    real(r8) hbps            ! bps/2 (ref. coeff. for lnps term in div. eq. / 2)
    real(r8) onepeps         ! decentering coefficient
    integer m                ! Fourier wavenumber index
    integer n                ! 2-d wavenumber index
    integer k,kk             ! level indices
    integer mr,mc            ! real and imaginary spectral indices
    integer ir,ii            ! real and imaginary spectral indices
!
!-----------------------------------------------------------------------
!
! Complete rhs of helmholtz eq.
!
    m  = locm(lm,iam)
    mr = nstart(m)
    mc = 2*mr
    onepeps = 1._r8 + epssld
!
    do k=1,plev
!
! Coefficients for diagonal terms
!
       hhref = onepeps*0.5_r8*href(k,k)
       hbps = onepeps*0.5_r8*bps(k)
!
! Loop along total wavenumber index (in spectral space)
! Add lnps and diagonal (vertical space) T terms to initialize hs
!
       do n=1,nlen(m)
          ir = mc + 2*n - 1
          ii = ir + 1
          hs(ir,k) = ztdtsq(n+m-1)*(hhref*t(ir,k) + hbps*alps(ir))
          hs(ii,k) = ztdtsq(n+m-1)*(hhref*t(ii,k) + hbps*alps(ii))
       end do
       if (k.lt.plev) then
          do kk=k+1,plev
!
! Add off-diagonal (vertical space) T terms to hs
!
             hhref = onepeps*0.5_r8*href(kk,k)
             do n=1,nlen(m)
                ir = mc + 2*n - 1
                ii = ir + 1
                hs(ir,k) = hs(ir,k) + ztdtsq(n+m-1)*hhref*t(ir,kk)
                hs(ii,k) = hs(ii,k) + ztdtsq(n+m-1)*hhref*t(ii,kk)
             end do
          end do
       end if
    end do                    ! k=1,plev (calculation level)
!
! Transform semi-implicit vectors to vertical normal mode space
!
    do k = 1,plev
       do n=1,nlen(m)
          ir = mc + 2*n - 1
          ii = ir + 1
          dsnm(ir,k) = 0._r8
          dsnm(ii,k) = 0._r8
          hsnm(ir,k) = 0._r8
          hsnm(ii,k) = 0._r8
          vznm(ir,k) = 0._r8
          vznm(ii,k) = 0._r8
          do kk = 1,plev
             dsnm(ir,k) = dsnm(ir,k) + bmi(kk,k)*d (ir,kk)
             dsnm(ii,k) = dsnm(ii,k) + bmi(kk,k)*d (ii,kk)
             hsnm(ir,k) = hsnm(ir,k) + bmi(kk,k)*hs(ir,kk)
             hsnm(ii,k) = hsnm(ii,k) + bmi(kk,k)*hs(ii,kk)
             vznm(ir,k) = vznm(ir,k) + bmi(kk,k)*vz(ir,kk)
             vznm(ii,k) = vznm(ii,k) + bmi(kk,k)*vz(ii,kk)
          end do
       end do
    end do
!
    return
  end subroutine tstep

