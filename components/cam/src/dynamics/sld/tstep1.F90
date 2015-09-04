  subroutine tstep1(lm      ,zdt     )
!-----------------------------------------------------------------------
!
! Purpose:
! Solution of the vertically coupled system of equations arising
! from the semi-impicit equations for each spectral element along
! two dimensional wavenumber n.  The inverse matrix depends
! only on two dimensional wavenumber and the reference atmosphere.
! It is precomputed and stored for use during the forecast. The routine
! overwrites the d,T and lnps coefficients with the new values.
!
!---------------------------Code history--------------------------------
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
    use hycoef, only : hypi, hypd
    implicit none

!------------------------------Arguments--------------------------------
!
    integer , intent(in)   :: lm  ! local Fourier wavenumber index
    real(r8), intent(in)   :: zdt ! timestep, dt (seconds)
!
!---------------------------Local workspace-----------------------------
!
    real(r8) z (2*pnmax,plev) ! workspace for computation of spectral array d
    real(r8) zz(2*pnmax,plev) ! workspace for computation of spectral array vz
    real(r8) ztemp            ! temporary workspace
    real(r8) onepeps          ! decentering coefficient
    integer m                 ! Fourier wavenumber
    integer n,j               ! 2-d wavenumber index
    integer k,kk              ! level indices
    integer mr,mc             ! real and imaginary spectral indices
    integer ir,ii             ! real and imaginary spectral indices
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
! Solution of helmholtz equation
! First: initialize temporary space for solution
!
    do k=1,plev
       do j=1,2*pnmax
          z (j,k) = 0._r8
          zz(j,k) = 0._r8
       end do
    end do
!
! Transform back from normal mode space
!
    do k=1,plev
       do kk=1,plev
          do n=1,nlen(m)
             ir = mc + 2*n - 1
             ii = ir + 1
             z (2*n-1,k) = z (2*n-1,k) + bm1(kk,k)*dnm (ir,kk)
             z (2*n  ,k) = z (2*n  ,k) + bm1(kk,k)*dnm (ii,kk)
             zz(2*n-1,k) = zz(2*n-1,k) + bm1(kk,k)*vznm(ir,kk)
             zz(2*n  ,k) = zz(2*n  ,k) + bm1(kk,k)*vznm(ii,kk)
          end do
       end do                  ! inner loop over levels
    end do                    ! outer loop over levels
!
! Move solution for divergence and vorticity to d and vz.
!
    do k=1,plev
       do n=1,nlen(m)
          ir = mc + 2*n - 1
          ii = ir + 1
          d (ir,k) = z (2*n-1,k)
          d (ii,k) = z (2*n  ,k)
          vz(ir,k) = zz(2*n-1,k)
          vz(ii,k) = zz(2*n  ,k)
       end do
    end do
!
! Complete ln(pstar) and T forecasts
! Add semi-implicit part to surface pressure (vector multiply)
!
    do k=1,plev
       ztemp = onepeps*zdt*hypd(k)/hypi(plevp)
       do n=1,nlen(m)
          ir = mc + 2*n - 1
          ii = ir + 1
          alps(ir) = alps(ir) - ztemp*d(ir,k)
          alps(ii) = alps(ii) - ztemp*d(ii,k)
       end do
    end do
!
! Add ln(Ps)star back in to get full ln(Ps)
!
    do n=1,nlen(m)
       ir = mc + 2*n - 1
       ii = ir + 1
       alps(ir) = alps(ir) + lnpstar(ir)
       alps(ii) = alps(ii) + lnpstar(ii)
    end do
!
! Add semi-implicit part to temperature (matrix multiply)
!
    do k=1,plev
       do kk=1,plev
          ztemp = onepeps*zdt*tau(kk,k)
          do n=1,nlen(m)
             ir = mc + 2*n - 1
             ii = ir + 1
             t(ir,k) = t(ir,k) - ztemp*d(ir,kk)
             t(ii,k) = t(ii,k) - ztemp*d(ii,kk)
          end do
       end do
    end do
!
    return
  end subroutine tstep1
