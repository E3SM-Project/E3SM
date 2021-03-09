
subroutine grcalcs (irow    ,ztodt   ,grts    ,grqs    ,grths   , &
                    grds    ,grus    ,gruhs   ,grvs    ,grvhs   , &
                    grpss   ,grdpss  ,grpms   ,grpls   ,grtms   , &
                    grtls   ,grqms   ,grqls   )
!-----------------------------------------------------------------------
!
! Purpose:
! Complete inverse legendre transforms from spectral to Fourier space at
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
! Thus "grts" contains the symmetric fourier coeffs of temperature and
! "grtha" contains the antisymmetric fourier coeffs of the temperature
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
    use rgrid
    use commap
    use physconst, only: ra
    use sld_control_mod
    use spmd_utils, only : iam
    implicit none

!
! Input arguments
!
    integer , intent(in)   :: irow              ! latitude pair index
    real(r8), intent(in)   :: ztodt             ! twice the timestep unless nstep = 0
!
! Output arguments: symmetric fourier coefficients
!
    real(r8), intent(out) :: grts(2*maxm,plev)  ! sum(n) of t(n,m)*P(n,m)
    real(r8), intent(out) :: grqs(2*maxm,plev)  ! sum(n) of q(n,m)*P(n,m)
    real(r8), intent(out) :: grths(2*maxm,plev) ! sum(n) of K(2i)*t(n,m)*P(n,m)
    real(r8), intent(out) :: grds(2*maxm,plev)  ! sum(n) of d(n,m)*P(n,m)
    real(r8), intent(out) :: grus(2*maxm,plev)  ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: gruhs(2*maxm,plev) ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grvs(2*maxm,plev)  ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grvhs(2*maxm,plev) ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grpss(2*maxm)      ! sum(n) of lnps(n,m)*P(n,m)
    real(r8), intent(out) :: grdpss(2*maxm)     ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt*lnps(n,m)*P(n,m)
    real(r8), intent(out) :: grpms(2*maxm)    ! sum(n) of lnps(n,m)*H(n,m)
    real(r8), intent(out) :: grpls(2*maxm)    ! sum(n) of lnps(n,m)*P(n,m)*m/a
    real(r8), intent(out) :: grtms (2*maxm,plev)
    real(r8), intent(out) :: grtls (2*maxm,plev)
    real(r8), intent(out) :: grqms (2*maxm,plev)
    real(r8), intent(out) :: grqls (2*maxm,plev)
!
!---------------------------Local workspace-----------------------------
!
    real(r8) gru1s (2*maxm)      ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) gruh1s(2*maxm)      ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) grv1s (2*maxm)      ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) grvh1s(2*maxm)      ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) alpn  (pspt)       ! (a*m/(n(n+1)))*Legendre functions (complex)
    real(r8) dalpn (pspt)       ! (a/(n(n+1)))*derivative of Legendre functions (complex)

    integer k                   ! level index
    integer lm, m               ! local and global Fourier wavenumber indices of spectral array
    integer mlength             ! number of local wavenumbers
    integer n                   ! meridional wavenumber index
    integer ir,ii               ! spectral indices
    integer mr,mc               ! spectral indices
    real(r8) tmp,raxm           ! temporary workspace
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
!
    mlength = numm(iam)
    do lm=1,mlength
       m = locm(lm,iam)
       mr = nstart(m)
       raxm = ra*xm(m)
       do n=1,nlen(m)
          alpn(mr+n) = alp(mr+n,irow)*rsq(m+n-1)*raxm
          dalpn(mr+n) = dalp(mr+n,irow)*rsq(m+n-1)*ra
       end do
    end do
!
! Initialize sums
!
    grts(:,:) = 0._r8
    grqs(:,:) = 0._r8
    grths(:,:) = 0._r8
    grds(:,:)  = 0._r8
    grus(:,:)  = 0._r8
    gruhs(:,:) = 0._r8
    grvs(:,:)  = 0._r8
    grvhs(:,:) = 0._r8
    grpss(:)   = 0._r8
    grdpss(:)   = 0._r8
    grpms(:)   = 0._r8
    grpls(:)   = 0._r8
    grtms(:,:)   = 0._r8
    grtls(:,:)   = 0._r8
    grqms(:,:)   = 0._r8
    grqls(:,:)   = 0._r8
!
!-----------------------------------------------------------------------
!
! Computation for multilevel variables
!
    do k=1,plev
!
! Initialize local sums
!
       gru1s(:) = 0._r8
       gruh1s(:) = 0._r8
       grv1s(:) = 0._r8
       grvh1s(:) = 0._r8
!
! Loop over n for t,q,d,and end of u and v
!
       do lm=1,mlength
          m = locm(lm,iam)
          mr = nstart(m)
          mc = 2*mr
          do n=1,nlen(m),2
             ir = mc + 2*n - 1
             ii = ir + 1

             grts (2*lm-1,k) = grts (2*lm-1,k) + t(ir,k)*alp(mr+n,irow)
             grts (2*lm  ,k) = grts (2*lm  ,k) + t(ii,k)*alp(mr+n,irow)
!
             grqs (2*lm-1,k) = grqs (2*lm-1,k) + q(ir,k)*alp(mr+n,irow)
             grqs (2*lm  ,k) = grqs (2*lm  ,k) + q(ii,k)*alp(mr+n,irow)
!
             tmp = alp(mr+n,irow)*hdiftq(n+m-1,k)
             grths(2*lm-1,k) = grths(2*lm-1,k) - t(ir,k)*tmp
             grths(2*lm  ,k) = grths(2*lm  ,k) - t(ii,k)*tmp
!
             grds(2*lm-1,k) = grds(2*lm-1,k) + d(ir,k)*alp(mr+n,irow)
             grds(2*lm  ,k) = grds(2*lm  ,k) + d(ii,k)*alp(mr+n,irow)
!
             gru1s (2*lm-1) = gru1s (2*lm-1) + d(ir,k)*alpn(mr+n)
             gru1s (2*lm  ) = gru1s (2*lm  ) + d(ii,k)*alpn(mr+n)
!
             tmp = alpn(mr+n)*hdifzd(n+m-1,k)
             gruh1s(2*lm-1) = gruh1s(2*lm-1) - d(ir,k)*tmp
             gruh1s(2*lm  ) = gruh1s(2*lm  ) - d(ii,k)*tmp
!
             grv1s (2*lm-1) = grv1s (2*lm-1) + vz(ir,k)*alpn(mr+n)
             grv1s (2*lm  ) = grv1s (2*lm  ) + vz(ii,k)*alpn(mr+n)
!
             grvh1s(2*lm-1) = grvh1s(2*lm-1) - vz(ir,k)*tmp
             grvh1s(2*lm  ) = grvh1s(2*lm  ) - vz(ii,k)*tmp
          end do
       end do
       do lm=1,mlength
          m = locm(lm,iam)
          mr = nstart(m)
          mc = 2*mr
          do n=2,nlen(m),2
             ir = mc + 2*n - 1
             ii = ir + 1
!
             grtms(2*lm-1,k) = grtms(2*lm-1,k) + t(ir,k)*dalp(mr+n,irow)*ra
             grtms(2*lm  ,k) = grtms(2*lm  ,k) + t(ii,k)*dalp(mr+n,irow)*ra
!
             grqms(2*lm-1,k) = grqms(2*lm-1,k) + q(ir,k)*dalp(mr+n,irow)*ra
             grqms(2*lm  ,k) = grqms(2*lm  ,k) + q(ii,k)*dalp(mr+n,irow)*ra
!
             grus (2*lm-1,k) = grus (2*lm-1,k) + vz(ir,k)*dalpn(mr+n)
             grus (2*lm  ,k) = grus (2*lm  ,k) + vz(ii,k)*dalpn(mr+n)
!
             tmp = dalpn(mr+n)*hdifzd(n+m-1,k)
             gruhs(2*lm-1,k) = gruhs(2*lm-1,k) - vz(ir,k)*tmp
             gruhs(2*lm  ,k) = gruhs(2*lm  ,k) - vz(ii,k)*tmp
!
             grvs (2*lm-1,k) = grvs (2*lm-1,k) - d(ir,k)*dalpn(mr+n)
             grvs (2*lm  ,k) = grvs (2*lm  ,k) - d(ii,k)*dalpn(mr+n)
!
             grvhs(2*lm-1,k) = grvhs(2*lm-1,k) + d(ir,k)*tmp
             grvhs(2*lm  ,k) = grvhs(2*lm  ,k) + d(ii,k)*tmp
          end do
       end do
!
! Combine the two parts of u(m) and v(m)
!
       do lm=1,mlength
          m = locm(lm,iam)
          grus (2*lm-1,k) = grus (2*lm-1,k) + gru1s (2*lm  )
          gruhs(2*lm-1,k) = gruhs(2*lm-1,k) + gruh1s(2*lm  )
          grus (2*lm  ,k) = grus (2*lm  ,k) - gru1s (2*lm-1)
          gruhs(2*lm  ,k) = gruhs(2*lm  ,k) - gruh1s(2*lm-1)
          grvs (2*lm-1,k) = grvs (2*lm-1,k) + grv1s (2*lm  )
          grvhs(2*lm-1,k) = grvhs(2*lm-1,k) + grvh1s(2*lm  )
          grvs (2*lm  ,k) = grvs (2*lm  ,k) - grv1s (2*lm-1)
          grvhs(2*lm  ,k) = grvhs(2*lm  ,k) - grvh1s(2*lm-1)
!
! Derivatives
!
          grtls(2*lm-1,k) = -grts(2*lm  ,k)*ra*xm(m)
          grtls(2*lm  ,k) =  grts(2*lm-1,k)*ra*xm(m)
          grqls(2*lm-1,k) = -grqs(2*lm  ,k)*ra*xm(m)
          grqls(2*lm  ,k) =  grqs(2*lm-1,k)*ra*xm(m)
       end do
    end do
!
!-----------------------------------------------------------------------
!
! Computation for 1-level variables (ln(p*) and derivatives).
!
    do lm=1,mlength
       m = locm(lm,iam)
       mr = nstart(m)
       mc = 2*mr
       do n=1,nlen(m),2
          ir = mc + 2*n - 1
          ii = ir + 1
!
          grpss(2*lm-1) = grpss(2*lm-1) + alps(ir)*alp(mr+n,irow)
          grpss(2*lm  ) = grpss(2*lm  ) + alps(ii)*alp(mr+n,irow)
!
          grdpss(2*lm-1) = grdpss(2*lm-1) + alps(ir)*alp(mr+n,irow)*hdfst4(m+n-1)*ztodt
          grdpss(2*lm  ) = grdpss(2*lm  ) + alps(ii)*alp(mr+n,irow)*hdfst4(m+n-1)*ztodt
!
       end do
    end do

    do lm=1,mlength
       m = locm(lm,iam)
       mr = nstart(m)
       mc = 2*mr
       do n=2,nlen(m),2
          ir = mc + 2*n - 1
          ii = ir + 1
!
          grpms(2*lm-1) = grpms(2*lm-1) + alps(ir)*dalp(mr+n,irow)*ra
          grpms(2*lm  ) = grpms(2*lm  ) + alps(ii)*dalp(mr+n,irow)*ra
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


subroutine grcalca (irow    ,ztodt   ,grta    ,grqa    ,grtha   , &
                    grda    ,grua    ,gruha   ,grva    ,grvha   , &
                    grpsa   ,grdpsa  ,grpma   ,grpla   ,grtma   , &
                    grtla   ,grqma   ,grqla   )

!-----------------------------------------------------------------------
!
! Purpose:
! Complete inverse legendre transforms from spectral to Fourier space at
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
! Thus "grts" contains the symmetric fourier coeffs of temperature and
! "grtha" contains the antisymmetric fourier coeffs of the temperature
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
    use rgrid
    use commap
    use physconst, only: ra
    use sld_control_mod
    use spmd_utils, only : iam
    implicit none

!
! Input arguments
!
    integer , intent(in)   :: irow              ! latitude pair index
    real(r8), intent(in)   :: ztodt             ! twice the timestep unless nstep = 0
!
! Output arguments: anti-symmetric fourier coefficients
!
    real(r8), intent(out) :: grta(2*maxm,plev)  ! sum(n) of t(n,m)*P(n,m)
    real(r8), intent(out) :: grqa(2*maxm,plev)  ! sum(n) of q(n,m)*P(n,m)
    real(r8), intent(out) :: grtha(2*maxm,plev) ! sum(n) of K(2i)*t(n,m)*P(n,m)
    real(r8), intent(out) :: grda(2*maxm,plev)  ! sum(n) of d(n,m)*P(n,m)
    real(r8), intent(out) :: grua(2*maxm,plev)  ! sum(n) of z(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: gruha(2*maxm,plev) ! sum(n) of K(2i)*z(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grva(2*maxm,plev)  ! sum(n) of d(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grvha(2*maxm,plev) ! sum(n) of K(2i)*d(n,m)*H(n,m)*a/(n(n+1))
    real(r8), intent(out) :: grpsa(2*maxm)      ! sum(n) of lnps(n,m)*P(n,m)
    real(r8), intent(out) :: grdpsa(2*maxm)     ! sum(n) of K(4)*(n(n+1)/a**2)**2*2dt
!                                               ! *lnps(n,m)*P(n,m)
    real(r8), intent(out) :: grpma(2*maxm)      ! sum(n) of lnps(n,m)*H(n,m)
    real(r8), intent(out) :: grpla(2*maxm)        ! sum(n) of lnps(n,m)*P(n,m)*m/a
    real(r8), intent(out) :: grtma (2*maxm,plev)
    real(r8), intent(out) :: grtla (2*maxm,plev)
    real(r8), intent(out) :: grqma (2*maxm,plev)
    real(r8), intent(out) :: grqla (2*maxm,plev)
!
!---------------------------Local workspace-----------------------------
!
    real(r8) gru1a (2*maxm)      ! sum(n) of d(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) gruh1a(2*maxm)      ! sum(n) of K(2i)*d(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) grv1a (2*maxm)      ! sum(n) of z(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) grvh1a(2*maxm)      ! sum(n) of K(2i)*z(n,m)*P(n,m)*m*a/(n(n+1))
    real(r8) alpn  (pspt)       ! (a*m/(n(n+1)))*Legendre functions (complex)
    real(r8) dalpn (pspt)       ! (a/(n(n+1)))*derivative of Legendre functions (complex)

    integer k                   ! level index
    integer lm, m               ! local and global Fourier wavenumber indices of spectral array
    integer mlength             ! number of local wavenumbers
    integer n                   ! meridional wavenumber index
    integer ir,ii               ! spectral indices
    integer mr,mc               ! spectral indices
    real(r8) tmp,raxm           ! temporary workspace
!
!-----------------------------------------------------------------------
!
! Compute alpn and dalpn
!
    mlength = numm(iam)
    do lm=1,mlength
       m = locm(lm,iam)
       mr = nstart(m)
       raxm = ra*xm(m)
       do n=1,nlen(m)
          alpn(mr+n) = alp(mr+n,irow)*rsq(m+n-1)*raxm
          dalpn(mr+n) = dalp(mr+n,irow)*rsq(m+n-1)*ra
       end do
    end do
!
! Initialize sums
!
    grta(:,:) = 0._r8
    grqa(:,:) = 0._r8
    grtha(:,:) = 0._r8
    grda(:,:)  = 0._r8
    grua(:,:)  = 0._r8
    gruha(:,:) = 0._r8
    grva(:,:)  = 0._r8
    grvha(:,:) = 0._r8
    grpsa(:)   = 0._r8
    grdpsa(:)   = 0._r8
    grpma(:)   = 0._r8
    grpla(:)   = 0._r8
    grtma(:,:)   = 0._r8
    grtla(:,:)   = 0._r8
    grqma(:,:)   = 0._r8
    grqla(:,:)   = 0._r8
!
!-----------------------------------------------------------------------
!
! Computation for multilevel variables
!
    do k=1,plev
!
! Initialize local sums
!
       gru1a(:) = 0._r8
       gruh1a(:) = 0._r8
       grv1a(:) = 0._r8
       grvh1a(:) = 0._r8
!
! Loop over n for t,q,d,and end of u and v
!
       do lm=1,mlength
          m = locm(lm,iam)
          mr = nstart(m)
          mc = 2*mr
          do n=1,nlen(m),2
             ir = mc + 2*n - 1
             ii = ir + 1
!
             grtma(2*lm-1,k) = grtma(2*lm-1,k) + t(ir,k)*dalp(mr+n,irow)*ra
             grtma(2*lm  ,k) = grtma(2*lm  ,k) + t(ii,k)*dalp(mr+n,irow)*ra
!
             grqma(2*lm-1,k) = grqma(2*lm-1,k) + q(ir,k)*dalp(mr+n,irow)*ra
             grqma(2*lm  ,k) = grqma(2*lm  ,k) + q(ii,k)*dalp(mr+n,irow)*ra
!
             grua (2*lm-1,k) = grua (2*lm-1,k) + vz(ir,k)*dalpn(mr+n)
             grua (2*lm  ,k) = grua (2*lm  ,k) + vz(ii,k)*dalpn(mr+n)
!
             tmp = dalpn(mr+n)*hdifzd(n+m-1,k)
             gruha(2*lm-1,k) = gruha(2*lm-1,k) - vz(ir,k)*tmp
             gruha(2*lm  ,k) = gruha(2*lm  ,k) - vz(ii,k)*tmp
!
             grva (2*lm-1,k) = grva (2*lm-1,k) - d(ir,k)*dalpn(mr+n)
             grva (2*lm  ,k) = grva (2*lm  ,k) - d(ii,k)*dalpn(mr+n)
!
             grvha(2*lm-1,k) = grvha(2*lm-1,k) + d(ir,k)*tmp
             grvha(2*lm  ,k) = grvha(2*lm  ,k) + d(ii,k)*tmp
          end do
       end do

       do lm=1,mlength
          m = locm(lm,iam)
          mr = nstart(m)
          mc = 2*mr
          do n=2,nlen(m),2
             ir = mc + 2*n - 1
             ii = ir + 1
             grta (2*lm-1,k) = grta (2*lm-1,k) + t(ir,k)*alp(mr+n,irow)
             grta (2*lm  ,k) = grta (2*lm  ,k) + t(ii,k)*alp(mr+n,irow)
!
             grqa (2*lm-1,k) = grqa (2*lm-1,k) + q(ir,k)*alp(mr+n,irow)
             grqa (2*lm  ,k) = grqa (2*lm  ,k) + q(ii,k)*alp(mr+n,irow)
!
             tmp = alp(mr+n,irow)*hdiftq(n+m-1,k)
             grtha(2*lm-1,k) = grtha(2*lm-1,k) - t(ir,k)*tmp
             grtha(2*lm  ,k) = grtha(2*lm  ,k) - t(ii,k)*tmp
!
             grda(2*lm-1,k) = grda(2*lm-1,k) + d(ir,k)*alp(mr+n,irow)
             grda(2*lm  ,k) = grda(2*lm  ,k) + d(ii,k)*alp(mr+n,irow)
!
             gru1a (2*lm-1) = gru1a (2*lm-1) + d(ir,k)*alpn(mr+n)
             gru1a (2*lm  ) = gru1a (2*lm  ) + d(ii,k)*alpn(mr+n)
!
             tmp = alpn(mr+n)*hdifzd(n+m-1,k)
             gruh1a(2*lm-1) = gruh1a(2*lm-1) - d(ir,k)*tmp
             gruh1a(2*lm  ) = gruh1a(2*lm  ) - d(ii,k)*tmp
!
             grv1a (2*lm-1) = grv1a (2*lm-1) + vz(ir,k)*alpn(mr+n)
             grv1a (2*lm  ) = grv1a (2*lm  ) + vz(ii,k)*alpn(mr+n)
!
             grvh1a(2*lm-1) = grvh1a(2*lm-1) - vz(ir,k)*tmp
             grvh1a(2*lm  ) = grvh1a(2*lm  ) - vz(ii,k)*tmp
          end do
       end do
!
! Combine the two parts of u(m) and v(m)
!
       do lm=1,mlength
          m = locm(lm,iam)
          grua (2*lm-1,k) = grua (2*lm-1,k) + gru1a (2*lm  )
          gruha(2*lm-1,k) = gruha(2*lm-1,k) + gruh1a(2*lm  )
          grua (2*lm  ,k) = grua (2*lm  ,k) - gru1a (2*lm-1)
          gruha(2*lm  ,k) = gruha(2*lm  ,k) - gruh1a(2*lm-1)
          grva (2*lm-1,k) = grva (2*lm-1,k) + grv1a (2*lm  )
          grvha(2*lm-1,k) = grvha(2*lm-1,k) + grvh1a(2*lm  )
          grva (2*lm  ,k) = grva (2*lm  ,k) - grv1a (2*lm-1)
          grvha(2*lm  ,k) = grvha(2*lm  ,k) - grvh1a(2*lm-1)
!
! Derivatives
!
          grtla(2*lm-1,k) = -grta(2*lm  ,k)*ra*xm(m)
          grtla(2*lm  ,k) =  grta(2*lm-1,k)*ra*xm(m)
          grqla(2*lm-1,k) = -grqa(2*lm  ,k)*ra*xm(m)
          grqla(2*lm  ,k) =  grqa(2*lm-1,k)*ra*xm(m)
       end do
    end do
!
!-----------------------------------------------------------------------
!
! Computation for 1-level variables (ln(p*) and derivatives).
!
    do lm=1,mlength
       m = locm(lm,iam)
       mr = nstart(m)
       mc = 2*mr
       do n=1,nlen(m),2
          ir = mc + 2*n - 1
          ii = ir + 1
!
          grpma(2*lm-1) = grpma(2*lm-1) + alps(ir)*dalp(mr+n,irow)*ra
          grpma(2*lm  ) = grpma(2*lm  ) + alps(ii)*dalp(mr+n,irow)*ra
       end do
    end do

    do lm=1,mlength
       m = locm(lm,iam)
       mr = nstart(m)
       mc = 2*mr
       do n=2,nlen(m),2
          ir = mc + 2*n - 1
          ii = ir + 1
!
          grpsa(2*lm-1) = grpsa(2*lm-1) + alps(ir)*alp(mr+n,irow)
          grpsa(2*lm  ) = grpsa(2*lm  ) + alps(ii)*alp(mr+n,irow)
!
          grdpsa(2*lm-1) = grdpsa(2*lm-1) + alps(ir)*alp(mr+n,irow)*hdfst4(m+n-1)*ztodt
          grdpsa(2*lm  ) = grdpsa(2*lm  ) + alps(ii)*alp(mr+n,irow)*hdfst4(m+n-1)*ztodt
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

