module commap

   use shr_kind_mod, only: r8 => shr_kind_r8
   use pmgrid,       only: plon, plev, plat
   use pspect,       only: pmmax, pnmax

   private

   public bps, sq, rsq, slat, w, cs, href, ecref, clat, clon
   public latdeg, zcr, bm1, bmi, tau, londeg, t0, xm

   real(r8) :: bps(plev)          ! coefficient for ln(ps) term in divergence eqn
   real(r8) :: sq(pnmax)          ! n(n+1)/a^2  (del^2 response function)
   real(r8) :: rsq(pnmax)         ! a^2/(n(n+1))
   real(r8) :: slat(plat/2)       ! |sine latitude| (hemisphere)
   real(r8), target :: w(plat)            ! gaussian weights (hemisphere)
   real(r8) :: cs(plat/2)         ! cosine squared latitude (hemisphere)
   real(r8) :: href(plev,plev)    ! reference hydrostatic equation matrix
   real(r8) :: ecref(plev,plev)   ! reference energy conversion matrix
   real(r8), target :: clat(plat)         ! model latitudes (radians)
   real(r8), target :: clon(plon,plat)    ! model longitudes (radians)
   real(r8), target :: latdeg(plat)       ! model latitudes (degrees)
   real(r8) :: zcr(pnmax,plev)    ! SLD only: real(r8) eigenvalues of semi-implicit matrix 
!                                 ! (divergence equation)
   real(r8) :: bm1(plev,plev)     ! transpose of right eigenvectors of semi-implicit matrix
   real(r8) :: bmi(plev,plev)     ! SLD only: transpose of inverse of the right eigenvectors
!                                 ! (left eigenvectors) of the semi-implicit matrix
   real(r8) :: tau(plev,plev )    ! matrix for reference d term in thermodynamic eqn
   real(r8), target :: londeg(plon,plat)  ! model longitudes (degrees)
   real(r8) :: t0(plev)           ! Reference temperature for t-prime computations
   real(r8) :: xm(pmmax)          ! m (longitudinal wave number)
end module commap
