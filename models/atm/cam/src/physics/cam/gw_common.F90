module gw_common

!
! This module contains code common to different gravity wave
! parameterizations.
!
use gw_utils, only: r8

implicit none
private
save

! Public interface.
public :: gw_common_init
public :: gw_prof
public :: momentum_energy_conservation
public :: gw_drag_prof

public :: pver, pgwv
public :: dc
public :: cref
public :: do_molec_diff
public :: nbot_molec
public :: west, east, north, south
public :: fcrit2
public :: kwv
public :: gravit
public :: rair

! This flag preserves answers for vanilla CAM by making a few changes (e.g.
! order of operations) when only orographic waves are on.
logical, public :: orographic_only = .false.

! Number of levels in the atmosphere.
integer, protected :: pver = 0

! Maximum number of waves allowed (i.e. wavenumbers are -pgwv:pgwv).
integer, protected :: pgwv = 0

! Bin width for spectrum.
real(r8), protected :: dc = huge(1._r8)

! Reference speeds for the spectrum.
real(r8), allocatable, protected :: cref(:)

! Whether or not molecular diffusion is being done, and bottom level where
! it is done.
logical, protected :: do_molec_diff = .false.
integer, protected :: nbot_molec = huge(1)

! Index the cardinal directions.
integer, parameter :: west = 1
integer, parameter :: east = 2
integer, parameter :: south = 3
integer, parameter :: north = 4

! Critical Froude number.
real(r8), protected :: fcrit2 = huge(1._r8)

! Effective horizontal wave number.
real(r8), protected :: kwv = huge(1._r8)

! Acceleration due to gravity.
real(r8), protected :: gravit = huge(1._r8)

! Gas constant for dry air.
real(r8), protected :: rair = huge(1._r8)

!
! Private variables
!

! Interface levels for gravity wave sources.
integer :: ktop = huge(1)
integer :: kbotbg = huge(1)

! Effective wavenumber.
real(r8) :: effkwv = huge(1._r8)

! Background diffusivity.
real(r8), parameter :: dback = 0.05_r8

! rair/gravit
real(r8) :: rog = huge(1._r8)

! Newtonian cooling coefficients.
real(r8), allocatable :: alpha(:)

!
! Limits to keep values reasonable.
!

! Minimum non-zero stress.
real(r8), parameter :: taumin = 1.e-10_r8
! Maximum wind tendency from stress divergence (before efficiency applied).
real(r8) :: tndmax = huge(1._r8)
! Maximum allowed change in u-c (before efficiency applied).
real(r8), parameter :: umcfac = 0.5_r8
! Minimum value of (u-c)**2.
real(r8), parameter :: ubmc2mn = 0.01_r8

contains

!==========================================================================

subroutine gw_common_init(pver_in, pgwv_in, dc_in, cref_in, &
     do_molec_diff_in, nbot_molec_in, ktop_in, kbotbg_in, fcrit2_in, &
     kwv_in, gravit_in, rair_in, alpha_in, errstring)

  integer,  intent(in) :: pver_in
  integer,  intent(in) :: pgwv_in
  real(r8), intent(in) :: dc_in
  real(r8), intent(in) :: cref_in(-pgwv_in:)
  logical,  intent(in) :: do_molec_diff_in
  integer,  intent(in) :: nbot_molec_in
  integer,  intent(in) :: ktop_in
  integer,  intent(in) :: kbotbg_in
  real(r8), intent(in) :: fcrit2_in
  real(r8), intent(in) :: kwv_in
  real(r8), intent(in) :: gravit_in
  real(r8), intent(in) :: rair_in
  real(r8), intent(in) :: alpha_in(0:)
  ! Report any errors from this routine.
  character(len=*), intent(out) :: errstring

  integer :: ierr

  errstring = ""

  pver = pver_in
  pgwv = pgwv_in
  dc = dc_in
  allocate(cref(-pgwv:pgwv), stat=ierr, errmsg=errstring)
  if (ierr /= 0) return
  cref = cref_in
  do_molec_diff = do_molec_diff_in
  nbot_molec = nbot_molec_in
  ktop = ktop_in
  kbotbg = kbotbg_in
  fcrit2 = fcrit2_in
  kwv = kwv_in
  gravit = gravit_in
  rair = rair_in
  allocate(alpha(0:pver), stat=ierr, errmsg=errstring)
  if (ierr /= 0) return
  alpha = alpha_in

  effkwv = kwv * fcrit2
  rog = rair/gravit

  if (.not. orographic_only) then
     ! 400 m/s/day
     tndmax  = 400._r8 / 86400._r8
  else
     ! 500 m/s/day
     tndmax  = 500._r8 / 86400._r8
  end if

end subroutine gw_common_init

!==========================================================================

subroutine gw_prof (ncol, cpair, t, pmid, pint, rhoi, ti, nm, ni)
  !-----------------------------------------------------------------------
  ! Compute profiles of background state quantities for the multiple
  ! gravity wave drag parameterization.
  !
  ! The parameterization is assumed to operate only where water vapor
  ! concentrations are negligible in determining the density.
  !-----------------------------------------------------------------------
  use gw_utils, only: midpoint_interp
  !------------------------------Arguments--------------------------------
  ! Column dimension.
  integer, intent(in) :: ncol

  ! Specific heat of dry air, constant pressure.
  real(r8), intent(in) :: cpair
  ! Midpoint temperatures.
  real(r8), intent(in) :: t(ncol,pver)
  ! Midpoint and interface pressures.
  real(r8), intent(in) :: pmid(ncol,pver), pint(ncol,0:pver)

  ! Interface density.
  real(r8), intent(out) :: rhoi(ncol,0:pver)
  ! Interface temperature.
  real(r8), intent(out) :: ti(ncol,0:pver)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(out) :: nm(ncol,pver), ni(ncol,0:pver)

  !---------------------------Local Storage-------------------------------
  ! Column and level indices.
  integer :: i,k

  ! dt/dp
  real(r8) :: dtdp
  ! Brunt-Vaisalla frequency squared.
  real(r8) :: n2

  ! Minimum value of Brunt-Vaisalla frequency squared.
  real(r8), parameter :: n2min = 1.e-8_r8

  !------------------------------------------------------------------------
  ! Determine the interface densities and Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------

  ! The top interface values are calculated assuming an isothermal
  ! atmosphere above the top level.
  k = 0
  do i = 1, ncol
     ti(i,k) = t(i,k+1)
     rhoi(i,k) = pint(i,k) / (rair*ti(i,k))
     ni(i,k) = sqrt(gravit*gravit / (cpair*ti(i,k)))
  end do

  ! Interior points use centered differences.
  ti(:,1:pver-1) = midpoint_interp(t)
  do k = 1, pver-1
     do i = 1, ncol
        rhoi(i,k) = pint(i,k) / (rair*ti(i,k))
        dtdp = (t(i,k+1)-t(i,k)) / (pmid(i,k+1)-pmid(i,k))
        n2 = gravit*gravit/ti(i,k) * (1._r8/cpair - rhoi(i,k)*dtdp)
        ni(i,k) = sqrt(max(n2min, n2))
     end do
  end do

  ! Bottom interface uses bottom level temperature, density; next interface
  ! B-V frequency.
  k = pver
  do i = 1, ncol
     ti(i,k) = t(i,k)
     rhoi(i,k) = pint(i,k) / (rair*ti(i,k))
     ni(i,k) = ni(i,k-1)
  end do

  !------------------------------------------------------------------------
  ! Determine the midpoint Brunt-Vaisala frequencies.
  !------------------------------------------------------------------------
  nm = midpoint_interp(ni)

end subroutine gw_prof

!==========================================================================

subroutine momentum_energy_conservation(ncol, tend_level, dt, taucd, &
     pint, pdel, u, v, dudt, dvdt, dsdt, utgw, vtgw, ttgw)

  ! C.-C. Chen, momentum & energy conservation

  integer, intent(in) :: ncol
  integer, intent(in) :: tend_level(ncol)
  real(r8), intent(in) :: dt
  real(r8), intent(in) :: taucd(ncol,0:pver,4)
  real(r8), intent(in) :: pint(ncol,pver)
  real(r8), intent(in) :: pdel(ncol,pver)
  real(r8), intent(in) :: u(ncol,pver)
  real(r8), intent(in) :: v(ncol,pver)
  ! The following are assumed-shape because, for now, CAM is passing them
  ! with size pcols.
  real(r8), intent(inout) :: dudt(:,:)
  real(r8), intent(inout) :: dvdt(:,:)
  real(r8), intent(inout) :: dsdt(:,:)
  real(r8), intent(inout) :: utgw(ncol,pver)
  real(r8), intent(inout) :: vtgw(ncol,pver)
  real(r8), intent(inout) :: ttgw(ncol,pver)
  ! Local variables.
  integer :: i, k
  real(r8) :: dz(ncol), dE(ncol), ut_dz(ncol), vt_dz(ncol)

  ! Total mass from ground to source level: rho*dz = dp/gravit
  dz = 0.0_r8
  do k = pver, minval(tend_level)+1, -1
     where (k > tend_level)
        dz = dz + pdel(:,k)/gravit
     end where
  end do

  ! Tendency for U & V below source level.
  do i = 1, ncol
     ut_dz(i) = -(taucd(i,tend_level(i), east) + &
          taucd(i,tend_level(i), west))/dz(i)
     vt_dz(i) = -(taucd(i,tend_level(i),north) + &
          taucd(i,tend_level(i),south))/dz(i)
  end do

  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        dudt(:ncol,k) = dudt(:ncol,k)+ut_dz
        dvdt(:ncol,k) = dvdt(:ncol,k)+vt_dz
        utgw(:,k) = utgw(:,k)+ut_dz
        vtgw(:,k) = vtgw(:,k)+vt_dz
     end where
  end do

  ! Net gain/loss of total energy in the column.
  dE = 0.0_r8
  do k = 1, pver
     dE = dE + pdel(:,k) * (dsdt(:ncol,k) + &
          dudt(:ncol,k)*(u(:,k)+dudt(:ncol,k)*0.5_r8*dt) + &
          dvdt(:ncol,k)*(v(:,k)+dvdt(:ncol,k)*0.5_r8*dt) )
  end do

  do i = 1, ncol
     dE(i) = dE(i)/(pint(i,pver)-pint(i,tend_level(i)))
  end do

  ! Subtract net gain/loss of total energy below source level.
  do k = minval(tend_level)+1, pver
     where (k > tend_level)
        dsdt(:ncol,k) = dsdt(:ncol,k)-dE
        ttgw(:,k) = ttgw(:,k)-dE
     end where
  end do

end subroutine momentum_energy_conservation

!==========================================================================

subroutine gw_drag_prof(ncol, ngwv, src_level, tend_level, do_taper, dt, &
     lat,           t,    ti,  pmid, pint, dpm,   rdpm, &
     piln, rhoi,    nm,   ni,  ubm,  ubi,  xv,    yv,   &
     effgw,      c, kvtt, q,   dse,  tau,  utgw,  vtgw, &
     ttgw, qtgw, taucd,   egwdffi,   gwut, dttdf, dttke)

  !-----------------------------------------------------------------------
  ! Solve for the drag profile from the multiple gravity wave drag
  ! parameterization.
  ! 1. scan up from the wave source to determine the stress profile
  ! 2. scan down the stress profile to determine the tendencies
  !     => apply bounds to the tendency
  !          a. from wkb solution
  !          b. from computational stability constraints
  !     => adjust stress on interface below to reflect actual bounded
  !        tendency
  !-----------------------------------------------------------------------

  use gw_diffusion, only: gw_ediff, gw_diff_tend
  use vdiff_lu_solver, only: lu_decomp

  !------------------------------Arguments--------------------------------
  ! Column and gravity wave spectrum dimensions.
  integer, intent(in) :: ncol, ngwv
  ! Level from which gravity waves are propagated upward.
  integer, intent(in) :: src_level(ncol)
  ! Lowest level where wind tendencies are calculated.
  integer, intent(in) :: tend_level(ncol)
  ! Using tend_level > src_level allows the orographic waves to prescribe
  ! wave propagation up to a certain level, but then allow wind tendencies
  ! and adjustments to tau below that level.

  ! Whether or not to apply the polar taper.
  logical, intent(in) :: do_taper
  ! Time step.
  real(r8), intent(in) :: dt

  ! Latitudes for each column.
  real(r8), intent(in) :: lat(ncol)
  ! Midpoint and interface temperatures.
  real(r8), intent(in) :: t(ncol,pver), ti(ncol,0:pver)
  ! Midpoint and interface pressures.
  real(r8), intent(in) :: pmid(ncol,pver), pint(ncol,0:pver)
  ! Layer thickness (delta pressure) and reciprocal of layer thickness.
  real(r8), intent(in) :: dpm(ncol,pver), rdpm(ncol,pver)
  ! Log of interface pressures.
  real(r8), intent(in) :: piln(ncol,0:pver)
  ! Interface densities.
  real(r8), intent(in) :: rhoi(ncol,0:pver)
  ! Midpoint and interface Brunt-Vaisalla frequencies.
  real(r8), intent(in) :: nm(ncol,pver), ni(ncol,0:pver)
  ! Projection of wind at midpoints and interfaces.
  real(r8), intent(in) :: ubm(ncol,pver), ubi(ncol,0:pver)
  ! Unit vectors of source wind (zonal and meridional components).
  real(r8), intent(in) :: xv(ncol), yv(ncol)
  ! Tendency efficiency.
  real(r8), intent(in) :: effgw
  ! Wave phase speeds for each column.
  real(r8), intent(in) :: c(ncol,-pgwv:pgwv)
  ! Molecular thermal diffusivity.
  real(r8), intent(in) :: kvtt(ncol,0:pver)
  ! Constituent array.
  real(r8), intent(in) :: q(:,:,:)
  ! Dry static energy.
  real(r8), intent(in) :: dse(ncol,pver)

  ! Wave Reynolds stress.
  real(r8), intent(inout) :: tau(ncol,-pgwv:pgwv,0:pver)
  ! Zonal/meridional wind tendencies.
  real(r8), intent(out) :: utgw(ncol,pver), vtgw(ncol,pver)
  ! Gravity wave heating tendency.
  real(r8), intent(out) :: ttgw(ncol,pver)
  ! Gravity wave constituent tendency.
  real(r8), intent(out) :: qtgw(:,:,:)

  ! Reynolds stress for waves propagating in each cardinal direction.
  real(r8), intent(out) :: taucd(ncol,0:pver,4)

  ! Effective gravity wave diffusivity at interfaces.
  real(r8), intent(out) :: egwdffi(ncol,0:pver)

  ! Gravity wave wind tendency for each wave.
  real(r8), intent(out) :: gwut(ncol,pver,-ngwv:ngwv)

  ! Temperature tendencies from diffusion and kinetic energy.
  real(r8), intent(out) :: dttdf(ncol,pver)
  real(r8), intent(out) :: dttke(ncol,pver)

  !---------------------------Local storage-------------------------------
  ! Column, level, wavenumber, and constituent loop indices.
  integer :: i, k, l, m

  ! "Total" and saturation diffusivity.
  real(r8) :: d(ncol), dsat(ncol)
  ! Fraction of dsat to use.
  real(r8) :: dscal(ncol)
  ! Imaginary part of vertical wavenumber.
  real(r8) :: mi(ncol)
  ! Stress after damping.
  real(r8) :: taudmp(ncol)
  ! Saturation stress.
  real(r8) :: tausat(ncol,-ngwv:ngwv)
  ! ubi at tend_level.
  real(r8) :: ubi_tend(ncol)
  ! (ub-c) and (ub-c)**2
  real(r8) :: ubmc(ncol,-ngwv:ngwv), ubmc2(ncol)
  ! Temporary ubar tendencies (overall, at wave l, and at saturation).
  real(r8) :: ubt(ncol,pver), ubtl(ncol), ubtlsat(ncol)
  real(r8) :: wrk(ncol)

  ! Signed wave Reynolds stress.
  real(r8) :: tausg(ncol)

  ! Reynolds stress for waves propagating behind and forward of the wind.
  real(r8) :: taub(ncol)
  real(r8) :: tauf(ncol)

  ! Polar taper.
  real(r8) :: ptaper(ncol)

  ! Recalculated rhoi to preserve answers.
  real(r8) :: rhoi_kludge(ncol,pver+1)

  ! LU decomposition.
  type(lu_decomp) :: decomp

  !------------------------------------------------------------------------

  ! Initialize gravity wave drag tendencies to zero.

  utgw = 0._r8
  vtgw = 0._r8

  taucd = 0._r8

  gwut = 0._r8

  dttke = 0._r8
  ttgw = 0._r8

  ! Workaround floating point exception issues on Intel by initializing
  ! everything that's first set in a where block.
  dsat = 0._r8
  dscal = 0._r8
  mi = 0._r8
  taudmp = 0._r8
  tausat = 0._r8
  tausg = 0._r8
  ubmc = 0._r8
  ubmc2 = 0._r8
  wrk = 0._r8

  !------------------------------------------------------------------------
  ! Compute the stress profiles and diffusivities
  !------------------------------------------------------------------------

  ! Loop from bottom to top to get stress profiles.
  do k = maxval(src_level)-1, ktop, -1

     ! Determine the absolute value of the saturation stress.
     ! Define critical levels where the sign of (u-c) changes between
     ! interfaces.

     do l = -ngwv, ngwv
        where (src_level > k)

           ubmc(:,l) = ubi(:,k) - c(:,l)

           ! Test to see if u-c has the same sign here as the level below.
           where (ubmc(:,l) * (ubi(:,k+1)-c(:,l)) > 0.0_r8)
              tausat(:,l) = abs(effkwv * rhoi(:,k) * ubmc(:,l)**3 / &
                   (2._r8*ni(:,k)))
              where (tausat(:,l) <= taumin) tausat(:,l) = 0._r8
           elsewhere
              tausat(:,l) = 0.0_r8
           end where

        end where
     end do

     ! Determine the diffusivity for each column.

     d = dback
     if (do_molec_diff) then
        d = d + kvtt(:,k)
     else
        do l = -ngwv, ngwv
           where (src_level > k)
              dsat = (ubmc(:,l) / ni(:,k))**2 * &
                   (effkwv * ubmc(:,l)**2 / (rog * ti(:,k) * ni(:,k)) - &
                   alpha(k))
              dscal = min(1.0_r8, tau(:,l,k+1) / (tausat(:,l)+taumin))
              d = max(d, dscal * dsat)
           end where
        end do
     end if

     ! Compute stress for each wave. The stress at this level is the min of
     ! the saturation stress and the stress at the level below reduced by
     ! damping. The sign of the stress must be the same as at the level
     ! below.

     ! If molecular diffusion is on, only do this in levels with molecular
     ! diffusion. Otherwise, do it everywhere.
     if (k <= nbot_molec .or. .not. do_molec_diff) then
        do l = -ngwv, ngwv
           where (src_level > k)

              ubmc2 = max(ubmc(:,l)**2, ubmc2mn)
              mi = ni(:,k) / (2._r8 * kwv * ubmc2) * &
                   (alpha(k) + ni(:,k)**2/ubmc2 * d)
              wrk = -2._r8*mi*rog*t(:,k+1)*(piln(:,k+1) - piln(:,k))

              where (wrk >= -150._r8 .or. .not. do_molec_diff)
                 taudmp = tau(:,l,k+1) * exp(wrk)
              elsewhere
                 taudmp = 0._r8
              end where
              where (taudmp <= taumin) taudmp = 0._r8
              tau(:,l,k) = min(taudmp, tausat(:,l))

           end where
        end do
     else
        do l = -ngwv, ngwv
           where (src_level > k)
              tau(:,l,k) = min(tau(:,l,k+1), tausat(:,l))
           end where
        end do
     end if

  end do


  ! Tau projected in the four cardinal directions, for the momentum
  ! conservation routine and for diagnostic output.
  if ( ngwv > 0) then

     ubi_tend = (/ (ubi(i,tend_level(i)), i = 1, ncol) /)

     do k = ktop, maxval(tend_level)

        taub = 0._r8
        tauf = 0._r8

        do l = -ngwv, ngwv
           where (k <= tend_level)

              tausg = sign(tau(:,l,k), c(:,l)-ubi(:,k))

              where ( c(:,l) < ubi_tend )
                 taub = taub + tausg
              elsewhere (  c(:,l) > ubi_tend )
                 tauf = tauf + tausg
              end where

           end where
        end do

        where (k <= tend_level)
           where (xv > 0._r8)
              taucd(:,k,east) = tauf * xv
              taucd(:,k,west) = taub * xv
           elsewhere ( xv < 0._r8)
              taucd(:,k,east) = taub * xv
              taucd(:,k,west) = tauf * xv
           end where

           where ( yv > 0._r8)
              taucd(:,k,north) = tauf * yv
              taucd(:,k,south) = taub * yv
           elsewhere ( yv < 0._r8)
              taucd(:,k,north) = taub * yv
              taucd(:,k,south) = tauf * yv
           end where
        end where

     end do

  end if

  !------------------------------------------------------------------------
  ! Compute the tendencies from the stress divergence.
  !------------------------------------------------------------------------

  if (do_taper) then    ! taper CM only
     ptaper = cos(lat)
  else                  ! do not taper other cases
     ptaper = 1._r8
  endif

  ! Loop over levels from top to bottom
  do k = ktop+1, maxval(tend_level)

     ! Accumulate the mean wind tendency over wavenumber.
     ubt(:,k) = 0.0_r8

     do l = -ngwv, ngwv    ! loop over wave

        ! Determine the wind tendency, including excess stress carried down
        ! from above.
        ubtl = gravit * (tau(:,l,k)-tau(:,l,k-1)) * rdpm(:,k)

        if (orographic_only) then
           ! Require that the tendency be no larger than the analytic
           ! solution for a saturated region [proportional to (u-c)^3].
           ubtlsat = effkwv * abs((c(:,l)-ubm(:,k))**3) / &
                (2._r8*rog*t(:,k)*nm(:,k))
           ubtl = min(ubtl, ubtlsat)
        end if

        ! Apply tendency limits to maintain numerical stability.
        ! 1. du/dt < |c-u|/dt  so u-c cannot change sign
        !    (u^n+1 = u^n + du/dt * dt)
        ! 2. du/dt < tndmax    so that ridicuously large tendencies are not
        !    permitted
        ubtl = min(ubtl, umcfac * abs(c(:,l)-ubm(:,k)) / dt)
        ubtl = min(ubtl, tndmax)

        where (k <= tend_level)

           ! Save tendency for each wave (for later computation of kzz),
           ! applying efficiency and taper:
           gwut(:,k,l) = sign(ubtl, c(:,l)-ubm(:,k)) * effgw * ptaper

        end where

        if (.not. orographic_only) then
           where (k <= tend_level)
              ubt(:,k) = ubt(:,k) + gwut(:,k,l)
           end where
        else
           where (k <= tend_level)
              ubt(:,k) = ubt(:,k) + sign(ubtl, c(:,l)-ubm(:,k))
           end where
        end if

        where (k <= tend_level)

           ! Redetermine the effective stress on the interface below from
           ! the wind tendency. If the wind tendency was limited above,
           ! then the new stress will be smaller than the old stress, 
           ! causing stress divergence in the next layer down. This
           ! smoothes large stress divergences downward while conserving
           ! total stress.
           tau(:,l,k) = tau(:,l,k-1) + ubtl * dpm(:,k) / gravit

        end where

     end do

     ! Project the mean wind tendency onto the components.
     if (.not. orographic_only) then
        where (k <= tend_level)
           utgw(:,k) = ubt(:,k) * xv
           vtgw(:,k) = ubt(:,k) * yv
        end where
     else
        where (k <= tend_level)
           utgw(:,k) = ubt(:,k) * xv * effgw * ptaper
           vtgw(:,k) = ubt(:,k) * yv * effgw * ptaper
        end where
     end if

     ! End of level loop.
  end do

  if (ngwv > 0) then

     ! Precalculate rhoi for the following routines. We have rhoi, but
     ! recalculate it here to preserve answers. (Ideal gas law.)

     ! Note: pint is from 0 to pver instead of 1 to pver+1.
     rhoi_kludge(:,1) = pint(:,0) / (rair * t(:,1))
     do k = 2, pver
        rhoi_kludge(:,k) = pint(:,k-1) * 2._r8 / &
             (rair * (t(:,k) + t(:,k-1)))
     end do
     rhoi_kludge(:,pver+1) = pint(:,pver) / (rair * t(:,pver))

     ! Calculate effective diffusivity and LU decomposition for the
     ! vertical diffusion solver.
     call gw_ediff (ncol, pver, ngwv, kbotbg, ktop, tend_level, &
          gwut, ubm, nm, rhoi_kludge, dt, gravit, pmid, rdpm, c, &
          egwdffi, decomp)

     ! Calculate tendency on each constituent.
     do m = 1, size(q,3)

        call gw_diff_tend(ncol, pver, kbotbg, ktop, q(:,:,m), dt, &
             decomp, qtgw(:,:,m))

     enddo

     ! Calculate tendency from diffusing dry static energy (dttdf).
     call gw_diff_tend(ncol, pver, kbotbg, ktop, dse, dt, decomp, dttdf)

     ! Evaluate second temperature tendency term: Conversion of kinetic
     ! energy into thermal.
     do l = -ngwv, ngwv
        do k = ktop+1, kbotbg
           dttke(:,k) = dttke(:,k) + c(:,l) * gwut(:,k,l)
        end do
     end do

     ttgw = dttke + dttdf

     ! Deallocate decomp.
     call decomp%finalize()

  else

     qtgw = 0._r8
     ttgw = 0._r8

  end if

end subroutine gw_drag_prof

end module gw_common
