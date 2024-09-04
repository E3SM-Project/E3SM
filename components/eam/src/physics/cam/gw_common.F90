module gw_common

!
! This module contains code common to different gravity wave
! parameterizations.
!
use gw_utils, only: r8 
!!======Jinbo Xie=======
use ppgrid,        only: nvar_dirOA,nvar_dirOL!pcols,pver,pverp,
use cam_logfile,   only: iulog
!!======Jinbo Xie=======

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

!!======Jinbo Xie=======
public :: gwdo_gsd,pblh_get_level_idx,grid_size
!!======Jinbo Xie=======

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

! Whether or not to enforce an upper boundary condition of tau = 0.
logical :: tau_0_ubc = .false.

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
     do_molec_diff_in, tau_0_ubc_in, nbot_molec_in, ktop_in, kbotbg_in, &
     fcrit2_in, kwv_in, gravit_in, rair_in, alpha_in, errstring)

  integer,  intent(in) :: pver_in
  integer,  intent(in) :: pgwv_in
  real(r8), intent(in) :: dc_in
  real(r8), intent(in) :: cref_in(-pgwv_in:)
  logical,  intent(in) :: do_molec_diff_in
  logical,  intent(in) :: tau_0_ubc_in
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
  tau_0_ubc = tau_0_ubc_in
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
  real(r8), intent(in) :: pint(ncol,pver+1)
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
     dE(i) = dE(i)/(pint(i,pver+1)-pint(i,tend_level(i)+1))
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

  ! Force tau at the top of the model to zero, if requested.
  if (tau_0_ubc) tau(:,:,0) = 0._r8

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




!===========Jinbo Xie===============================
function pblh_get_level_idx(height_array ,pblheight)
implicit none
real(8),intent(in),dimension(30) :: height_array
real(8),intent(in) :: pblheight
integer :: pblh_get_level_idx
       
!local 
integer :: i
logical :: found

pblh_get_level_idx = -1
found=.False.
          
do i = 1, pver
        if((pblheight >= height_array(i+1).and.pblheight <height_array(i)))then
                pblh_get_level_idx =  pver+1-i           
                found=.True.
                return
        endif
enddo
end function
!================================Jinbo Xie====================
subroutine grid_size(state, grid_dx, grid_dy)
  ! Determine the size of the grid for each of the columns in state

  use phys_grid,       only: get_area_p
  use shr_const_mod,   only: shr_const_pi
  use physics_types,   only: physics_state
  use ppgrid,          only: pver, pverp, pcols

  type(physics_state), intent(in) :: state
  real(r8), intent(out)           :: grid_dx(pcols), grid_dy(pcols)   ! E3SM grid [m]

  real(r8), parameter :: earth_ellipsoid1 = 111132.92_r8 ! World Geodetic System 1984 (WGS84) 
                                                         ! first coefficient, meters per degree longitude at equator
  real(r8), parameter :: earth_ellipsoid2 = 559.82_r8 ! second expansion coefficient for WGS84 ellipsoid
  real(r8), parameter :: earth_ellipsoid3 = 1.175_r8 ! third expansion coefficient for WGS84 ellipsoid

  real(r8) :: mpdeglat, column_area, degree, lat_in_rad
  integer  :: i

  do i=1,state%ncol
      ! determine the column area in radians
      column_area = get_area_p(state%lchnk,i)
      ! convert to degrees
      degree = sqrt(column_area)*(180._r8/shr_const_pi)

      ! convert latitude to radians
      lat_in_rad = state%lat(i)*(shr_const_pi/180._r8)

      ! Now find meters per degree latitude
      ! Below equation finds distance between two points on an ellipsoid, derived from expansion
      !  taking into account ellipsoid using World Geodetic System (WGS84) reference 
      mpdeglat = earth_ellipsoid1 - earth_ellipsoid2 * cos(2._r8*lat_in_rad) + earth_ellipsoid3 * cos(4._r8*lat_in_rad)
      grid_dx(i) = mpdeglat * degree
      grid_dy(i) = grid_dx(i) ! Assume these are the same
  enddo
end subroutine grid_size
!================================Jinbo Xie====================
!-------------------------------------------------------------------------------
   subroutine gwdo_gsd(u3d,v3d,t3d,qv3d,p3d,p3di,pi3d,z,                       &
                  rublten,rvblten,rthblten,                                    &
                  dtaux3d_ls,dtauy3d_ls,dtaux3d_bl,dtauy3d_bl,                 &
                  dtaux3d_ss,dtauy3d_ss,dtaux3d_fd,dtauy3d_fd,                 &
                  dusfcg_ls,dvsfcg_ls,dusfcg_bl,dvsfcg_bl,dusfcg_ss,dvsfcg_ss, &
                  dusfcg_fd,dvsfcg_fd,xland,br,                                &
                  var2d,oc12d,oa2d,ol2d,znu,znw,p_top,dz,pblh,                 &
                  cp,g,rd,rv,ep1,pi,bnvbg,                                     &
                  dt,dx,dy,kpbl2d,itimestep,gwd_opt,                           &
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte,                                 &
                    gwd_ls,gwd_bl,gwd_ss,gwd_fd)!Jinbo Xie added
!Jinbo Xie add dy, since global model is not dx=dy
!=================================
!   Jinbo Xie0 Modification
!changed the index used in WRF
!to chunk in CAM
!==============================
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
!                                                                       
!-- u3d         3d u-velocity interpolated to theta points (m/s)
!-- v3d         3d v-velocity interpolated to theta points (m/s)
!-- t3d         temperature (k)
!-- qv3d        3d water vapor mixing ratio (kg/kg)
!-- p3d         3d pressure (pa)
!-- p3di        3d pressure (pa) at interface level
!-- pi3d        3d exner function (dimensionless)
!-- rublten     u tendency due to pbl parameterization (m/s/s) 
!-- rvblten     v tendency due to pbl parameterization (m/s/s)
!-- rthblten    theta tendency due to pbl parameterization (K/s)
!-- znu         eta values (sigma values)
!-- cp          heat capacity at constant pressure for dry air (j/kg/k)
!-- g           acceleration due to gravity (m/s^2)
!-- rd          gas constant for dry air (j/kg/k)
!-- z           height above sea level (m)
!-- rv          gas constant for water vapor (j/kg/k)
!-- dt          time step (s)
!-- dx          model grid interval (m)
!-- dz          height of model layers (m)
!-- xland       land mask (1 for land, 2 for water)
!-- br          bulk richardson number in surface layer
!-- pblh        planetary boundary layer height (m)
!-- ep1         constant for virtual temperature (r_v/r_d - 1) (dimensionless)
!-- ids         start index for i in domain
!-- ide         end index for i in domain
!-- jds         start index for j in domain
!-- jde         end index for j in domain
!-- kds         start index for k in domain
!-- kde         end index for k in domain
!-- ims         start index for i in memory
!-- ime         end index for i in memory
!-- jms         start index for j in memory
!-- jme         end index for j in memory
!-- kms         start index for k in memory
!-- kme         end index for k in memory
!-- its         start index for i in tile
!-- ite         end index for i in tile
!-- jts         start index for j in tile
!-- jte         end index for j in tile
!-- kts         start index for k in tile
!-- kte         end index for k in tile
!-------------------------------------------------------------------------------
  integer,  intent(in   )   ::      ids,ide, jds,jde, kds,kde,                 &
                                     ims,ime, jms,jme, kms,kme,                &
                                     its,ite, jts,jte, kts,kte
  integer,  intent(in   )   ::      itimestep,gwd_opt
        real(r8),     intent(in   )   ::      cp,g,rd,rv,ep1,pi
!======Jinbo Xie=========
        real(r8),     intent(in), optional   ::  dt
        real(r8),     intent(in), dimension( ims:ime, kms:kme ),optional   ::  bnvbg
!======Jinbo Xie=========
        real(r8),    intent(in)   ::     dx(:)
        real(r8),    intent(in)   ::     dy(:)
!==========================
!
  real(r8),     dimension( ims:ime, kms:kme )             ,  &

            intent(in)   ::                      qv3d, &
                                                              p3d, &
                                                             pi3d, &
                                                              t3d, &
                                                                z, &
                                                               dz
  real(r8),     dimension( ims:ime, kms:kme )                    , &
     intent(in   )   ::                                 p3di
  real(r8),     dimension( ims:ime, kms:kme )                    , &
           intent(inout)   ::                   rublten, &
                                                          rvblten, &
                                                          rthblten
  real(r8),     dimension( ims:ime, kms:kme ), optional                 , &
            intent(inout)   ::  dtaux3d_ls,dtauy3d_ls,dtaux3d_bl,dtauy3d_bl,   &
                                dtaux3d_ss,dtauy3d_ss,dtaux3d_fd,dtauy3d_fd
!
  real(r8),     dimension( ims:ime, kms:kme)   ::                                    &
                                  dtaux2d_ls,dtauy2d_ls,dtaux2d_bl,dtauy2d_bl, &
                                  dtaux2d_ss,dtauy2d_ss,dtaux2d_fd,dtauy2d_fd

  real(r8),      dimension( ims:ime, kms:kme )                          , &

             intent(in   )   ::                        u3d, &
                                                                v3d
!
  integer,   dimension( ims:ime )                                   , &
             intent(in  )   ::                             kpbl2d

  real(r8),   dimension( ims:ime )                                      , &
        intent(in  )   ::                                   pblh, &
                                                                 br, &
                                                                 xland

  real(r8),   dimension( ims:ime ), optional                            , &
             intent(inout  )   ::  dusfcg_ls,dvsfcg_ls,dusfcg_bl,dvsfcg_bl,    &
                                   dusfcg_ss,dvsfcg_ss,dusfcg_fd,dvsfcg_fd

  real(r8),   dimension( ims:ime ) ::  dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,        &
                                       dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd
!Jinbo Xie
        real(r8),   dimension( ims:ime ),     optional               , &
             intent(in  )   ::                                  var2d, &
                                                                oc12d
       real(r8),dimension(ims:ime,nvar_dirOL),optional, intent(in) :: ol2d
       real(r8),dimension(ims:ime,nvar_dirOA),optional, intent(in) :: oa2d
!Jinbo Xie
            real(r8),   optional                                                         , &
            intent(in)   ::                                             znu(:), &
                                                                        znw(:)
!
  real(r8),     optional, intent(in)         ::                           p_top
!
!local
!
  real(r8),   dimension( its:ite, kts:kte )  ::                           delprsi, &
                                                                          pdh
  real(r8),   dimension( its:ite, kts:kte+1 )   ::     pdhi
  real(r8),   dimension( its:ite, nvar_dirOA )  ::     oa4
  real(r8),   dimension( its:ite, nvar_dirOL )  ::     ol4
  integer ::  i,j,k,kdt,kpblmax
!
        !===========Jinbo Xie=========
        integer , intent(in) :: gwd_ls,gwd_bl,gwd_ss,gwd_fd!Jinbo Xie 
        !===========Jinbo Xie=========
   do k = kts,kte
     if(znu(k).gt.0.6_r8) kpblmax = k + 1
   enddo
!
      do k = kts,kte+1
         do i = its,ite
            if(k.le.kte)pdh(i,k) = p3d(i,k)
             pdhi(i,k) = p3di(i,k)
         enddo
      enddo
!
      do k = kts,kte
        do i = its,ite
          delprsi(i,k) = pdhi(i,k)-pdhi(i,k+1)
        enddo
      enddo
!
!=======Jinbo Xie=================
!no need when there is no large drag
IF ( (gwd_ls .EQ. 1).and.(gwd_bl .EQ. 1)) then

        do i = its,ite
            oa4(i,:) = oa2d(i,:)
            ol4(i,:) = ol2d(i,:)
        enddo
ENDIF
!========Jinbo Xie================
!=================================================================
        ! Jinbo Xie1  changed all ,j out, turn 3d into 2d, for cam formation
              !=================================================================
      call gwdo2d(dudt=rublten(ims,kms),dvdt=rvblten(ims,kms)                  &
             ,dthdt=rthblten(ims,kms)                                          &
              ,dtaux2d_ls=dtaux2d_ls,dtauy2d_ls=dtauy2d_ls                     &
              ,dtaux2d_bl=dtaux2d_bl,dtauy2d_bl=dtauy2d_bl                     &
              ,dtaux2d_ss=dtaux2d_ss,dtauy2d_ss=dtauy2d_ss                     &
              ,dtaux2d_fd=dtaux2d_fd,dtauy2d_fd=dtauy2d_fd                     &
              ,u1=u3d(ims,kms),v1=v3d(ims,kms)                                 &
              ,t1=t3d(ims,kms)                                                 &
              ,q1=qv3d(ims,kms)                                                &
              ,del=delprsi(its,kts)                                            &
              ,prsi=pdhi(its,kts)                                              &
              ,prsl=pdh(its,kts),prslk=pi3d(ims,kms)                           &
              ,zl=z(ims,kms),rcl=1.0_r8                                        &
              ,xland1=xland(ims),br1=br(ims),hpbl=pblh(ims)                    &
              ,bnv_in=bnvbg(ims,kms)                                           &
              ,dz2=dz(ims,kms)                                                 &
              ,kpblmax=kpblmax                                                 &
              ,dusfc_ls=dusfc_ls,dvsfc_ls=dvsfc_ls                             &
              ,dusfc_bl=dusfc_bl,dvsfc_bl=dvsfc_bl                             &
              ,dusfc_ss=dusfc_ss,dvsfc_ss=dvsfc_ss                             &
              ,dusfc_fd=dusfc_fd,dvsfc_fd=dvsfc_fd                             &
              ,var=var2d(ims),oc1=oc12d(ims)                                   &
              ,oa4=oa4,ol4=ol4                                                 &
              ,g=g,cp=cp,rd=rd,rv=rv,fv=ep1,pi=pi                              &
              ,dxmeter=dx,dymeter=dy,deltim=dt                                 &
              ,kpbl=kpbl2d(ims),kdt=itimestep,lat=j                            &
              ,ids=ids,ide=ide, jds=jds,jde=jde, kds=kds,kde=kde               &
              ,ims=ims,ime=ime, jms=jms,jme=jme, kms=kms,kme=kme               &
              ,its=its,ite=ite, jts=jts,jte=jte, kts=kts,kte=kte               &
              ,gsd_gwd_ls=gwd_ls,gsd_gwd_bl=gwd_bl,gsd_gwd_ss=gwd_ss,gsd_gwd_fd=gwd_fd)
                !=============Jinbo Xie==============
                !Jinbo Xie
                !added dymeter in here
                !IF (gwd_opt == 33) then !research mode
!1
IF ( (gwd_ls .EQ. 1).and.(gwd_bl .EQ. 1)) then
                do i = its,ite
                dusfcg_ls(i)=dusfc_ls(i)
                dvsfcg_ls(i)=dvsfc_ls(i)
                dusfcg_bl(i)=dusfc_bl(i)
                dvsfcg_bl(i)=dvsfc_bl(i)
                enddo

             dtaux3d_ls=dtaux2d_ls
             dtaux3d_bl=dtaux2d_bl
             dtauy3d_ls=dtauy2d_ls
             dtauy3d_bl=dtauy2d_bl
!Jinbo Xie temporarily for base flux
!IF (gwd_ss .EQ. 1) then
                do i = its,ite
                dusfcg_ss(i)=dusfc_ss(i)
                dvsfcg_ss(i)=dvsfc_ss(i)
                end do
                    
                dtaux3d_ss=dtaux2d_ss
                dtauy3d_ss=dtauy2d_ss
!ENDIF
ENDIF
IF (gwd_fd .EQ. 1) then

                do i = its,ite
                dusfcg_fd(i)=dusfc_fd(i)
                dvsfcg_fd(i)=dvsfc_fd(i)
                enddo
                dtaux3d_fd=dtaux2d_fd
                dtauy3d_fd=dtauy2d_fd
ENDIF
!#endif
!
   end subroutine gwdo_gsd
!
!-------------------------------------------------------------------------------
!
!-------------------------------------------------------------------------------
   subroutine gwdo2d(dudt,dvdt,dthdt,dtaux2d_ls,dtauy2d_ls,                    &
                    dtaux2d_bl,dtauy2d_bl,dtaux2d_ss,dtauy2d_ss,               &
                    dtaux2d_fd,dtauy2d_fd,u1,v1,t1,q1,                         &
                    del,                                                       &
                    prsi,prsl,prslk,zl,rcl,                                    &
                    xland1,br1,hpbl,bnv_in,dz2,                               &
                    kpblmax,dusfc_ls,dvsfc_ls,dusfc_bl,dvsfc_bl,               &
                    dusfc_ss,dvsfc_ss,dusfc_fd,dvsfc_fd,var,oc1,oa4,ol4,&
                    g,cp,rd,rv,fv,pi,dxmeter,dymeter,deltim,kpbl,kdt,lat,      &
                    ids,ide, jds,jde, kds,kde,                                 &
                    ims,ime, jms,jme, kms,kme,                                 &
                    its,ite, jts,jte, kts,kte,                                 &
                    gsd_gwd_ls,gsd_gwd_bl,gsd_gwd_ss,gsd_gwd_fd)!Jinbo Xie 
!===============================
! Jinbo Xie add another dymeter
!===============================
!=====Jinbo Xie=====
use sub_xjb,only:OLgrid,dxygrid
!=====Jinbo Xie=====
!  this code handles the time tendencies of u v due to the effect of mountain 
!  induced gravity wave drag from sub-grid scale orography. this routine 
!  not only treats the traditional upper-level wave breaking due to mountain 
!  variance (alpert 1988), but also the enhanced lower-tropospheric wave 
!  breaking due to mountain convexity and asymmetry (kim and arakawa 1995). 
!  thus, in addition to the terrain height data in a model grid gox, 
!  additional 10-2d topographic statistics files are needed, including 
!  orographic standard  deviation (var), convexity (oc1), asymmetry (oa4) 
!  and ol (ol4). these data sets are prepared based on the 30 sec usgs orography
!  hong (1999). the current scheme was implmented as in hong et al.(2008)
!
!  coded by song-you hong and young-joon kim and implemented by song-you hong
!
!  The Code is further added with flow-blocking (Kim et al. 2005), turbulent orographic form drag
!  (TOFD) (Beljaars et al. 2004), and small-scale drag (Tsiringakis et al.
!  2017).
!  Xie et al.,(2020) recently extended the orographic anisotropy to all flow-directions.
!
!  program history log:
!
!           Activation of each component is done by specifying the integer-parameters
!           (defined below) to 0: inactive or 1: active
!                    gsd_gwd_ls = 0 or 1: large-scale
!                    gsd_gwd_bl = 0 or 1: blocking drag 
!                    gsd_gwd_ss = 0 or 1: small-scale gravity wave drag
!                    gsd_gwd_fd = 0 or 1: topographic form drag
!
!  references:
!        hong et al. (2008), wea. and forecasting
!        kim and doyle (2005), Q. J. R. Meteor. Soc.
!        kim and arakawa (1995), j. atmos. sci.
!        alpet et al. (1988), NWP conference.
!        hong (1999), NCEP office note 424.
!        steeneveld et al (2008), JAMC
!        Tsiringakis et al. (2017), Q. J. R. Meteor. Soc.
!        Xie et al. (2020), JAMES
!
!  notice : comparible or lower resolution orography files than model resolution
!           are desirable in preprocess (wps) to prevent weakening of the drag
!-------------------------------------------------------------------------------
!
!  input                                                                
!        dudt (ims:ime,kms:kme)  non-lin tendency for u wind component
!        dvdt (ims:ime,kms:kme)  non-lin tendency for v wind component
!        u1(ims:ime,kms:kme) zonal wind / sqrt(rcl)  m/sec  at t0-dt
!        v1(ims:ime,kms:kme) meridional wind / sqrt(rcl) m/sec at t0-dt
!        t1(ims:ime,kms:kme) temperature deg k at t0-dt
!        q1(ims:ime,kms:kme) specific humidity at t0-dt
!
!        rcl     a scaling factor = reciprocal of square of cos(lat)
!                for gmp.  rcl=1 if u1 and v1 are wind components.
!        deltim  time step    secs                                       
!        del(kts:kte)  positive increment of pressure across layer (pa)
!                                                                       
!  output
!        dudt, dvdt    wind tendency due to gwdo
!
!-------------------------------------------------------------------------------
   implicit none
!-------------------------------------------------------------------------------
   integer              ::  kdt,lat,latd,lond,kpblmax,                         &
                            ids,ide, jds,jde, kds,kde,                         &
                            ims,ime, jms,jme, kms,kme,                         &
                            its,ite, jts,jte, kts,kte
!
   real(r8)                 ::  g,rd,rv,fv,cp,pi,deltim,rcl!dxmeter,deltim,rcl
!=======================
!!!!Jinbo Xie, add dymeter
        real(r8),dimension(:)    :: dymeter
        real(r8),dimension(:)    :: dxmeter
!======================
   real(r8)                ::  dudt(ims:ime,kms:kme),dvdt(ims:ime,kms:kme),          &
                                dthdt(ims:ime,kms:kme),&
                            dtaux2d_ls(ims:ime,kms:kme),dtauy2d_ls(ims:ime,kms:kme), &
                            dtaux2d_bl(ims:ime,kms:kme),dtauy2d_bl(ims:ime,kms:kme), &
                            dtaux2d_ss(ims:ime,kms:kme),dtauy2d_ss(ims:ime,kms:kme), &
                            dtaux2d_fd(ims:ime,kms:kme),dtauy2d_fd(ims:ime,kms:kme), &
                            u1(ims:ime,kms:kme),v1(ims:ime,kms:kme),           &
                            t1(ims:ime,kms:kme),q1(ims:ime,kms:kme),           &
                            zl(ims:ime,kms:kme),prsl(its:ite,kts:kte),         &
                            prslk(ims:ime,kms:kme)
   real(r8),intent(in)                ::  prsi(its:ite,kts:kte+1),del(its:ite,kts:kte)
!=========Jinbo Xie=======
   real(r8),intent(in),optional    ::  oa4(its:ite,nvar_dirOA)
   real(r8),intent(in),optional    ::  ol4(its:ite,nvar_dirOL)
!=========Jinbo Xie=======
!
! GSD surface drag options to regulate specific components
! Each component is tapered off automatically as a function of dx, so best to 
! keep them activated (=1).
!integer, parameter ::                                                          &
!   gsd_gwd_ls      = 1,       & ! large-scale gravity wave drag
!   gsd_gwd_bl      = 1,       & ! blocking drag 
!   gsd_gwd_ss      = 1,       & ! small-scale gravity wave drag (Steeneveld et al. 2008)
!   gsd_gwd_fd      = 0,       & ! form drag (Beljaars et al. 2004, QJRMS)

!
! added for small-scale orographic wave drag
   real(r8), dimension(its:ite,kts:kte)     :: utendwave,vtendwave,thx,thvx,za
   real(r8), dimension(ims:ime), intent(in) :: br1,hpbl,xland1
   real(r8), dimension(its:ite)             :: govrth
   real(r8), dimension(ims:ime,kms:kme), intent(in) :: dz2
   real(r8), dimension(its:ite,kts:kte+1)   :: zq
   real(r8)                 :: tauwavex0,tauwavey0,XNBV,density,tvcon,hpbl2
   integer              :: kpbl2,kvar
   real(r8), parameter      :: varmax = 200._r8
!
   integer              ::  kpbl(ims:ime)
   real(r8)                 ::  var(ims:ime),oc1(ims:ime),                         &
                            dusfc_ls(ims:ime),dvsfc_ls(ims:ime),               &
                            dusfc_bl(ims:ime),dvsfc_bl(ims:ime),               &
                            dusfc_ss(ims:ime),dvsfc_ss(ims:ime),               &
                            dusfc_fd(ims:ime),dvsfc_fd(ims:ime)
! Variables for scale-awareness:
! Small-scale GWD + turbulent form drag
   real(r8), parameter   :: dxmin_ss = 1000._r8, dxmax_ss = 12000._r8  ! min,max range of tapering (m)
! Large-scale GWD
   real(r8), parameter   :: dxmin_ls = 3000._r8, dxmax_ls = 13000._r8  ! min,max range of tapering (m)
!===========================
!!!!!Jinbo Xie
!!!!!Add y axis for taper consider
   real(r8), parameter   :: dymin_ls = 3000._r8, dymax_ls = 13000._r8  ! min,maxrange of tapering (m)
   real(r8), parameter   :: dymin_ss = 3000._r8, dymax_ss = 13000._r8  ! min,maxrange of tapering (m)
!==========================
   real(r8)              :: ss_taper, ls_taper  ! small- and large-scale tapering factors (-)
!
! added Beljaars orographic form drag
   real(r8), dimension(its:ite,kts:kte)     :: utendform,vtendform
   real(r8)                 :: a1,a2,wsp
! critical richardson number for wave breaking : ! larger drag with larger value
!!==========Jinbo Xie=========
   !real(r8),parameter       ::  ric     = 0.25_r8 
real(r8),parameter       ::  ric     = 1._r8
real(r8),parameter       ::  ric_rig  = 0.25_r8
!==========Jinbo Xie==========

   real(r8),parameter       ::  dw2min  = 1._r8
   real(r8),parameter       ::  rimin   = -100._r8
   real(r8),parameter       ::  bnv2min = 1.0e-5_r8
   real(r8),parameter       ::  efmin   = 0.0_r8
   real(r8),parameter       ::  efmax   = 10.0_r8
   real(r8),parameter       ::  xl      = 4.0e4_r8
   real(r8),parameter       ::  critac  = 1.0e-5_r8
   real(r8),parameter       ::  gmax    = 1._r8
   real(r8),parameter       ::  veleps  = 1.0_r8
   real(r8),parameter       ::  factop  = 0.5_r8
   real(r8),parameter       ::  frc     = 1.0_r8
   real(r8),parameter       ::  ce      = 0.8_r8
   real(r8),parameter       ::  cg      = 0.5_r8
   integer,parameter    ::  kpblmin = 2
!
!  local variables
!
   integer              ::  j,i,k,lcap,lcapp1,nwd,idir,                          &
                            klcap,kp1,ikount,kk,nwd1!added nwd1 Jinbo Xie
!
   real(r8)                 ::  rcs,rclcs,csg,fdir,cleff,cs,rcsks,                 &
                            wdir,ti,rdz,temp,tem2,dw2,shr2,bvf2,rdelks,        &
                            wtkbj,tem,gfobnv,hd,fro,rim,temc,tem1,efact,       &
                            temv,dtaux,dtauy,eng0,eng1,theta,rad,wdir1!Jinbo Xie added theta,rad,wdir1
!=====Jinbo Xie=====
real(r8),dimension(its:ite,kts:kte),intent(in), optional :: bnv_in
!=====Jinbo Xie=====
!
   logical              ::  ldrag(its:ite),icrilv(its:ite),                    &
                            flag(its:ite),kloop1(its:ite)
!                                                                       
   real(r8)                 ::  taub(its:ite),taup(its:ite,kts:kte+1),             &
                            xn(its:ite),yn(its:ite),                           &
                            ubar(its:ite),vbar(its:ite),                       &
                            fr(its:ite),ulow(its:ite),                         &
                            rulow(its:ite),bnv(its:ite),                       &
                            oa1(its:ite),ol(its:ite),                          &
                            roll(its:ite),dtfac(its:ite),                      &
                            brvf(its:ite),xlinv(its:ite),                      &
                            delks(its:ite),delks1(its:ite),                    &
                            bnv2(its:ite,kts:kte),usqj(its:ite,kts:kte),       &
                            taud_ls(its:ite,kts:kte),taud_bl(its:ite,kts:kte), &
                            ro(its:ite,kts:kte),                               &
                            vtk(its:ite,kts:kte),vtj(its:ite,kts:kte),         &
                            zlowtop(its:ite),velco(its:ite,kts:kte-1),         &
                            coefm(its:ite)
!
   integer              ::  kbl(its:ite),klowtop(its:ite)
!
   logical :: iope
   integer,parameter    ::  mdir=2*nvar_dirOL!the number of directions
!integer              ::  nwdir(mdir)
!data nwdir/6,7,5,8,2,3,1,4/
!  variables for flow-blocking drag
!
   real(r8),parameter       :: frmax  = 10._r8
   real(r8),parameter       :: olmin  = 1.0e-5_r8
   real(r8),parameter       :: odmin  = 0.1_r8
   real(r8),parameter       :: odmax  = 10._r8
   real(r8),parameter       :: erad   = 6371.315e+3_r8
   integer              :: komax(its:ite)
   integer              :: kblk
   real(r8)                 :: cd
   real(r8)                 :: zblk,tautem

!Jinbo Xie
real(r8) :: zblk_col(its:ite)
real(r8) :: taub_xjb(its:ite)
real(r8) :: taufb_xjb(its:ite)
real(r8) :: wdir1_xjb(its:ite)
integer  :: kblk_xjb(its:ite)

!Jinbo Xie
real(r8)                 :: pe,ke
!================================
   real(r8)                 :: dely(its:ite),dxy4(its:ite,nvar_dirOL),&
                     delx(its:ite),dxy4p(its:ite,nvar_dirOL)
!=====Jinbo Xie==================
   real(r8)                 :: dxy(its:ite),dxyp(its:ite)
   real(r8)                 ::  olp(its:ite),&
                                 od(its:ite)
   real(r8)                 :: taufb(its:ite,kts:kte+1)
        !===========Jinbo Xie=========
        integer , intent(in) :: gsd_gwd_ls,gsd_gwd_bl,gsd_gwd_ss,gsd_gwd_fd        !Jinbo Xie 
        !===========Jinbo Xie=========


!=====Jinbo Xie=====
integer :: wdir_add_xjb(mdir,its:ite)
real(r8):: ncleff  !!tunable parameter
real(r8):: ncd
!=====Jinbo Xie=====


 !===================================
 !Jinbo Xie readdata
 !real(r8),allocatable :: need3(:,:)
 !real(r8) :: oc11(ims:ime)
 !real(r8) :: wind_xjb(kts:kte)
 !real(r8) :: shr2_xjb(its:ite,kts:kte)
 real(r8) :: l1,l2,S!,shrrok1,shrrok0,gamma1
 !
 logical  :: iint
 real(r8) :: zl_hint(its:ite)
 !===================================
!
!---- constants                                                         
!                                                                       
   rcs    = sqrt(rcl)
   cs     = 1._r8 / sqrt(rcl)
   csg    = cs * g
   lcap   = kte
   lcapp1 = lcap + 1
   fdir   = mdir / (2.0_r8*pi)
!
!--- calculate scale-aware tapering factors
!
!=========================================
!=========================================
!!Jinbo Xie  add criteria for dymeter
!Taper for small GWD only, currently assumes equal length in both direction
!Taper matters not much
#if 0
if ( dxmeter .ge. dxmax_ls .and. dymeter .ge. dymax_ls) then
!=========================================
   ls_taper = 1.
else
   if ( dxmeter .le. dxmin_ls) then
      ls_taper = 0.
   else
      ls_taper = 0.5 * ( SIN(pi*(dxmeter-0.5*(dxmax_ls+dxmin_ls))/    &
                                (dxmax_ls-dxmin_ls)) + 1. )
   end if
end if
if ( dxmeter .ge. dxmax_ss ) then
   ss_taper = 1.
else
   if ( dxmeter .le. dxmin_ss) then
      ss_taper = 0.
   else
      ss_taper = dxmax_ss * (1. - dxmin_ss/dxmeter)/(dxmax_ss-dxmin_ss)
   end if
end if
#endif
!Jinbo Xie currently use smoothed topo, taper sets to none-taper
!Jinbo Xie maybe only used when using directly derived data from 30s
ls_taper=1._r8
ss_taper=1._r8
!
!--- calculate length of grid for flow-blocking drag
!
   delx=dxmeter
!============
! Jinbo Xie2
!============
   dely=dymeter !Jinbo Xie, add dy, since global model dx/=dy
!============
! Jinbo Xie2
!============
!Jinbo Xie 
!varied delx,so everything needs add another dim
!
!
!-----initialize arrays
   dtaux = 0.0_r8
   dtauy = 0.0_r8
   do i = its,ite
     klowtop(i)    = 0
     kbl(i)        = 0
   enddo
!
   do i = its,ite
     xn(i)         = 0.0_r8
     yn(i)         = 0.0_r8
     ubar (i)      = 0.0_r8
     vbar (i)      = 0.0_r8
     roll (i)      = 0.0_r8
     taub (i)      = 0.0_r8
     oa1(i)        = 0.0_r8
     ol(i)         = 0.0_r8
     ulow (i)      = 0.0_r8
     dtfac(i)      = 1.0_r8
     ldrag(i)      = .false.
     icrilv(i)     = .false.
     flag(i)       = .true.
   enddo
!

   do k = kts,kte
     do i = its,ite
       usqj(i,k) = 0.0_r8
       bnv2(i,k) = 0.0_r8
       vtj(i,k)  = 0.0_r8
       vtk(i,k)  = 0.0_r8
       taup(i,k) = 0.0_r8
       taud_ls(i,k) = 0.0_r8
       taud_bl(i,k) = 0.0_r8
       dtaux2d_ls(i,k)= 0.0_r8
       dtauy2d_ls(i,k)= 0.0_r8
       dtaux2d_bl(i,k)= 0.0_r8
       dtauy2d_bl(i,k)= 0.0_r8
       dtaux2d_ss(i,k)= 0.0_r8
       dtauy2d_ss(i,k)= 0.0_r8
       dtaux2d_fd(i,k)= 0.0_r8
       dtauy2d_fd(i,k)= 0.0_r8
     enddo
   enddo
!
   do i = its,ite
     dusfc_ls(i) = 0.0_r8
     dvsfc_ls(i) = 0.0_r8
     dusfc_bl(i) = 0.0_r8
     dvsfc_bl(i) = 0.0_r8
     dusfc_ss(i) = 0.0_r8
     dvsfc_ss(i) = 0.0_r8
     dusfc_fd(i) = 0.0_r8
     dvsfc_fd(i) = 0.0_r8

!temp for base flux xjb
taub_xjb(i)=0.0_r8
taufb_xjb(i)=0.0_r8
zblk_col(i)=0.0_r8
wdir1_xjb(i)=0.0_r8
kblk_xjb(i)=0.0_r8
   enddo

!
   do i = its,ite
     taup(i,kte+1) = 0.0_r8
     xlinv(i)     = 1.0_r8/xl
   enddo
!
!  initialize array for flow-blocking drag
!
   taufb(its:ite,kts:kte+1) = 0.0_r8
   komax(its:ite) = 0
!
   do k = kts,kte
     do i = its,ite
       vtj(i,k)  = t1(i,k)  * (1._r8+fv*q1(i,k))
       vtk(i,k)  = vtj(i,k) / prslk(i,k)
       ro(i,k)   = 1._r8/rd * prsl(i,k) / vtj(i,k) ! density kg/m**3
     enddo
   enddo
!
!  determine reference level: maximum of 2*var and pbl heights
!
   do i = its,ite
     zlowtop(i) = 2._r8 * var(i)
   enddo
!
   do i = its,ite
     kloop1(i) = .true.
   enddo
!
   do k = kts+1,kte
     do i = its,ite
         if(kloop1(i).and.zl(i,k)-zl(i,1).ge.zlowtop(i)) then
         klowtop(i) = k+1
         kloop1(i)  = .false.
         endif
     enddo
   enddo
!
   do i = its,ite
     kbl(i)   = max(kpbl(i), klowtop(i))
     kbl(i)   = max(min(kbl(i),kpblmax),kpblmin)
   enddo
!
!  determine the level of maximum orographic height
!
   komax(:) = kbl(:)
!
   do i = its,ite
     delks(i)  = 1.0_r8 / (prsi(i,1) - prsi(i,kbl(i)))
     delks1(i) = 1.0_r8 / (prsl(i,1) - prsl(i,kbl(i)))
   enddo
!
!  compute low level averages within pbl
!
   do k = kts,kpblmax
     do i = its,ite
       if (k.lt.kbl(i)) then
         rcsks   = rcs     * del(i,k) * delks(i)
         rdelks  = del(i,k)  * delks(i)
         ubar(i) = ubar(i) + rcsks  * u1(i,k)      ! pbl u  mean
         vbar(i) = vbar(i) + rcsks  * v1(i,k)      ! pbl v  mean
         roll(i) = roll(i) + rdelks * ro(i,k)      ! ro mean
       endif
     enddo
   enddo
!
!=======Jinbo Xie=======
!For ls and bl only
IF  ((gsd_gwd_ls .EQ. 1).or.(gsd_gwd_bl .EQ. 1)) then
!     figure out low-level horizontal wind direction 
!====Jinbo Xie order into a counterclockwise index instead====
!no more 1-8 index
   do i = its,ite
     wdir   = atan2(vbar(i),ubar(i)) + pi!changed into y/x Jinbo Xie
     !idir   = MOD(nint(fdir*wdir),mdir) + 1!starts from pi already
      !nwd   = idir
        wdir1=wdir-pi
        if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
        nwd  = MOD(nint(fdir*wdir1),mdir) + 1
        else!(-pi,0)
        nwd  = MOD(nint(fdir*(wdir1+2._r8*pi)),mdir) + 1
        endif
        !turn backwords because start is pi
        !need turning
        rad=4.0_r8*atan(1.0_r8)/180.0_r8
        theta=(real(nwd,kind=r8)-1._r8)*(360._r8/real(mdir,kind=r8))
        !
wdir1_xjb(i)=wdir1/rad
        !
	oa1(i)= oa4(i,1)*cos(theta*rad)+oa4(i,2)*sin(theta*rad)
        !select OL
        ol(i)  = ol4(i,MOD(nwd-1,int(mdir/2))+1)
        !calculate dxygrid, not so slow
        call dxygrid(dxmeter(i),dymeter(i),theta,dxy(i))
	!
	!----- compute orographic width along (ol) and perpendicular (olp)
!----- the direction of wind
!
!====Jinbo Xie====
!put wdir inside the (0,2*pi) section
!changing pi/2 either way is perpendicular
                !wdir1=wdir-pi
                if (wdir1.ge.0._r8.and.wdir1.lt.pi) then
                nwd1  = MOD(nint(fdir*(wdir1+pi/2._r8)),mdir) + 1
                olp(i)=ol4(i,MOD(nwd1-1,int(mdir/2))+1)
                else!(-pi,0)
                nwd1  = MOD(nint(fdir*(wdir1-pi/2._r8+2._r8*pi)),mdir) + 1
                olp(i)=ol4(i,MOD(nwd1-1,int(mdir/2))+1)
                endif
                theta=(real(nwd1,kind=r8)-1._r8)*(360._r8/real(mdir,kind=r8))
                call dxygrid(dxmeter(i),dymeter(i),theta,dxyp(i))
		!
		!====Jinbo Xie====
	!
!----- compute orographic direction (horizontal orographic aspect ratio)
!
     od(i) = olp(i)/max(ol(i),olmin)
     od(i) = min(od(i),odmax)
     od(i) = max(od(i),odmin)
!
!----- compute length of grid in the along(dxy) and cross(dxyp) wind directions
!
!==========================================
!Jinbo Xie
   enddo
!===Jinbo Xie===
ENDIF
!===Jinbo Xie===

!Jinbo Xie Since variable grid,change dxy4 to larger
!==========================================
!
! END INITIALIZATION; BEGIN GWD CALCULATIONS:
!
IF ( ((gsd_gwd_ls .EQ. 1).or.(gsd_gwd_bl .EQ. 1)).and.   &
               (ls_taper .GT. 1.E-02) ) THEN   !====

!                                                                       
!---  saving richardson number in usqj for migwdi                       
!
   do k = kts,kte-1
     do i = its,ite
       ti        = 2.0_r8 / (t1(i,k)+t1(i,k+1))
       rdz       = 1._r8/(zl(i,k+1) - zl(i,k))
       tem1      = u1(i,k) - u1(i,k+1)
       tem2      = v1(i,k) - v1(i,k+1)
       dw2       = rcl*(tem1*tem1 + tem2*tem2)
       shr2      = max(dw2,dw2min) * rdz * rdz
       bvf2      = g*(g/cp+rdz*(vtj(i,k+1)-vtj(i,k))) * ti
       usqj(i,k) = max(bvf2/shr2,rimin)
       !bnv2(i,k) = 2.0_r8*g*rdz*(vtk(i,k+1)-vtk(i,k))/(vtk(i,k+1)+vtk(i,k))
       !bnv2(i,k) = max( bnv2(i,k), bnv2min )
       bnv2(i,k) = max(bnv_in(i,k)**2,bnv2min )
     enddo
   enddo
!
!----compute the "low level" or 1/3 wind magnitude (m/s)                
!                                                                       
   do i = its,ite
     ulow(i) = max(sqrt(ubar(i)*ubar(i) + vbar(i)*vbar(i)), 1.0_r8)
     rulow(i) = 1._r8/ulow(i)
   enddo
!
   do k = kts,kte-1
     do i = its,ite
       velco(i,k)  = (0.5_r8*rcs) * ((u1(i,k)+u1(i,k+1)) * ubar(i)                &
                                   + (v1(i,k)+v1(i,k+1)) * vbar(i))
       velco(i,k)  = velco(i,k) * rulow(i)
       if ((velco(i,k).lt.veleps) .and. (velco(i,k).gt.0._r8)) then
         velco(i,k) = veleps
       endif
     enddo
   enddo
!                                                                       
!  no drag when critical level in the base layer                        
!                                                                       
   do i = its,ite
     ldrag(i) = velco(i,1).le.0._r8
   enddo
!
!  no drag when velco.lt.0                                               
! 
   do k = kpblmin,kpblmax
     do i = its,ite
       if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. velco(i,k).le.0._r8
     enddo
   enddo
!                                                                       
!  no drag when bnv2.lt.0                                               
!                                                                       
   do k = kts,kpblmax
     do i = its,ite
       if (k .lt. kbl(i)) ldrag(i) = ldrag(i).or. bnv2(i,k).lt.0._r8
     enddo
   enddo

!                                                                       
!-----the low level weighted average ri is stored in usqj(1,1; im)      
!-----the low level weighted average n**2 is stored in bnv2(1,1; im)    
!---- this is called bnvl2 in phys_gwd_alpert_sub not bnv2                           
!---- rdelks (del(k)/delks) vert ave factor so we can * instead of /    
!                                                                       
   do i = its,ite
     wtkbj     = (prsl(i,1)-prsl(i,2)) * delks1(i)
     bnv2(i,1) = wtkbj * bnv2(i,1)
     usqj(i,1) = wtkbj * usqj(i,1)
   enddo
!
   do k = kpblmin,kpblmax
     do i = its,ite
       if (k .lt. kbl(i)) then
         rdelks    = (prsl(i,k)-prsl(i,k+1)) * delks1(i)
         bnv2(i,1) = bnv2(i,1) + bnv2(i,k) * rdelks
         usqj(i,1) = usqj(i,1) + usqj(i,k) * rdelks
       endif
     enddo
   enddo
!                                                                       
   do i = its,ite
     ldrag(i) = ldrag(i) .or. bnv2(i,1).le.0.0_r8
     ldrag(i) = ldrag(i) .or. ulow(i).eq.1.0_r8
     ldrag(i) = ldrag(i) .or. var(i) .le. 0.0_r8
   enddo

!                                                                       
!  set all ri low level values to the low level value          
!                                                                       
   do k = kpblmin,kpblmax
     do i = its,ite
       if (k .lt. kbl(i)) usqj(i,k) = usqj(i,1)
     enddo
   enddo
!
   do i = its,ite
     if (.not.ldrag(i))   then
       bnv(i) = sqrt( bnv2(i,1) )
       fr(i) = bnv(i)  * rulow(i) * 2._r8 * var(i) * od(i)
       fr(i) = min(fr(i),frmax)
       xn(i)  = ubar(i) * rulow(i)
       yn(i)  = vbar(i) * rulow(i)
     endif
   enddo

!
!  compute the base level stress and store it in taub
!  calculate enhancement factor, number of mountains & aspect        
!  ratio const. use simplified relationship between standard            
!  deviation & critical hgt                                          
!

   do i = its,ite
     if (.not. ldrag(i))   then
       efact    = (oa1(i) + 2._r8) ** (ce*fr(i)/frc)
       efact    = min( max(efact,efmin), efmax )
!!!!!!! cleff (effective grid length) is highly tunable parameter
!!!!!!! the bigger (smaller) value produce weaker (stronger) wave drag
       cleff    = sqrt(dxy(i)**2._r8 + dxyp(i)**2._r8)
!==============Jinbo Xie=============================================
       !cleff    = 3._r8 * max(dxmeter(i),cleff)!turned dxmeter to array
        ncleff   = 3._r8
        cleff    = 3._r8/ncleff * max(dxmax_ls,cleff)
        !cleff    = 2._r8 * max(dxmax_ls,cleff)
!==============Jinbo Xie=============================================
       coefm(i) = (1._r8 + ol(i)) ** (oa1(i)+1._r8)
       xlinv(i) = coefm(i) / cleff
       tem      = fr(i) * fr(i) * oc1(i)
       gfobnv   = gmax * tem / ((tem + cg)*bnv(i))

       if ( gsd_gwd_ls .NE. 0 ) then
          taub(i)  = xlinv(i) * roll(i) * ulow(i) * ulow(i)                       &
                   * ulow(i) * gfobnv * efact
       else     ! We've gotten what we need for the blocking scheme
          taub(i) = 0.0_r8
       end if
!Jinbo Xie for base flux
taub_xjb(i)=taub(i)
!Jinbo Xie for base flux
     else
       taub(i) = 0.0_r8
       xn(i)   = 0.0_r8
       yn(i)   = 0.0_r8
     endif
   enddo

ENDIF   ! (gsd_gwd_ls .EQ. 1).or.(gsd_gwd_bl .EQ. 1)
!=========================================================
! add small-scale wavedrag for stable boundary layer
!=========================================================
  XNBV=0._r8
  tauwavex0=0._r8
  tauwavey0=0._r8
  density=1.2_r8
  utendwave=0._r8
  vtendwave=0._r8
  zq=0._r8
!
  IF ( (gsd_gwd_ss .EQ. 1).and.(ss_taper.GT.1.E-02) ) THEN
!
! declaring potential temperature
!
    do k = kts,kte
      do i = its,ite
        thx(i,k) = t1(i,k)/prslk(i,k)
      enddo
    enddo
!
    do k = kts,kte
      do i = its,ite
        tvcon = (1._r8+fv*q1(i,k))
        thvx(i,k) = thx(i,k)*tvcon
      enddo
    enddo
!
! Defining layer height
!
    do k = kts,kte
      do i = its,ite
        zq(i,k+1) = dz2(i,k)+zq(i,k)
      enddo
    enddo
!
    do k = kts,kte
      do i = its,ite
        za(i,k) = 0.5_r8*(zq(i,k)+zq(i,k+1))
      enddo
    enddo

    do i=its,ite
       hpbl2 = hpbl(i)+10._r8
       kpbl2 = kpbl(i)
       kvar = 1
       do k=kts+1,MAX(kpbl(i),kts+1)
          IF (za(i,k)>300._r8) then
             kpbl2 = k
             IF (k == kpbl(i)) then
                hpbl2 = hpbl(i)+10._r8
             ELSE
                hpbl2 = za(i,k)+10._r8
             ENDIF
             exit
          ENDIF
       enddo

        if(xland1(i).gt.0._r8 .and. 2._r8*var(i).le.hpbl(i))then
          if(br1(i).gt.0._r8 .and. thvx(i,kpbl2)-thvx(i,kts) > 0._r8)then
            cleff    = sqrt(dxy(i)**2_r8 + dxyp(i)**2_r8)
            cleff    = 2.0_r8 * max(dxmax_ss,cleff)
            coefm(i) = (1._r8 + ol(i)) ** (oa1(i)+1._r8)
            xlinv(i) = coefm(i) / cleff
            govrth(i)=g/(0.5_r8*(thvx(i,kpbl2)+thvx(i,kts)))
            XNBV=sqrt(govrth(i)*(thvx(i,kpbl2)-thvx(i,kts))/hpbl2)
!
            if(abs(XNBV/u1(i,kpbl2)).gt.xlinv(i))then
              tauwavex0=0.5_r8*XNBV*xlinv(i)*(2._r8*MIN(var(i),varmax))**2_r8*ro(i,kvar)*u1(i,kvar)
              tauwavex0=tauwavex0*ss_taper   ! "Scale-awareness"
            else
              tauwavex0=0._r8
            endif
!
            if(abs(XNBV/v1(i,kpbl2)).gt.xlinv(i))then
              tauwavey0=0.5_r8*XNBV*xlinv(i)*(2._r8*MIN(var(i),varmax))**2._r8*ro(i,kvar)*v1(i,kvar)
              tauwavey0=tauwavey0*ss_taper   ! "Scale-awareness"
            else
              tauwavey0=0._r8
            endif
!

            do k=kts,kpbl(i) !MIN(kpbl2+1,kte-1)
              utendwave(i,k)=-1._r8*tauwavex0*2._r8*max((1._r8-za(i,k)/hpbl2),0._r8)/hpbl2
              vtendwave(i,k)=-1._r8*tauwavey0*2._r8*max((1._r8-za(i,k)/hpbl2),0._r8)/hpbl2
            enddo
          endif
       endif
    enddo ! end i loop

    do k = kts,kte
       do i = its,ite
         dudt(i,k)  = dudt(i,k) + utendwave(i,k)
         dvdt(i,k)  = dvdt(i,k) + vtendwave(i,k)
         dtaux2d_ss(i,k) = utendwave(i,k)
         dtauy2d_ss(i,k) = vtendwave(i,k)
         dusfc_ss(i) = dusfc_ss(i) + utendwave(i,k) * del(i,k)
         dvsfc_ss(i) = dvsfc_ss(i) + vtendwave(i,k) * del(i,k)
       enddo
    enddo

ENDIF  ! end if gsd_gwd_ss == 1
!================================================================
!add Beljaars et al. (2004, QJRMS, equ. 16) form drag:
!================================================================
IF ( (gsd_gwd_fd .EQ. 1).and.(ss_taper.GT.1.E-02) ) THEN

   utendform=0._r8
   vtendform=0._r8
   zq=0._r8

   IF ( (gsd_gwd_ss .NE. 1).and.(ss_taper.GT.1.E-02) ) THEN
      ! Defining layer height. This is already done above is small-scale GWD is used
      do k = kts,kte
        do i = its,ite
          zq(i,k+1) = dz2(i,k)+zq(i,k)
        enddo
      enddo

      do k = kts,kte
        do i = its,ite
          za(i,k) = 0.5_r8*(zq(i,k)+zq(i,k+1))
        enddo
      enddo
   ENDIF

   DO i=its,ite
      !IF ((xland1(i)-1.5) .le. 0.) then
       !IF (xland1(i) .gt. 0..and.2._r8*var(i).gt.0) then
        IF (xland1(i) .gt. 0.) then
          a1=0.00026615161_r8*var(i)**2_r8
          a2=a1*0.005363_r8
         DO k=kts,kte
            wsp=SQRT(u1(i,k)**2_r8 + v1(i,k)**2_r8)
            ! alpha*beta*Cmd*Ccorr*2.109 = 12.*1.*0.005*0.6*2.109 = 0.0759 
            utendform(i,k)=-0.0759_r8*wsp*u1(i,k)* &
                           EXP(-(za(i,k)/1500._r8)**1.5_r8)*a2*za(i,k)**(-1.2_r8)*ss_taper
            vtendform(i,k)=-0.0759_r8*wsp*v1(i,k)* &
                           EXP(-(za(i,k)/1500._r8)**1.5_r8)*a2*za(i,k)**(-1.2_r8)*ss_taper
            !IF(za(i,k) > 4000.) exit
            !!
            !write(iulog,*) "Jinbo Xie var(i),utendform(i,k),vtendform(i,k)",var(i),utendform(i,k),vtendform(i,k)
            !!
         ENDDO
      ENDIF
   ENDDO
   !!
   do k = kts,kte
      do i = its,ite
         dudt(i,k)  = dudt(i,k) + utendform(i,k)
         dvdt(i,k)  = dvdt(i,k) + vtendform(i,k)
         dtaux2d_fd(i,k) = utendform(i,k)
         dtauy2d_fd(i,k) = vtendform(i,k)
         dusfc_fd(i) = dusfc_fd(i) + utendform(i,k) * del(i,k)
         dvsfc_fd(i) = dvsfc_fd(i) + vtendform(i,k) * del(i,k)
      enddo
   enddo
   ENDIF  ! end if gsd_gwd_fd == 1
!=======================================================
! More for the large-scale gwd component
IF ( (gsd_gwd_ls .EQ. 1).and.(ls_taper.GT.1.E-02) ) THEN
!                                                                       
!   now compute vertical structure of the stress.
!
   do k = kts,kpblmax
      do i = its,ite
         if (k .le. kbl(i)) taup(i,k) = taub(i)
      enddo
   enddo
!
!================================
!Jinbo Xie
!determination of the interface height
do i=its,ite
iint=.false.
        do k=kpblmin,kte-1
        if (k.gt.kbl(i).and.usqj(1,k)-usqj(1,k-1).lt.0.and.(.not.iint)) then
        iint=.true.
        zl_hint(i)=zl(i,k+1)
        endif
        enddo
enddo
!print*,"zl_hint",zl_hint
!!stop
!Jinbo Xie
!================================
   do k = kpblmin, kte-1                   ! vertical level k loop!
      kp1 = k + 1
      do i = its,ite
!
!   unstablelayer if ri < ric
!   unstable layer if upper air vel comp along surf vel <=0 (crit lay)
!   at (u-c)=0. crit layer exists and bit vector should be set (.le.)
!
         if (k .ge. kbl(i)) then
           !icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric)                  &
           !                      .or. (velco(i,k) .le. 0.0_r8)
!============================
!Jinbo Xie
!we modify the criteria for unstable layer
!that the lv is critical under 0.25
!while we keep wave breaking ric for
!other larger lv
           icrilv(i) = icrilv(i) .or. ( usqj(i,k) .lt. ric_rig)&
                                 .or. (velco(i,k) .le. 0.0_r8)
!Jinbo Xie
!============================
           brvf(i)  = max(bnv2(i,k),bnv2min) ! brunt-vaisala frequency squared
           brvf(i)  = sqrt(brvf(i))          ! brunt-vaisala frequency
         endif
      enddo
!
      do i = its,ite
        if (k .ge. kbl(i) .and. (.not. ldrag(i)))   then
          if (.not.icrilv(i) .and. taup(i,k) .gt. 0.0_r8 ) then
            temv = 1.0_r8 / velco(i,k)
            tem1 = coefm(i)/(dxy(i)/ncleff)*(ro(i,kp1)+ro(i,k))*brvf(i)*velco(i,k)*0.5_r8
            hd   = sqrt(taup(i,k) / tem1)
            fro  = brvf(i) * hd * temv

!
!  rim is the minimum-richardson number by shutts (1985)
!
            tem2   = sqrt(usqj(i,k))
            tem    = 1._r8 + tem2 * fro
            rim    = usqj(i,k) * (1._r8-fro) / (tem * tem)

!
!  check stability to employ the 'saturation hypothesis'
!  of lindzen (1981) except at tropospheric downstream regions
!
            if (rim .le. ric) then  ! saturation hypothesis!
              if ((oa1(i) .le. 0._r8).or.(kp1 .ge. kpblmin )) then
                temc = 2.0_r8 + 1.0_r8 / tem2
                hd   = velco(i,k) * (2.0_r8*sqrt(temc)-temc) / brvf(i)
                taup(i,kp1) = tem1 * hd * hd
!==============================================
!taup is restricted to monotoncally decrease
!to avoid unexpected high taup with taup cal
taup(i,kp1)=min(tem1*hd*hd,taup(i,k))
!add vertical decrease at low level below hint (Kim and Doyle 2005)
!where Ri first decreases
!#if 0
if (k.gt.klowtop(i).and.zl(i,k).le.zl_hint(i)) then
l1=(9.81_r8*bnv2(i,kp1)/velco(i,kp1)**2)!-(shr2_xjb(i,kp1)/velco(i,kp1))
l2=(9.81_r8*bnv2(i,k)/velco(i,k)**2)!-(shr2_xjb(i,k)/velco(i,k))
!print*,"l1,l2,l1/l2",l1,l2,l1/l2
taup(i,kp1)=min(taup(i,k),taup(i,k)*(l1/l2),tem1*hd*hd)
!taup(i,kp1)=max(0.2*taup(i,k),min(taup(i,k),taup(i,k)*(l1/l2),tem1*hd*hd))
!taup(i,k)*(l1/l2)
!print*,"taup(i,kp1)",taup(i,kp1)
!print*,"k",k
endif
!#endif
!==============================================
              endif
            else                    ! no wavebreaking!
              taup(i,kp1) = taup(i,k)
            endif
          endif
        endif
      enddo
   enddo
!


   if(lcap.lt.kte) then
      do klcap = lcapp1,kte
         do i = its,ite
           taup(i,klcap) = prsi(i,klcap) / prsi(i,lcap) * taup(i,lcap)
         enddo
      enddo
   endif

ENDIF !END LARGE-SCALE TAU CALCULATION
!===============================================================
!COMPUTE BLOCKING COMPONENT 
!===============================================================
IF ( (gsd_gwd_bl .EQ. 1) .and. (ls_taper .GT. 1.E-02) ) THEN

   do i = its,ite
      if(.not.ldrag(i)) then
!
!------- determine the height of flow-blocking layer
!
        kblk = 0
        pe = 0.0_r8

        do k = kte, kpblmin, -1
          if(kblk.eq.0 .and. k.le.komax(i)) then
!Jinbo Xie
!flow block appears within the reference level
!compare potential energy and kinetic energy
!divided by g*ro is to turn del(pa) into height
            pe = pe + bnv2(i,k)*(zl(i,komax(i))-zl(i,k))*del(i,k)/g/ro(i,k)
            ke = 0.5_r8*((rcs*u1(i,k))**2._r8+(rcs*v1(i,k))**2._r8)
!
!---------- apply flow-blocking drag when pe >= ke 
!
            if(pe.ge.ke) then
              kblk = k
              kblk = min(kblk,kbl(i))
              zblk = zl(i,kblk)-zl(i,kts)
zblk_col(i)=zblk
!write(iulog,*)"Jinbo Xie kblk,k",kblk,k,zblk
            endif
          endif
        enddo
        if(kblk.ne.0) then
!
!--------- compute flow-blocking stress
!

!Jinbo Xie the max(dxmax_ls,dxy(i))**2
!Jinbo Xie here is a crude estimate since the cam is uneven 0.9*1.25deg
!dxmax_ls is different than the usual one
!because the taper is very different
!Jinbo Xie dxy is a length scale mostly in the direction of the flow to the ridge
!so it is good and not needed for an uneven grid area
!ref Lott and Miller (1997) original scheme
          cd = max(2.0_r8-1.0_r8/od(i),0.0_r8)

ncd=3._r8
cd=ncd*cd

          taufb(i,kts) = 0.5_r8 * roll(i) * coefm(i) / max(dxmax_ls,dxy(i))**2 * cd * dxyp(i)   &
                         * olp(i) * zblk * ulow(i)**2
!Jinbo Xie for base flux
taufb_xjb(i)=taufb(i,kts)
!Jinbo Xie for base flux

        !changed grid box area into dy*dy
          tautem = taufb(i,kts)/float(kblk-kts)
          do k = kts+1, kblk
            taufb(i,k) = taufb(i,k-1) - tautem
          enddo

!
!----------sum orographic GW stress and flow-blocking stress
!
           !taup(i,:) = taup(i,:) + taufb(i,:)   ! Keep taup and taufb separate for now
        endif
!!Jinbo Xie
!kblk_xjb(i)=kblk
!!Jinbo Xie
      endif
   enddo

ENDIF   ! end blocking drag
!===========================================================
IF ( (gsd_gwd_ls .EQ. 1 .OR. gsd_gwd_bl .EQ. 1) .and. (ls_taper .GT. 1.E-02) ) THEN

!                                                                       
!  calculate - (g)*d(tau)/d(pressure) and deceleration terms dtaux, dtauy
!
   do k = kts,kte
     do i = its,ite
       taud_ls(i,k) = 1._r8 * (taup(i,k+1) - taup(i,k)) * csg / del(i,k)
       taud_bl(i,k) = 1._r8 * (taufb(i,k+1) - taufb(i,k)) * csg / del(i,k)
     enddo
   enddo
!                                                                       
!  limit de-acceleration (momentum deposition ) at top to 1/2 value 
!  the idea is some stuff must go out the 'top'                     
!                                                                       

   do klcap = lcap,kte
     do i = its,ite
       taud_ls(i,klcap) = taud_ls(i,klcap) * factop
       taud_bl(i,klcap) = taud_bl(i,klcap) * factop
     enddo
   enddo

!                                                                       
!  if the gravity wave drag would force a critical line             
!  in the lower ksmm1 layers during the next deltim timestep,     
!  then only apply drag until that critical line is reached.        
!                                                                       
   do k = kts,kpblmax-1
      do i = its,ite
         if (k .le. kbl(i)) then
           if((taud_ls(i,k)+taud_bl(i,k)).ne.0._r8)                      &
              dtfac(i) = min(dtfac(i),abs(velco(i,k)                     &
                   /(deltim*rcs*(taud_ls(i,k)+taud_bl(i,k)))))
         endif
      enddo
   enddo
!

   do k = kts,kte
      do i = its,ite
         taud_ls(i,k)  = taud_ls(i,k) * dtfac(i) * ls_taper
         taud_bl(i,k)  = taud_bl(i,k) * dtfac(i) * ls_taper
         dtaux2d_ls(i,k) = taud_ls(i,k) * xn(i)
         dtauy2d_ls(i,k) = taud_ls(i,k) * yn(i)
         dtaux2d_bl(i,k) = taud_bl(i,k) * xn(i)
         dtauy2d_bl(i,k) = taud_bl(i,k) * yn(i)
         dudt(i,k)  = dtaux2d_ls(i,k) + dtaux2d_bl(i,k) + dudt(i,k)
         dvdt(i,k)  = dtauy2d_ls(i,k) + dtauy2d_bl(i,k) + dvdt(i,k)
         dusfc_ls(i)  = dusfc_ls(i) + dtaux2d_ls(i,k) * del(i,k)
         dvsfc_ls(i)  = dvsfc_ls(i) + dtauy2d_ls(i,k) * del(i,k)
         dusfc_bl(i)  = dusfc_bl(i) + dtaux2d_bl(i,k) * del(i,k)
         dvsfc_bl(i)  = dvsfc_bl(i) + dtauy2d_bl(i,k) * del(i,k)
      enddo
   enddo

ENDIF


!  Finalize dusfc and dvsfc diagnoses
do i = its,ite
   dusfc_ls(i) = (-1._r8/g*rcs) * dusfc_ls(i)
   dvsfc_ls(i) = (-1._r8/g*rcs) * dvsfc_ls(i)
   dusfc_bl(i) = (-1._r8/g*rcs) * dusfc_bl(i)
   dvsfc_bl(i) = (-1._r8/g*rcs) * dvsfc_bl(i)
   dusfc_ss(i) = (-1._r8/g*rcs) * dusfc_ss(i)
   dvsfc_ss(i) = (-1._r8/g*rcs) * dvsfc_ss(i)
   dusfc_fd(i) = (-1._r8/g*rcs) * dusfc_fd(i)
   dvsfc_fd(i) = (-1._r8/g*rcs) * dvsfc_fd(i)
enddo


!#Jinbo get base flux
#if 0
do i = its,ite
!dusfc_ss(i)=wdir1_xjb(i)
!dvsfc_ss(i)=ol(i)
dtaux2d_ss(i,1)=oa1(i)
dtaux2d_ss(i,2)=ol(i)
dtaux2d_ss(i,3)=olp(i)
dtaux2d_ss(i,4)=dxy(i)
dtaux2d_ss(i,5)=dxyp(i)
dtaux2d_ss(i,6)=var(i)
dtaux2d_ss(i,7)=klowtop(i)!pe_xjb(i)!bnv2(i,2)!vtk(i,k+1)-vtk(i,k)
dtaux2d_ss(i,8)=kpbl(i)
dtaux2d_ss(i,9)=zblk_col(i)
dtaux2d_ss(i,10)=oa4(i,1)
dtaux2d_ss(i,11)=oa4(i,2)
dtaux2d_ss(i,12)=oa4(i,3)
dtaux2d_ss(i,13)=taub_xjb(i)
dtaux2d_ss(i,14)=taufb_xjb(i)
dtaux2d_ss(i,15)=xn(i)
dtaux2d_ss(i,16)=yn(i)
dtaux2d_ss(i,17)=dtfac(i)
dtaux2d_ss(i,18)=ulow(i)
dtaux2d_ss(i,19)=kblk_xjb(i)
dtaux2d_ss(i,20)=br1(i)
enddo
#endif
!stop
#if 0
do i = its,ite
dusfc_ss(i)=taub_xjb(i)
dvsfc_ss(i)=taufb_xjb(i)
enddo
#endif

   return
   end subroutine gwdo2d
!-------------------------------------------------------------------

end module gw_common
