module micro_p3_utils

    use physconst, only: pi, cpair, gravit, rair, rh2o, mwh2o, mwdry, rhoh2o, cpliq, &
                         latvap, latice, tmelt
    use shr_kind_mod,   only: r8=>shr_kind_r8, i8=>shr_kind_i8

    implicit none
    private
    save

    public :: get_latent_heat, get_precip_fraction, micro_p3_utils_init, size_dist_param_liq, &
              size_dist_param_basic,avg_diameter, rising_factorial


    ! 8 byte real and integer
!    integer, parameter, public :: i8 = selected_int_kind(18)

    ! Signaling NaN bit pattern that represents a limiter that's turned off.
    integer(i8), parameter :: limiter_off = int(Z'7FF1111111111111', i8)

    real(r8), public, parameter :: qsmall = 1.e-14
    real(r8), public, parameter :: nsmall = 1.e-16

    real(r8),public :: rhosur,rhosui,ar,br,f1r,f2r,ecr,rhow,kr,kc,bimm,aimm,rin,mi0,nccnst,  &
       eci,eri,bcn,cpw,cons1,cons2,cons3,cons4,cons5,cons6,cons7,         &
       inv_rhow,bsmall,cp,g,rd,rv,ep_2,inv_cp,mw,osm,   &
       vi,epsm,rhoa,map,ma,rr,bact,inv_rm1,inv_rm2,sig1,nanew1,f11,f21,sig2, &
       nanew2,f12,f22,thrd,sxth,piov3,piov6,rho_rimeMin,     &
       rho_rimeMax,inv_rho_rimeMax,max_total_Ni,dbrk,nmltratio,clbfact_sub,  &
       clbfact_dep
    real(r8),dimension(16), public :: dnu

    real(r8),public :: zerodegc  ! Temperature at zero degree celcius ~K
    real(r8),public :: rainfrze  ! Contact and immersion freexing temp, -4C  ~K
    real(r8),public :: homogfrze ! Homogeneous freezing temperature, -40C  ~K
    real(r8),public :: icenuct   ! Ice nucleation temperature, -5C ~K

    real(r8),public,parameter :: pi_e3sm = pi
    ! ice microphysics lookup table array dimensions
    integer, public,parameter :: isize        = 50
    integer, public,parameter :: densize      =  5
    integer, public,parameter :: rimsize      =  4
    integer, public,parameter :: rcollsize    = 30
    integer, public,parameter :: tabsize      = 12  ! number of quantities used from lookup table
    integer, public,parameter :: colltabsize  =  2  ! number of ice-rain collection  quantities used from lookup table
    ! switch for warm-rain parameterization
    ! = 1 Seifert and Beheng 2001
    ! = 2 Beheng 1994
    ! = 3 Khairoutdinov and Kogan 2000
    integer, public,parameter :: iparam = 3

    real(r8), parameter, public :: mincld=0.0001
    real(r8), parameter, public :: rhosn = 250._r8  ! bulk density snow
    real(r8), parameter, public :: rhoi = 500._r8   ! bulk density ice
    real(r8), parameter, public :: rhows = 917._r8  ! bulk density water solid


public :: MGHydrometeorProps

type :: MGHydrometeorProps
   ! Density (kg/m^3)
   real(r8) :: rho
   ! Information for size calculations.
   ! Basic calculation of mean size is:
   !     lambda = (shape_coef*nic/qic)^(1/eff_dim)
   ! Then lambda is constrained by bounds.
   real(r8) :: eff_dim
   real(r8) :: shape_coef
   real(r8) :: lambda_bounds(2)
   ! Minimum average particle mass (kg).
   ! Limit is applied at the beginning of the size distribution calculations.
   real(r8) :: min_mean_mass
end type MGHydrometeorProps

interface MGHydrometeorProps
   module procedure NewMGHydrometeorProps
end interface

type(MGHydrometeorProps), public :: mg_liq_props
type(MGHydrometeorProps), public :: mg_ice_props
type(MGHydrometeorProps), public :: mg_rain_props
type(MGHydrometeorProps), public :: mg_snow_props

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d
! exponent
real(r8), parameter :: dsph = 3._r8

! Bounds for mean diameter for different constituents.
real(r8), parameter :: lam_bnd_rain(2) = 1._r8/[500.e-6_r8, 20.e-6_r8]
real(r8), parameter :: lam_bnd_snow(2) = 1._r8/[2000.e-6_r8, 10.e-6_r8]

! Minimum average mass of particles.
real(r8), parameter :: min_mean_mass_liq = 1.e-20_r8
real(r8), parameter :: min_mean_mass_ice = 1.e-20_r8

!=========================================================
! Utilities that are cheaper if the compiler knows that
! some argument is an integer.
!=========================================================

interface rising_factorial
   module procedure rising_factorial_r8
   module procedure rising_factorial_integer
end interface rising_factorial

interface var_coef
   module procedure var_coef_r8
   module procedure var_coef_integer
end interface var_coef

    contains
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine micro_p3_utils_init()


    real(r8) :: ice_lambda_bounds(2)

    ! mathematical/optimization constants
    thrd  = 1./3.
    sxth  = 1./6.
    piov3 = pi*thrd
    piov6 = pi*sxth

    ! maximum total ice concentration (sum of all categories)
     max_total_Ni = 500.e+3  !(m)

    ! droplet concentration (m-3)
    nccnst = 200.e+6

    ! parameters for Seifert and Beheng (2001) autoconversion/accretion
    kc     = 9.44e+9
    kr     = 5.78e+3

    ! Temperature parameters
    zerodegc  = tmelt 
    homogfrze = tmelt-40.
    icenuct   = tmelt-15.
    rainfrze  = tmelt-4.

    ! physical constants
    cp     = cpair ! specific heat of dry air (J/K/kg) !1005.
    inv_cp = 1./cp ! inverse of cp
    g      = gravit ! Gravity (m/s^2) !9.816
    rd     = rair ! Dry air gas constant     ~ J/K/kg     !287.15
    rv     = rh2o ! Water vapor gas constant ~ J/K/kg     !461.51
    ep_2   = mwh2o/mwdry  ! ratio of molecular mass of water to the molecular mass of dry air !0.622
    rhosur = 100000./(rd*zerodegc) ! density of air at surface
    rhosui = 60000./(rd*253.15)
    ar     = 841.99667 
    br     = 0.8
    f1r    = 0.78
    f2r    = 0.32
    ecr    = 1.
    rhow   = rhoh2o ! Density of liquid water (STP) !997.
    cpw    = cpliq  ! specific heat of fresh h2o (J/K/kg) !4218.
    inv_rhow = 1./rhow  !inverse of (max.) density of liquid water

    ! limits for rime density [kg m-3]
    rho_rimeMin     =  50.
    rho_rimeMax     = 900.
    inv_rho_rimeMax = 1./rho_rimeMax

    ! minium allowable prognostic variables
    bsmall = qsmall*inv_rho_rimeMax

    ! Bigg (1953)
    !bimm   = 100.
    !aimm   = 0.66
    ! Barklie and Gokhale (1959)
    bimm   = 2.
    aimm   = 0.65
    rin    = 0.1e-6
    mi0    = 4.*piov3*900.*1.e-18

    eci    = 0.5
    eri    = 1.
    bcn    = 2.

    ! mean size for soft lambda_r limiter [microns]
    dbrk   = 600.e-6
    ! ratio of rain number produced to ice number loss from melting
    nmltratio = 0.2


    cons1 = piov6*rhow
    cons2 = 4.*piov3*rhow
    cons3 = 1./(cons2*(25.e-6)**3)
    cons4 = 1./(dbrk**3*pi*rhow)
    cons5 = piov6*bimm
    cons6 = piov6**2*rhow*bimm
    cons7 = 4.*piov3*rhow*(1.e-6)**3

    ! aerosol/droplet activation parameters
    mw     = 0.018
    osm    = 1.
    vi     = 3.
    epsm   = 0.9
    rhoa   = 1777.
    map    = 0.132
    ma     = 0.0284
    rr     = 8.3187
    bact   = vi*osm*epsm*mw*rhoa/(map*rhow)
    ! inv_bact = (map*rhow)/(vi*osm*epsm*mw*rhoa)    *** to replace /bact **

    ! mode 1
    inv_rm1 = 2.e+7           ! inverse aerosol mean size (m-1)
    sig1    = 2.0             ! aerosol standard deviation
    nanew1  = 300.e6          ! aerosol number mixing ratio (kg-1)
    f11     = 0.5*exp(2.5*(log(sig1))**2)
    f21     = 1. + 0.25*log(sig1)

    ! note: currently only set for a single mode, droplet activation code needs to
    !       be modified to include the second mode
    ! mode 2
    inv_rm2 = 7.6923076e+5    ! inverse aerosol mean size (m-1)
    sig2    = 2.5             ! aerosol standard deviation
    nanew2  = 0.              ! aerosol number mixing ratio (kg-1)
    f12     = 0.5*exp(2.5*(log(sig2))**2)
    f22     = 1. + 0.25*log(sig2)

    ! droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
    ! warm rain autoconversion/accretion option only (iparam = 1)
!    allocate(dnu(16))
    dnu(1)  =  0.000
    dnu(2)  = -0.557
    dnu(3)  = -0.430
    dnu(4)  = -0.307
    dnu(5)  = -0.186
    dnu(6)  = -0.067
    dnu(7)  = -0.050
    dnu(8)  = -0.167
    dnu(9)  = -0.282
    dnu(10) = -0.397
    dnu(11) = -0.512
    dnu(12) = -0.626
    dnu(13) = -0.739
    dnu(14) = -0.853
    dnu(15) = -0.966
    dnu(16) = -0.966

    ! calibration factors for ice deposition and sublimation
    !   These are adjustable ad hoc factors used to increase or decrease deposition and/or
    !   sublimation rates.  The representation of the ice capacitances are highly simplified
    !   and the appropriate values in the diffusional growth equation are uncertain.
    clbfact_dep = 1.
    clbfact_sub = 1.

    ! Don't specify lambda bounds for cloud liquid, as they are determined by
    ! pgam dynamically.
    mg_liq_props = MGHydrometeorProps(rhow, dsph, &
         min_mean_mass=min_mean_mass_liq)
  
    ! Mean ice diameter can not grow bigger than twice the autoconversion
    ! threshold for snow.
    ice_lambda_bounds = 1._r8/[2._r8*400.e-6, 10.e-6_r8] !! dcs 400.e-6
    mg_ice_props = MGHydrometeorProps(rhoi, dsph, &
         ice_lambda_bounds, min_mean_mass_ice)
  
    mg_rain_props = MGHydrometeorProps(rhow, dsph, lam_bnd_rain)
    mg_snow_props = MGHydrometeorProps(rhosn, dsph, lam_bnd_snow)
    
    return
    end subroutine micro_p3_utils_init
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine get_latent_heat(its,ite,kts,kte,v,s,f)

       integer,intent(in) :: its,ite,kts,kte
       real(r8),dimension(its:ite,kts:kte),intent(out) :: v,s,f

!       integer i,k

       v(:,:) = latvap           ! latent heat of vaporization
       s(:,:) = latvap + latice  ! latent heat of sublimation
       f(:,:) = latice           ! latent heat of fusion
 
! Original P3 definition of latent heats:   
!       do i = its,ite
!          do k = kts,kte
!          xxlv(i,k)    = 3.1484e6-2370.*t(i,k)
!          xxls(i,k)    = xxlv(i,k)+0.3337e6
!          xlf(i,k)     = xxls(i,k)-xxlv(i,k)
!          end do
!       end do
       return
    end subroutine get_latent_heat

!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine get_precip_fraction(its,ite,kts,kte,kbot,ktop,kdir,ast,qc,qr,qitot,method, &
                cldm,icldm,lcldm,rcldm)
      
       integer,intent(in)                              :: its,ite,kts,kte,kbot,ktop,kdir
       real(r8),dimension(its:ite,kts:kte),intent(in)  :: ast, qc, qr, qitot
       character(len=16),intent(in)                    :: method
       real(r8),dimension(its:ite,kts:kte),intent(out) :: cldm, icldm, lcldm, rcldm

       integer  :: i,k

       cldm(:,:)  = mincld
       icldm(:,:) = mincld
       lcldm(:,:) = mincld
       do k = kbot,ktop,kdir
          do i=its,ite
             cldm(i,k)  = max(ast(i,k), mincld)
             icldm(i,k) = max(ast(i,k), mincld)
             lcldm(i,k) = max(ast(i,k), mincld)
          end do
       end do

       DO k = ktop,kbot,-kdir
          DO i=its,ite
       !! 
       !! precipitation fraction 
       !! 
          rcldm(i,k) = cldm(i,k)
          IF (trim(method) == 'in_cloud') THEN
             IF (k /= ktop) THEN
                IF (qc(i,k) .lt. qsmall .and. qitot(i,k) .lt. qsmall) THEN
                   rcldm(i,k) = rcldm(i,k+kdir)
                END IF
             END IF
          ELSE IF (trim(method) == 'max_overlap') THEN
          ! calculate precip fraction based on maximum overlap assumption

          ! IF rain or snow mix ratios are smaller than threshold,
          ! then leave rcldm as cloud fraction at current level
             IF (k /= ktop) THEN
                IF (qr(i,k+kdir) .ge. qsmall .or. qitot(i,k+kdir) .ge. qsmall) THEN
                   rcldm(i,k) = max(cldm(i,k+kdir),rcldm(i,k))
                END IF
             END IF
          END IF
          END DO ! i
       END DO    ! k


       return
    end subroutine get_precip_fraction
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
! get cloud droplet size distribution parameters
elemental subroutine size_dist_param_liq(props, qcic, ncic, rho, pgam, lamc)
  type(MGHydrometeorProps), intent(in) :: props
  real(r8), intent(in) :: qcic
  real(r8), intent(inout) :: ncic
  real(r8), intent(in) :: rho

  real(r8), intent(out) :: pgam
  real(r8), intent(out) :: lamc

  type(MGHydrometeorProps) :: props_loc

  if (qcic > qsmall) then

     ! Local copy of properties that can be modified.
     ! (Elemental routines that operate on arrays can't modify scalar
     ! arguments.)
     props_loc = props

     ! Get pgam from fit to observations of martin et al. 1994
#if ! defined(CLUBB_BFB_S2) && ! defined(CLUBB_BFB_ALL)
     pgam = 0.0005714_r8*(ncic/1.e6_r8*rho) + 0.2714_r8
     pgam = 1._r8/(pgam**2) - 1._r8
     pgam = max(pgam, 2._r8)
     pgam = min(pgam, 15._r8)

     ! Set coefficient for use in size_dist_param_basic.
     props_loc%shape_coef = pi * props_loc%rho / 6._r8 * &
          rising_factorial(pgam+1._r8, props_loc%eff_dim)
#else
     pgam = 0.0005714_r8*1.e-6_r8*ncic*rho + 0.2714_r8
     pgam = 1._r8/(pgam**2) - 1._r8
     pgam = max(pgam, 2._r8)

     ! Set coefficient for use in size_dist_param_basic.
     ! The 3D case is so common and optimizable that we specialize it:
     if (props_loc%eff_dim == 3._r8) then
        props_loc%shape_coef = pi / 6._r8 * props_loc%rho * &
             rising_factorial(pgam+1._r8, 3)
     else
        props_loc%shape_coef = pi / 6._r8 * props_loc%rho * &
             rising_factorial(pgam+1._r8, props_loc%eff_dim)
     end if
#endif

     ! Limit to between 2 and 50 microns mean size.
     props_loc%lambda_bounds = (pgam+1._r8)*1._r8/[50.e-6_r8, 2.e-6_r8]

     call size_dist_param_basic(props_loc, qcic, ncic, lamc)

  else
     ! pgam not calculated in this case, so set it to a value likely to
     ! cause an error if it is accidentally used
     ! (gamma function undefined for negative integers)
     pgam = -100._r8
     lamc = 0._r8
  end if

end subroutine size_dist_param_liq
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
! Basic routine for getting size distribution parameters.
elemental subroutine size_dist_param_basic(props, qic, nic, lam, n0)
  type(MGHydrometeorProps), intent(in) :: props
  real(r8), intent(in) :: qic
  real(r8), intent(inout) :: nic

  real(r8), intent(out) :: lam
  real(r8), intent(out), optional :: n0

  if (qic > qsmall) then

     ! add upper limit to in-cloud number concentration to prevent
     ! numerical error
     if (limiter_is_on(props%min_mean_mass)) then
        nic = min(nic, qic / props%min_mean_mass)
     end if

     ! lambda = (c n/q)^(1/d)
     lam = (props%shape_coef * nic/qic)**(1._r8/props%eff_dim)

     ! check for slope
     ! adjust vars
     if (lam < props%lambda_bounds(1)) then
        lam = props%lambda_bounds(1)
        nic = lam**(props%eff_dim) * qic/props%shape_coef
     else if (lam > props%lambda_bounds(2)) then
        lam = props%lambda_bounds(2)
        nic = lam**(props%eff_dim) * qic/props%shape_coef
     end if

  else
     lam = 0._r8
  end if

  if (present(n0)) n0 = nic * lam

end subroutine size_dist_param_basic
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!

real(r8) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  real(r8), intent(in) :: q         ! mass mixing ratio
  real(r8), intent(in) :: n         ! number concentration (per volume)
  real(r8), intent(in) :: rho_air   ! local density of the air
  real(r8), intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (pi * rho_sub * n/(q*rho_air))**(-1._r8/3._r8)

end function avg_diameter
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
! Constructor for a constituent property object.
function NewMGHydrometeorProps(rho, eff_dim, lambda_bounds, min_mean_mass) &
     result(res)
  real(r8), intent(in) :: rho, eff_dim
  real(r8), intent(in), optional :: lambda_bounds(2), min_mean_mass
  type(MGHydrometeorProps) :: res

  res%rho = rho
  res%eff_dim = eff_dim
  if (present(lambda_bounds)) then
     res%lambda_bounds = lambda_bounds
  else
     res%lambda_bounds = no_limiter()
  end if
  if (present(min_mean_mass)) then
     res%min_mean_mass = min_mean_mass
  else
     res%min_mean_mass = no_limiter()
  end if

  res%shape_coef = rho*pi*gamma(eff_dim+1._r8)/6._r8

end function NewMGHydrometeorProps
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
!========================================================================
!UTILITIES
!========================================================================

pure function no_limiter()
  real(r8) :: no_limiter

  no_limiter = transfer(limiter_off, no_limiter)

end function no_limiter

pure function limiter_is_on(lim)
  real(r8), intent(in) :: lim
  logical :: limiter_is_on

  limiter_is_on = transfer(lim, limiter_off) /= limiter_off

end function limiter_is_on
!========================================================================
!FORMULAS
!========================================================================

! Use gamma function to implement rising factorial extended to the reals.
pure function rising_factorial_r8(x, n) result(res)
  real(r8), intent(in) :: x, n
  real(r8) :: res

  res = gamma(x+n)/gamma(x)

end function rising_factorial_r8

! Rising factorial can be performed much cheaper if n is a small integer.
pure function rising_factorial_integer(x, n) result(res)
  real(r8), intent(in) :: x
  integer, intent(in) :: n
  real(r8) :: res

  integer :: i
  real(r8) :: factor

  res = 1._r8
  factor = x

  do i = 1, n
     res = res * factor
     factor = factor + 1._r8
  end do

end function rising_factorial_integer

elemental function var_coef_r8(relvar, a) result(res)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(r8), intent(in) :: relvar
  real(r8), intent(in) :: a
  real(r8) :: res

  res = rising_factorial(relvar, a) / relvar**a

end function var_coef_r8

elemental function var_coef_integer(relvar, a) result(res)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(r8), intent(in) :: relvar
  integer, intent(in) :: a
  real(r8) :: res

  res = rising_factorial(relvar, a) / relvar**a

end function var_coef_integer




end module micro_p3_utils
