module micro_p3_utils

#ifdef SCREAM_CONFIG_IS_CMAKE
    use iso_c_binding, only: c_double, c_float
#else
    use shr_kind_mod,   only: rtype=>shr_kind_r8, itype=>shr_kind_i8
#endif

    implicit none
    private
    save

#ifdef SCREAM_CONFIG_IS_CMAKE
#include "scream_config.f"

#ifdef SCREAM_DOUBLE_PRECISION
  integer,parameter,public :: rtype = c_double ! 8 byte real, compatible with c type double
#else
  integer,parameter,public :: rtype = c_float ! 4 byte real, compatible with c type float
#endif
  integer,parameter :: itype = selected_int_kind (13) ! 8 byte integer
#else
    public :: rtype
#endif

    public :: get_latent_heat, micro_p3_utils_init, size_dist_param_liq, &
              size_dist_param_basic,avg_diameter, rising_factorial, calculate_incloud_mixingratios

    integer, public :: iulog_e3sm
    logical, public :: masterproc_e3sm

    ! Signaling NaN bit pattern that represents a limiter that's turned off.
    integer(itype), parameter :: limiter_off = int(Z'7FF1111111111111', itype)

    real(rtype), public, parameter :: qsmall = 1.e-14_rtype
    real(rtype), public, parameter :: nsmall = 1.e-16_rtype

    real(rtype) :: xxlv, xxls, xlf

    real(rtype),public :: rhosur,rhosui,ar,br,f1r,f2r,ecr,rhow,kr,kc,aimm,bimm,rin,mi0,nccnst,  &
       eci,eri,bcn,cpw,cons1,cons2,cons3,cons4,cons5,cons6,cons7,         &
       inv_rhow,cp,g,rd,rv,ep_2,inv_cp,   &
       thrd,sxth,piov3,piov6,rho_rimeMin,     &
       rho_rimeMax,inv_rho_rimeMax,max_total_Ni,dbrk,nmltratio,clbfact_sub,  &
       clbfact_dep
    real(rtype),dimension(16), public :: dnu

    real(rtype), public, parameter :: mu_r_constant = 1.0_rtype
    real(rtype), public, parameter :: lookup_table_1a_dum1_c = 1.0_rtype/(0.1_rtype*log10(261.7_rtype))

    real(rtype),public :: zerodegc  ! Temperature at zero degree celcius ~K
    real(rtype),public :: rainfrze  ! Contact and immersion freexing temp, -4C  ~K
    real(rtype),public :: homogfrze ! Homogeneous freezing temperature, -40C  ~K
    real(rtype),public :: icenuct   ! Ice nucleation temperature, -5C ~K

    real(rtype),public :: pi_e3sm
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

    real(rtype), parameter, public :: mincld=0.0001_rtype
    real(rtype), parameter, public :: rhosn = 250._rtype  ! bulk density snow
    real(rtype), parameter, public :: rhoi = 500._rtype   ! bulk density ice
    real(rtype), parameter, public :: rhows = 917._rtype  ! bulk density water solid


public :: MicroHydrometeorProps

type :: MicroHydrometeorProps
   ! Density (kg/m^3)
   real(rtype) :: rho
   ! Information for size calculations.
   ! Basic calculation of mean size is:
   !     lambda = (shape_coef*nic/qic)^(1/eff_dim)
   ! Then lambda is constrained by bounds.
   real(rtype) :: eff_dim
   real(rtype) :: shape_coef
   real(rtype) :: lambda_bounds(2)
   ! Minimum average particle mass (kg).
   ! Limit is applied at the beginning of the size distribution calculations.
   real(rtype) :: min_mean_mass
end type MicroHydrometeorProps

interface MicroHydrometeorProps
   module procedure NewMicroHydrometeorProps
end interface

type(MicroHydrometeorProps), public :: micro_liq_props
type(MicroHydrometeorProps), public :: micro_ice_props
type(MicroHydrometeorProps), public :: micro_rain_props
type(MicroHydrometeorProps), public :: micro_snow_props

! particle mass-diameter relationship
! currently we assume spherical particles for cloud ice/snow
! m = cD^d
! exponent
real(rtype), parameter :: dsph = 3._rtype

! Bounds for mean diameter for different constituents.
real(rtype), parameter :: lam_bnd_rain(2) = 1._rtype/[500.e-6_rtype, 20.e-6_rtype]
real(rtype), parameter :: lam_bnd_snow(2) = 1._rtype/[2000.e-6_rtype, 10.e-6_rtype]

! Minimum average mass of particles.
real(rtype), parameter :: min_mean_mass_liq = 1.e-20_rtype
real(rtype), parameter :: min_mean_mass_ice = 1.e-20_rtype

! in-cloud values
REAL(rtype), PARAMETER :: cldm_min   = 1.e-20_rtype !! threshold min value for cloud fraction
real(rtype), parameter :: incloud_limit = 5.1E-3
real(rtype), parameter :: precip_limit  = 1.0E-2
!=========================================================
! Utilities that are cheaper if the compiler knows that
! some argument is an integer.
!=========================================================

interface rising_factorial
   module procedure rising_factorial_rtype
   module procedure rising_factorial_integer
end interface rising_factorial

interface var_coef
   module procedure var_coef_rtype
   module procedure var_coef_integer
end interface var_coef

    contains
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine micro_p3_utils_init(cpair,rair,rh2o,rhoh2o,mwh2o,mwdry,gravit,latvap,latice, &
                   cpliq,tmelt,pi,iulog,masterproc)

    real(rtype), intent(in) :: cpair
    real(rtype), intent(in) :: rair
    real(rtype), intent(in) :: rh2o
    real(rtype), intent(in) :: rhoh2o
    real(rtype), intent(in) :: mwh2o
    real(rtype), intent(in) :: mwdry
    real(rtype), intent(in) :: gravit
    real(rtype), intent(in) :: latvap
    real(rtype), intent(in) :: latice
    real(rtype), intent(in) :: cpliq
    real(rtype), intent(in) :: tmelt
    real(rtype), intent(in) :: pi
    integer, intent(in)     :: iulog
    logical, intent(in)     :: masterproc

    real(rtype) :: ice_lambda_bounds(2)

    ! logfile info
    iulog_e3sm      = iulog
    masterproc_e3sm = masterproc

    ! mathematical/optimization constants
    thrd  = 1._rtype/3._rtype
    sxth  = 1._rtype/6._rtype 
    pi_e3sm = pi
    piov3 = pi*thrd
    piov6 = pi*sxth

    ! maximum total ice concentration (sum of all categories)
     max_total_Ni = 500.e+3_rtype  !(m)

    ! droplet concentration (m-3)
    nccnst = 200.e+6_rtype

    ! parameters for Seifert and Beheng (2001) autoconversion/accretion
    kc     = 9.44e+9_rtype
    kr     = 5.78e+3_rtype

    ! Temperature parameters
    zerodegc  = tmelt 
    homogfrze = tmelt-40._rtype
    icenuct   = tmelt-15._rtype
    rainfrze  = tmelt-4._rtype

    ! physical constants
    cp     = cpair ! specific heat of dry air (J/K/kg) !1005.
    inv_cp = 1._rtype/cp ! inverse of cp
    g      = gravit ! Gravity (m/s^2) !9.816
    rd     = rair ! Dry air gas constant     ~ J/K/kg     !287.15
    rv     = rh2o ! Water vapor gas constant ~ J/K/kg     !461.51
    ep_2   = mwh2o/mwdry  ! ratio of molecular mass of water to the molecular mass of dry air !0.622
    rhosur = 100000._rtype/(rd*zerodegc) ! density of air at surface
    rhosui = 60000._rtype/(rd*253.15_rtype)
    ar     = 841.99667_rtype 
    br     = 0.8_rtype
    f1r    = 0.78_rtype
    f2r    = 0.32_rtype
    ecr    = 1._rtype
    rhow   = rhoh2o ! Density of liquid water (STP) !997.
    cpw    = cpliq  ! specific heat of fresh h2o (J/K/kg) !4218.
    inv_rhow = 1._rtype/rhow  !inverse of (max.) density of liquid water

    xxlv = latvap           ! latent heat of vaporization
    xxls = latvap + latice  ! latent heat of sublimation
    xlf  = latice           ! latent heat of fusion

    ! limits for rime density [kg m-3]
    rho_rimeMin     =  50._rtype
    rho_rimeMax     = 900._rtype
    inv_rho_rimeMax =   1._rtype/rho_rimeMax

    ! Bigg (1953)
    !bimm   = 100.
    !aimm   = 0.66
    ! Barklie and Gokhale (1959)
    bimm   = 2._rtype
    aimm   = 0.65_rtype
    rin    = 0.1e-6_rtype
    mi0    = 4._rtype*piov3*900._rtype*1.e-18_rtype

    eci    = 0.5_rtype
    eri    = 1._rtype
    bcn    = 2._rtype

    ! mean size for soft lambda_r limiter [microns]
    dbrk   = 600.e-6_rtype
    ! ratio of rain number produced to ice number loss from melting
    nmltratio = 0.2_rtype

    cons1 = piov6*rhow
    cons2 = 4._rtype*piov3*rhow
    cons3 = 1._rtype/(cons2*(25.e-6_rtype)**3)
    cons4 = 1._rtype/(dbrk**3*pi*rhow)
    cons5 = piov6*bimm
    cons6 = piov6**2*rhow*bimm
    cons7 = 4._rtype*piov3*rhow*(1.e-6_rtype)**3

    ! droplet spectral shape parameter for mass spectra, used for Seifert and Beheng (2001)
    ! warm rain autoconversion/accretion option only (iparam = 1)
!    allocate(dnu(16))
    dnu(1)  =  0.000_rtype
    dnu(2)  = -0.557_rtype
    dnu(3)  = -0.430_rtype
    dnu(4)  = -0.307_rtype
    dnu(5)  = -0.186_rtype
    dnu(6)  = -0.067_rtype
    dnu(7)  = -0.050_rtype
    dnu(8)  = -0.167_rtype
    dnu(9)  = -0.282_rtype
    dnu(10) = -0.397_rtype
    dnu(11) = -0.512_rtype
    dnu(12) = -0.626_rtype
    dnu(13) = -0.739_rtype
    dnu(14) = -0.853_rtype
    dnu(15) = -0.966_rtype
    dnu(16) = -0.966_rtype

    ! calibration factors for ice deposition and sublimation
    !   These are adjustable ad hoc factors used to increase or decrease deposition and/or
    !   sublimation rates.  The representation of the ice capacitances are highly simplified
    !   and the appropriate values in the diffusional growth equation are uncertain.
    clbfact_dep = 1._rtype
    clbfact_sub = 1._rtype

    ! Don't specify lambda bounds for cloud liquid, as they are determined by
    ! pgam dynamically.
    micro_liq_props = MicroHydrometeorProps(rhow, dsph, &
         min_mean_mass=min_mean_mass_liq)
  
    ! Mean ice diameter can not grow bigger than twice the autoconversion
    ! threshold for snow.
    ice_lambda_bounds = 1._rtype/[2._rtype*400.e-6_rtype, 10.e-6_rtype] !! dcs 400.e-6
    micro_ice_props = MicroHydrometeorProps(rhoi, dsph, &
         ice_lambda_bounds, min_mean_mass_ice)
  
    micro_rain_props = MicroHydrometeorProps(rhow, dsph, lam_bnd_rain)
    micro_snow_props = MicroHydrometeorProps(rhosn, dsph, lam_bnd_snow)
    
    return
    end subroutine micro_p3_utils_init
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
    subroutine get_latent_heat(its,ite,kts,kte,v,s,f)

       integer,intent(in) :: its,ite,kts,kte
       real(rtype),dimension(its:ite,kts:kte),intent(out) :: v,s,f

!       integer i,k

       v(:,:) = xxlv !latvap           ! latent heat of vaporization
       s(:,:) = xxls !latvap + latice  ! latent heat of sublimation
       f(:,:) = xlf  !latice           ! latent heat of fusion
 
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
    subroutine calculate_incloud_mixingratios(qc,qr,qitot,qirim,nc,nr,nitot,birim, &
          inv_lcldm,inv_icldm,inv_rcldm, &
          qc_incld,qr_incld,qitot_incld,qirim_incld,nc_incld,nr_incld,nitot_incld,birim_incld)

       real(rtype),intent(in)   :: qc, qr, qitot, qirim
       real(rtype),intent(in)   :: nc, nr, nitot, birim
       real(rtype),intent(in)   :: inv_lcldm, inv_icldm, inv_rcldm
       real(rtype),intent(out)  :: qc_incld, qr_incld, qitot_incld, qirim_incld
       real(rtype),intent(out)  :: nc_incld, nr_incld, nitot_incld, birim_incld


       if (qc.ge.qsmall) then
          qc_incld = qc*inv_lcldm
          nc_incld = max(nc*inv_lcldm,0._rtype)
          !AaronDonahue, kai has something about if nccons then nc=ncnst/rho
       else
          qc_incld = 0._rtype
          nc_incld = 0._rtype
       end if 
       if (qitot.ge.qsmall) then
          qitot_incld = qitot*inv_icldm
          nitot_incld = max(nitot*inv_icldm,0._rtype)
          !AaronDonahue, kai has something about if nicons then ni=ninst/rho
       else
          qitot_incld = 0._rtype
          nitot_incld = 0._rtype
       end if 
       if (qirim.ge.qsmall.and.qitot.ge.qsmall) then
          qirim_incld = qirim*inv_icldm
          birim_incld = max(birim*inv_lcldm,0._rtype)
       else
          qirim_incld = 0._rtype
          birim_incld = 0._rtype
       end if 
       if (qr.ge.qsmall) then
          qr_incld = qr*inv_rcldm
          nr_incld = max(nr*inv_rcldm,0._rtype)
          !AaronDonahue, kai has something about if nccons then nc=ncnst/rho
       else
          qr_incld = 0._rtype
          nr_incld = 0._rtype
       end if
       if (qc_incld.gt.incloud_limit .or.qitot_incld.gt.incloud_limit .or. qr_incld.gt.precip_limit .or.birim_incld.gt.incloud_limit) then
!          write(errmsg,'(a3,i4,3(a5,1x,e16.8,1x))') 'k: ', k, ', qc:',qc_incld, ', qi:',qitot_incld,', qr:',qr_incld
          qc_incld    = max(qc_incld,incloud_limit)
          qitot_incld = max(qitot_incld,incloud_limit)
          birim_incld = max(birim_incld,incloud_limit)
          qr_incld    = max(qr_incld,precip_limit)
!          if (masterproc) write(iulog,*)  errmsg

!          call handle_errmsg('Micro-P3 (Init)',subname='In-cloud mixing
!          ratio too large',extra_msg=errmsg)
       end if
    end subroutine calculate_incloud_mixingratios
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
! get cloud droplet size distribution parameters
elemental subroutine size_dist_param_liq(props, qcic, ncic, rho, pgam, lamc)
  type(MicroHydrometeorProps), intent(in) :: props
  real(rtype), intent(in) :: qcic
  real(rtype), intent(inout) :: ncic
  real(rtype), intent(in) :: rho

  real(rtype), intent(out) :: pgam
  real(rtype), intent(out) :: lamc

  type(MicroHydrometeorProps) :: props_loc

  if (qcic > qsmall) then

     ! Local copy of properties that can be modified.
     ! (Elemental routines that operate on arrays can't modify scalar
     ! arguments.)
     props_loc = props

     ! Get pgam from fit to observations of martin et al. 1994
#if ! defined(CLUBB_BFB_S2) && ! defined(CLUBB_BFB_ALL)
     pgam = 0.0005714_rtype*(ncic/1.e6_rtype*rho) + 0.2714_rtype
     pgam = 1._rtype/(pgam**2) - 1._rtype
     pgam = max(pgam, 2._rtype)
     pgam = min(pgam, 15._rtype)

     ! Set coefficient for use in size_dist_param_basic.
     props_loc%shape_coef = pi_e3sm * props_loc%rho / 6._rtype * &
          rising_factorial(pgam+1._rtype, props_loc%eff_dim)
#else
     pgam = 0.0005714_rtype*1.e-6_rtype*ncic*rho + 0.2714_rtype
     pgam = 1._rtype/(pgam**2) - 1._rtype
     pgam = max(pgam, 2._rtype)

     ! Set coefficient for use in size_dist_param_basic.
     ! The 3D case is so common and optimizable that we specialize it:
     if (props_loc%eff_dim == 3._rtype) then
        props_loc%shape_coef = pi_e3sm / 6._rtype * props_loc%rho * &
             rising_factorial(pgam+1._rtype, 3)
     else
        props_loc%shape_coef = pi_e3sm / 6._rtype * props_loc%rho * &
             rising_factorial(pgam+1._rtype, props_loc%eff_dim)
     end if
#endif

     ! Limit to between 2 and 50 microns mean size.
     props_loc%lambda_bounds = (pgam+1._rtype)*1._rtype/[50.e-6_rtype, 2.e-6_rtype]

     call size_dist_param_basic(props_loc, qcic, ncic, lamc)

  else
     ! pgam not calculated in this case, so set it to a value likely to
     ! cause an error if it is accidentally used
     ! (gamma function undefined for negative integers)
     pgam = -100._rtype
     lamc = 0._rtype
  end if

end subroutine size_dist_param_liq
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
! Basic routine for getting size distribution parameters.
elemental subroutine size_dist_param_basic(props, qic, nic, lam, n0)
  type(MicroHydrometeorProps), intent(in) :: props
  real(rtype), intent(in) :: qic
  real(rtype), intent(inout) :: nic

  real(rtype), intent(out) :: lam
  real(rtype), intent(out), optional :: n0

  if (qic > qsmall) then

     ! add upper limit to in-cloud number concentration to prevent
     ! numerical error
     if (limiter_is_on(props%min_mean_mass)) then
        nic = min(nic, qic / props%min_mean_mass)
     end if

     ! lambda = (c n/q)^(1/d)
     lam = (props%shape_coef * nic/qic)**(1._rtype/props%eff_dim)

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
     lam = 0._rtype
  end if

  if (present(n0)) n0 = nic * lam

end subroutine size_dist_param_basic
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!

real(rtype) elemental function avg_diameter(q, n, rho_air, rho_sub)
  ! Finds the average diameter of particles given their density, and
  ! mass/number concentrations in the air.
  ! Assumes that diameter follows an exponential distribution.
  real(rtype), intent(in) :: q         ! mass mixing ratio
  real(rtype), intent(in) :: n         ! number concentration (per volume)
  real(rtype), intent(in) :: rho_air   ! local density of the air
  real(rtype), intent(in) :: rho_sub   ! density of the particle substance

  avg_diameter = (pi_e3sm * rho_sub * n/(q*rho_air))**(-1._rtype/3._rtype)

end function avg_diameter
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
! Constructor for a constituent property object.
function NewMicroHydrometeorProps(rho, eff_dim, lambda_bounds, min_mean_mass) &
     result(res)
  real(rtype), intent(in) :: rho, eff_dim
  real(rtype), intent(in), optional :: lambda_bounds(2), min_mean_mass
  type(MicroHydrometeorProps) :: res

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

  res%shape_coef = rho*pi_e3sm*gamma(eff_dim+1._rtype)/6._rtype

end function NewMicroHydrometeorProps
!__________________________________________________________________________________________!
!                                                                                          !
!__________________________________________________________________________________________!
!========================================================================
!UTILITIES
!========================================================================

pure function no_limiter()
  real(rtype) :: no_limiter

  no_limiter = transfer(limiter_off, no_limiter)

end function no_limiter

pure function limiter_is_on(lim)
  real(rtype), intent(in) :: lim
  logical :: limiter_is_on

  limiter_is_on = transfer(lim, limiter_off) /= limiter_off

end function limiter_is_on
!========================================================================
!FORMULAS
!========================================================================

! Use gamma function to implement rising factorial extended to the reals.
pure function rising_factorial_rtype(x, n) result(res)
  real(rtype), intent(in) :: x, n
  real(rtype) :: res

  res = gamma(x+n)/gamma(x)

end function rising_factorial_rtype

! Rising factorial can be performed much cheaper if n is a small integer.
pure function rising_factorial_integer(x, n) result(res)
  real(rtype), intent(in) :: x
  integer, intent(in) :: n
  real(rtype) :: res

  integer :: i
  real(rtype) :: factor

  res = 1._rtype
  factor = x

  do i = 1, n
     res = res * factor
     factor = factor + 1._rtype
  end do

end function rising_factorial_integer

elemental function var_coef_rtype(relvar, a) result(res)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(rtype), intent(in) :: relvar
  real(rtype), intent(in) :: a
  real(rtype) :: res

  res = rising_factorial(relvar, a) / relvar**a

end function var_coef_rtype

elemental function var_coef_integer(relvar, a) result(res)
  ! Finds a coefficient for process rates based on the relative variance
  ! of cloud water.
  real(rtype), intent(in) :: relvar
  integer, intent(in) :: a
  real(rtype) :: res

  res = rising_factorial(relvar, a) / relvar**a

end function var_coef_integer




end module micro_p3_utils
