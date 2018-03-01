
module EOSWaterMod

#ifdef USE_PETSC_LIB

#include <petsc/finclude/petsc.h>

  ! !USES:
  use petscsys
  use mpp_varctl                , only : iulog
  use mpp_abortutils            , only : endrun
  use mpp_shr_log_mod           , only : errMsg => shr_log_errMsg
  use MultiPhysicsProbConstants , only : FMWH2O
  !
  ! !PUBLIC TYPES:
  implicit none
  private

  PetscInt, parameter, public :: DENSITY_CONSTANT         = 1
  PetscInt, parameter, public :: DENSITY_TGDPB01          = 2
  PetscInt, parameter, public :: DENSITY_IFC67            = 3

  PetscInt, parameter, public :: INT_ENERGY_ENTHALPY_CONSTANT  = 1
  PetscInt, parameter, public :: INT_ENERGY_ENTHALPY_IFC67     = 2

  ! For IFC-67
  PetscReal, parameter        :: H2O_CRITICAL_TEMPERATURE = 647.3d0  ! K
  PetscReal, parameter        :: H2O_CRITICAL_PRESSURE    = 22.064d6 ! Pa

  public :: Density
  public :: Viscosity
  public :: InternalEnergyAndEnthalpy

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine Density(p, t_K, density_itype, den, dden_dp, dden_dT)
    !
    ! !DESCRIPTION:
    ! Given pressure, temperature and type of density formulation, compute:
    ! - density of water, and
    ! - first derivative of density w.r.t pressure
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal , intent(in)  :: p
    PetscReal , intent(in)  :: t_K
    PetscInt  , intent(in)  :: density_itype
    PetscReal , intent(out) :: den
    PetscReal , intent(out) :: dden_dp
    PetscReal , intent(out) :: dden_dT
    !
    PetscReal               :: t_C
    PetscReal               :: den_kg
    PetscBool               :: cal_deriv
    PetscErrorCode          :: ierr

    select case(density_itype)
    case (DENSITY_CONSTANT)
       call DensityConstant(den, dden_dp, dden_dT)
    case (DENSITY_TGDPB01)
       call DensityTGDPB01 (p , t_K, den, dden_dp,dden_dT )
    case (DENSITY_IFC67)
       t_C = t_K - 273.15d0
       cal_deriv = PETSC_TRUE
       call DensityIFC67(t_C, p, cal_deriv, den_kg, den, &
            dden_dp,dden_dT, ierr)
    case default
       write(iulog,*)'Density: Unknown denity_itype. '
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine Density

  !------------------------------------------------------------------------
  subroutine DensityConstant(den, dden_dp, dden_dT)
    !
    ! !DESCRIPTION:
    ! Return constant density of water
    !
    ! !USES:
    use mpp_varcon                , only : denh2o
    !
    implicit none
    !
    ! !ARGUMENTS    
    PetscReal, intent(out) :: den       ! [kmol m^{-3}]
    PetscReal, intent(out) :: dden_dp   ! [kmol m^{-3} Pa^{-1}]
    PetscReal, intent(out) :: dden_dT   ! [kmol m^{-3} K^{-1}]

    den     = denh2o/FMWH2O ! [kmol m^{-3}]
    dden_dp = 0.d0
    dden_dT = 0.d0

  end subroutine DensityConstant

  !------------------------------------------------------------------------
  subroutine DensityTGDPB01(p, t_K, den, dden_dp, dden_dT)
    !
    ! !DESCRIPTION:
    ! Return density and deriv. w.r.t. pressure based on Tanaka et al. (2001)
    !
    ! Reference:
    ! Tanaka M. , G. Girard, R. Davis, A. Peuto, and N. Bignell. 2001.
    ! Recommended table for the density of water between 0 °C
    ! and 40 °C based on recent experimental reports. Metrologia,
    ! 38:301-309 [doi:10.1088/0026-1394/38/4/3].
    !
    ! !USES:
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: p                ! [Pa]
    PetscReal, intent(in)  :: t_K              ! [K]
    PetscReal, intent(out) :: den              ! [kmol m^{-3}]
    PetscReal, intent(out) :: dden_dp          ! [kmol m^{-3} Pa^{-1}]
    PetscReal, intent(out) :: dden_dT          ! [kmol m^{-3} K^{-1}]

    !
    PetscReal,parameter    :: a1 = -3.983035d0     ! [degC]
    PetscReal,parameter    :: a2 = 301.797d0       ! [degC]
    PetscReal,parameter    :: a3 = 522528.9d0      ! [degC^{2}]
    PetscReal,parameter    :: a4 = 69.34881d0      ! [degC]
    PetscReal,parameter    :: a5 = 999.974950d0    ! [kg m^{-3}]
    PetscReal,parameter    :: k0 = 50.74d-11       ! [Pa^{-1}]
    PetscReal,parameter    :: k1 = -0.326d-11      ! [Pa^{-1} degC^{-1}]
    PetscReal,parameter    :: k2 = 0.00416d-11     ! [Pa^{-1} degC^{-2}]
    PetscReal,parameter    :: p0 = 101325.d0       ! [Pa]
    PetscReal              :: t_c
    PetscReal              :: dent
    PetscReal              :: kappa
    PetscReal              :: ddent_dt
    PetscReal              :: ddent_dt_1
    PetscReal              :: ddent_dt_2
    PetscReal              :: ddent_dt_3
    PetscReal              :: ddent_dp
    PetscReal              :: dkappa_dp
    PetscReal              :: dkappa_dt

    t_c = t_K - 273.15d0

    ! Density of water as function of t_K
    dent = a5*(1.d0 - ((t_c + a1)**2.d0)*(t_c + a2)/a3/(t_c + a4))

    ! Compressibility of water
    if (p > p0) then
       kappa = (1.d0 + (k0 + k1*t_c + k2*t_c**2.d0)*(p - p0))
    else
       kappa = 1.d0
    endif

    ! Density of water
    den = dent*kappa/FMWH2O ! [kmol m^{-3}]

    ! Derivative
    ddent_dp   = 0.d0
    ddent_dt_1 = -((t_c + a1)**2.d0)/a3/(t_c + a4)
    ddent_dt_2 = -2.d0*(t_c + a1)*(t_c + a2)/a3/(t_c + a4)
    ddent_dt_3 =  ((t_c + a1)**2.d0)*(t_c + a2)/a3/((t_c + a4)**2.d0)
    ddent_dt   = a5*(ddent_dt_1 + ddent_dt_2 + ddent_dt_3)

    if (p > p0) then
       dkappa_dp  = (k0 + k1*t_c + k2*t_c**2.d0)
       dkappa_dt  = (k1 + 2.d0*k2*t_c)*(p - p0)
    else
       dkappa_dp  = 0.d0
       dkappa_dt  = 0.d0
    endif

    dden_dT    = (ddent_dt*kappa + dent*dkappa_dt)/FMWH2O
    dden_dp    = (ddent_dp*kappa + dent*dkappa_dp)/FMWH2O

  end subroutine DensityTGDPB01

  !------------------------------------------------------------------------
  subroutine DensityIFC67(t,p,calculate_derivatives,dw,dwmol, &
       dwp,dwt,ierr)

    !  This subroutine calculates water and steam-gas mixture properties.
    !  The water and steam properties are valid in the range of:
    !
    !            0 < p < 165.4 * 10^5 pascals (165.4 bars)
    !            0 < t < 350 centigrade (623.15 Kelvin)
    !
    !  The properties cover densities, enthalpies, internal energies,
    !  and partial derivatives of these quanties with respect to
    !  pressure and temperature.
    !
    !  For saturated fluid, it will also calculate water saturation
    !  temperature for the specified pressure using Newton-Raphson and
    !  the derivative dts/dp (=tsp) or Ps for a given temperature.
    !
    !  Ref.: International Formulation Committee of the Sixth International
    !       Conference on Properties of Steam (1967).
    !
    implicit none

    PetscReal      , intent(in)  :: t   ! Temperature in centigrade
    PetscReal      , intent(in)  :: p   ! Pressure in Pascals
    PetscBool      , intent(in)  :: calculate_derivatives
    PetscReal      , intent(out) :: dw,dwmol,dwp,dwt
    PetscErrorCode , intent(out) :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt             :: i

    PetscReal, save      :: aa(0:22)
    PetscReal, save      :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12

    PetscReal            :: beta,beta2x,beta4,theta,utheta,theta2x,theta18,theta20
    PetscReal            :: xx,yy,zz
    PetscReal            :: u0,u1,u2,u3,u4,u5,u6,u7,u8,u9
    PetscReal            :: tempreal
    PetscReal            :: vr,ypt,yptt,zpt,zpp,vrpt,vrpp,cnv
    PetscReal            :: tc1,pc1,vc1,utc1,upc1,vc1mol
    PetscReal            :: d2z_dp2    ! 2nd derivative of z w.r.t. pressure
    PetscReal            :: d2vr_dp2   ! 2nd derivative of vr w.r.t. pressure
    PetscReal            :: d2wmol_dp2 ! 2nd derivative of density w.r.t. pressure
    PetscReal, parameter :: zero = 0.d0
    PetscReal, parameter :: one = 1.d0
    PetscReal, parameter :: two = 2.d0
    PetscReal, parameter :: three = 3.d0
    PetscReal, parameter :: four = 4.d0
    PetscReal, parameter :: five = 5.d0
    PetscReal, parameter :: six = 6.d0
    PetscReal, parameter :: seven = 7.d0
    PetscReal, parameter :: eight = 8.d0
    PetscReal, parameter :: nine = 9.d0
    PetscReal, parameter :: ten = 10.d0

    data aa/ &
         !-----data aa0,aa1,aa2,aa3/
         6.824687741d03,-5.422063673d02,-2.096666205d04, 3.941286787d04, &
         !-----data aa4,aa5,aa6,aa7/
         -6.733277739d04, 9.902381028d04,-1.093911774d05, 8.590841667d04, &
         !-----data aa8,aa9,aa10,aa11/
         -4.511168742d04, 1.418138926d04,-2.017271113d03, 7.982692717d00, &
         !-----data aa12,aa13,aa14,aa15/
         -2.616571843d-2, 1.522411790d-3, 2.284279054d-2, 2.421647003d02, &
         !-----data aa16,aa17,aa18,aa19/
         1.269716088d-10,2.074838328d-7, 2.174020350d-8, 1.105710498d-9, &
         !-----data aa20,aa21,aa22/
         1.293441934d01, 1.308119072d-5, 6.047626338d-14/

    data a1,a2,a3,a4/ &
         8.438375405d-1, 5.362162162d-4, 1.720000000d00, 7.342278489d-2/
    data a5,a6,a7,a8/ &
         4.975858870d-2, 6.537154300d-1, 1.150000000d-6, 1.510800000d-5/
    data a9,a10,a11,a12/ &
         1.418800000d-1, 7.002753165d00, 2.995284926d-4, 2.040000000d-1/

    ierr = 0

    tc1     = H2O_CRITICAL_TEMPERATURE ! K
    pc1     = H2O_CRITICAL_PRESSURE    ! Pa
    vc1     = 0.00317d0                ! m^3/kg
    utc1    = one/tc1                  ! 1/C
    upc1    = one/pc1                  ! 1/Pa
    vc1mol  = vc1*FMWH2O               ! m^3/kmol

    theta   = (t + 273.15d0)*utc1
    theta2x = theta*theta
    theta18 = theta**18.d0
    theta20 = theta18*theta2x

    beta    = p*upc1
    beta2x  = beta*beta
    beta4   = beta2x*beta2x

    yy      = one-a1*theta2x-a2*theta**(-6.d0)
    xx      = a3*yy*yy-two*(a4*theta-a5*beta)

    !   Note: xx may become negative near the critical point-pcl.
    if (xx.gt.zero) then
       xx = sqrt(xx)
    else
       write(iulog,*)'Negative term in density ','t= ',t,' p= ',p,' xx= ',xx
       call endrun(msg=errMsg(__FILE__, __LINE__))
       ierr = 1
       xx = 1.e-6  !set arbitrarily
    end if
    zz = yy + xx
    u0 = -five/17.d0
    u1 = aa(11)*a5*zz**u0
    u2 = one/(a8+theta**11.d0)
    u3 = aa(17)+(two*aa(18)+three*aa(19)*beta)*beta
    u4 = one/(a7+theta18*theta)
    u5 = (a10+beta)**(-4.d0)
    u6 = a11-three*u5
    u7 = aa(20)*theta18*(a9+theta2x)
    u8 = aa(15)*(a6-theta)**9.d0

    vr = u1+aa(12)+theta*(aa(13)+aa(14)*theta)+u8*(a6-theta) &
         +aa(16)*u4-u2*u3-u6*u7+(three*aa(21)*(a12-theta) &
         +four*aa(22)*beta/theta20)*beta2x

    dwmol = one/(vr*vc1mol) ! kmol/m^3
    dw = one/(vr*vc1) ! kg/m^3

    ! ypt used for enthalpy even if derivative not calculated
    ypt = six*a2*theta**(-7.d0)-two*a1*theta

    !---calculate derivatives for water density
    if (calculate_derivatives) then
       zpt = ypt+(a3*yy*ypt-a4)/xx
       zpp = a5/xx
       u9 = u0*u1/zz
       vrpt = u9*zpt+aa(13)+two*aa(14)*theta-ten*u8 &
            -19.d0*aa(16)*u4*u4*theta18+11.d0*u2*u2*u3*theta**10.d0 &
            -aa(20)*u6*(18.d0*a9*theta18+20.d0*theta20)/theta &
            -(three*aa(21)+80.d0*aa(22)*beta/(theta20*theta))*beta2x

       vrpp = u9*zpp-u2*(two*aa(18)+six*aa(19)*beta)-12.d0*u7*u5/ &
            (a10+beta)+(six*aa(21)*(a12-theta)+12.d0*aa(22)*beta/ &
            theta20)*beta

       cnv = -one/(vc1mol*vr*vr)
       dwt = cnv*vrpt*utc1 ! kmol/m^3/C
       dwp = cnv*vrpp*upc1 ! kmol/m^3/Pa

       ! 2nd derivative w.r.t pressure
       d2z_dp2 = -a5*a5/xx/xx/xx

       d2vr_dp2 = u9*(u0-one)/zz*zpp*zpp + &
            u0*u1/zz*d2z_dp2 + &
            six*u2*aa(19) + &
            60.d0*u7*u5/(a10+beta)/(a10+beta) + &
            six*(a12 - theta) + &
            24.d0*aa(22)*beta/theta20

       d2wmol_dp2 = -cnv*upc1*upc1*(2.d0/vr*vrpp*vrpp -  d2vr_dp2) ! kmol/m^3/Pa^2

    else
       dwt = 0.d0
       dwp = 0.d0
       d2wmol_dp2 = 0.d0
    endif

  end subroutine DensityIFC67

  !------------------------------------------------------------------------
  subroutine EnthalpyIFC67(t,p,calculate_derivatives,hw, &
       hwp,hwt,ierr)

    !  This subroutine calculates water and steam-gas mixture properties.
    !  The water and steam properties are valid in the range of:
    !
    !            0 < p < 165.4 * 10^5 pascals (165.4 bars)
    !            0 < t < 350 centigrade (623.15 Kelvin)
    !
    !  The properties cover densities, enthalpies, internal energies,
    !  and partial derivatives of these quanties with respect to
    !  pressure and temperature.
    !
    !  For saturated fluid, it will also calculate water saturation
    !  temperature for the specified pressure using Newton-Raphson and
    !  the derivative dts/dp (=tsp) or Ps for a given temperature.
    !
    !  Ref.: International Formulation Committee of the Sixth International
    !       Conference on Properties of Steam (1967).

    implicit none
    !
    ! !ARGUMENTS
    PetscReal      , intent(in)  :: t   ! Temperature in centigrade
    PetscReal      , intent(in)  :: p   ! Pressure in Pascals
    PetscBool      , intent(in)  :: calculate_derivatives
    PetscReal      , intent(out) :: hw,hwp,hwt
    PetscErrorCode , intent(out) :: ierr
    !
    ! !LOCAL VARIABLES:
    PetscInt :: i

    PetscReal, save      :: aa(0:22)
    PetscReal, save      :: a1,a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12

    PetscReal            :: beta,beta2x,beta4,theta,utheta,theta2x,theta18,theta20
    PetscReal            :: xx,yy,zz
    PetscReal            :: u0,u1
    PetscReal            :: tempreal
    PetscReal            :: v0_1, v1_1, v2_1, v3_1, v4_1
    PetscReal            :: v1_2, v2_2, v3_2, v4_2, v20_2, v40_2
    PetscReal            :: v1_3, v2_3, v3_3, v4_3
    PetscReal            :: v1_4, v2_4, v3_4
    PetscReal            :: v1_5, v2_5
    PetscReal            :: v1_6
    PetscReal            :: term1,term2,term2t,term3,term3t,term3p,term4,term4t,term4p, &
                            term5,term5t,term5p,term6,term6t,term6p,term7,term7t,term7p
    PetscReal            :: dv2t,dv2p,dv3t
    PetscReal            :: ypt,yptt,zpt,zpp
    PetscReal            :: tc1,pc1,vc1,utc1,upc1,vc1mol
    PetscReal, parameter :: zero = 0.d0
    PetscReal, parameter :: one = 1.d0
    PetscReal, parameter :: two = 2.d0
    PetscReal, parameter :: three = 3.d0
    PetscReal, parameter :: four = 4.d0
    PetscReal, parameter :: five = 5.d0
    PetscReal, parameter :: six = 6.d0
    PetscReal, parameter :: seven = 7.d0
    PetscReal, parameter :: eight = 8.d0
    PetscReal, parameter :: nine = 9.d0
    PetscReal, parameter :: ten = 10.d0

    data aa/ &
         !-----data aa0,aa1,aa2,aa3/
         6.824687741d03,-5.422063673d02,-2.096666205d04, 3.941286787d04, &
         !-----data aa4,aa5,aa6,aa7/
         -6.733277739d04, 9.902381028d04,-1.093911774d05, 8.590841667d04, &
         !-----data aa8,aa9,aa10,aa11/
         -4.511168742d04, 1.418138926d04,-2.017271113d03, 7.982692717d00, &
         !-----data aa12,aa13,aa14,aa15/
         -2.616571843d-2, 1.522411790d-3, 2.284279054d-2, 2.421647003d02, &
         !-----data aa16,aa17,aa18,aa19/
         1.269716088d-10,2.074838328d-7, 2.174020350d-8, 1.105710498d-9, &
         !-----data aa20,aa21,aa22/
         1.293441934d01, 1.308119072d-5, 6.047626338d-14/

    data a1,a2,a3,a4/ &
         8.438375405d-1, 5.362162162d-4, 1.720000000d00, 7.342278489d-2/
    data a5,a6,a7,a8/ &
         4.975858870d-2, 6.537154300d-1, 1.150000000d-6, 1.510800000d-5/
    data a9,a10,a11,a12/ &
         1.418800000d-1, 7.002753165d00, 2.995284926d-4, 2.040000000d-1/

    ierr = 0

    tc1     = H2O_CRITICAL_TEMPERATURE ! K
    pc1     = H2O_CRITICAL_PRESSURE    ! Pa
    vc1     = 0.00317d0                ! m^3/kg
    utc1    = one/tc1                  ! 1/C
    upc1    = one/pc1                  ! 1/Pa
    vc1mol  = vc1*FMWH2O               ! m^3/kmol

    theta   = (t+273.15d0)*utc1
    theta2x = theta*theta
    theta18 = theta**18.d0
    theta20 = theta18*theta2x

    beta    = p*upc1
    beta2x  = beta*beta
    beta4   = beta2x*beta2x

    yy = one-a1*theta2x-a2*theta**(-6.d0)
    xx = a3*yy*yy-two*(a4*theta-a5*beta)
    !   Note: xx may become negative near the critical point-pcl.
    if (xx.gt.zero) then
       xx = sqrt(xx)
    else
       write(iulog,*)'Negative term in density ','t= ',t,' p= ',p,' xx= ',xx
       call endrun(msg=errMsg(__FILE__, __LINE__))
       ierr = 1
       xx = 1.e-6               !set arbitrarily
    end if
    zz = yy + xx
    u0 = -five/17.d0
    u1 = aa(11)*a5*zz**u0

    ! ypt used for enthalpy even if derivative not calculated
    ypt = six*a2*theta**(-7.d0)-two*a1*theta


    !---compute enthalpy internal energy and derivatives for water
    utheta = one/theta
    term1 = aa(0)*theta
    term2 = -aa(1)
    ! term2t is part of the derivative calc., but left here to avoid
    ! recomputing the expensive do loop
    term2t = zero
    do i = 3,10
       tempreal = dfloat(i-2)*aa(i)*theta**(i-1)
       term2t = term2t+tempreal*utheta*dfloat(i-1)
       term2 = term2+tempreal
    end do

    ! "v" section 1
    v0_1 = u1/a5
    v2_1 = 17.d0*(zz/29.d0-yy/12.d0)+five*theta*ypt/12.d0
    v3_1 = a4*theta-(a3-one)*theta*yy*ypt
    v1_1 = zz*v2_1+v3_1
    term3 = v0_1*v1_1

    ! block 1 removed from here

    ! "v" section 2
    v1_2 = nine*theta+a6
    v20_2 = (a6-theta)
    v2_2 = v20_2**9.d0
    v3_2 = a7+20.d0*theta**19.d0
    v40_2 = a7+theta**19.d0
    v4_2 = one/(v40_2*v40_2)
    ! term4p is a derivative, but left due to dependency in term4
    term4p = aa(12)-aa(14)*theta2x+aa(15)*v1_2*v2_2+aa(16)*v3_2*v4_2
    term4 = term4p*beta

    ! block 2 removed from here

    ! "v" section 3
    v1_3 = beta*(aa(17)+aa(18)*beta+aa(19)*beta2x)
    v2_3 = 12.d0*theta**11.d0+a8
    v4_3 = one/(a8+theta**11.d0)
    v3_3 = v4_3*v4_3
    term5 = v1_3*v2_3*v3_3

    ! block 3 removed from here

    ! "v" section 4
    v1_4 = (a10+beta)**(-3.d0)+a11*beta
    v3_4 = (17.d0*a9+19.d0*theta2x)
    v2_4 = aa(20)*theta18*v3_4
    term6 = v1_4*v2_4

    ! block 4 removed from here

    ! "v" section 5
    v1_5 = 21.d0*aa(22)/theta20*beta4
    v2_5 = aa(21)*a12*beta2x*beta
    term7 = v1_5+v2_5

    ! "v" section 6
    v1_6 = pc1*vc1mol
    hw = (term1-term2+term3+term4-term5+term6+term7)*v1_6

    if (calculate_derivatives) then

       zpt = ypt+(a3*yy*ypt-a4)/xx
       zpp = a5/xx

       ! block 1
       yptt = -two*a1-42.d0*a2/theta**8.d0
       dv2t = 17.d0*(zpt/29.d0-ypt/12.d0)+five/12.d0*(ypt+theta*yptt)
       dv3t = a4-(a3-one)*(theta*yy*yptt+yy*ypt+theta*ypt*ypt)
       dv2p = 17.d0*zpp/29.d0
       v4_1 = five*v1_1/(17.d0*zz)
       term3t = v0_1*(zz*dv2t+(v2_1-v4_1)*zpt+dv3t)
       term3p = v0_1*(zz*dv2p+(v2_1-v4_1)*zpp)

       ! block 2
       term4t = (-two*aa(14)*theta+nine*aa(15)*(v2_2-v1_2*v2_2/v20_2) &
            +38.d0*theta18*aa(16)*(ten*v4_2-v3_2*v4_2/v40_2))*beta

       ! block 3
       term5p = v3_3*v2_3*(aa(17)+two*aa(18)*beta+three*aa(19)*beta2x)
       term5t = v1_3*(132.d0*v3_3*theta**10.d0-22.d0*v2_3*v3_3*v4_3*theta**10.d0)

       ! block 4
       term6p = v2_4*(a11-three*(a10+beta)**(-4.d0))
       term6t = v1_4*aa(20)*theta18*(18.d0*v3_4*utheta+38.d0*theta)

       ! block 5
       term7p = beta2x*(three*aa(21)*a12+84.d0*aa(22)*beta/theta20)
       term7t = -420.d0*aa(22)*beta4/(theta20*theta)

       hwp = (term3p+term4p-term5p+term6p+term7p)*vc1mol
       hwt = (aa(0)-term2t+term3t+term4t-term5t+term6t+term7t)*v1_6*utc1
    else
       hwp = 0.d0
       hwt = 0.d0
    endif

  end subroutine EnthalpyIFC67

  !------------------------------------------------------------------------
  subroutine Viscosity(p, t_K, vis, dvis_dp, dvis_dT)
    !
    ! !DESCRIPTION:
    ! Return viscosity of water
    !
    implicit none
    !
    ! !ARGUMENTS    
    PetscReal, intent(in)  :: p
    PetscReal, intent(in)  :: t_K
    PetscReal, intent(out) :: vis
    PetscReal, intent(out) :: dvis_dp
    PetscReal, intent(out) :: dvis_dT

    vis     = 8.904156d-4 ! [Pa s]
    dvis_dp = 0.d0
    dvis_dT = 0.d0

  end subroutine Viscosity

  !------------------------------------------------------------------------
  subroutine InternalEnergyAndEnthalpy(P, t_K, itype, den, dden_dT, dden_dP, &
       U, H, dU_dT, dH_dT, dU_dP, dH_dP)
    !
    ! !DESCRIPTION:
    ! Return internal energy of water
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: P        ! [Pa]
    PetscReal, intent(in)  :: T_K      ! [K]
    PetscInt               :: itype    ! [-]
    PetscReal, intent(in)  :: den      ! [kg m^{-3}]
    PetscReal, intent(in)  :: dden_dT  ! [kg m^{-3} K^{-1}]
    PetscReal, intent(in)  :: dden_dP  ! [kg m^{-3} Pa^{-1}]
    PetscReal, intent(out) :: U        ! [J kmol^{-1}]
    PetscReal, intent(out) :: H        ! [J kmol^{-1}]
    PetscReal, intent(out) :: dU_dT    ! [J kmol^{-1} K^{-1}]
    PetscReal, intent(out) :: dH_dT    ! [J kmol^{-1} K^{-1}]
    PetscReal, intent(out) :: dU_dP    ! [J kmol^{-1} Pa^{-1}]
    PetscReal, intent(out) :: dH_dP    ! [J kmol^{-1} Pa^{-1}]

    select case(itype)
    case (INT_ENERGY_ENTHALPY_CONSTANT)
       call InternalEnergyAndEnthalpyConstant(P, t_K, den, dden_dT, dden_dP, &
            U, H, dU_dT, dH_dT, dU_dP, dH_dP)

    case (INT_ENERGY_ENTHALPY_IFC67)
       call InternalEnergyAndEnthalpyIFC67(P, t_K, den, dden_dT, dden_dP, &
            U, H, dU_dT, dH_dT, dU_dP, dH_dP)

    case default
       write(iulog,*)'Density: Unknown denity_itype. '
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine InternalEnergyAndEnthalpy

  !------------------------------------------------------------------------
  subroutine InternalEnergyAndEnthalpyConstant(P, t_K, den, dden_dT, dden_dP, &
       U, H, dU_dT, dH_dT, dU_dP, dH_dP)
    !
    ! !DESCRIPTION:
    ! Return internal energy of water
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: P                      ! [Pa]
    PetscReal, intent(in)  :: T_K                    ! [K]
    PetscReal, intent(in)  :: den                    ! [kg m^{-3}]
    PetscReal, intent(in)  :: dden_dT                ! [kg m^{-3} K^{-1}]
    PetscReal, intent(in)  :: dden_dP                ! [kg m^{-3} Pa^{-1}]
    PetscReal, intent(out) :: U                      ! [J kmol^{-1}]
    PetscReal, intent(out) :: H                      ! [J kmol^{-1}]
    PetscReal, intent(out) :: dU_dT                  ! [J kmol^{-1} K^{-1}]
    PetscReal, intent(out) :: dH_dT                  ! [J kmol^{-1} K^{-1}]
    PetscReal, intent(out) :: dU_dP                  ! [J kmol^{-1} Pa^{-1}]
    PetscReal, intent(out) :: dH_dP                  ! [J kmol^{-1} Pa^{-1}]
                                                     !
    PetscReal, parameter   :: u0 = 4.217 * 1.d3      ! [J kg^{-1} K^{-1}]
    PetscReal              :: T_C
    PetscBool              :: cal_deriv
    PetscErrorCode         :: ierr

    U     = u0 * (t_K - 273.15d0)                    ! [J kg^{-1}]
    dU_dT = u0                                       ! [J kg^{-1} K^{-1}]
    dU_dP = 0.d0                                     ! [J kg^{-1} Pa^{-1}]

    H     = U + P/den                                ! [J kg^{-1}]
    dH_dT = dU_dT - P/(den**2.d0)*dden_dT            ! [J kg^{-1} K^{-1}]
    dH_dP = dU_dP + 1.d0/den - P/(den**2.d0)*dden_dP ! [J kg^{-1} Pa^{-1}]

    U     = U * FMWH2O                               ! [J kmol^{-1}]
    H     = H * FMWH2O                               ! [J kmol^{-1}]
    dU_dT = dU_dT * FMWH2O                           ! [J kmol^{-1} K^{-1}]
    dH_dT = dH_dT * FMWH2O                           ! [J kmol^{-1} K^{-1}]
    dH_dP = dH_dP * FMWH2O                           ! [J kmol^{-1} Pa^{-1}]

  end subroutine InternalEnergyAndEnthalpyConstant

  !------------------------------------------------------------------------
  subroutine InternalEnergyAndEnthalpyIFC67(P, t_K, den, dden_dT, dden_dP, &
       U, H, dU_dT, dH_dT, dU_dP, dH_dP)
    !
    ! !DESCRIPTION:
    ! Return internal energy of water
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: P        ! [Pa]
    PetscReal, intent(in)  :: T_K      ! [K]
    PetscReal, intent(in)  :: den      ! [kg m^{-3}]
    PetscReal, intent(in)  :: dden_dT  ! [kg m^{-3} K^{-1}]
    PetscReal, intent(in)  :: dden_dP  ! [kg m^{-3} Pa^{-1}]
    PetscReal, intent(out) :: U        ! [J kmol^{-1}]
    PetscReal, intent(out) :: H        ! [J kmol^{-1}]
    PetscReal, intent(out) :: dU_dT    ! [J kmol^{-1} K^{-1}]
    PetscReal, intent(out) :: dH_dT    ! [J kmol^{-1} K^{-1}]
    PetscReal, intent(out) :: dU_dP    ! [J kmol^{-1} Pa^{-1}]
    PetscReal, intent(out) :: dH_dP    ! [J kmol^{-1} Pa^{-1}]
    !
    PetscReal, parameter   :: u0 = 4.217 * 1.d3 ! [J kg^{-1} K^{-1}]
    PetscReal              :: T_C
    PetscBool              :: cal_deriv
    PetscErrorCode         :: ierr

    T_C = T_K - 273.15d0
    cal_deriv = PETSC_TRUE

    call EnthalpyIFC67(T_C, P, cal_deriv, H, dH_dP, dH_dT, ierr)

    U = H - P / (den/FMWH2O)

    dU_dT = dH_dT + P/((den/FMWH2O)**2.d0)*(dden_dT/FMWH2O)
    dU_dP = dH_dP - 1.d0/(den/FMWH2O) + P/((den/FMWH2O)**2.d0)*(dden_dP/FMWH2O)

  end subroutine InternalEnergyAndEnthalpyIFC67

#endif

end module EOSWaterMod
