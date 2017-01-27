#ifdef USE_PETSC_LIB


module EOSWaterMod


  ! !USES:
  use clm_varctl         , only : iulog
  use abortutils         , only : endrun
  use shr_log_mod        , only : errMsg => shr_log_errMsg
  !
  ! !PUBLIC TYPES:
  implicit none
  private

#include "finclude/petscsys.h"

  PetscInt, parameter, public :: DENSITY_CONSTANT  = 1
  PetscInt, parameter, public :: DENSITY_TGDPB01   = 2

  public :: Density
  public :: Viscosity

  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine Density(p, t_K, denity_itype, den, dden_dp)
    !
    ! !DESCRIPTION:
    ! Given pressure, temperature and type of density formulation, compute:
    ! - density of water, and
    ! - first derivative of density w.r.t pressure
    !
    ! !USES:
    use MultiPhysicsProbConstants, only : FMWH2O
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)   :: p
    PetscReal, intent(in)   :: t_K
    PetscInt, intent(in)    :: denity_itype
    PetscReal, intent(out)  :: den
    PetscReal, intent(out)  :: dden_dp

    select case(denity_itype)
    case (DENSITY_CONSTANT)
       call DensityConstant(p, t_K, den, dden_dp)
    case (DENSITY_TGDPB01)
       call DensityTGDPB01 (p , t_K, den, dden_dp)
    case default
       write(iulog,*)'Density: Unknown denity_itype. '
       call endrun(msg=errMsg(__FILE__, __LINE__))
    end select

  end subroutine Density

  !------------------------------------------------------------------------
  subroutine DensityConstant(p, t_K, den, dden_dp)
    !
    ! !DESCRIPTION:
    ! Return constant density of water
    !
    ! !USES:
    use MultiPhysicsProbConstants , only : FMWH2O
    use clm_varcon                , only : denh2o
    !
    implicit none
    !
    ! !ARGUMENTS    
    PetscReal, intent(in)  :: p         ! [Pa]
    PetscReal, intent(in)  :: t_K       ! [K]
    PetscReal, intent(out) :: den       ! [kmol m^{-3}]
    PetscReal, intent(out) :: dden_dp   ! [kmol m^{-3} Pa^{-1}]

    den     = denh2o/FMWH2O ! [kmol m^{-3}]
    dden_dp = 0.d0

  end subroutine DensityConstant

  !------------------------------------------------------------------------
  subroutine DensityTGDPB01(p, t_K, den, dden_dp)
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
    use MultiPhysicsProbConstants, only : FMWH2O
    !
    implicit none
    !
    ! !ARGUMENTS
    PetscReal, intent(in)  :: p                ! [Pa]
    PetscReal, intent(in)  :: t_K              ! [K]
    PetscReal, intent(out) :: den              ! [kmol m^{-3}]
    PetscReal, intent(out) :: dden_dp          ! [kmol m^{-3} Pa^{-1}]

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
    PetscReal              :: dden_dt
    PetscReal              :: dp
    PetscBool              :: unsaturated

    t_c = t_K - 273.15d0

    ! Density of water as function of t_K
    dent = a5*(1.d0 - ((t_c + a1)**2.d0)*(t_c + a2)/a3/(t_c + a4))

    unsaturated = PETSC_FALSE
    dp = (p - p0)

    if (p <= p0) then
       unsaturated = PETSC_TRUE
       dp = 0.d0
    endif

    ! Compressibility of water
    kappa = (1.d0 + (k0 + k1*t_c + k2*t_c**2.d0)*dp)

    ! Density of water
    den = dent*kappa/FMWH2O ! [kmol m^{-3}]

    ! Derivative
    ddent_dp   = 0.d0
    ddent_dt_1 = -((t_c + a1)**2.d0)/a3/(t_c + a4)
    ddent_dt_2 = -2.d0*(t_c + a1)*(t_c + a2)/a3/(t_c + a4)
    ddent_dt_3 =  ((t_c + a1)**2.d0)*(t_c + a2)/a3/((t_c + a4)**2.d0)
    ddent_dt   = a5*(ddent_dt_1 + ddent_dt_2 + ddent_dt_3)

    dkappa_dp  = (k0 + k1*t_c + k2*t_c**2.d0)
    dkappa_dt  = (k1 + 2.d0*k2*t_c)*dp

    dden_dt    = (ddent_dt*kappa + dent*dkappa_dt)/FMWH2O
    dden_dp    = (ddent_dp*kappa + dent*dkappa_dp)/FMWH2O

    if (unsaturated) then
       dden_dp = 0.d0
    endif

  end subroutine DensityTGDPB01

  !------------------------------------------------------------------------
  subroutine Viscosity(p, t_K, vis, dvis_dp)
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

    vis     = 8.904156d-4 ! [Pa s]
    dvis_dp = 0.d0

  end subroutine Viscosity

end module EOSWaterMod
#endif
