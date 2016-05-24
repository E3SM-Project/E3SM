module SoilWaterRetentionCurveClappHornberg1978Mod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Implementation of soil_water_retention_curve_type using the Clapp-Hornberg 1978
  ! parameterizations.
  !
  ! !USES:
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use SoilWaterRetentionCurveMod, only : soil_water_retention_curve_type
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: soil_water_retention_curve_clapp_hornberg_1978_type
  
  type, extends(soil_water_retention_curve_type) :: &
       soil_water_retention_curve_clapp_hornberg_1978_type
     private
   contains
     procedure :: soil_hk              ! compute hydraulic conductivity
     procedure :: soil_suction         ! compute soil suction potential
     procedure :: soil_suction_inverse ! compute relative saturation at which soil suction is equal to a target value
  end type soil_water_retention_curve_clapp_hornberg_1978_type

  interface soil_water_retention_curve_clapp_hornberg_1978_type
     ! initialize a new soil_water_retention_curve_clapp_hornberg_1978_type object
     module procedure constructor  
  end interface soil_water_retention_curve_clapp_hornberg_1978_type

contains

  !-----------------------------------------------------------------------
  type(soil_water_retention_curve_clapp_hornberg_1978_type) function constructor()
    !
    ! !DESCRIPTION:
    ! Creates an object of type soil_water_retention_curve_clapp_hornberg_1978_type.
    ! For now, this is simply a place-holder.
    !-----------------------------------------------------------------------

  end function constructor

  !-----------------------------------------------------------------------
  subroutine soil_hk(this, hksat, imped, s, bsw, hk, dhkds)
    !
    ! !DESCRIPTION:
    ! Compute hydraulic conductivity
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_clapp_hornberg_1978_type), intent(in) :: this
    real(r8), intent(in) :: hksat    !saturated hydraulic conductivity [mm/s]
    real(r8), intent(in) :: imped    !ice impedance
    real(r8), intent(in) :: s        !reletive saturation, [0, 1]
    real(r8), intent(in) :: bsw      !shape parameter
    real(r8), intent(out):: hk       !hydraulic conductivity [mm/s]
    real(r8), optional, intent(out):: dhkds    !d[hk]/ds   [mm/s]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_hk'
    !-----------------------------------------------------------------------
    
    !compute hydraulic conductivity
    hk=imped*hksat*s**(2._r8*bsw+3._r8)

    !compute the derivative
    if(present(dhkds))then
       dhkds=(2._r8*bsw+3._r8)*hk/s
    endif

  end subroutine soil_hk

  !-----------------------------------------------------------------------
  subroutine soil_suction(this, smpsat, s, bsw, smp, dsmpds)
    !
    ! !DESCRIPTION:
    ! Compute soil suction potential
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_clapp_hornberg_1978_type), intent(in) :: this
    real(r8), intent(in)            :: smpsat   !minimum soil suction, positive [mm]
    real(r8), intent(in)            :: s        !reletive saturation, [0, 1]
    real(r8), intent(in)            :: bsw      !shape parameter
    real(r8), intent(out)           :: smp      !soil suction, negative, [mm]
    real(r8), optional, intent(out) :: dsmpds   !d[smp]/ds, [mm]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_suction'
    !-----------------------------------------------------------------------
    
    !compute soil suction potential, negative
    smp = -smpsat*s**(-bsw)

    !compute derivative
    if(present(dsmpds))then
       dsmpds=-bsw*smp/s
    endif

  end subroutine soil_suction

  !-----------------------------------------------------------------------
  subroutine soil_suction_inverse(this, smp_target, smpsat, bsw, s_target)
    !
    ! !DESCRIPTION:
    ! Compute relative saturation at which soil suction is equal to a target value.
    ! This is done by inverting the soil_suction equation to solve for s.
    !
    ! !USES:
    !
    ! !ARGUMENTS:
    class(soil_water_retention_curve_clapp_hornberg_1978_type), intent(in) :: this
    real(r8) , intent(in)  :: smp_target ! target soil suction, negative [mm]
    real(r8) , intent(in)  :: smpsat     ! minimum soil suction, positive [mm]
    real(r8) , intent(in)  :: bsw        ! shape parameter
    real(r8) , intent(out) :: s_target   ! relative saturation at which smp = smp_target [0,1]
    !
    ! !LOCAL VARIABLES:
    
    character(len=*), parameter :: subname = 'soil_suction_inverse'
    !-----------------------------------------------------------------------
    
    s_target = (-smp_target/smpsat)**(-1/bsw)

  end subroutine soil_suction_inverse


end module SoilWaterRetentionCurveClappHornberg1978Mod

