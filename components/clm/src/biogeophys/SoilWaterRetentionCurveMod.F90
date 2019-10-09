module SoilWaterRetentionCurveMod

  !---------------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Abstract base class for functions to compute soil water retention curve
  !
  ! !USES:
  implicit none
  save
  private
  !
  ! !PUBLIC TYPES:
  public :: soil_water_retention_curve_type

  type, abstract :: soil_water_retention_curve_type
     private
   contains
     ! compute hydraulic conductivity
     procedure(soil_hk_interface), deferred :: soil_hk

     ! compute soil suction potential
     procedure(soil_suction_interface), deferred :: soil_suction

     ! compute relative saturation at which soil suction is equal to a target value
     procedure(soil_suction_inverse_interface), deferred :: soil_suction_inverse
  end type soil_water_retention_curve_type

  abstract interface

     ! Note: The following interfaces are set up based on the arguments needed for the
     ! clapphornberg1978 implementations. It's likely that these interfaces are not
     ! totally general for all desired implementations. In that case, we'll need to think
     ! about how to support different interfaces. Some possible solutions are:
     !
     ! - Make the interfaces contain all possible inputs that are needed by any
     !   implementation; each implementation will then ignore the inputs it doesn't need.
     !
     ! - For inputs that are needed only by particular implementations - and particularly
     !   for inputs that are constant in time (e.g., this is the case for bsw, I think):
     !   pass these into the constructor, and save pointers to these inputs as components
     !   of the child type that needs them. Then they aren't needed as inputs to the
     !   individual routines, allowing the interfaces for these routines to remain more
     !   consistent between different implementations.

     subroutine soil_hk_interface(this, hksat, imped, s, bsw, hk, dhkds)
       ! !DESCRIPTION:
       ! Compute hydraulic conductivity
       !
       ! !USES:
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: soil_water_retention_curve_type
       !
       ! !ARGUMENTS:
       class(soil_water_retention_curve_type), intent(in) :: this
       real(r8), intent(in) :: hksat    !saturated hydraulic conductivity [mm/s]
       real(r8), intent(in) :: imped    !ice impedance
       real(r8), intent(in) :: s        !reletive saturation, [0, 1]
       real(r8), intent(in) :: bsw      !shape parameter
       real(r8), intent(out):: hk       !hydraulic conductivity [mm/s]
       real(r8), optional, intent(out):: dhkds    !d[hk]/ds   [mm/s]
     end subroutine soil_hk_interface


     subroutine soil_suction_interface(this, smpsat, s, bsw, smp, dsmpds)
       ! !DESCRIPTION:
       ! Compute soil suction potential
       !
       ! !USES:
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: soil_water_retention_curve_type
       !
       ! !ARGUMENTS:
       class(soil_water_retention_curve_type), intent(in) :: this
       real(r8), intent(in)            :: smpsat   !minimum soil suction, positive [mm]
       real(r8), intent(in)            :: s        !reletive saturation, [0, 1]
       real(r8), intent(in)            :: bsw      !shape parameter
       real(r8), intent(out)           :: smp      !soil suction, negative, [mm]
       real(r8), optional, intent(out) :: dsmpds   !d[smp]/ds, [mm]
     end subroutine soil_suction_interface

     subroutine soil_suction_inverse_interface(this, smp_target, smpsat, bsw, s_target)
       ! !DESCRIPTION:
       ! Compute relative saturation at which soil suction is equal to a target value.
       ! This is done by inverting the soil_suction equation to solve for s.
       !
       ! !USES:
       use shr_kind_mod   , only : r8 => shr_kind_r8
       import :: soil_water_retention_curve_type
       !
       ! !ARGUMENTS:
       class(soil_water_retention_curve_type), intent(in) :: this
       real(r8) , intent(in)  :: smp_target ! target soil suction, negative [mm]
       real(r8) , intent(in)  :: smpsat     ! minimum soil suction, positive [mm]
       real(r8) , intent(in)  :: bsw        ! shape parameter
       real(r8) , intent(out) :: s_target   ! relative saturation at which smp = smp_target [0,1]
     end subroutine soil_suction_inverse_interface

  end interface

end module SoilWaterRetentionCurveMod
