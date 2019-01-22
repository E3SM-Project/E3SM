module shr_precip_mod

  ! This module contains methods for manipulating precipitation quantities

  use shr_kind_mod, only : r8 => SHR_KIND_R8

  implicit none
  private
  save

  ! determine a rain-snow partitioning using a ramp method based on temperature
  public :: shr_precip_partition_rain_snow_ramp

contains

  !-----------------------------------------------------------------------
  subroutine shr_precip_partition_rain_snow_ramp(temperature, frac_rain)
    !
    ! !DESCRIPTION:
    ! Determine a rain-snow partitioning using a ramp method based on temperature.
    !
    ! Returns fractional mass of precipitation falling as rain. The rest (1 - frac_rain)
    ! falls as snow.
    !
    ! This is meant to be used for precipitation at the surface, e.g., to force CLM.
    !
    ! !USES:
    use shr_const_mod, only : SHR_CONST_TKFRZ
    !
    ! !ARGUMENTS:
    real(r8), intent(in)  :: temperature  ! temperature (K)
    real(r8), intent(out) :: frac_rain    ! fraction of precipitation falling as rain
    !
    ! !LOCAL VARIABLES:

    character(len=*), parameter :: subname = 'shr_precip_partition_rain_snow_ramp'
    !-----------------------------------------------------------------------

    ! ramp near freezing
    frac_rain = (temperature - SHR_CONST_TKFRZ) * 0.5_r8

    ! bound in [0,1]
    frac_rain = min(1.0_r8,max(0.0_r8,frac_rain))

  end subroutine shr_precip_partition_rain_snow_ramp

end module shr_precip_mod
