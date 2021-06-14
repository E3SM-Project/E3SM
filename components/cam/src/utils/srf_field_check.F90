module srf_field_check

  ! Utility functions called from the coupling layer to allow CAM to query
  ! whether specific fields are being provided by the coupler.  This provides
  ! a coupler independent way for CAM to make these queries.

  implicit none
  private

  ! Input to atm
  logical, public, protected :: active_Sl_ram1      = .false.
  logical, public, protected :: active_Sl_fv        = .false.
  logical, public, protected :: active_Sl_soilw     = .false.
  logical, public, protected :: active_Fall_flxdst1 = .false.
  logical, public, protected :: active_Fall_flxvoc  = .false.
  logical, public, protected :: active_Fall_flxfire = .false.
  logical, public, protected :: active_Fall_fco2_lnd = .false.
  logical, public, protected :: active_Faoo_fco2_ocn = .false.

  ! output from atm
  logical, public, protected :: active_Faxa_nhx = .false.
  logical, public, protected :: active_Faxa_noy = .false.

  public :: set_active_Sl_ram1
  public :: set_active_Sl_fv
  public :: set_active_Sl_soilw
  public :: set_active_Fall_flxdst1
  public :: set_active_Fall_flxvoc
  public :: set_active_Fall_flxfire
  public :: set_active_Fall_fco2_lnd
  public :: set_active_Faoo_fco2_ocn
  public :: set_active_Faxa_nhx
  public :: set_active_Faxa_noy

!===============================================================================
contains
!===============================================================================

  subroutine set_active_Sl_ram1(is_active)
    logical, intent(in) :: is_active
    active_Sl_ram1 = is_active
  end subroutine set_active_Sl_ram1

  subroutine set_active_Sl_fv(is_active)
    logical, intent(in) :: is_active
    active_Sl_fv = is_active
  end subroutine set_active_Sl_fv

  subroutine set_active_Sl_soilw(is_active)
    logical, intent(in) :: is_active
    active_Sl_soilw = is_active
  end subroutine set_active_Sl_soilw

  subroutine set_active_Fall_flxdst1(is_active)
    logical, intent(in) :: is_active
    active_Fall_flxdst1 = is_active
  end subroutine set_active_Fall_flxdst1

  subroutine set_active_Fall_flxvoc(is_active)
    logical, intent(in) :: is_active
    active_Fall_flxvoc = is_active
  end subroutine set_active_Fall_flxvoc

  subroutine set_active_Fall_flxfire(is_active)
    logical, intent(in) :: is_active
    active_Fall_flxfire = is_active
  end subroutine set_active_Fall_flxfire

  subroutine set_active_Fall_fco2_lnd(is_active)
    logical, intent(in) :: is_active
    active_Fall_fco2_lnd = is_active
  end subroutine set_active_Fall_fco2_lnd

  subroutine set_active_Faoo_fco2_ocn(is_active)
    logical, intent(in) :: is_active
    active_Faoo_fco2_ocn = is_active
  end subroutine set_active_Faoo_fco2_ocn

  subroutine set_active_Faxa_nhx(is_active)
    logical, intent(in) :: is_active
    active_Faxa_nhx = is_active
  end subroutine set_active_Faxa_nhx

  subroutine set_active_Faxa_noy(is_active)
    logical, intent(in) :: is_active
    active_Faxa_noy = is_active
  end subroutine set_active_Faxa_noy

end module srf_field_check
