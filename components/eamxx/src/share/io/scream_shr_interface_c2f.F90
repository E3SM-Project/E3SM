module eamxx_shr_interface_mod
  use iso_c_binding, only: c_int
  implicit none

!
! This file contains bridges from EAMxx to E3SM shr module
!
contains

!=====================================================================!
  function shr_get_iosysid_c2f(atm_id) result (iosysid) bind(c)
    use pio, only: iosystem_desc_t
    use shr_pio_mod, only: shr_pio_getiosys
    integer(kind=c_int), value, intent(in) :: atm_id

    integer (kind=c_int) :: iosysid
    type(iosystem_desc_t), pointer :: iosystem

    iosystem => shr_pio_getiosys(atm_id)
    iosysid = iosystem%iosysid
  end function shr_get_iosysid_c2f

!=====================================================================!
  function shr_get_iotype_c2f(atm_id) result (iotype) bind(c)
    use shr_pio_mod, only: shr_pio_getiotype
    integer(kind=c_int), value, intent(in) :: atm_id

    integer (kind=c_int) :: iotype

    iotype = shr_pio_getiotype(atm_id)
  end function shr_get_iotype_c2f

!=====================================================================!
  function shr_get_rearranger_c2f(atm_id) result (rearr) bind(c)
    use shr_pio_mod, only: shr_pio_getrearranger
    integer(kind=c_int), value, intent(in) :: atm_id

    integer (kind=c_int) :: rearr

    rearr = shr_pio_getrearranger(atm_id)
  end function shr_get_rearranger_c2f

!=====================================================================!
  function shr_get_ioformat_c2f(atm_id) result (ioformat) bind(c)
    use shr_pio_mod, only: shr_pio_getioformat
    integer(kind=c_int), value, intent(in) :: atm_id

    integer (kind=c_int) :: ioformat

    ioformat = shr_pio_getioformat(atm_id)
  end function shr_get_ioformat_c2f

end module eamxx_shr_interface_mod
