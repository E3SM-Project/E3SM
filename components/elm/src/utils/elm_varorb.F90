
module elm_varorb

  use shr_kind_mod , only: r8 => shr_kind_r8
  use seq_infodata_mod, only: seq_infodata_type, seq_infodata_SolarDecl
  use seq_infodata_mod, only: seq_infodata_GetData, seq_infodata_PutData, seq_infodata_UseEquationOfTime
  use shr_sys_mod, only: shr_sys_abort
  implicit none
  
  ! Orbital information needed as input to orbit_parms

  real(r8) :: eccen    ! Earth's eccentricity factor (unitless) (typically 0 to 0.1)

  ! Orbital information after processed by orbit_params

  real(r8) :: obliqr  ! Earth's obliquity in radians
  real(r8) :: lambm0  ! Mean longitude of perihelion at the vernal equinox (radians)
  real(r8) :: mvelpp  ! Earth's moving vernal equinox longitude of perihelion plus pi (radians)

  type(seq_infodata_type), private, save :: solar_context
  logical, private, save :: solar_context_initialized = .false.

contains

  subroutine elm_set_orbital_context(infodata)
    ! Refresh before forcing import, independently of EOT, so copied orbital policy cannot become stale.
    type(seq_infodata_type), intent(in) :: infodata
    solar_context = infodata
    solar_context_initialized = .true.
  end subroutine elm_set_orbital_context

  subroutine elm_set_nextsw_solar_time(ymd, tod, valid, calday)
    integer, intent(in) :: ymd
    real(r8), intent(in) :: tod
    logical, intent(in) :: valid
        real(r8), intent(in) :: calday
        call seq_infodata_PutData(solar_context, nextsw_ymd=ymd, nextsw_tod=tod, &
          nextsw_valid=valid, nextsw_cday=calday)
  end subroutine elm_set_nextsw_solar_time

  function elm_nextsw_eot_seconds() result(eot_seconds)
    use seq_timemgr_mod, only: seq_timemgr_GetCalendar
    real(r8) :: eot_seconds, tod, calday, delta, eccf
    integer :: ymd
    logical :: valid

    eot_seconds = 0.0_r8
    if (.not. solar_context_initialized) return
    if (.not. seq_infodata_UseEquationOfTime(solar_context)) return
    call seq_infodata_GetData(solar_context, nextsw_cday=calday, nextsw_ymd=ymd, nextsw_tod=tod, nextsw_valid=valid)
    if (calday < -0.5_r8) return
    if (.not. valid) call shr_sys_abort('ELM EOT geometry requires a full next-shortwave timestamp')
    call seq_infodata_SolarDecl(solar_context, ymd, tod, seq_timemgr_GetCalendar(), delta, eccf, eot_seconds)
  end function elm_nextsw_eot_seconds

  subroutine elm_orbital_decl(ymd, tod, calendar, delta, eccf, eot_seconds)
    integer, intent(in) :: ymd
    real(r8), intent(in) :: tod
    character(len=*), intent(in) :: calendar
    real(r8), intent(out) :: delta, eccf
    real(r8), intent(out), optional :: eot_seconds

    if (.not. solar_context_initialized) then
      call shr_sys_abort('ELM Gregorian solar geometry requires initialized orbital context')
    end if
    call seq_infodata_SolarDecl(solar_context, ymd, tod, calendar, delta, eccf, eot_seconds)
  end subroutine elm_orbital_decl

  function elm_solar_cosine(ymd, tod, calendar, latitude, longitude) result(cosine)
    ! Instantaneous cosine: full civil date, fractional seconds, and coordinates in radians.
    ! Requires elm_set_orbital_context; scheduled albedo and bypass sampling choose their own valid times.
    use shr_cal_mod, only: shr_cal_date2julian
    use shr_orb_mod, only: shr_orb_cosz
    integer, intent(in) :: ymd
    real(r8), intent(in) :: tod, latitude, longitude
    character(len=*), intent(in) :: calendar
    real(r8) :: cosine, calday, delta, eccf, eot_seconds

    call elm_orbital_decl(ymd, tod, calendar, delta, eccf, eot_seconds)
    call shr_cal_date2julian(ymd, int(tod), calday, calendar)
    calday = calday + (tod - real(int(tod), r8)) / 86400._r8
    cosine = shr_orb_cosz(calday, latitude, longitude, delta, eot_seconds=eot_seconds)
  end function elm_solar_cosine

end module elm_varorb
