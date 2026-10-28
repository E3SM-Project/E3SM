program test_input_calendar
  use iso_fortran_env, only : real64
  use iso_c_binding, only : c_char, c_int, c_null_char
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use shr_input_calendar_mod
  implicit none
  character(64) :: calendar, provenance
  character(1024) :: error, rule
  character(64) :: requested(2), values(2)
  character(c_char) :: c_calendar(64), c_origin(64), c_flavor(64), c_error(1024)
  integer(c_int) :: status

  call require(input_calendar_normalize(' NO_LEAP ') == 'NO_LEAP', 'noleap alias')
  call require(input_calendar_normalize('365day') == 'NO_LEAP', '365day alias')
  call require(input_calendar_normalize('standard') == 'GREGORIAN', 'Gregorian alias')
  call require(input_calendar_normalize('no leap') == '', 'internal whitespace rejected')
  call input_calendar_resolve('noleap', .true., '', calendar, provenance, error)
  call require(error == '' .and. calendar == 'NO_LEAP' .and. provenance == 'metadata', 'metadata')
  call input_calendar_resolve('', .false., '', calendar, provenance, error)
  call require(index(error,'explicit time calendar') > 0, 'missing metadata')
  call input_calendar_resolve('', .true., '', calendar, provenance, error)
  call require(index(error,'empty calendar') > 0, 'empty metadata')
  call input_calendar_resolve('', .false., '365_day', calendar, provenance, error)
  call require(error == '' .and. calendar == 'NO_LEAP' .and. provenance == 'override', 'explicit override')
  call input_calendar_resolve('gregorian', .true., 'noleap', calendar, provenance, error)
  call require(index(error,'conflicts') > 0, 'metadata conflict')
  call input_calendar_resolve('360_day', .true., 'noleap', calendar, provenance, error)
  call require(index(error,'unsupported calendar') > 0, 'unsupported metadata cannot be overridden')
  call input_calendar_check_times([0._real64, 0.25_real64, 1._real64], error)
  call require(error == '', 'finite ordered coordinates')
  call input_calendar_check_times([1._real64, 1._real64], error)
  call require(index(error,'strictly increasing') > 0, 'duplicate coordinates')
  call input_calendar_check_times([1._real64, 0._real64], error)
  call require(index(error,'strictly increasing') > 0, 'reversed coordinates')
  call input_calendar_check_times([ieee_value(0._real64, ieee_quiet_nan)], error)
  call require(index(error,'non-finite') > 0, 'nonfinite coordinate')
  call input_calendar_check_date(20000229, 0._real64, 'noleap', error)
  call require(error /= '', 'impossible noleap February 29')
  call input_calendar_check_date(20000229, 0._real64, 'gregorian', error)
  call require(error == '', 'valid Gregorian leap day')
  call require(input_calendar_flavor('standard') == 'STANDARD', 'mixed Gregorian flavor')
  call require(input_calendar_flavor('proleptic_gregorian') == 'PROLEPTIC_GREGORIAN', 'proleptic flavor')
  call input_calendar_resolve('standard', .true., 'proleptic_gregorian', calendar, provenance, error)
  call require(index(error,'conflicts') > 0, 'override cannot relabel mixed metadata')
  call input_calendar_check_range('standard',10101,0._real64,'days',[730181._real64],error)
  call require(error /= '', 'year-one mixed epoch rejected even for modern records')
  call input_calendar_check_range('proleptic_gregorian',10101,0._real64,'days',[730181._real64],error)
  call require(error == '', 'historical proleptic epoch retained')
  call input_calendar_check_range('gregorian',20000101,0._real64,'days',[0._real64,60._real64],error)
  call require(error == '', 'modern mixed Gregorian accepted')
  call input_calendar_check_range('standard',15821015,0._real64,'seconds',[-1._real64,1._real64],error)
  call require(error /= '', 'negative offset crosses reform cutoff')
  call input_calendar_check_range('standard',15821015,0._real64,'days',[0._real64,1._real64],error)
  call require(error == '', 'reform cutoff is inclusive')
  call input_calendar_check_range('standard',15821010,0._real64,'days',[10._real64],error)
  call require(error /= '', 'reform gap epoch rejected')
  call input_calendar_check_date(19000229, 0._real64, 'gregorian', error)
  call require(error /= '', 'invalid Gregorian century leap day')
  call input_calendar_check_date(20000101, 86400._real64, 'noleap', error)
  call require(error /= '', 'invalid seconds of day')
  requested = ['file:/forcing/a.nc  ', 'stream:/case/forcing']
  values = ['noleap  ', '365_day ']
  call input_calendar_set_overrides(requested, values, error)
  call require(error == '', 'register overrides')
  status = input_calendar_resolve_c('/forcing/a.nc',13_c_int,'',0_c_int,0_c_int,c_calendar,c_origin,c_flavor,c_error)
  call require(status == 0 .and. c_calendar(1) == 'N' .and. c_origin(1) == 'o', 'C bridge exact override')
  call require(any(c_calendar == c_null_char) .and. any(c_origin == c_null_char), 'C bridge termination')
    status = input_calendar_resolve_c('/forcing/a.nc',13_c_int,'gregorian',9_c_int,1_c_int,c_calendar,c_origin,c_flavor,c_error)
  call require(status /= 0 .and. c_error(1) == 'N', 'C bridge metadata conflict')
    status = input_calendar_resolve_c('/forcing/b.nc',13_c_int,'proleptic_gregorian',19_c_int,1_c_int, &
      c_calendar,c_origin,c_flavor,c_error)
    call require(status == 0 .and. c_flavor(1) == 'P', 'C bridge preserves proleptic flavor')
  call input_calendar_get_override('/forcing/a.nc', '', calendar, rule, error)
  call require(error == '' .and. calendar == 'NO_LEAP' .and. rule == 'file:/forcing/a.nc', 'file match')
  call input_calendar_get_override('/forcing/b.nc', '/case/forcing', calendar, rule, error)
  call require(error == '' .and. calendar == 'NO_LEAP', 'stream match')
  call input_calendar_get_override('/elsewhere/a.nc', '', calendar, rule, error)
  call require(error == '' .and. calendar == '', 'no basename guessing')
  call input_calendar_set_overrides(requested(2:1:-1), values(2:1:-1), error)
  call require(error == '', 'entry order does not change policy')
  values(1) = 'gregorian'
  call input_calendar_set_overrides(requested, values, error)
  call require(index(error,'cannot change') > 0, 'immutable overrides')
  call require(index(input_calendar_signature(), 'file:/forcing/a.nc=NO_LEAP') > 0, 'provenance signature')
  print *, 'Input calendar resolver: PASS'

contains

  subroutine require(condition, label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label
    if (.not. condition) then
      print *, 'FAIL: '//label
      stop 1
    end if
  end subroutine require
end program test_input_calendar