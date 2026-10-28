module shr_input_calendar_mod
  use iso_fortran_env, only : real64
  use iso_c_binding, only : c_char, c_int, c_null_char
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  private

  integer, parameter, public :: input_calendar_max_overrides = 128
  integer, parameter, public :: input_calendar_selector_length = 1024
  integer, parameter, public :: input_calendar_signature_length = &
       input_calendar_max_overrides * (input_calendar_selector_length + 66)
  public :: input_calendar_normalize, input_calendar_resolve
  public :: input_calendar_check_times, input_calendar_check_date
  public :: input_calendar_flavor, input_calendar_check_range
  public :: input_calendar_set_overrides, input_calendar_get_override, input_calendar_signature
  public :: input_calendar_resolve_c

  character(input_calendar_selector_length), save :: selectors(input_calendar_max_overrides) = ''
  character(64), save :: calendars(input_calendar_max_overrides) = ''
  logical, save :: initialized = .false.

contains

  function input_calendar_normalize(text) result(calendar)
    character(*), intent(in) :: text
    character(64) :: calendar
    character(len(text)) :: token
    integer :: position, code, last

    token = text
    do position = 1, len(token)
      code = iachar(token(position:position))
      if (code >= iachar('A') .and. code <= iachar('Z')) token(position:position) = achar(code + 32)
      if (index(achar(9)//achar(10)//achar(13)//achar(0), token(position:position)) > 0) &
           token(position:position) = ' '
    end do
    token = adjustl(token)
    last = len_trim(token)
    calendar = ''
    select case (token(:last))
    case ('noleap', 'no_leap', '365_day', '365day')
      calendar = 'NO_LEAP'
    case ('gregorian', 'standard', 'proleptic_gregorian')
      calendar = 'GREGORIAN'
    end select
  end function input_calendar_normalize

  function input_calendar_flavor(text) result(flavor)
    character(*), intent(in) :: text
    character(64) :: flavor
    character(len(text)) :: token
    integer :: position, code

    flavor = input_calendar_normalize(text)
    if (flavor /= 'GREGORIAN') return
    token = adjustl(text)
    do position = 1, len(token)
      code = iachar(token(position:position))
      if (code >= iachar('A') .and. code <= iachar('Z')) token(position:position) = achar(code + 32)
      if (index(achar(9)//achar(10)//achar(13)//achar(0), token(position:position)) > 0) &
           token(position:position) = ' '
    end do
    token = adjustl(token)
    flavor = 'STANDARD'
    if (trim(token) == 'proleptic_gregorian') flavor = 'PROLEPTIC_GREGORIAN'
  end function input_calendar_flavor

  subroutine input_calendar_resolve(declared, metadata_present, override, calendar, provenance, error, flavor)
    character(*), intent(in) :: declared, override
    logical, intent(in) :: metadata_present
    character(*), intent(out) :: calendar, provenance, error
    character(*), optional, intent(out) :: flavor
    character(64) :: metadata_calendar, override_calendar

    calendar = ''
    provenance = ''
    error = ''
    if (present(flavor)) flavor = ''
    metadata_calendar = input_calendar_normalize(declared)
    if (.not. metadata_present) metadata_calendar = ''
    override_calendar = input_calendar_normalize(override)
    if (len_trim(override) > 0 .and. override_calendar == '') then
      error = 'NOLEAP_INPUT: unsupported override calendar='//trim(override)
      return
    end if
      if (metadata_present .and. verify(declared, ' '//achar(9)//achar(10)//achar(13)//achar(0)) > 0 .and. &
          metadata_calendar == '') then
      error = 'NOLEAP_INPUT: unsupported calendar='//trim(declared)
      return
    end if
    if (metadata_calendar /= '' .and. override_calendar /= '') then
      if (input_calendar_flavor(declared) /= input_calendar_flavor(override)) then
        error = 'NOLEAP_INPUT: calendar override conflicts with file metadata'
        return
      end if
    end if
    if (metadata_present .and. metadata_calendar /= '') then
      calendar = metadata_calendar
      provenance = 'metadata'
      if (override_calendar /= '') provenance = 'metadata+override'
      if (present(flavor)) flavor = input_calendar_flavor(declared)
    else if (override_calendar /= '') then
      calendar = override_calendar
      provenance = 'override'
      if (present(flavor)) flavor = input_calendar_flavor(override)
    else if (metadata_present) then
      error = 'NOLEAP_INPUT: empty calendar metadata'
    else
      error = 'NOLEAP_INPUT: explicit time calendar metadata is required'
    end if
  end subroutine input_calendar_resolve

  subroutine input_calendar_check_range(flavor, epoch, seconds, units, values, error)
    character(*), intent(in) :: flavor, units
    integer, intent(in) :: epoch
    real(real64), intent(in) :: seconds, values(:)
    character(*), intent(out) :: error
    real(real64) :: scale, elapsed
    integer :: year, month, day, ordinal, cutoff
    integer, parameter :: preceding(12) = [0,31,59,90,120,151,181,212,243,273,304,334]

    error = ''
    if (input_calendar_flavor(flavor) /= 'STANDARD') return
    call input_calendar_check_date(epoch, seconds, 'proleptic_gregorian', error)
    if (error /= '') return
    if (size(values) == 0) then
      error = 'INPUT_CALENDAR: empty mixed Gregorian source range'
      return
    end if
    if (.not. all(ieee_is_finite(values))) then
      error = 'INPUT_CALENDAR: non-finite mixed Gregorian source range'
      return
    end if
    select case (trim(units))
    case ('days', 'day')
      scale = 86400._real64
    case ('hours', 'hour')
      scale = 3600._real64
    case ('minutes', 'minute')
      scale = 60._real64
    case ('seconds', 'second', 'secs', 'sec')
      scale = 1._real64
    case default
      error = 'INPUT_CALENDAR: unsupported mixed Gregorian time units='//trim(units)
      return
    end select
    year = epoch / 10000
    month = mod(epoch / 100,100)
    day = mod(epoch,100)
    ordinal = 365*(year-1)+(year-1)/4-(year-1)/100+(year-1)/400+preceding(month)+day
    if (month > 2 .and. mod(year,4) == 0 .and. &
         (mod(year,100) /= 0 .or. mod(year,400) == 0)) ordinal = ordinal + 1
    cutoff = 577736
    elapsed = real(ordinal-cutoff,real64)*86400._real64 + seconds + minval(values)*scale
    if (epoch < 15821015 .or. elapsed < 0._real64 .or. .not. ieee_is_finite(elapsed)) then
      write(error,'(a,a,a,i8.8,a,es24.16,a,es24.16,a)') &
           'INPUT_CALENDAR: unsupported historical calendar=',trim(flavor),'; epoch=',epoch, &
           '; range=[',minval(values),',',maxval(values),']; require epoch and range >= 1582-10-15'
    end if
  end subroutine input_calendar_check_range

  subroutine input_calendar_check_times(values, error)
    real(real64), intent(in) :: values(:)
    character(*), intent(out) :: error
    integer :: record

    error = ''
    if (size(values) == 0) then
      error = 'NOLEAP_INPUT: empty time coordinate'
      return
    end if
    do record = 1, size(values)
      if (.not. ieee_is_finite(values(record))) then
        error = 'NOLEAP_INPUT: non-finite source time coordinate'
        return
      end if
      if (record > 1) then
        if (values(record) <= values(record-1)) then
          error = 'NOLEAP_INPUT: source time coordinate must be strictly increasing'
          return
        end if
      end if
    end do
  end subroutine input_calendar_check_times

  subroutine input_calendar_check_date(date, seconds, calendar, error)
    integer, intent(in) :: date
    real(real64), intent(in) :: seconds
    character(*), intent(in) :: calendar
    character(*), intent(out) :: error
    integer :: year, month, day, days(12)
    character(64) :: canonical

    error = 'NOLEAP_INPUT: invalid source date for declared calendar'
    canonical = input_calendar_normalize(calendar)
    if (canonical == '' .or. date < 0 .or. .not. ieee_is_finite(seconds)) return
    if (seconds < 0 .or. seconds >= 86400) return
    year = date / 10000
    month = mod(date / 100, 100)
    day = mod(date, 100)
    if (month < 1 .or. month > 12) return
    days = [31,28,31,30,31,30,31,31,30,31,30,31]
    if (canonical == 'GREGORIAN') then
      if (mod(year,4) == 0 .and. (mod(year,100) /= 0 .or. mod(year,400) == 0)) days(2) = 29
    end if
    if (day < 1 .or. day > days(month)) return
    error = ''
  end subroutine input_calendar_check_date

  subroutine input_calendar_set_overrides(requested_selectors, requested_calendars, error)
    character(*), intent(in) :: requested_selectors(:), requested_calendars(:)
    character(*), intent(out) :: error
    character(input_calendar_selector_length), allocatable :: next_selectors(:)
    character(input_calendar_selector_length) :: saved_selector
    character(64) :: next_calendars(input_calendar_max_overrides), saved_calendar
    integer :: entry, previous, count

    error = ''
    allocate(next_selectors(input_calendar_max_overrides))
    next_selectors = ''
    next_calendars = ''
    if (size(requested_selectors) /= size(requested_calendars) .or. &
         size(requested_selectors) > input_calendar_max_overrides) then
      error = 'NOLEAP_INPUT: override selector/calendar arrays have incompatible sizes'
      return
    end if
    count = 0
    do entry = 1, size(requested_selectors)
      if (len_trim(requested_selectors(entry)) == 0 .and. len_trim(requested_calendars(entry)) == 0) cycle
      if (len_trim(requested_selectors(entry)) > input_calendar_selector_length .or. &
           index(requested_selectors(entry), achar(10)) > 0 .or. &
           index(requested_selectors(entry), achar(13)) > 0) then
        error = 'NOLEAP_INPUT: invalid override selector'
        return
      end if
      saved_selector = trim(adjustl(requested_selectors(entry)))
      if (.not. ((index(saved_selector,'file:') == 1 .and. len_trim(saved_selector) > 5) .or. &
           (index(saved_selector,'stream:') == 1 .and. len_trim(saved_selector) > 7))) then
        error = 'NOLEAP_INPUT: override selector requires an exact file: or stream: identifier'
        return
      end if
      saved_calendar = input_calendar_flavor(requested_calendars(entry))
      if (saved_calendar == '') then
        error = 'NOLEAP_INPUT: unsupported override calendar='//trim(requested_calendars(entry))
        return
      end if
      do previous = 1, count
        if (next_selectors(previous) == saved_selector) then
          error = 'NOLEAP_INPUT: duplicate calendar override selector'
          return
        end if
      end do
      count = count + 1
      previous = count
      do while (previous > 1)
        if (next_selectors(previous-1) < saved_selector) exit
        next_selectors(previous) = next_selectors(previous-1)
        next_calendars(previous) = next_calendars(previous-1)
        previous = previous - 1
      end do
      next_selectors(previous) = saved_selector
      next_calendars(previous) = saved_calendar
    end do
    if (initialized) then
      if (any(selectors /= next_selectors) .or. any(calendars /= next_calendars)) &
           error = 'NOLEAP_INPUT: calendar overrides cannot change after initialization'
      return
    end if
    selectors = next_selectors
    calendars = next_calendars
    initialized = .true.
  end subroutine input_calendar_set_overrides

  subroutine input_calendar_get_override(filename, stream, calendar, rule, error)
    character(*), intent(in) :: filename, stream
    character(*), intent(out) :: calendar, rule, error
    integer :: entry

    calendar = ''
    rule = ''
    error = ''
    do entry = 1, input_calendar_max_overrides
      if (selectors(entry) /= 'file:'//trim(filename)) then
        if (len_trim(stream) == 0) cycle
        if (selectors(entry) /= 'stream:'//trim(stream)) cycle
      end if
      if (calendar /= '' .and. calendar /= calendars(entry)) then
        error = 'NOLEAP_INPUT: conflicting file and stream calendar overrides'
        return
      end if
      calendar = calendars(entry)
      rule = selectors(entry)
    end do
  end subroutine input_calendar_get_override

  function input_calendar_signature() result(signature)
    character(input_calendar_signature_length) :: signature
    integer :: entry, length

    signature = ''
    length = 0
    do entry = 1, input_calendar_max_overrides
      if (selectors(entry) == '') cycle
      signature = signature(:length)//trim(selectors(entry))//'='//trim(calendars(entry))//achar(10)
      length = len_trim(signature)
    end do
  end function input_calendar_signature

  integer(c_int) function input_calendar_resolve_c(filename, filename_length, declared, declared_length, &
      present, resolved, origin, flavor, diagnostic) bind(c)
    character(c_char), intent(in) :: filename(*), declared(*)
    integer(c_int), value :: filename_length, declared_length, present
    character(c_char), intent(out) :: resolved(*), origin(*), flavor(*), diagnostic(*)
    character(input_calendar_selector_length) :: path, metadata, rule, error
    character(64) :: override, calendar, provenance, source_flavor
    integer :: position

    error = ''
    calendar = ''
    provenance = ''
    source_flavor = ''
    if (filename_length < 0 .or. filename_length > len(path) .or. &
         declared_length < 0 .or. declared_length > len(metadata)) then
      error = 'NOLEAP_INPUT: calendar resolver string exceeds supported length'
    else
      path = ''
      metadata = ''
      do position = 1, filename_length
        path(position:position) = filename(position)
      end do
      do position = 1, declared_length
        metadata(position:position) = declared(position)
      end do
      call input_calendar_get_override(trim(path), '', override, rule, error)
      if (error == '') call input_calendar_resolve(metadata, present /= 0, override, calendar, provenance, error, source_flavor)
    end if
    call copy_c_string(calendar, resolved, 64)
    call copy_c_string(provenance, origin, 64)
    call copy_c_string(source_flavor, flavor, 64)
    call copy_c_string(error, diagnostic, 1024)
    input_calendar_resolve_c = 0
    if (error /= '') input_calendar_resolve_c = 1
  end function input_calendar_resolve_c

  subroutine copy_c_string(text, buffer, capacity)
    character(*), intent(in) :: text
    character(c_char), intent(out) :: buffer(*)
    integer, intent(in) :: capacity
    integer :: position

    buffer(1:capacity) = c_null_char
    do position = 1, min(len_trim(text), capacity-1)
      buffer(position) = text(position:position)
    end do
  end subroutine copy_c_string
end module shr_input_calendar_mod