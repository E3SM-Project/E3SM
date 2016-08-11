module shr_wv_sat_mod

! This portable module contains all CAM methods for estimating the saturation
! vapor pressure of water.
!
! wv_saturation provides CAM-specific interfaces and utilities based on these
! formulae.
!
! Typical usage of this module:
!
! Init:
! call shr_wv_sat_init(<constants>, errstring)
!
! Get scheme index from a name string:
! scheme_idx = shr_wv_sat_get_scheme_idx("GoffGratch")
! if (.not. shr_wv_sat_valid_idx(scheme_idx)) <throw some error>
!
! Get pressures:
! es = shr_wv_sat_svp_liquid(t, scheme_idx)
! es = shr_wv_sat_svp_ice(t, scheme_idx)
!
! Use ice/water transition range:
! es = shr_wv_sat_svp_mixed(t, scheme_idx)
!
! Note that elemental functions cannot be pointed to, nor passed as
! arguments. If you need to do either, it is recommended to wrap the function so
! that it can be given an explicit (non-elemental) interface.
!
! Since usually only one scheme is used at a time, the scheme index is an
! optional argument. If omitted a default scheme will be used, which is
! initially the Goff & Gratch scheme. To change the default, you can make a call
! like this:
!
! call shr_wv_sat_set_default("MurphyKoop")
!
! This module has the ability to store lookup tables to cache values. To do so,
! create a ShrWVSatTableSpec instance for water and ice with the desired size
! and range, then call shr_wv_sat_make_tables like so:
!
! type(ShrWVSatTableSpec) :: liquid_spec, ice_spec, mixed_spec
!
! liquid_spec = ShrWVSatTableSpec(151, tmelt-50._rk, 1._rk)
! ice_spec = ShrWVSatTableSpec(106, tmelt-100._rk, 1._rk)
! mixed_spec = ShrWVSatTableSpec(21, tmelt-20._rk, 1._rk)
! call shr_wv_sat_make_tables(liquid_spec, ice_spec, mixed_spec)
!
! Currently, this module only supports making tables for the default scheme, and
! shr_wv_sat_make_tables must be invoked again to produce new tables if the default
! is changed.
!
! Once tables are produced, all uses of the default scheme will attempt to
! linearly interpolate values from the lookup tables. If a temperature is
! outside the table bounds, the original scheme will be invoked as if no table
! was present. That is, the tables will *not* be extrapolated.
!
! Threading note: The physical constants, the current default, and the lookup
! tables are stored in global, thread-shared variables. Therefore the following
! procedures are not thread-safe:
!
!  - shr_wv_sat_init
!  - shr_wv_sat_final
!  - shr_wv_sat_set_default
!  - shr_wv_sat_make_tables
!
! If a multi-threaded application calls any of these procedures, the user is
! responsible for ensuring that each call is perfomed by only one thread, and
! that synchronization is performed to before and after each of these calls.
!
! All other public routines are thread-safe.

#ifdef SHR_WV_SAT_USE_CLUBB_KIND
! If running within CLUBB, use the same precision as CLUBB.
use clubb_precision, only: core_rkind
#endif

implicit none
private
save

! Initialize/finalize module and pick schemes
public shr_wv_sat_init
public shr_wv_sat_final
public shr_wv_sat_get_scheme_idx
public shr_wv_sat_valid_idx
public shr_wv_sat_set_default

! Manage lookup tables
public ShrWVSatTableSpec
public shr_wv_sat_make_tables

! Basic SVP calculations
public shr_wv_sat_svp_liquid
public shr_wv_sat_svp_ice
public shr_wv_sat_svp_mixed

! pressure -> humidity conversion
public shr_wv_sat_svp_to_qsat

! pressure -> mass mixing ratio conversion
public shr_wv_sat_svp_to_qmmr

! Combined qsat operations
public shr_wv_sat_qsat_liquid
public shr_wv_sat_qsat_ice
public shr_wv_sat_qsat_mixed

#ifdef SHR_WV_SAT_USE_CLUBB_KIND
integer, parameter :: rk = core_rkind
#else
! Double precision
integer, parameter :: rk = selected_real_kind(12)
#endif

real(rk) :: tmelt   ! Melting point of water at 1 atm (K)
real(rk) :: h2otrip ! Triple point temperature of water (K)

real(rk) :: ttrice  ! Ice-water transition range

real(rk) :: epsilo  ! Ice-water transition range
real(rk) :: omeps   ! 1._rk - epsilo

! Indices representing individual schemes
integer, parameter :: Invalid_idx = -1
! Skip 0; this was the old implementation of Goff & Gratch.
integer, parameter :: GoffGratch_idx = 1
integer, parameter :: MurphyKoop_idx = 2
integer, parameter :: Bolton_idx = 3
integer, parameter :: Flatau_idx = 4

! Index representing the current default scheme.
integer, parameter :: initial_default_idx = GoffGratch_idx
integer :: default_idx = initial_default_idx

! Type to represent a table specification.
type ShrWVSatTableSpec
   integer :: table_size
   real(rk) :: minimum
   real(rk) :: spacing
end type ShrWVSatTableSpec

type(ShrWVSatTableSpec) :: liquid_table_spec
real(rk), allocatable :: liquid_table(:)

type(ShrWVSatTableSpec) :: ice_table_spec
real(rk), allocatable :: ice_table(:)

type(ShrWVSatTableSpec) :: mixed_table_spec
real(rk), allocatable :: mixed_table(:)

interface shr_wv_sat_svp_liquid
   module procedure shr_wv_sat_svp_liquid
   module procedure shr_wv_sat_svp_liquid_vec
end interface shr_wv_sat_svp_liquid

interface shr_wv_sat_svp_ice
   module procedure shr_wv_sat_svp_ice
   module procedure shr_wv_sat_svp_ice_vec
end interface shr_wv_sat_svp_ice

interface shr_wv_sat_svp_mixed
   module procedure shr_wv_sat_svp_mixed
   module procedure shr_wv_sat_svp_mixed_vec
end interface shr_wv_sat_svp_mixed

interface shr_wv_sat_svp_to_qsat
   module procedure shr_wv_sat_svp_to_qsat
   module procedure shr_wv_sat_svp_to_qsat_vec
end interface shr_wv_sat_svp_to_qsat

interface shr_wv_sat_svp_to_qmmr
   module procedure shr_wv_sat_svp_to_qmmr
   module procedure shr_wv_sat_svp_to_qmmr_vec
end interface shr_wv_sat_svp_to_qmmr

interface shr_wv_sat_qsat_liquid
   module procedure shr_wv_sat_qsat_liquid
   module procedure shr_wv_sat_qsat_liquid_vec
end interface shr_wv_sat_qsat_liquid

interface shr_wv_sat_qsat_ice
   module procedure shr_wv_sat_qsat_ice
   module procedure shr_wv_sat_qsat_ice_vec
end interface shr_wv_sat_qsat_ice

interface shr_wv_sat_qsat_mixed
   module procedure shr_wv_sat_qsat_mixed
   module procedure shr_wv_sat_qsat_mixed_vec
end interface shr_wv_sat_qsat_mixed

contains

!---------------------------------------------------------------------
! ADMINISTRATIVE FUNCTIONS
!---------------------------------------------------------------------

! Get physical constants
subroutine shr_wv_sat_init(tmelt_in, h2otrip_in, ttrice_in, epsilo_in, &
     errstring)
  real(rk), intent(in) :: tmelt_in
  real(rk), intent(in) :: h2otrip_in
  real(rk), intent(in) :: ttrice_in
  real(rk), intent(in) :: epsilo_in
  character(len=*), intent(out)  :: errstring

  errstring = ' '

  if (ttrice_in < 0._rk) then
     write(errstring,*) 'shr_wv_sat_init: ERROR: ', &
          ttrice_in,' was input for ttrice, but negative range is invalid.'
     return
  end if

  tmelt = tmelt_in
  h2otrip = h2otrip_in
  ttrice = ttrice_in
  epsilo = epsilo_in

  omeps = 1._rk - epsilo

end subroutine shr_wv_sat_init

! Reset module data to the original state (primarily for testing purposes).
! It doesn't seem worthwhile to reset the constants, so just deal with options
! and dynamic memory here.
subroutine shr_wv_sat_final()

  default_idx = initial_default_idx

  if (allocated(liquid_table)) deallocate(liquid_table)
  if (allocated(ice_table)) deallocate(ice_table)
  if (allocated(mixed_table)) deallocate(mixed_table)

end subroutine shr_wv_sat_final

! Look up index by name.
pure function shr_wv_sat_get_scheme_idx(name) result(idx)
  character(len=*), intent(in) :: name
  integer :: idx

  ! Several names are given to most methods in order to support the names that
  ! CLUBB accepts.
  select case (name)
  case("GoffGratch", "gfdl", "GFDL")
     idx = GoffGratch_idx
  case("MurphyKoop")
     idx = MurphyKoop_idx
  case("Bolton", "bolton")
     idx = Bolton_idx
  case("Flatau", "flatau")
     idx = Flatau_idx
  case default
     idx = Invalid_idx
  end select

end function shr_wv_sat_get_scheme_idx

! Check validity of an index from the above routine.
pure function shr_wv_sat_valid_idx(idx) result(status)
  integer, intent(in) :: idx
  logical :: status

  status = (idx /= Invalid_idx)

end function shr_wv_sat_valid_idx

! Set default scheme (otherwise, Goff & Gratch is default)
! Returns a logical representing success (.true.) or
! failure (.false.).
function shr_wv_sat_set_default(name) result(status)
  character(len=*), intent(in) :: name
  logical :: status

  ! Don't want to overwrite valid default with invalid,
  ! so assign to temporary and check it first.
  integer :: tmp_idx

  tmp_idx = shr_wv_sat_get_scheme_idx(name)

  status = shr_wv_sat_valid_idx(tmp_idx)

  ! If we have changed the default, deallocated the tables as well as setting
  ! the new default.
  if (status .and. tmp_idx /= default_idx) then
     if (allocated(liquid_table)) deallocate(liquid_table)
     if (allocated(ice_table)) deallocate(ice_table)
     if (allocated(mixed_table)) deallocate(mixed_table)
     default_idx = tmp_idx
  end if

end function shr_wv_sat_set_default

subroutine shr_wv_sat_make_tables(liquid_spec_in, ice_spec_in, mixed_spec_in)
  type(ShrWVSatTableSpec), intent(in), optional :: liquid_spec_in
  type(ShrWVSatTableSpec), intent(in), optional :: ice_spec_in
  type(ShrWVSatTableSpec), intent(in), optional :: mixed_spec_in

  if (present(liquid_spec_in)) then
     liquid_table_spec = liquid_spec_in
     call shr_wv_sat_make_one_table(liquid_table_spec, shr_wv_sat_svp_liquid_no_table, &
          liquid_table)
  end if

  if (present(ice_spec_in)) then
     ice_table_spec = ice_spec_in
     call shr_wv_sat_make_one_table(ice_table_spec, shr_wv_sat_svp_ice_no_table, &
          ice_table)
  end if

  if (present(mixed_spec_in)) then
     mixed_table_spec = mixed_spec_in
     call shr_wv_sat_make_one_table(mixed_table_spec, shr_wv_sat_svp_mixed_no_table, &
          mixed_table)
  end if

end subroutine shr_wv_sat_make_tables

! Table-generating generic function (would be simpler with an object-oriented
! design, but we want this code to be runnable on compilers with poor Fortran
! 2003 support). This means we can't attach function pointers or methods to
! derived types.
subroutine shr_wv_sat_make_one_table(table_spec, svp_function, table)
  type(ShrWVSatTableSpec), intent(in) :: table_spec
  interface
     function svp_function(t, idx) result(es)
       import :: rk
       real(rk), intent(in) :: t
       integer, intent(in) :: idx
       real(rk) :: es
     end function svp_function
  end interface
  real(rk), intent(out), allocatable :: table(:)

  integer :: i

  allocate(table(table_spec%table_size))

  do i = 1, table_spec%table_size
     table(i) = svp_function( &
          table_spec%minimum + table_spec%spacing*real(i-1, rk), &
          default_idx)
  end do

end subroutine shr_wv_sat_make_one_table

!---------------------------------------------------------------------
! UTILITIES
!---------------------------------------------------------------------

! Get saturation specific humidity given pressure and SVP.
! Specific humidity is limited to the range 0-1.
elemental function shr_wv_sat_svp_to_qsat(es, p) result(qs)

  real(rk), intent(in) :: es  ! SVP
  real(rk), intent(in) :: p   ! Current pressure.
  real(rk) :: qs

  ! If pressure is less than SVP, set qs to maximum of 1.
  if ( (p - es) <= 0._rk ) then
     qs = 1.0_rk
  else
     qs = epsilo*es / (p - omeps*es)
  end if

end function shr_wv_sat_svp_to_qsat

pure function shr_wv_sat_svp_to_qsat_vec(n, es, p) result(qs)

  integer,  intent(in) :: n      ! Size of input arrays
  real(rk), intent(in) :: es(n)  ! SVP
  real(rk), intent(in) :: p(n)   ! Current pressure
  real(rk) :: qs(n)

  integer :: i

  ! If pressure is less than SVP, set qs to maximum of 1.
  do i = 1, n
     if ( (p(i) - es(i)) <= 0._rk ) then
        qs(i) = 1.0_rk
     else
        qs(i) = epsilo*es(i) / (p(i) - omeps*es(i))
     end if
  end do

end function shr_wv_sat_svp_to_qsat_vec

! Get saturation mass mixing ratio (over dry air) given pressure and SVP.
! Output is limited to the range 0-epsilo.
elemental function shr_wv_sat_svp_to_qmmr(es, p) result(qs)

  real(rk), intent(in) :: es  ! SVP
  real(rk), intent(in) :: p   ! Current pressure
  real(rk) :: qs

  ! If pressure is less than SVP, set qs to maximum of 1.
  if ( (p - es) <= es ) then
     qs = epsilo
  else
     qs = epsilo*es / (p - es)
  end if

end function shr_wv_sat_svp_to_qmmr

pure function shr_wv_sat_svp_to_qmmr_vec(n, es, p) result(qs)

  integer,  intent(in) :: n      ! Size of input arrays
  real(rk), intent(in) :: es(n)  ! SVP
  real(rk), intent(in) :: p(n)   ! Current pressure
  real(rk) :: qs(n)

  integer :: i

  ! If pressure is less than SVP, set qs to maximum of 1.
  do i = 1, n
     if ( (p(i) - es(i)) <= es(i) ) then
        qs(i) = epsilo
     else
        qs(i) = epsilo*es(i) / (p(i) - es(i))
     end if
  end do

end function shr_wv_sat_svp_to_qmmr_vec

elemental subroutine shr_wv_sat_qsat_liquid(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(rk), intent(in) :: t    ! Temperature
  real(rk), intent(in) :: p    ! Pressure
  ! Outputs
  real(rk), intent(out) :: es  ! Saturation vapor pressure
  real(rk), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = shr_wv_sat_svp_liquid(t, idx)

  qs = shr_wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine shr_wv_sat_qsat_liquid

pure subroutine shr_wv_sat_qsat_liquid_vec(n, t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over water at a given temperature, and then      !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  integer,  intent(in) :: n       ! Size of input arrays
  real(rk), intent(in) :: t(n)    ! Temperature
  real(rk), intent(in) :: p(n)    ! Pressure
  ! Outputs
  real(rk), intent(out) :: es(n)  ! Saturation vapor pressure
  real(rk), intent(out) :: qs(n)  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = shr_wv_sat_svp_liquid(n, t, idx)

  qs = shr_wv_sat_svp_to_qsat(n, es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine shr_wv_sat_qsat_liquid_vec

elemental subroutine shr_wv_sat_qsat_ice(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(rk), intent(in) :: t    ! Temperature
  real(rk), intent(in) :: p    ! Pressure
  ! Outputs
  real(rk), intent(out) :: es  ! Saturation vapor pressure
  real(rk), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = shr_wv_sat_svp_ice(t, idx)

  qs = shr_wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine shr_wv_sat_qsat_ice

pure subroutine shr_wv_sat_qsat_ice_vec(n, t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  integer,  intent(in) :: n       ! Size of input arrays
  real(rk), intent(in) :: t(n)    ! Temperature
  real(rk), intent(in) :: p(n)    ! Pressure
  ! Outputs
  real(rk), intent(out) :: es(n)  ! Saturation vapor pressure
  real(rk), intent(out) :: qs(n)  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = shr_wv_sat_svp_ice(n, t, idx)

  qs = shr_wv_sat_svp_to_qsat(n, es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine shr_wv_sat_qsat_ice_vec

elemental subroutine shr_wv_sat_qsat_mixed(t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  real(rk), intent(in) :: t    ! Temperature
  real(rk), intent(in) :: p    ! Pressure
  ! Outputs
  real(rk), intent(out) :: es  ! Saturation vapor pressure
  real(rk), intent(out) :: qs  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = shr_wv_sat_svp_mixed(t, idx)

  qs = shr_wv_sat_svp_to_qsat(es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine shr_wv_sat_qsat_mixed

pure subroutine shr_wv_sat_qsat_mixed_vec(n, t, p, es, qs, idx)
  !------------------------------------------------------------------!
  ! Purpose:                                                         !
  !   Calculate SVP over ice at a given temperature, and then        !
  !   calculate and return saturation specific humidity.             !
  !------------------------------------------------------------------!

  ! Inputs
  integer,  intent(in) :: n       ! Size of input arrays
  real(rk), intent(in) :: t(n)    ! Temperature
  real(rk), intent(in) :: p(n)    ! Pressure
  ! Outputs
  real(rk), intent(out) :: es(n)  ! Saturation vapor pressure
  real(rk), intent(out) :: qs(n)  ! Saturation specific humidity

  integer,  intent(in), optional :: idx ! Scheme index

  es = shr_wv_sat_svp_mixed(n, t, idx)

  qs = shr_wv_sat_svp_to_qsat(n, es, p)

  ! Ensures returned es is consistent with limiters on qs.
  es = min(es, p)

end subroutine shr_wv_sat_qsat_mixed_vec

!---------------------------------------------------------------------
! SVP INTERFACE FUNCTIONS
!---------------------------------------------------------------------

elemental function shr_wv_sat_svp_liquid(t, idx) result(es)
  real(rk), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(rk) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  if (use_idx == default_idx .and. allocated(liquid_table)) then
     es = lookup_svp_in_table(t, liquid_table_spec, liquid_table, &
          shr_wv_sat_svp_liquid_no_table)
  else
     es = shr_wv_sat_svp_liquid_no_table(t, use_idx)
  end if

end function shr_wv_sat_svp_liquid

pure function shr_wv_sat_svp_liquid_vec(n, t, idx) result(es)
  integer,  intent(in) :: n
  real(rk), intent(in) :: t(n)
  integer,  intent(in), optional :: idx
  real(rk) :: es(n)

  integer :: use_idx
  integer :: i

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  if (use_idx == default_idx .and. allocated(liquid_table)) then
     do i = 1, n
        es(i) = lookup_svp_in_table(t(i), liquid_table_spec, liquid_table, &
             shr_wv_sat_svp_liquid_no_table)
     end do
  else
     do i = 1, n
        es(i) = shr_wv_sat_svp_liquid_no_table(t(i), use_idx)
     end do
  end if

end function shr_wv_sat_svp_liquid_vec

pure function shr_wv_sat_svp_liquid_no_table(t, idx) result(es)
  real(rk), intent(in) :: t
  integer,  intent(in) :: idx
  real(rk) :: es

  select case (idx)
  case(GoffGratch_idx)
     es = GoffGratch_svp_liquid(t)
  case(MurphyKoop_idx)
     es = MurphyKoop_svp_liquid(t)
  case(Bolton_idx)
     es = Bolton_svp_liquid(t)
  case (Flatau_idx)
     es = Flatau_svp_liquid(t)
  case default
     ! Providing a correct index is an important precondition for these
     ! functions. Since we don't have a way of signaling an error, produce an
     ! obviously unreasonable answer.
     es = -huge(1._rk)
  end select

end function shr_wv_sat_svp_liquid_no_table

elemental function shr_wv_sat_svp_ice(t, idx) result(es)
  real(rk), intent(in) :: t
  integer,  intent(in), optional :: idx
  real(rk) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  if (use_idx == default_idx .and. allocated(ice_table)) then
     es = lookup_svp_in_table(t, ice_table_spec, ice_table, &
          shr_wv_sat_svp_ice_no_table)
  else
     es = shr_wv_sat_svp_ice_no_table(t, use_idx)
  end if

end function shr_wv_sat_svp_ice

pure function shr_wv_sat_svp_ice_vec(n, t, idx) result(es)
  integer,  intent(in) :: n
  real(rk), intent(in) :: t(n)
  integer,  intent(in), optional :: idx
  real(rk) :: es(n)

  integer :: use_idx
  integer :: i

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  if (use_idx == default_idx .and. allocated(ice_table)) then
     do i = 1, n
        es(i) = lookup_svp_in_table(t(i), ice_table_spec, ice_table, &
             shr_wv_sat_svp_ice_no_table)
     end do
  else
     do i = 1, n
        es(i) = shr_wv_sat_svp_ice_no_table(t(i), use_idx)
     end do
  end if

end function shr_wv_sat_svp_ice_vec

pure function shr_wv_sat_svp_ice_no_table(t, idx) result(es)
  real(rk), intent(in) :: t
  integer,  intent(in) :: idx
  real(rk) :: es

  select case (idx)
  case(GoffGratch_idx)
     es = GoffGratch_svp_ice(t)
  case(MurphyKoop_idx)
     es = MurphyKoop_svp_ice(t)
  case(Bolton_idx)
     es = Bolton_svp_ice(t)
  case (Flatau_idx)
     es = Flatau_svp_ice(t)
  case default
     ! Providing a correct index is an important precondition for these
     ! functions. Since we don't have a way of signaling an error, produce an
     ! obviously unreasonable answer.
     es = -huge(1._rk)
  end select

end function shr_wv_sat_svp_ice_no_table

elemental function shr_wv_sat_svp_mixed(t, idx) result (es)

  real(rk), intent(in) :: t
  integer,  intent(in), optional :: idx

  real(rk) :: es

  integer :: use_idx

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  if (use_idx == default_idx .and. allocated(mixed_table)) then
     es = lookup_svp_in_table(t, mixed_table_spec, mixed_table, &
          shr_wv_sat_svp_mixed_no_table)
  else
     es = shr_wv_sat_svp_mixed_no_table(t, use_idx)
  end if

end function shr_wv_sat_svp_mixed

pure function shr_wv_sat_svp_mixed_vec(n, t, idx) result (es)
  integer,  intent(in) :: n
  real(rk), intent(in) :: t(n)
  integer,  intent(in), optional :: idx

  real(rk) :: es(n)

  integer :: use_idx
  integer :: i

  if (present(idx)) then
     use_idx = idx
  else
     use_idx = default_idx
  end if

  if (use_idx == default_idx .and. allocated(mixed_table)) then
     do i = 1, n
        es(i) = lookup_svp_in_table(t(i), mixed_table_spec, mixed_table, &
             shr_wv_sat_svp_mixed_no_table)
     end do
  else
     es = shr_wv_sat_svp_mixed_no_table_vec(n, t, use_idx)
  end if

end function shr_wv_sat_svp_mixed_vec

pure function shr_wv_sat_svp_mixed_no_table(t, idx) result (es)

  real(rk), intent(in) :: t
  integer,  intent(in) :: idx

  real(rk) :: es

  real(rk) :: esice      ! Saturation vapor pressure over ice
  real(rk) :: weight     ! Intermediate scratch variable for es transition

!
! Water
!
  if (t >= (tmelt - ttrice)) then
     es = shr_wv_sat_svp_liquid(t,idx)
  else
     es = 0.0_rk
  end if

!
! Ice
!
  if (t < tmelt) then

     esice = shr_wv_sat_svp_ice(t,idx)

     if ( (tmelt - t) > ttrice ) then
        weight = 1.0_rk
     else
        weight = (tmelt - t)/ttrice
     end if

     es = weight*esice + (1.0_rk - weight)*es
  end if

end function shr_wv_sat_svp_mixed_no_table

! The reason why we need a separate vector version for this function (but
! not the other "no_table" functions) is that even if the functions called
! are all inlined, we may still have a table lookup for the liquid/ice
! tables (rather than the mixed-phase one). That table lookup is not
! vectorizable, meaning that the above function can't vectorize.
!
! Cripes this is ugly, though.
pure function shr_wv_sat_svp_mixed_no_table_vec(n, t, idx) result (es)

  integer,  intent(in) :: n
  real(rk), intent(in) :: t(n)
  integer,  intent(in) :: idx

  real(rk) :: es(n)

  real(rk) :: esice      ! Saturation vapor pressure over ice
  real(rk) :: weight     ! Intermediate scratch variable for es transition

  integer :: i

!
! Water
!
  if (idx == default_idx .and. allocated(liquid_table)) then
     do i = 1, n
        if (t(i) >= (tmelt - ttrice)) then
           es(i) = lookup_svp_in_table(t(i), liquid_table_spec, liquid_table, &
                shr_wv_sat_svp_liquid_no_table)
        else
           es(i) = 0.0_rk
        end if
     end do
  else
     do i = 1, n
        if (t(i) >= (tmelt - ttrice)) then
           es(i) = shr_wv_sat_svp_liquid_no_table(t(i), idx)
        else
           es(i) = 0.0_rk
        end if
     end do
  end if

!
! Ice
!
  if (idx == default_idx .and. allocated(ice_table)) then
     do i = 1, n
        if (t(i) < tmelt) then
           esice = lookup_svp_in_table(t(i), ice_table_spec, ice_table, &
                shr_wv_sat_svp_ice_no_table)

           if ( (tmelt - t(i)) > ttrice ) then
              weight = 1.0_rk
           else
              weight = (tmelt - t(i))/ttrice
           end if

           es(i) = weight*esice + (1.0_rk - weight)*es(i)
        end if
     end do
  else
     do i = 1, n
        if (t(i) < tmelt) then
           esice = shr_wv_sat_svp_ice_no_table(t(i), idx)

           if ( (tmelt - t(i)) > ttrice ) then
              weight = 1.0_rk
           else
              weight = (tmelt - t(i))/ttrice
           end if

           es(i) = weight*esice + (1.0_rk - weight)*es(i)
        end if
     end do
  end if

end function shr_wv_sat_svp_mixed_no_table_vec

!---------------------------------------------------------------------
! SVP METHODS
!---------------------------------------------------------------------

! Use the lookup table
! Note that a function that takes a procedure argument can't be elemental, but
! it must be pure to be called from the elemental interfaces above.
recursive pure function lookup_svp_in_table(t, table_spec, table, fallback) &
     result(es)
  ! Temperature in Kelvin
  real(rk), intent(in) :: t
  ! Table range specification
  type(ShrWVSatTableSpec), intent(in) :: table_spec
  ! The table itself
  real(rk), intent(in) :: table(table_spec%table_size)
  ! Fallback used when we're outside the table bounds
  interface
     pure function fallback(t, idx) result(es)
       import :: rk
       real(rk), intent(in) :: t
       integer, intent(in) :: idx
       real(rk) :: es
     end function fallback
  end interface
  ! SVP in Pa
  real(rk) :: es

  integer :: idx
  real(rk) :: residual

  residual = (t - table_spec%minimum)/table_spec%spacing + 1._rk
  idx = int(residual)
  ! Deal with the case where we're outside the table bounds.
  if (idx < 1 .or. idx+1 > table_spec%table_size) then
     es = fallback(t, default_idx)
     return
  end if

  ! Now we want the "residual" to be how far this temperature is past the
  ! current index.
  residual = residual - real(idx, rk)

  ! Just use linear interpolation. Some iterative methods might do better if we
  ! used a cubic.
  es = table(idx) * (1._rk-residual) + table(idx+1) * residual

end function lookup_svp_in_table

! Goff & Gratch (1946)

pure function GoffGratch_svp_liquid(t) result(es)
  real(rk), intent(in) :: t  ! Temperature in Kelvin
  real(rk) :: es             ! SVP in Pa

  ! Boiling point of water at 1 atm (K)
  ! This is not really the most accurate value, but it is the one that the
  ! Goff & Gratch scheme was designed to use, so we hard-code it.
  real(rk), parameter :: tboil = 373.16_rk

  ! uncertain below -70 C
  es = 10._rk**(-7.90298_rk*(tboil/t-1._rk)+ &
       5.02808_rk*log10(tboil/t)- &
       1.3816e-7_rk*(10._rk**(11.344_rk*(1._rk-t/tboil))-1._rk)+ &
       8.1328e-3_rk*(10._rk**(-3.49149_rk*(tboil/t-1._rk))-1._rk)+ &
       log10(1013.246_rk))*100._rk

end function GoffGratch_svp_liquid

pure function GoffGratch_svp_ice(t) result(es)
  real(rk), intent(in) :: t  ! Temperature in Kelvin
  real(rk) :: es             ! SVP in Pa

  ! good down to -100 C
  es = 10._rk**(-9.09718_rk*(h2otrip/t-1._rk)-3.56654_rk* &
       log10(h2otrip/t)+0.876793_rk*(1._rk-t/h2otrip)+ &
       log10(6.1071_rk))*100._rk

end function GoffGratch_svp_ice

! Murphy & Koop (2005)

pure function MurphyKoop_svp_liquid(t) result(es)
  real(rk), intent(in) :: t  ! Temperature in Kelvin
  real(rk) :: es             ! SVP in Pa

  ! (good for 123 < T < 332 K)
  es = exp(54.842763_rk - (6763.22_rk / t) - (4.210_rk * log(t)) + &
       (0.000367_rk * t) + (tanh(0.0415_rk * (t - 218.8_rk)) * &
       (53.878_rk - (1331.22_rk / t) - (9.44523_rk * log(t)) + &
       0.014025_rk * t)))

end function MurphyKoop_svp_liquid

pure function MurphyKoop_svp_ice(t) result(es)
  real(rk), intent(in) :: t  ! Temperature in Kelvin
  real(rk) :: es             ! SVP in Pa

  ! (good down to 110 K)
  es = exp(9.550426_rk - (5723.265_rk / t) + (3.53068_rk * log(t)) &
       - (0.00728332_rk * t))

end function MurphyKoop_svp_ice

! Taken from CLUBB, based on Bolton (1980).

pure function Bolton_svp_liquid(t) result(es)
  real(rk), parameter :: c1 = 611.2_rk
  real(rk), parameter :: c2 = 17.67_rk
  real(rk), parameter :: c3 = 29.65_rk

  real(rk), intent(in) :: t  ! Temperature in Kelvin
  real(rk) :: es             ! SVP in Pa

  es = c1*exp( (c2*(t - tmelt))/(t - c3) )

end function Bolton_svp_liquid

pure function Bolton_svp_ice(t) result(es)

  real(rk), intent(in) :: t  ! Temperature in Kelvin
  real(rk) :: es             ! SVP in Pa

  es = exp( 27.93603_rk - (6111.72784_rk/t) + (0.15215_rk*log(t)) )

end function Bolton_svp_ice

! "Flatau" scheme modified from CLUBB, based on:
!
! ``Polynomial Fits to Saturation Vapor Pressure'' Flatau, Walko,
!   and Cotton.  (1992)  Journal of Applied Meteorology, Vol. 31,
!   pp. 1507--1513

pure function Flatau_svp_liquid(t) result(es)

  real(rk), intent(in) :: t  ! Temperature in Kelvin
  real(rk) :: es             ! SVP in Pa

  real(rk), dimension(9), parameter :: coef = [ &
       6.11583699E+02_rk, 4.44606896E+01_rk, 1.43177157E+00_rk, &
       2.64224321E-02_rk, 2.99291081E-04_rk, 2.03154182E-06_rk, &
       7.02620698E-09_rk, 3.79534310E-12_rk,-3.21582393E-14_rk ]

  real(rk), parameter :: min_T_in_C = -85._rk

  real(rk) :: T_in_C

  T_in_C = max(t - tmelt, min_T_in_C)

  es = coef(1) + T_in_C * (coef(2) + T_in_C * (coef(3) + T_in_C * &
       (coef(4) + T_in_C * (coef(5) + T_in_C * (coef(6) + T_in_C * &
       (coef(7) + T_in_C * (coef(8) + T_in_C * coef(9))))))))

end function Flatau_svp_liquid

pure function Flatau_svp_ice(t) result(es)

  real(rk), intent(in) :: t  ! Temperature in Kelvin
  real(rk) :: es             ! SVP in Pa

  real(rk), dimension(9), parameter :: coef = [ &
       6.09868993E+02_rk, 4.99320233E+01_rk, 1.84672631E+00_rk, &
       4.02737184E-02_rk, 5.65392987E-04_rk, 5.21693933E-06_rk, &
       3.07839583E-08_rk, 1.05785160E-10_rk, 1.61444444E-13_rk ]

  real(rk), parameter :: min_T_in_C = -90._rk

  real(rk) :: T_in_C

  T_in_C = max(t - tmelt, min_T_in_C)

  es = coef(1) + T_in_C * (coef(2) + T_in_C * (coef(3) + T_in_C * &
       (coef(4) + T_in_C * (coef(5) + T_in_C * (coef(6) + T_in_C * &
       (coef(7) + T_in_C * (coef(8) + T_in_C * coef(9))))))))

end function Flatau_svp_ice

end module shr_wv_sat_mod
