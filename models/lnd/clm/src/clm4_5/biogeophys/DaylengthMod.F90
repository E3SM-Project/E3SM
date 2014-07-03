module DaylengthMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Computes daylength
  !
  use shr_kind_mod   , only : r8 => shr_kind_r8
  use decompMod      , only : bounds_type

  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: daylength         ! function to compute daylength
  public :: InitDaylength     ! initialize daylength for all grid cells
  public :: UpdateDaylength   ! update daylength for all grid cells
  !
  ! !PRIVATE DATA MEMBERS:
  logical :: first_step       ! is this the first step since initialization?
  !
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  elemental real(r8) function daylength(lat, decl)
    !
    ! !DESCRIPTION:
    ! Computes daylength (in seconds)
    !
    ! Latitude and solar declination angle should both be specified in radians. decl must
    ! be strictly less than pi/2; lat must be less than pi/2 within a small tolerance.
    !
    ! !USES:
    use shr_infnan_mod, only : nan => shr_infnan_nan, &
                               assignment(=)
    use shr_const_mod , only : SHR_CONST_PI
    !
    ! !ARGUMENTS:
    real(r8), intent(in) :: lat    ! latitude (radians)
    real(r8), intent(in) :: decl   ! solar declination angle (radians)
    !
    ! !LOCAL VARIABLES:
    real(r8) :: my_lat             ! local version of lat, possibly adjusted slightly
    real(r8) :: temp               ! temporary variable

    ! number of seconds per radian of hour-angle
    real(r8), parameter :: secs_per_radian = 13750.9871_r8

    ! epsilon for defining latitudes "near" the pole
    real(r8), parameter :: lat_epsilon = 10._r8 * epsilon(1._r8)

    ! Define an offset pole as slightly less than pi/2 to avoid problems with cos(lat) being negative
    real(r8), parameter :: pole = SHR_CONST_PI/2.0_r8
    real(r8), parameter :: offset_pole = pole - lat_epsilon
    !-----------------------------------------------------------------------

    ! Can't SHR_ASSERT in an elemental function; instead, return a bad value if any
    ! preconditions are violated

    ! lat must be less than pi/2 within a small tolerance
    if (abs(lat) >= (pole + lat_epsilon)) then
       daylength = nan

    ! decl must be strictly less than pi/2
    else if (abs(decl) >= pole) then
       daylength = nan

    ! normal case
    else    
       ! Ensure that latitude isn't too close to pole, to avoid problems with cos(lat) being negative
       my_lat = min(offset_pole, max(-1._r8 * offset_pole, lat))

       temp = -(sin(my_lat)*sin(decl))/(cos(my_lat) * cos(decl))
       temp = min(1._r8,max(-1._r8,temp))
       daylength = 2.0_r8 * secs_per_radian * acos(temp) 
    end if

  end function daylength


  !-----------------------------------------------------------------------
  subroutine InitDaylength(bounds, declin, declinm1)
    !
    ! !DESCRIPTION:
    ! Initialize daylength for all grid cells, and initialize previous daylength.
    !
    ! This should be called with declin set at the value for the first model time step,
    ! and declinm1 at the value for the previous time step
    !
    ! !USES:
    use clmtype
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in) :: declin              ! solar declination angle for the first model time step (radians)
    real(r8), intent(in) :: declinm1            ! solar declination angle for the previous time step (radians)
    !
    !-----------------------------------------------------------------------

    associate(&
    lat       => grc%lat,       & ! Input:   [real(r8) (:)] latitude (radians)
    dayl      => gps%dayl,      & ! Output:  [real(r8) (:)] day length (s)
    prev_dayl => gps%prev_dayl, & ! Output:  [real(r8) (:)] day length from previous time step (s)

    begg      => bounds%begg  , & ! beginning grid cell index
    endg      => bounds%endg    & ! ending grid cell index
    )

    prev_dayl(begg:endg) = daylength(lat(begg:endg), declinm1)
    dayl(begg:endg) = daylength(lat(begg:endg), declin)

    first_step = .true.

    end associate

  end subroutine InitDaylength


  !-----------------------------------------------------------------------
  subroutine UpdateDaylength(bounds, declin)
    !
    ! !DESCRIPTION:
    ! Update daylength for all grid cells, and set previous daylength. This should be
    ! called exactly once per time step.
    !
    ! Assumes that InitDaylength has been called in initialization. This Update routine
    ! should NOT be called in initialization.
    !
    ! !USES:
    use clmtype
    !
    ! !ARGUMENTS:
    type(bounds_type), intent(in) :: bounds
    real(r8), intent(in) :: declin            ! solar declination angle (radians)
    !
    !-----------------------------------------------------------------------

    associate(&
    lat       => grc%lat,       & ! Input:  [real(r8) (:)] latitude (radians)
    dayl      => gps%dayl,      & ! InOut:  [real(r8) (:)] day length (s)
    prev_dayl => gps%prev_dayl, & ! Output: [real(r8) (:)] day length from previous time step (s)

    begg      => bounds%begg  , & ! beginning grid cell index
    endg      => bounds%endg    & ! ending grid cell index
    )

    ! In the first time step, we simply use dayl & prev_dayl that were set in
    ! initialization. (We do NOT want to run the normal code in that case, because that
    ! would incorrectly set prev_dayl to be the same as the current dayl in the first
    ! time step, because of the way prev_dayl is initialized.)
    if (first_step) then
       first_step = .false.
    else 
       prev_dayl(begg:endg) = dayl(begg:endg)
       dayl(begg:endg) = daylength(lat(begg:endg), declin)
    end if

    end associate

  end subroutine UpdateDaylength

end module DaylengthMod
