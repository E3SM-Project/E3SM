module SphericalCoordinates_mod
   implicit none
   private

   public :: SphericalCoordinates

   type SphericalCoordinates
      private
      real :: latitude
      real :: longitude
   contains
      procedure :: restrict
      procedure :: getLat
      procedure :: getLon
   end type SphericalCoordinates

   interface SphericalCoordinates
      module procedure newSphericalCoordinates
   end interface SphericalCoordinates

   integer, parameter :: MAX_LAT = +90
   integer, parameter :: MIN_LAT = -90
   integer, parameter :: MAX_LON = 360
   integer, parameter :: MIN_LON =   0

contains

   function newSphericalCoordinates(latitude, longitude) result(coordinates)
      type (SphericalCoordinates) :: coordinates
      real, intent(in) :: latitude
      real, intent(in) :: longitude

      coordinates%latitude = latitude
      coordinates%longitude = longitude

   end function newSphericalCoordinates

   subroutine restrict(this)
      class (SphericalCoordinates), intent(inout) :: this

      integer :: numPoleCrossing
      real :: absLat
      real :: absLon
      real :: adjustLat
      integer :: quadrant

      associate(lat => this%latitude, lon => this%longitude)

        absLat = abs(lat)
        adjustLat = absLat - (360 * floor(absLat / 360))
        quadrant = mod(floor(adjustLat/90), 4)

        select case(quadrant)
        case (0)
           lat = sign(adjustLat, lat)
        case (1)
           lat = sign(180 - adjustLat, lat)
           lon = lon + 180
        case (2)
           lat = - sign(180 - adjustLat, lat)
           lon = lon + 180
        case (3)
           lat = - sign(adjustLat-360, lat)
        end select

        absLon = abs(lon)
        lon = sign(absLon - 360 * floor(absLon/360), lon)
        if (lon < 0) lon = lon + 360

      end associate

   end subroutine restrict

   real function getLat(this) result(lat)
      class (SphericalCoordinates), intent(in) :: this
      lat = this%latitude
   end function getLat

   real function getLon(this) result(lon)
      class (SphericalCoordinates), intent(in) :: this
      lon = this%longitude
   end function getLon

end module SphericalCoordinates_mod
