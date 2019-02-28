!------------------------------------------------------------------------------
module m_constants
!------------------------------------------------------------------------------
!
! physical constants
!
real grav    ! gravitational acceleration
real sqrtg   ! square root of grav
real gsq     ! square of grav
real nu      ! kinematic viscosity of water
!
real d_water ! density of water
real d_air   ! density of air
!
! mathematical constants
!
real pi     ! circular constant, 3.1415...
real pi2    ! 2*pi
real pih    ! pi/2
real dera   ! conversion from degrees to radians
real rade   ! conversion from radians to degrees
real expmin ! min argument for exp. function to avoid underflow
real expmax ! max argument for exp. function to avoid overflow
real sqrt2  ! square root of 2 ~ 1.41
!
contains
!
!------------------------------------------------------------------------------
subroutine init_constants
!------------------------------------------------------------------------------
!
pi   = 4.*atan(1.)
pi2  = 2.*pi
pih  = 0.5*pi
dera = pi/180.
rade = 180./pi
!
expmin = -20.
expmax =  20.
!
!  physical constants
!
grav    = 9.81
sqrtg   = sqrt(grav)
gsq     = grav*grav
nu      = 1.e-6
d_air   = 1.2
d_water = 1027
!
end subroutine
!
end module
