
      module mo_params

      use shr_kind_mod, only : r8 => shr_kind_r8
      use physconst,     only : rearth

      implicit none

      integer, parameter :: kz = 100                  ! altitudes
      integer, parameter :: kw = 650                  ! wavelengths
!----------------------------------------------------------------------------
! 	... number of weighting functions
!           wavelength dependent
!----------------------------------------------------------------------------
      integer, parameter :: ks = 60
!----------------------------------------------------------------------------
!  	... wavelength and altitude dependent
!----------------------------------------------------------------------------
      integer, parameter :: kj = 70

!----------------------------------------------------------------------------
!  	... number of photorates to use from tuv
!----------------------------------------------------------------------------
      integer, parameter :: tuv_jmax = 31

!----------------------------------------------------------------------------
! 	... delta for adding points at beginning or end of data grids
!----------------------------------------------------------------------------
      real(r8), parameter :: deltax = 1.e-4_r8

      real(r8) :: largest
      real(r8) :: smallest

      end module mo_params
