! Copyright (c) 2015,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_constants
!
!> \brief MPAS land ice specific constants
!> \author Matthew Hoffman
!> \date   17 Feb. 2015
!> \details
!>  This module contains constants specific to the land ice model.
!
!-----------------------------------------------------------------------

module li_constants

   use mpas_grid_types
   use mpas_kind_types

#ifdef MPAS_CESM
   use shr_const_mod, only: &
          cp_ice => SHR_CONST_CPICE,&
          latent_heat_ice => SHR_CONST_LATICE,&
          triple_point => SHR_CONST_TKTRIP
   implicit none
   save

#else

   implicit none
   save

   ! physical constants
   real (kind=RKIND), parameter, public :: cp_ice = 2009.0_RKIND       !< heat capacity of ice (J/kg/K)
   real (kind=RKIND), parameter, public :: latent_heat_ice = 335.0d3   !< Latent heat of melting of ice (J/kg)
   real (kind=RKIND), parameter, public :: triple_point = 273.16_RKIND  !< Triple point of water (K)

#endif

   real (kind=RKIND), parameter, public :: idealGasConstant = 8.314_RKIND  !< ideal gas constant (J mol^-1 K^-1)

   !  conversion factors
   real (kind=RKIND), parameter, public :: kelvin_to_celsius = 273.15_RKIND    !< factor to convert Kelvin to Celsius


!***********************************************************************


!***********************************************************************

end module li_constants
