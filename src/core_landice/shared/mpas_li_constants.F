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

   use mpas_derived_types
   use mpas_kind_types

#ifdef MPAS_ESM_SHR_CONST
   use shr_const_mod, only: &
          cp_ice => SHR_CONST_CPICE,&
          latent_heat_ice => SHR_CONST_LATICE,&
          triple_point => SHR_CONST_TKTRIP,&
          rho_water => SHR_CONST_RHOFW, &
          cp_freshwater => SHR_CONST_CPFW
   implicit none
   save

#else

   implicit none
   save

   ! physical constants
   real (kind=RKIND), parameter, public :: cp_ice = 2009.0_RKIND        !< heat capacity of ice (J/kg/K)
   real (kind=RKIND), parameter, public :: latent_heat_ice = 335.0d3    !< Latent heat of melting of ice (J/kg)
   real (kind=RKIND), parameter, public :: triple_point = 273.16_RKIND  !< Triple point of water (K)
   real (kind=RKIND), parameter, public :: rho_water = 1000.0_RKIND     !< Density of fresh water (kg m^-3)
   real (kind=RKIND), parameter, public :: cp_freshwater = 4.188e3_RKIND        !< heat capacity of freshwater (J/kg/K)

#endif

   real (kind=RKIND), parameter, public :: idealGasConstant = 8.314_RKIND  !< ideal gas constant (J mol^-1 K^-1)
   real (kind=RKIND), parameter, public :: iceConductivity = 2.1_RKIND  !< thermal conductivity of ice (W m^-1 K^-1)

   real (kind=RKIND), parameter, public :: &
        oceanFreezingTempSurface = -1.92_RKIND,  &         !< Freezing temperature of seawater (deg C) at surface pressure, given S = 35 PSU
        oceanFreezingTempDepthDependence = -7.53e-4_RKIND  !< Rate of change of freezing temperature of seawater with depth (deg m^-1), given S = 35 PSU
                                                           !< These values are from the Ocean Water Freezing Point Calculator,
                                                           !< http://www.csgnetwork.com/h2ofreezecalc.html (25 Nov. 2014)

   real (kind=RKIND), parameter, public :: &
        iceMeltingPointPressureDependence = 9.7456e-8_RKIND  ! Dependence of ice melting point on pressure (K Pa^-1)

   !  conversion factors
   real (kind=RKIND), parameter, public :: kelvin_to_celsius = 273.15_RKIND    !< factor to convert Kelvin to Celsius
   real (kind=RKIND), parameter, public :: scyr = 31536000.0_RKIND       !< seconds in a 365-day year; used for diagnostics
   real (kind=RKIND), parameter, public :: seconds_per_day = 86400.0_RKIND  !< Conversion from days to seconds

!***********************************************************************


!***********************************************************************

end module li_constants
