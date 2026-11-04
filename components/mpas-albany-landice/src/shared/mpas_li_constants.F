! Copyright (c) 2015-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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

   ! Shared physical constants come from shr_const_mod in all builds, so E3SM and
   ! standalone runs use the same values. shr_const_mod takes its values from the
   ! E3SM Physical Constants Dictionary (share/util/pcd_const.F90). Standalone
   ! builds compile both from share/util (see shared/Makefile).
   use shr_const_mod, only: &
          cp_ice => SHR_CONST_CPICE, &            !< heat capacity of ice (J/kg/K)
          latent_heat_ice => SHR_CONST_LATICE, &  !< latent heat of melting of ice (J/kg)
          rho_water => SHR_CONST_RHOFW, &         !< density of fresh water (kg m^-3)
          cp_freshwater => SHR_CONST_CPFW, &      !< heat capacity of freshwater (J/kg/K)
          cp_seawater => SHR_CONST_CPSW, &        !< heat capacity of seawater (J/kg/K)
          pii => SHR_CONST_PI, &                  !< pi
          gravity => SHR_CONST_G, &               !< acceleration due to gravity (m s^-2)
          iceConductivity => SHR_CONST_CONDICE, & !< thermal conductivity of ice (W m^-1 K^-1)
          kelvin_to_celsius => SHR_CONST_TKFRZ, & !< offset to convert Kelvin to Celsius (K)
          seconds_per_day => SHR_CONST_CDAY, &    !< seconds in a calendar day (s)
          rho_ice => SHR_CONST_RHOICE, &          !< default density of ice (kg m^-3)
          shr_gas_constant => SHR_CONST_RGAS      !< universal gas constant (J kmol^-1 K^-1)

   implicit none
   save

   private :: shr_gas_constant

   real (kind=RKIND), parameter, public :: idealGasConstant = shr_gas_constant * 1.0e-3_RKIND
      !< ideal gas constant (J mol^-1 K^-1), converted from the per-kmol SHR_CONST_RGAS

   ! The constants below are empirical or specific to MALI's formulations and stay local.
   real (kind=RKIND), parameter, public :: kin_visc_water = 1.793e-6_RKIND  !< kinematic viscosity of freshwater at 0 deg C (m^2 s^-1)
   !> haline contraction coefficient of seawater near the freezing point (PSU^-1),
   !> used by the ISMIP7 local quadratic ice-shelf melt parameterization
   !> (Burgard et al. 2022, following Lazeroms et al. 2019)
   real (kind=RKIND), parameter, public :: haline_contraction_coeff = 7.86e-4_RKIND

   real (kind=RKIND), parameter, public :: &
        oceanFreezingTempSurface = -1.92_RKIND,  &  !< Freezing temperature of seawater (deg C) at surface pressure,
                                                    !< given S = 35 PSU
        oceanFreezingTempDepthDependence = -7.53e-4_RKIND  !< Rate of change of freezing temperature of seawater
                                                           !< with depth (deg m^-1), given S = 35 PSU
                                                           !< These values are from the Ocean Water Freezing Point Calculator,
                                                           !< http://www.csgnetwork.com/h2ofreezecalc.html (25 Nov. 2014)

   real (kind=RKIND), parameter, public :: iceMeltingPointPressureDependence = 7.9e-8_RKIND
      !< Clausius-Clapeyron constant: Dependence of ice melting point on pressure (K Pa^-1)

   !  conversion factors
   real (kind=RKIND), parameter, public :: scyr = 365.0_RKIND * seconds_per_day  !< seconds in a 365-day year

!***********************************************************************


!***********************************************************************

end module li_constants
