! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_equation_of_state
!
!> \brief MPAS ocean equation of state driver
!> \author Mark Petersen
!> \date   September 2011
!> \details
!>  This module contains the main driver routine for calling
!>  the equation of state.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_equation_of_state_linear
   use ocn_equation_of_state_jm
   use ocn_constants
   use mpas_io_units

   implicit none
   private
   save

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_equation_of_state_density, &
             ocn_equation_of_state_init, &
             ocn_freezing_temperature, &
             ocn_freezing_temperature_salinity_deriv

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: linearEos, jmEos


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state
!
!> \brief   Calls equation of state
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine calls the equation of state to update the density
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_density(statePool, diagnosticsPool, meshPool, scratchPool, nCells, k_displaced, & !{{{
                                            displacement_type, density, err, thermalExpansionCoeff, &
                                            salineContractionCoeff, timeLevelIn)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !  This module contains routines necessary for computing the density
   !  from model temperature and salinity using an equation of state.
   !
   ! Input: mesh - mesh metadata
   !        s - state: activeTracers
   !        k_displaced
   !
   !  If k_displaced==0, density is returned with no displacement
   !
   !  If k_displaced~=0, density is returned, and is for
   !  a parcel adiabatically displaced from its original level to level
   !  k_displaced.  When using the linear EOS, state % displacedDensity is
   !  still filled, but depth (i.e. pressure) does not modify the output.
   !
   ! Output: s - state: computed density
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      implicit none

      type (mpas_pool_type), intent(in) :: statePool
      type (mpas_pool_type), intent(inout) :: diagnosticsPool
      type (mpas_pool_type), intent(in) :: meshPool
      type (mpas_pool_type), intent(in) :: scratchPool !< Input/Output: Scratch structure
      integer, intent(in) :: nCells
      integer, intent(in), optional :: timeLevelIn
      type (mpas_pool_type), pointer :: tracersPool
      integer :: k_displaced
      character(len=*), intent(in) :: displacement_type
      real (kind=RKIND), dimension(:,:), intent(out) :: density
      integer, intent(out) :: err
      real (kind=RKIND), dimension(:,:), intent(out), optional :: &
         thermalExpansionCoeff,  &! Thermal expansion coefficient (alpha), defined as $-1/\rho d\rho/dT$ (note negative sign)
         salineContractionCoeff   ! Saline contraction coefficient (beta), defined as $1/\rho d\rho/dS$

      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:,:), pointer :: tracersSurfaceValue
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer :: iCell, k
      integer, pointer :: indexT, indexS
      type (dm_info) :: dminfo
      integer :: timeLevel

      err = 0

      call mpas_timer_start("equation of state")

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_array(diagnosticsPool, 'tracersSurfaceValue', tracersSurfaceValue)
      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
      call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', indexT)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', indexS)

      if (linearEos) then

         call ocn_equation_of_state_linear_density(meshPool, nCells, k_displaced, displacement_type, indexT, indexS, &
                                                   activeTracers, density, err, tracersSurfaceValue, &
                                                   thermalExpansionCoeff, salineContractionCoeff)

      elseif (jmEos) then

         call ocn_equation_of_state_jm_density(meshPool, scratchPool, nCells, k_displaced, displacement_type, indexT, indexS, &
                                               activeTracers, density, err, tracersSurfaceValue, thermalExpansionCoeff, &
                                               salineContractionCoeff)

      endif

      call mpas_timer_stop("equation of state")

   end subroutine ocn_equation_of_state_density!}}}

!***********************************************************************
!
!  routine ocn_equation_of_stateInit
!
!> \brief   Initializes ocean momentum horizontal mixing quantities
!> \author  Mark Petersen
!> \date    September 2011
!> \details
!>  This routine initializes a variety of quantities related to
!>  horizontal velocity mixing in the ocean. Since a variety of
!>  parameterizations are available, this routine primarily calls the
!>  individual init routines for each parameterization.
!
!----------------------------------------------------------------------

   subroutine ocn_equation_of_state_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      character (len=StrKIND), pointer :: config_eos_type

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_eos_type', config_eos_type)

      linearEos = .false.
      jmEos = .false.

      if (config_eos_type.eq.'linear') then
         linearEos = .true.
         call ocn_equation_of_state_linear_init(err)
      elseif (config_eos_type.eq.'jm') then
         jmEos = .true.
         call ocn_equation_of_state_jm_init(err)
      else
         write (stderrUnit,*) 'Invalid choice for config_eos_type.'
         write (stderrUnit,*) '  Choices are: linear, jm'
         err = 1
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_equation_of_state_init!}}}


!***********************************************************************
!
!  function ocn_freezing_temperature
!
!> \brief   Computes the freezing temperature of the ocean.
!> \author  Xylar Asay-Davis
!> \date    11/16/2016
!> \details
!>  This routine computes the freezing temperature of the ocean at a given
!>  salinity and pressure.  Different coefficients are used in the open ocean
!>  (and under sea ice) than in land-ice cavities.
!
!-----------------------------------------------------------------------
    real (kind=RKIND) function ocn_freezing_temperature(salinity, pressure, inLandIceCavity)!{{{
      real (kind=RKIND), intent(in) :: salinity !< Input: Salinity value of water for freezing temperature
      real (kind=RKIND), intent(in) :: pressure !< Input: Pressure value for freezing temperature
      logical, intent(in) :: inLandIceCavity !< Input: flag indicating if the freezing temperature is computed
                                             !         in land ice cavities or in open ocean

      real (kind=RKIND), pointer :: coeff_0
      real (kind=RKIND), pointer :: coeff_S
      real (kind=RKIND), pointer :: coeff_p
      real (kind=RKIND), pointer :: coeff_pS
      real (kind=RKIND), pointer :: reference_pressure
      real (kind=RKIND) :: pressureOffset

      if(inLandIceCavity) then
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_cavity_freezing_temperature_coeff_0', &
                                   coeff_0)
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_cavity_freezing_temperature_coeff_S', &
                                   coeff_S)
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_cavity_freezing_temperature_coeff_p', &
                                   coeff_p)
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_cavity_freezing_temperature_coeff_pS', &
                                   coeff_pS)
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_cavity_freezing_temperature_reference_pressure', &
                                   reference_pressure)
      else
         call mpas_pool_get_config(ocnConfigs, 'config_open_ocean_freezing_temperature_coeff_0', &
                                   coeff_0)
         call mpas_pool_get_config(ocnConfigs, 'config_open_ocean_freezing_temperature_coeff_S', &
                                   coeff_S)
         call mpas_pool_get_config(ocnConfigs, 'config_open_ocean_freezing_temperature_coeff_p', &
                                   coeff_p)
         call mpas_pool_get_config(ocnConfigs, 'config_open_ocean_freezing_temperature_coeff_pS', &
                                   coeff_pS)
         call mpas_pool_get_config(ocnConfigs, 'config_open_ocean_freezing_temperature_reference_pressure', &
                                   reference_pressure)
      end if


      pressureOffset = max(pressure - reference_pressure, 0.0_RKIND)

      ocn_freezing_temperature = coeff_0 &
         + coeff_S * salinity &
         + coeff_p * pressureOffset &
         + coeff_pS * pressureOffset * salinity

    end function ocn_freezing_temperature!}}}

!***********************************************************************
!
!  function ocn_freezing_temperature_salinity_deriv
!
!> \brief   Computes the freezing-temperature salinity derivative
!> \author  Xylar Asay-Davis
!> \date    11/16/2016
!> \details
!>  This routine computes the derivative of the freezing temperature of the ocean with
!>  respect to salinity at a given salinity and pressure. Different coefficients are
!>  used in the open ocean (and under sea ice) than in land-ice cavities.
!
!-----------------------------------------------------------------------
    real (kind=RKIND) function ocn_freezing_temperature_salinity_deriv(salinity, pressure, inLandIceCavity)!{{{
      real (kind=RKIND), intent(in) :: salinity !< Input: Salinity value of water for freezing temperature
      real (kind=RKIND), intent(in) :: pressure !< Input: Pressure value for freezing temperature
      logical, intent(in) :: inLandIceCavity !< Input: flag indicating if the freezing temperature is computed
                                             !         in land ice cavities or in open ocean

      real (kind=RKIND), pointer :: coeff_S
      real (kind=RKIND), pointer :: coeff_pS
      real (kind=RKIND), pointer :: reference_pressure
      real (kind=RKIND) :: pressureOffset

      if(inLandIceCavity) then
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_cavity_freezing_temperature_coeff_S', &
                                   coeff_S)
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_cavity_freezing_temperature_coeff_pS', &
                                   coeff_pS)
         call mpas_pool_get_config(ocnConfigs, 'config_land_ice_cavity_freezing_temperature_reference_pressure', &
                                   reference_pressure)
      else
         call mpas_pool_get_config(ocnConfigs, 'config_open_ocean_freezing_temperature_coeff_S', &
                                   coeff_S)
         call mpas_pool_get_config(ocnConfigs, 'config_open_ocean_freezing_temperature_coeff_pS', &
                                   coeff_pS)
         call mpas_pool_get_config(ocnConfigs, 'config_open_ocean_freezing_temperature_reference_pressure', &
                                   reference_pressure)
      end if


      pressureOffset = max(pressure - reference_pressure, 0.0_RKIND)

      ocn_freezing_temperature_salinity_deriv = coeff_S + coeff_pS * pressureOffset

    end function ocn_freezing_temperature_salinity_deriv!}}}

!***********************************************************************

end module ocn_equation_of_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
