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
!> \author Mark Petersen, modified Phil Jones
!> \date   September 2011, modified May 2018
!> \details
!>  This module contains the main driver routines for calling
!>  the equation of state and related quantities.
!
!-----------------------------------------------------------------------

module ocn_equation_of_state

   use mpas_timer
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use ocn_equation_of_state_linear
   use ocn_equation_of_state_jm
   use ocn_equation_of_state_wright
   use ocn_constants
   use ocn_config
   use ocn_diagnostics_variables
   use mpas_log

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

   interface ocn_equation_of_state_density
      module procedure ocn_equation_of_state_density_openacc
      module procedure ocn_equation_of_state_density_original
   end interface

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer ::                   &
      ocnEqStateChoice           ! user choice of eq state method

   integer, parameter ::        &! supported equation of state methods
      ocnEqStateTypeLinear = 1, &! linear equation of state
      ocnEqStateTypeJM     = 2, &! Jackett-McDougall eqn of state
      ocnEqStateTypeWright = 3   ! Wright eqn of state

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_equation_of_state
!
!> \brief   Computes density using an ocean equation of state
!> \author  Mark Petersen, modified by Phil Jones
!> \date    September 2011, modified May 2018
!> \details
!>  This routine calls the equation of state to compute a density
!>  from model temperature and salinity.
!>  If kDisplaced equals 0, density for the current vertical level is
!>           returned
!>  If kDisplaced is not 0, density is returned for a parcel
!>           adiabatically displaced from its original level to level
!>           kDisplaced.
!>  When using the linear EOS, these options are ignored since the
!>           density is independent of pressure/depth.
!
!-----------------------------------------------------------------------

   subroutine ocn_equation_of_state_density_openacc(statePool, meshPool, activeTracers, tracersSurfaceValue, &
                           nCells, kDisplaced, displacement_type, density, err,       &
                           thermalExpansionCoeff, salineContractionCoeff, timeLevelIn)
      !{{{
      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horiz mesh cells
         kDisplaced           ! target vert level for adiab displacement

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         activeTracers

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tracersSurfaceValue

      integer, intent(in), optional :: &
         timeLevelIn          ! time level for state variables

      character(len=*), intent(in) :: &
         displacement_type    ! choice for adiabatic displacement

      type (mpas_pool_type), intent(in) :: &
         statePool,           &! pool containing state variables
         meshPool              ! pool containing mesh quantities

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: &
         err                   ! returned error flag

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         density               ! computed density of sea water

      ! optional expansion coefficients
      ! Thermal expansion, defined as $-1/\rho d\rho/dT$ (note neg sign)
      ! Saline contraction, defined as $1/\rho d\rho/dS$
      real (kind=RKIND), dimension(:,:), intent(out), optional :: &
         thermalExpansionCoeff,  &! Thermal expansion coefficient (alpha)
         salineContractionCoeff   ! Saline contraction coefficient (beta)

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: tracersPool

      integer, dimension(:), pointer :: maxLevelCell
      integer :: nVertLevels, indexT, indexS
      integer, pointer :: indexTptr, indexSptr
      integer :: timeLevel

      !-----------------------------------------------------------------

      !*** initialize error flag, start timer
      err = 0

      call mpas_timer_start("equation of state")

      !*** initialize or extract relevant fields

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_subpool  (statePool,   'tracers', &
                                                 tracersPool)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                 indexTptr)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                 indexSptr)
      call mpas_pool_get_array    (meshPool,    'maxLevelCell', &
                                                 maxLevelCell)

      indexT      = indexTptr
      indexS      = indexSptr
      nVertLevels = size(density,1)

      !*** call specific equation of state based on user choice

      select case (ocnEqStateChoice)

      case (ocnEqStateTypeLinear)

         if (present(thermalExpansionCoeff) .or. &
             present(salineContractionCoeff)) then
            call ocn_equation_of_state_linear_density(nVertLevels,  &
                     nCells, kDisplaced, displacement_type,         &
                     indexT, indexS, activeTracers, density, err,   &
                     thermalExpansionCoeff, salineContractionCoeff, &
                     tracersSurfaceValue)
         else
            call ocn_equation_of_state_linear_density(nVertLevels, &
                     nCells, kDisplaced, displacement_type,        &
                     indexT, indexS, activeTracers, density, err,  &
                     tracersSurfaceValue)
         endif

      case (ocnEqStateTypeJM)

         if (present(thermalExpansionCoeff) .or. &
             present(salineContractionCoeff)) then
            call ocn_equation_of_state_jm_density(nVertLevels,        &
                     nCells, kDisplaced, displacement_type,           &
                     indexT, indexS, activeTracers, density, err,     &
                     thermalExpansionCoeff, salineContractionCoeff,   &
                     tracersSurfaceValue)

         else
            call ocn_equation_of_state_jm_density(nVertLevels,        &
                     nCells, kDisplaced, displacement_type,           &
                     indexT, indexS, activeTracers, density, err,     &
                     tracersSurfaceValue)
         endif

      case (ocnEqStateTypeWright)

         if (present(thermalExpansionCoeff) .or. &
             present(salineContractionCoeff)) then
            call ocn_equation_of_state_wright_density(nVertLevels,    &
                     nCells, kDisplaced, displacement_type,           &
                     indexT, indexS, activeTracers, zMid,             &
                     maxLevelCell, density, err,                      &
                     thermalExpansionCoeff, salineContractionCoeff,   &
                     tracersSurfaceValue)

         else
            call ocn_equation_of_state_wright_density(nVertLevels,    &
                     nCells, kDisplaced, displacement_type,           &
                     indexT, indexS, activeTracers, zMid,             &
                     maxLevelCell, density, err, tracersSurfaceValue)
         endif

      case default
         !*** error handled during init
      end select

      !*** stop timer and exit

      call mpas_timer_stop("equation of state")

      !-----------------------------------------------------------------

   end subroutine ocn_equation_of_state_density_openacc!}}}

   subroutine ocn_equation_of_state_density_original(statePool, meshPool, tracersSurfaceValue, &
                           nCells, kDisplaced, displacement_type, density, err,       &
                           thermalExpansionCoeff, salineContractionCoeff, timeLevelIn)
      !{{{
      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      integer, intent(in) :: &
         nCells,             &! number of horiz mesh cells
         kDisplaced           ! target vert level for adiab displacement

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         tracersSurfaceValue

      integer, intent(in), optional :: &
         timeLevelIn          ! time level for state variables

      character(len=*), intent(in) :: &
         displacement_type    ! choice for adiabatic displacement

      type (mpas_pool_type), intent(in) :: &
         statePool,           &! pool containing state variables
         meshPool              ! pool containing mesh quantities

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: &
         err                   ! returned error flag

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         density               ! computed density of sea water

      ! optional expansion coefficients
      ! Thermal expansion, defined as $-1/\rho d\rho/dT$ (note neg sign)
      ! Saline contraction, defined as $1/\rho d\rho/dS$
      real (kind=RKIND), dimension(:,:), intent(out), optional :: &
         thermalExpansionCoeff,  &! Thermal expansion coefficient (alpha)
         salineContractionCoeff   ! Saline contraction coefficient (beta)

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), pointer :: tracersPool

      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers
      integer, dimension(:), pointer :: maxLevelCell
      integer :: nVertLevels, indexT, indexS
      integer, pointer :: indexTptr, indexSptr
      integer :: timeLevel

      !-----------------------------------------------------------------

      !*** initialize error flag, start timer
      err = 0

      call mpas_timer_start("equation of state")

      !*** initialize or extract relevant fields

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      call mpas_pool_get_subpool  (statePool,   'tracers', &
                                                 tracersPool)
      call mpas_pool_get_array    (tracersPool, 'activeTracers', &
                                                 activeTracers, timeLevel)
      call mpas_pool_get_dimension(tracersPool, 'index_temperature', &
                                                 indexTptr)
      call mpas_pool_get_dimension(tracersPool, 'index_salinity', &
                                                 indexSptr)
      call mpas_pool_get_array    (meshPool,    'maxLevelCell', &
                                                 maxLevelCell)

      indexT      = indexTptr
      indexS      = indexSptr
      nVertLevels = size(density,1)

      !*** call specific equation of state based on user choice

      select case (ocnEqStateChoice)

      case (ocnEqStateTypeLinear)

         if (present(thermalExpansionCoeff) .or. &
             present(salineContractionCoeff)) then
            call ocn_equation_of_state_linear_density(nVertLevels,  &
                     nCells, kDisplaced, displacement_type,         &
                     indexT, indexS, activeTracers, density, err,   &
                     thermalExpansionCoeff, salineContractionCoeff, &
                     tracersSurfaceValue)
         else
            call ocn_equation_of_state_linear_density(nVertLevels, &
                     nCells, kDisplaced, displacement_type,        &
                     indexT, indexS, activeTracers, density, err,  &
                     tracersSurfaceValue)
         endif

      case (ocnEqStateTypeJM)

         if (present(thermalExpansionCoeff) .or. &
             present(salineContractionCoeff)) then
            call ocn_equation_of_state_jm_density(nVertLevels,        &
                     nCells, kDisplaced, displacement_type,           &
                     indexT, indexS, activeTracers, density, err,     &
                     thermalExpansionCoeff, salineContractionCoeff,   &
                     tracersSurfaceValue)

         else
            call ocn_equation_of_state_jm_density(nVertLevels,        &
                     nCells, kDisplaced, displacement_type,           &
                     indexT, indexS, activeTracers, density, err,     &
                     tracersSurfaceValue)
         endif

      case (ocnEqStateTypeWright)

         if (present(thermalExpansionCoeff) .or. &
             present(salineContractionCoeff)) then
            call ocn_equation_of_state_wright_density(nVertLevels,    &
                     nCells, kDisplaced, displacement_type,           &
                     indexT, indexS, activeTracers, zMid,             &
                     maxLevelCell, density, err,                      &
                     thermalExpansionCoeff, salineContractionCoeff,   &
                     tracersSurfaceValue)

         else
            call ocn_equation_of_state_wright_density(nVertLevels,    &
                     nCells, kDisplaced, displacement_type,           &
                     indexT, indexS, activeTracers, zMid,             &
                     maxLevelCell, density, err, tracersSurfaceValue)
         endif

      case default
         !*** error handled during init
      end select

      !*** stop timer and exit

      call mpas_timer_stop("equation of state")

      !-----------------------------------------------------------------

   end subroutine ocn_equation_of_state_density_original!}}}


!***********************************************************************
!
!  routine ocn_equation_of_state_init
!
!> \brief   Initializes ocean equation of state quantities
!> \author  Mark Petersen, modified by Phil Jones
!> \date    September 2011, modified May 2018
!> \details
!>  This routine initializes a variety of quantities related to
!>  the ocean equation of state for computing the density of sea water.
!
!----------------------------------------------------------------------

   subroutine ocn_equation_of_state_init(domain, err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      type (domain_type), intent(in) :: &
         domain            ! domain information primarily for vert grid

      !-----------------------------------------------------------------
      ! Output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: err

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------

      !*** initialize error flag

      err = 0

      !*** set use choice based on config and call appropriate
      !***  init routine

      select case (trim(config_eos_type))

      case ('linear','Linear','LINEAR')

         ocnEqStateChoice  = ocnEqStateTypeLinear
         call ocn_equation_of_state_linear_init(err)

      case ('jm','JM')

         ocnEqStateChoice  = ocnEqStateTypeJM
         call ocn_equation_of_state_jm_init(domain, err)

      case ('wright', 'Wright', 'WRIGHT')

         ocnEqStateChoice  = ocnEqStateTypeWright
         call ocn_equation_of_state_wright_init(domain, err)

      case default

         call mpas_log_write(&
         'Invalid choice for config_eos_type. Choices are: linear, jm, wright')
         err = 1

      end select

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

      real (kind=RKIND) :: coeff_0, coeff_S, coeff_p, coeff_pS, az1_liq
      real (kind=RKIND) :: coeff_mushy

      if(inLandIceCavity) then
         coeff_0 = varinp_land_ice_cavity_freezing_temperature_coeff_0
         coeff_S = varinp_land_ice_cavity_freezing_temperature_coeff_S
         coeff_p = varinp_land_ice_cavity_freezing_temperature_coeff_p
         coeff_pS = varinp_land_ice_cavity_freezing_temperature_coeff_pS
         coeff_mushy = 0.0_RKIND
      else
         coeff_0 = varinp_open_ocean_freezing_temperature_coeff_0
         coeff_S = varinp_open_ocean_freezing_temperature_coeff_S
         coeff_p = varinp_open_ocean_freezing_temperature_coeff_p
         coeff_pS = varinp_open_ocean_freezing_temperature_coeff_pS
         az1_liq = varinp_open_ocean_freezing_temperature_coeff_mushy_az1_liq
         coeff_mushy = 1.0_RKIND / az1_liq
      end if

      ocn_freezing_temperature = coeff_0 &
         + coeff_S * salinity &
         + coeff_p * pressure &
         + coeff_pS * pressure * salinity &
         + coeff_mushy * salinity / (1.0_RKIND - salinity / 1e3_RKIND)

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

      real (kind=RKIND) :: coeff_S, coeff_pS, az1_liq
      real (kind=RKIND) :: coeff_mushy


      if(inLandIceCavity) then
         coeff_S = varinp_land_ice_cavity_freezing_temperature_coeff_S
         coeff_pS = varinp_land_ice_cavity_freezing_temperature_coeff_pS
         coeff_mushy = 0.0_RKIND
      else
         coeff_S = varinp_open_ocean_freezing_temperature_coeff_S
         coeff_pS = varinp_open_ocean_freezing_temperature_coeff_pS
         az1_liq = varinp_open_ocean_freezing_temperature_coeff_mushy_az1_liq
         coeff_mushy = 1.0_RKIND / az1_liq
      end if

      ocn_freezing_temperature_salinity_deriv = coeff_S + coeff_pS * pressure &
         + coeff_mushy / (1.0_RKIND - salinity/1e3_RKIND)**2

    end function ocn_freezing_temperature_salinity_deriv!}}}

!***********************************************************************

end module ocn_equation_of_state

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
