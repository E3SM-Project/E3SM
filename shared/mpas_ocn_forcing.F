! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_forcing
!
!> \brief MPAS ocean forcing
!> \author Doug Jacobsen
!> \date   04/25/12
!> \details
!>  This module contains driver routines for building the forcing arrays.
!
!-----------------------------------------------------------------------

module ocn_forcing

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_timekeeping
   use mpas_io_units
   use mpas_dmpar
   use mpas_abort, only : mpas_dmpar_global_abort
   use ocn_forcing_bulk
   use ocn_forcing_restoring
   use ocn_constants

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

   public :: ocn_forcing_build_arrays, &
             ocn_forcing_init, &
             ocn_forcing_build_fraction_absorbed_array, &
             ocn_forcing_transmission

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) :: attenuationCoefficient

   logical :: restoringOn, bulkOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_forcing_build_arrays
!
!> \brief   Determines the forcing arrays.
!> \author  Doug Jacobsen
!> \date    12/13/12
!> \details 
!>  This routine computes the forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_build_arrays(meshPool, statePool, forcingPool, err, timeLevelIn)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         statePool, &  !< Input: State information
         meshPool      !< Input: mesh information

      integer, intent(in), optional :: timeLevelIn

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: forcingPool !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: timeLevel
      integer, pointer :: indexTemperature, indexSalinity
      integer, pointer :: indexSurfaceTemperatureFlux, indexSurfaceSalinityFlux

      real (kind=RKIND), dimension(:), pointer :: temperatureRestore, salinityRestore
      real (kind=RKIND), dimension(:,:), pointer :: surfaceTracerFlux
      real (kind=RKIND), dimension(:,:,:), pointer :: tracers

      if (present(timeLevelIn)) then
         timeLevel = timeLevelIn
      else
         timeLevel = 1
      end if

      if ( bulkOn ) then
         call ocn_forcing_bulk_build_arrays(meshPool, forcingPool, err)
      end if

      if ( restoringOn ) then
         call mpas_pool_get_dimension(statePool, 'index_temperature', indexTemperature)
         call mpas_pool_get_dimension(statePool, 'index_salinity', indexSalinity)

         call mpas_pool_get_dimension(forcingPool, 'index_surfaceTemperatureFlux', indexSurfaceTemperatureFlux)
         call mpas_pool_get_dimension(forcingPool, 'index_surfaceSalinityFlux', indexSurfaceSalinityFlux)

         call mpas_pool_get_array(statePool, 'tracers', tracers, timeLevel)

         call mpas_pool_get_array(meshPool, 'temperatureRestore', temperatureRestore)
         call mpas_pool_get_array(meshPool, 'salinityRestore', salinityRestore)

         call mpas_pool_get_array(forcingPool, 'surfaceTracerFlux', surfaceTracerFlux)

         call ocn_forcing_restoring_build_arrays(meshPool, indexTemperature, indexSalinity, &
              indexSurfaceTemperatureFlux, indexSurfaceSalinityFlux, &
              tracers, temperatureRestore, salinityRestore, &
              surfaceTracerFlux, err)
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_forcing_build_arrays!}}}

!***********************************************************************
!
!  routine ocn_forcing_init
!
!> \brief   Initializes forcing module
!> \author  Doug Jacobsen
!> \date    12/13/12
!> \details 
!>  This routine initializes the forcing modules.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      integer :: err1

      character (len=StrKIND), pointer :: config_forcing_type
      real (kind=RKIND), pointer :: config_flux_attenuation_coefficient

      err = 0
      err1 = 0

      call mpas_pool_get_config(ocnConfigs, 'config_flux_attenuation_coefficient', config_flux_attenuation_coefficient)
      call mpas_pool_get_config(ocnConfigs, 'config_forcing_type', config_forcing_type)

      attenuationCoefficient = config_flux_attenuation_coefficient

      if ( config_forcing_type == trim('bulk') ) then
         call ocn_forcing_bulk_init(err1)
         bulkOn = .true.
         restoringOn = .false.
      else if ( config_forcing_type == trim('restoring') ) then
         call ocn_forcing_restoring_init(err1)
         restoringOn = .true.
         bulkOn = .false.
      else if ( config_forcing_type == trim('off') ) then
         restoringOn = .false.
         bulkOn = .false.
      else
         write(stderrUnit, *) "ERROR: config_forcing_type not one of 'bulk' 'restoring', or 'off'."
         err = 1
         call mpas_dmpar_global_abort("ERROR: config_forcing_type not one of 'bulk', 'restoring', or 'off'.")
      end if

      err = ior(err,err1)

   end subroutine ocn_forcing_init!}}}

!***********************************************************************
!
!  routine ocn_forcing_build_fraction_absorbed_array
!
!> \brief   fraction absorbed coefficient array for surface forcing.
!> \author  Doug Jacobsen
!> \date    10/03/2013
!> \details 
!>  This subroutine builds the fractionAbsorbed coefficient array for use in
!>  applying surface fluxes deeper than the surface layer.
!
!-----------------------------------------------------------------------

    subroutine ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, forcingPool, err, timeLevelIn)!{{{
        type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
        type (mpas_pool_type), intent(in) :: statePool !< Input: State information
        type (mpas_pool_type), intent(inout) :: forcingPool !< Input/Output: Forcing information
        integer, intent(out) :: err !< Output: Error code
        integer, intent(in), optional :: timeLevelIn

        !************************************************
        !
        ! Local Variables
        !
        !************************************************

        real (kind=RKIND) :: zTop, zBot, transmissionCoeffTop, transmissionCoeffBot

        real (kind=RKIND), dimension(:,:), pointer :: layerThickness, fractionAbsorbed

        integer :: iCell, k, timeLevel
        integer, pointer :: nCells

        integer, dimension(:), pointer :: maxLevelCell

        err = 0

        if (present(timeLevelIn)) then
           timeLevel = timeLevelIn
        else
           timeLevel = 1
        end if

        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

        call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', fractionAbsorbed)

        do iCell = 1, nCells
           zTop = 0.0_RKIND
           transmissionCoeffTop = ocn_forcing_transmission(zTop)
           do k = 1, maxLevelCell(iCell)
              zBot = zTop - layerThickness(k,iCell)
              transmissionCoeffBot = ocn_forcing_transmission(zBot)

              fractionAbsorbed(k, iCell) = transmissionCoeffTop - transmissionCoeffBot

              zTop = zBot
              transmissionCoeffTop = transmissionCoeffBot
           end do
        end do

    end subroutine ocn_forcing_build_fraction_absorbed_array!}}}

!***********************************************************************
!
!  real function ocn_forcing_transmission
!
!> \brief   Transmission coefficient for surface forcing.
!> \author  Doug Jacobsen
!> \date    05/03/2013
!> \details 
!>  This function computes and returns the transmission coefficient for surface
!>  forcing based on depth. It uses an exponential decay function to determine the
!>  coefficients.
!
!-----------------------------------------------------------------------

   real (kind=RKIND) function ocn_forcing_transmission(z)!{{{
      real (kind=RKIND), intent(in) :: z 

      ocn_forcing_transmission = exp( z / attenuationCoefficient )

   end function ocn_forcing_transmission!}}}

!***********************************************************************

end module ocn_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
