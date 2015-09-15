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
   use ocn_constants

   implicit none
   private
   save

   ! TRACER-CLEAN-UP
   ! Need to figure out what to do with absorption coefficient computation.

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

   public :: ocn_forcing_init, &
             ocn_forcing_build_fraction_absorbed_array, &
             ocn_forcing_transmission

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) :: attenuationCoefficient

!***********************************************************************

contains

!***********************************************************************

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

      real (kind=RKIND), pointer :: config_flux_attenuation_coefficient

      call mpas_pool_get_config(ocnConfigs, 'config_flux_attenuation_coefficient', config_flux_attenuation_coefficient)

      attenuationCoefficient = config_flux_attenuation_coefficient

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
