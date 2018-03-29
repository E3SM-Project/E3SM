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
   use mpas_log
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

      err = 0

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

    subroutine ocn_forcing_build_fraction_absorbed_array(meshPool, statePool, diagnosticsPool, forcingPool, err, timeLevelIn)!{{{
        type (mpas_pool_type), intent(in) :: meshPool !< Input: Mesh information
        type (mpas_pool_type), intent(in) :: statePool !< Input: State information
        type (mpas_pool_type), intent(in) :: diagnosticsPool !< Input: Diagnostics information
        type (mpas_pool_type), intent(inout) :: forcingPool !< Input/Output: Forcing information
        integer, intent(out) :: err !< Output: Error code
        integer, intent(in), optional :: timeLevelIn

        !************************************************
        !
        ! Local Variables
        !
        !************************************************

        real (kind=RKIND) :: zTop, zBot, transmissionCoeffTop, transmissionCoeffBot

        real (kind=RKIND), dimension(:), pointer :: surfaceFluxAttenuationCoefficient
        real (kind=RKIND), dimension(:), pointer :: surfaceFluxAttenuationCoefficientRunoff
        real (kind=RKIND), dimension(:,:), pointer :: layerThickness, fractionAbsorbed, fractionAbsorbedRunoff

        integer :: iCell, k, timeLevel, nCells

        integer, dimension(:), pointer :: maxLevelCell, nCellsArray

        err = 0

        if (present(timeLevelIn)) then
           timeLevel = timeLevelIn
        else
           timeLevel = 1
        end if

        call mpas_pool_get_dimension(meshPool, 'nCellsArray', nCellsArray)

        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)

        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, timeLevel)

        call mpas_pool_get_array(diagnosticsPool, 'surfaceFluxAttenuationCoefficient', surfaceFluxAttenuationCoefficient)
        call mpas_pool_get_array(diagnosticsPool, 'surfaceFluxAttenuationCoefficientRunoff',  &
                                                   surfaceFluxAttenuationCoefficientRunoff)

        call mpas_pool_get_array(forcingPool, 'fractionAbsorbed', fractionAbsorbed)
        call mpas_pool_get_array(forcingPool, 'fractionAbsorbedRunoff', fractionAbsorbedRunoff)

        nCells = nCellsArray( 2 )

        do iCell = 1, nCells
           zTop = 0.0_RKIND
           transmissionCoeffTop = ocn_forcing_transmission(zTop, surfaceFluxAttenuationCoefficient(iCell))
           do k = 1, maxLevelCell(iCell)
              zBot = zTop - layerThickness(k,iCell)
              transmissionCoeffBot = ocn_forcing_transmission(zBot, surfaceFluxAttenuationCoefficient(iCell))

              fractionAbsorbed(k, iCell) = transmissionCoeffTop - transmissionCoeffBot

              zTop = zBot
              transmissionCoeffTop = transmissionCoeffBot
           end do
        end do

!  now do river runoff separately

        do iCell = 1, nCells
           zTop = 0.0_RKIND
           transmissionCoeffTop = ocn_forcing_transmission(zTop, surfaceFluxAttenuationCoefficientRunoff(iCell))
           do k = 1, maxLevelCell(iCell)
              zBot = zTop - layerThickness(k,iCell)
              transmissionCoeffBot = ocn_forcing_transmission(zBot, surfaceFluxAttenuationCoefficientRunoff(iCell))

              fractionAbsorbedRunoff(k, iCell) = transmissionCoeffTop - transmissionCoeffBot

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

   real (kind=RKIND) function ocn_forcing_transmission(z, attenuationCoefficient)!{{{
      real (kind=RKIND), intent(in) :: z, attenuationCoefficient

      ocn_forcing_transmission = exp( max(z / attenuationCoefficient, -100.0_RKIND) )

   end function ocn_forcing_transmission!}}}

!***********************************************************************

end module ocn_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
