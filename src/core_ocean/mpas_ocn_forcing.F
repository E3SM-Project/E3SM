! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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
!> \version SVN:$Id:$
!> \details
!>  This module contains driver routines for building the forcing arrays.
!
!-----------------------------------------------------------------------

module ocn_forcing

   use mpas_kind_types
   use mpas_grid_types
   use mpas_configure
   use mpas_timekeeping
   use ocn_forcing_coupled

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
             ocn_forcing_build_transmission_array, &
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
!
!  routine ocn_forcing_build_arrays
!
!> \brief   Determines the forcing arrays.
!> \author  Doug Jacobsen
!> \date    12/13/12
!> \version SVN:$Id$
!> \details 
!>  This routine computes the forcing arrays used later in MPAS.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_build_arrays(timeStamp, mesh, forcing, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type(MPAS_Time_type), intent(in) :: timeStamp

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(inout) :: &
         mesh          !< Input: mesh information

      type (forcing_type), intent(inout) :: forcing !< Input: Forcing information

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

      call ocn_forcing_coupled_build_arrays(timeStamp, mesh, forcing, err)

   !--------------------------------------------------------------------

   end subroutine ocn_forcing_build_arrays!}}}

!***********************************************************************
!
!  routine ocn_forcing_init
!
!> \brief   Initializes forcing module
!> \author  Doug Jacobsen
!> \date    12/13/12
!> \version SVN:$Id$
!> \details 
!>  This routine initializes the forcing modules.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      integer :: err1

      err = 0

      attenuationCoefficient = config_flux_attenuation_coefficient

      call ocn_forcing_coupled_init(err1)

      err = ior(err,err1)

   end subroutine ocn_forcing_init!}}}

!***********************************************************************
!
!  routine ocn_forcing_build_transmission_array
!
!> \brief   Transmission coefficient array for surface forcing.
!> \author  Doug Jacobsen
!> \date    10/03/2013
!> \version SVN:$Id$
!> \details 
!>  This subroutine builds the transmission coefficient array for use in
!>  applying surface fluxes deeper than the surface layer.
!
!-----------------------------------------------------------------------

    subroutine ocn_forcing_build_transmission_array(mesh, state, forcing, err)!{{{
        type (mesh_type), intent(in) :: mesh !< Input: Mesh information
        type (state_type), intent(in) :: state !< Input: State information
        type (forcing_type), intent(inout) :: forcing !< Input/Output: Forcing information
        integer, intent(out) :: err !< Output: Error code

        !************************************************
        !
        ! Local Variables
        !
        !************************************************

        real (kind=RKIND) :: zTop, zBot, transmissionCoeffTop, transmissionCoeffBot

        real (kind=RKIND), dimension(:,:), pointer :: layerThickness, transmissionCoefficients

        integer :: iCell, k
        integer :: nCells

        integer, dimension(:), pointer :: maxLevelCell

        err = 0

        nCells = mesh % nCells

        maxLevelCell             => mesh % maxLevelCell % array

        layerThickness           => state % layerThickness % array

        transmissionCoefficients => forcing % transmissionCoefficients % array

        do iCell = 1, nCells
           zTop = 0.0_RKIND
           transmissionCoeffTop = ocn_forcing_transmission(zTop)
           do k = 1, maxLevelCell(iCell)
              zBot = zTop - layerThickness(k,iCell)
              transmissionCoeffBot = ocn_forcing_transmission(zBot)

              transmissionCoefficients(k, iCell) = transmissionCoeffTop - transmissionCoeffBot

              zTop = zBot
              transmissionCoeffTop = transmissionCoeffBot
           end do
        end do

    end subroutine ocn_forcing_build_transmission_array!}}}

!***********************************************************************
!
!  real function ocn_forcing_transmission
!
!> \brief   Transmission coefficient for surface forcing.
!> \author  Doug Jacobsen
!> \date    05/03/2013
!> \version SVN:$Id$
!> \details 
!>  This function computes and returns the transmission coefficient for surface
!>  forcing based on depth. It uses an exponential decay function to determine the
!>  coefficients.
!
!-----------------------------------------------------------------------

   real (kind=RKIND) function ocn_forcing_transmission(z)!{{{
      real (kind=RKIND), intent(in) :: z 


      ocn_forcing_transmission = exp(z/attenuationCoefficient)

   end function ocn_forcing_transmission!}}}

!***********************************************************************

end module ocn_forcing

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
