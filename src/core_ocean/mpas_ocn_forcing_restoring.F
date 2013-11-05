! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_forcing_restoring
!
!> \brief MPAS ocean restoring
!> \author Doug Jacobsen
!> \date   10/28/2013
!> \details
!>  This module contains routines for building surface flux arrays based on restoring.
!
!-----------------------------------------------------------------------

module ocn_forcing_restoring

   use mpas_grid_types
   use mpas_configure

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

   public :: ocn_forcing_restoring_build_arrays, &
             ocn_forcing_restoring_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   real (kind=RKIND) :: temperatureTimeScale, salinityTimeScale !< restoring timescales
   real (kind=RKIND) :: temperatureLengthScale, salinityLengthScale !< restoring timescales


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_forcing_restoring_build_arrays
!
!> \brief   Builds the forcing array for restoring
!> \author  Doug Jacobsen
!> \date    10/29/2013
!> \details 
!>  This routine builds the forcing array based on surface restoring.
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_restoring_build_arrays(mesh, indexT, indexS, indexTFlux, indexSFlux, tracers, temperatureRestoring, salinityRestoring, surfaceFluxes, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
        tracers !< Input: tracer quantities

      real (kind=RKIND), dimension(:), intent(in) :: &
        temperatureRestoring, & !< Input: Restoring values for temperature
        salinityRestoring !< Input: Restoring values for salinity

      integer, intent(in) :: indexT !< Input: index for temperature
      integer, intent(in) :: indexS !< Input: index for salinity
      integer, intent(in) :: indexTFlux !< Input: index for temperature flux
      integer, intent(in) :: indexSFlux !< Input: index for salinity flux

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
        surfaceFluxes !< Input: tracer quantities

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

      integer :: iCell, nCellsSolve, k

      real (kind=RKIND) :: invTemp, invSalinity

      err = 0

      nCellsSolve = mesh % nCellsSolve

      invTemp = 1.0 / (temperatureTimeScale * 86400.0)
      invSalinity = 1.0 / (salinityTimeScale * 86400.0)

      k = 1  ! restoring only in top layer
      do iCell=1,nCellsSolve
        surfaceFluxes(indexTFlux, iCell) = - temperatureLengthScale * (tracers(indexT, k, iCell) - temperatureRestoring(iCell)) * invTemp
        surfaceFluxes(indexSFlux, iCell) = - salinityLengthScale * (tracers(indexS, k, iCell) - salinityRestoring(iCell)) * invSalinity
      enddo

   !--------------------------------------------------------------------

   end subroutine ocn_forcing_restoring_build_arrays!}}}

!***********************************************************************
!
!  routine ocn_forcing_restoring_init
!
!> \brief   Initializes ocean surface restoring
!> \author  Doug Jacobsen
!> \date    10/29/2013
!> \details 
!>  This routine initializes a variety of quantities related to 
!>  restoring in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_forcing_restoring_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

      temperatureTimeScale = config_restoreT_timescale
      salinityTimeScale = config_restoreS_timescale
      temperatureLengthScale = config_restoreT_lengthscale
      salinityLengthScale = config_restoreS_lengthscale

   !--------------------------------------------------------------------

   end subroutine ocn_forcing_restoring_init!}}}

!***********************************************************************

end module ocn_forcing_restoring

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
