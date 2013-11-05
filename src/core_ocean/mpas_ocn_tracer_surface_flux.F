! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_surface_flux
!
!> \brief MPAS ocean tracer surface flux
!> \author Doug Jacobsen
!> \date   12/17/12
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  surface flux tendencies.  
!
!-----------------------------------------------------------------------

module ocn_tracer_surface_flux

   use mpas_grid_types
   use mpas_configure

   use ocn_forcing

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

   public :: ocn_tracer_surface_flux_tend, &
             ocn_tracer_surface_flux_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_surface_flux_tend
!
!> \brief   Computes tendency term for surface fluxes
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \version SVN:$Id$
!> \details 
!>  This routine computes the tendency for tracers based on surface fluxes.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_flux_tend(mesh, transmissionCoefficients, layerThickness, surfaceTracerFlux, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        layerThickness !< Input: Layer thickness

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        surfaceTracerFlux !< Input: surface tracer fluxes

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        transmissionCoefficients !< Input: Coefficients for the application of surface fluxes

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend          !< Input/Output: velocity tendency

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, nCells, k, iTracer, nTracers, nVertLevels
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellMask

      real (kind=RKIND) :: remainingFlux


      err = 0

      nCells = mesh % nCells
      nVertLevels = mesh % nVertLevels
      nTracers = size(tend, dim=1)

      maxLevelCell => mesh % maxLevelCell % array
      cellMask => mesh % cellMask % array

      do iCell = 1, nCells
        remainingFlux = 1.0_RKIND
        do k = 1, maxLevelCell(iCell)
          remainingFlux = remainingFlux - transmissionCoefficients(k, iCell)

          do iTracer = 1, nTracers
            tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + cellMask(k, icell) * surfaceTracerFlux(iTracer, iCell) * transmissionCoefficients(k, iCell)
          end do
        end do

        if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
          do iTracer = 1, nTracers
            tend(iTracer, maxLevelCell(iCell), iCell) = tend(iTracer, maxLevelCell(iCell), iCell) + cellMask(k, iCell) * surfaceTracerFlux(iTracer, iCell) * remainingFlux
          end do
        end if
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_surface_flux_tend!}}}

!***********************************************************************
!
!  routine ocn_tracer_surface_flux_init
!
!> \brief   Initializes ocean tracer surface flux quantities
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to surface fluxes in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_surface_flux_init(err)!{{{

   !--------------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

   end subroutine ocn_tracer_surface_flux_init!}}}

!***********************************************************************

end module ocn_tracer_surface_flux

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
