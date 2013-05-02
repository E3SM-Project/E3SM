! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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

   logical :: useShortWave

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

   subroutine ocn_tracer_surface_flux_tend(grid, surfaceTracerFlux, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
        surfaceTracerFlux !< Input: surface tracer fluxes

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

      integer :: iCell, nCells, k, iTracer, nTracers

      integer, dimension(:), pointer :: maxLevelCell

      if(.not.useShortWave) return

      err = 0

      nCells = grid % nCells
      nTracers = size(tend, dim=1)

      maxLevelCell => grid % maxLevelCell % array

      do iCell = 1, nCells
        k = max(1, maxLevelCell(iCell))

        do iTracer = 1, nTracers
           tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + surfaceTracerFlux(iTracer, iCell)
        end do
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
      useShortWave = .false.

      if(config_use_coupled_forcing) then
         useShortWave = .true.
      end if

   end subroutine ocn_tracer_surface_flux_init!}}}

!***********************************************************************

end module ocn_tracer_surface_flux

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
