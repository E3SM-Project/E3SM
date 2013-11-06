! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_surface_flux
!
!> \brief MPAS ocean surface fluxes for thickness
!> \author Doug Jacobsen
!> \date   12/17/12
!> \version SVN:$Id:$
!> \details
!>  This module contains the routine for computing 
!>  tendencies for thickness from surface fluxes
!
!-----------------------------------------------------------------------

module ocn_thick_surface_flux

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

   public :: ocn_thick_surface_flux_tend, &
             ocn_thick_surface_flux_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: surfaceFluxOn


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_thick_surface_flux_tend
!
!> \brief   Computes tendency term from horizontal advection of thickness
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \version SVN:$Id$
!> \details 
!>  This routine computes the horizontal advection tendency for
!>  thicknes based on current state and user choices of forcings.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_surface_flux_tend(mesh, transmissionCoefficients, layerThickness, surfaceMassFlux, tend, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh              !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         transmissionCoefficients !< Input: Coefficients for the transmission of surface fluxes

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThickness   !< Input: Layer thickness

      real (kind=RKIND), dimension(:), intent(in) :: &
         surfaceMassFlux   !< Input: surface flux of mass


      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend          !< Input/Output: thickness tendency

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

      integer :: iCell, nCells, k, nVertLevels
      integer, dimension(:), pointer :: maxLevelCell
      integer, dimension(:,:), pointer :: cellMask

      real (kind=RKIND) :: remainingFlux

      err = 0

      if (.not. surfaceFluxOn) return

      maxLevelCell => mesh % maxLevelCell % array
      cellMask => mesh % cellMask % array

      nCells = mesh % nCells

      do iCell = 1, nCells
        remainingFlux = 1.0_RKIND
        do k = 1, maxLevelCell(iCell)
          remainingFlux = remainingFlux - transmissionCoefficients(k, iCell)

          tend(k, iCell) = tend(k, iCell) + cellMask(k, iCell) * (surfaceMassFlux(iCell)/config_density0) * transmissionCoefficients(k, iCell)
        end do

        if(maxLevelCell(iCell) > 0 .and. remainingFlux > 0.0_RKIND) then
          tend(maxLevelCell(iCell), iCell) = tend(maxLevelCell(iCell), iCell) + cellMask(maxLevelCell(iCell), iCell) * remainingFlux*surfaceMassFlux(iCell)/config_density0
        end if
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_thick_surface_flux_tend!}}}

!***********************************************************************
!
!  routine ocn_thick_surface_flux_init
!
!> \brief   Initializes ocean horizontal thickness surface fluxes
!> \author  Doug Jacobsen
!> \date    12/17/12
!> \version SVN:$Id$
!> \details 
!>  This routine initializes quantities related to thickness 
!>  surface fluxes in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_surface_flux_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      surfaceFluxOn = .true.

      if (config_disable_thick_sflux) then
         surfaceFluxOn = .false.
      end if

      if (config_forcing_type == trim('off')) then
         surfaceFluxOn = .false.
      end if


   !--------------------------------------------------------------------

   end subroutine ocn_thick_surface_flux_init!}}}

!***********************************************************************

end module ocn_thick_surface_flux

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
