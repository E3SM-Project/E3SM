! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_hadv
!
!> \brief MPAS ocean horizontal advection for thickness
!> \author Doug Jacobsen
!> \date   16 September 2011
!> \details
!>  This module contains the routine for computing 
!>  tendencies for thickness from horizontal advection
!
!-----------------------------------------------------------------------

module ocn_thick_hadv

   use mpas_grid_types
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

   public :: ocn_thick_hadv_tend, &
             ocn_thick_hadv_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------
   
   logical :: thickHadvOn

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_thick_hadv_tend
!
!> \brief   Computes tendency term from horizontal advection of thickness
!> \author  Doug Jacobsen
!> \date    15 September 2011
!> \details 
!>  This routine computes the horizontal advection tendency for
!>  thicknes based on current state and user choices of forcings.
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_hadv_tend(meshPool, normalVelocity, layerThicknessEdge, tend, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         normalVelocity    !< Input: velocity

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessEdge     !< Input: thickness at edge

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
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

      integer :: iEdge, cell1, cell2, k, i, iCell 
      integer, pointer :: nCells, nEdges, nVertLevels

      integer, dimension(:), pointer :: maxLevelEdgeBot, MaxLevelCell, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: flux, invAreaCell, invAreaCell1, invAreaCell2
      real (kind=RKIND), dimension(:), pointer :: dvEdge, areaCell

      !-----------------------------------------------------------------
      !
      ! call relevant routines for computing tendencies
      ! note that the user can choose multiple options and the 
      !   tendencies will be added together
      !
      !-----------------------------------------------------------------

      err = 0

      if(.not.thickHadvOn) return

      call mpas_pool_get_dimension(meshPool, 'nEdges', nEdges)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_array(meshPool, 'maxLevelEdgeBot', maxLevelEdgeBot)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      do iCell = 1, nCells
        invAreaCell = 1.0 / areaCell(iCell)
        do i = 1, nEdgesOnCell(iCell)
          iEdge = edgesOnCell(i, iCell)
          do k = 1, maxLevelEdgeBot(iEdge)
            flux = normalVelocity(k, iEdge) * dvEdge(iEdge) * layerThicknessEdge(k, iEdge)
            tend(k, iCell) = tend(k, iCell) + edgeSignOnCell(i, iCell) * flux * invAreaCell
          end do
        end do
      end do

   !--------------------------------------------------------------------

   end subroutine ocn_thick_hadv_tend!}}}

!***********************************************************************
!
!  routine ocn_thick_hadv_init
!
!> \brief   Initializes ocean horizontal thickness advection
!> \author  Doug Jacobsen
!> \date    16 September 2011
!> \details 
!>  This routine initializes quantities related to horizontal thickness 
!>  advection in the ocean. 
!
!-----------------------------------------------------------------------

   subroutine ocn_thick_hadv_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      logical, pointer :: config_disable_thick_hadv

      thickHadvOn = .true.

      call mpas_pool_get_config(ocnConfigs, 'config_disable_thick_hadv', config_disable_thick_hadv)

      if(config_disable_thick_hadv) thickHadvOn = .false.

      err = 0

   !--------------------------------------------------------------------

   end subroutine ocn_thick_hadv_init!}}}

!***********************************************************************

end module ocn_thick_hadv

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

! vim: foldmethod=marker
