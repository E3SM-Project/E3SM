! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_cell_markers
!
!> \brief MPAS ocean cell marker
!> \author Doug Jacobsen
!> \date   03/20/2015
!> \details
!>  This module contains the routines for marking
!>  cells for removing
!
!-----------------------------------------------------------------------
module ocn_init_cell_markers

   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_timer

   use ocn_constants

   implicit none
   private

   public :: ocn_mark_north_boundary, ocn_mark_south_boundary
   public :: ocn_mark_east_boundary, ocn_mark_west_boundary
   public :: ocn_mark_maxlevelcell

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

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

   !***********************************************************************
   !
   !  routine ocn_mark_north_boundary
   !
   !> \brief   North boundary marker
   !> \author  Doug Jacobsen
   !> \date    03/30/2015
   !> \details
   !>  This routine marks cells along the north boundary of a domain for removal.
   !>  It can only be applied to a planar mesh. North-south is defined as the y direction.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_mark_north_boundary(meshPool, yMax, edgeMin, iErr)!{{{
      implicit none

      type (mpas_pool_type), intent(in) :: meshPool
      real (kind=RKIND), intent(in) :: yMax
      real (kind=RKIND), intent(in) :: edgeMin
      integer, intent(out) :: iErr

      real (kind=RKIND), dimension(:), pointer :: yCell
      integer, dimension(:), pointer :: cullCell

      logical, pointer :: on_a_sphere
      integer, pointer :: nCells

      integer :: iCell

      iErr = 0

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write( 'WARNING: Can only mark north boundaries of planar meshes. ' &
                                           // 'Skipping marking of cells...')

      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'cullCell', cullCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      if ( associated(cullCell) ) then
         do iCell = 1, nCells
            if ( yCell(iCell) > yMax - 0.8_RKIND * edgeMin ) then
               cullCell(iCell) = 1
            end if
         end do
      end if

   end subroutine ocn_mark_north_boundary!}}}

   !***********************************************************************
   !
   !  routine ocn_mark_south_boundary
   !
   !> \brief   south boundary marker
   !> \author  Doug Jacobsen
   !> \date    03/30/2015
   !> \details
   !>  This routine marks cells along the south boundary of a domain for removal.
   !>  It can only be applied to a planar mesh. north-south is defined as the y direction.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_mark_south_boundary(meshPool, yMin, edgeMin, iErr)!{{{
      implicit none

      type (mpas_pool_type), intent(in) :: meshPool
      real (kind=RKIND), intent(in) :: yMin
      real (kind=RKIND), intent(in) :: edgeMin
      integer, intent(out) :: iErr

      real (kind=RKIND), dimension(:), pointer :: yCell
      integer, dimension(:), pointer :: cullCell

      logical, pointer :: on_a_sphere
      integer, pointer :: nCells

      integer :: iCell

      iErr = 0

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write( 'WARNING: Can only mark north boundaries of planar meshes. ' &
                                           // 'Skipping marking of cells...')

      call mpas_pool_get_array(meshPool, 'yCell', yCell)
      call mpas_pool_get_array(meshPool, 'cullCell', cullCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      if ( associated(cullCell) ) then
         do iCell = 1, nCells
            if ( yCell(iCell) < yMin + 0.8_RKIND * edgeMin ) then
               cullCell(iCell) = 1
            end if
         end do
      end if

   end subroutine ocn_mark_south_boundary!}}}

   !***********************************************************************
   !
   !  routine ocn_mark_east_boundary
   !
   !> \brief   East boundary marker
   !> \author  Doug Jacobsen
   !> \date    03/30/2015
   !> \details
   !>  This routine marks cells along the east boundary of a domain for removal.
   !>  It can only be applied to a planar mesh. west-east is defined as the x direction.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_mark_east_boundary(meshPool, xMax, edgeMin, iErr)!{{{
      implicit none

      type (mpas_pool_type), intent(in) :: meshPool
      real (kind=RKIND), intent(in) :: xMax
      real (kind=RKIND), intent(in) :: edgeMin
      integer, intent(out) :: iErr

      real (kind=RKIND), dimension(:), pointer :: xCell
      integer, dimension(:), pointer :: cullCell

      logical, pointer :: on_a_sphere
      integer, pointer :: nCells

      integer :: iCell

      iErr = 0

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write( 'WARNING: Can only mark north boundaries of planar meshes. ' &
                                           // 'Skipping marking of cells...')

      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'cullCell', cullCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      if ( associated(cullCell) ) then
         do iCell = 1, nCells
            if ( xCell(iCell) > xMax - 0.8_RKIND * edgeMin ) then
               cullCell(iCell) = 1
            end if
         end do
      end if

   end subroutine ocn_mark_east_boundary!}}}

   !***********************************************************************
   !
   !  routine ocn_mark_west_boundary
   !
   !> \brief   West boundary marker
   !> \author  Doug Jacobsen
   !> \date    03/30/2015
   !> \details
   !>  This routine marks cells along the west boundary of a domain for removal.
   !>  It can only be applied to a planar mesh. west-east is defined as the x direction.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_mark_west_boundary(meshPool, xMin, edgeMin, iErr)!{{{
      implicit none

      type (mpas_pool_type), intent(in) :: meshPool
      real (kind=RKIND), intent(in) :: xMin
      real (kind=RKIND), intent(in) :: edgeMin
      integer, intent(out) :: iErr

      real (kind=RKIND), dimension(:), pointer :: xCell
      integer, dimension(:), pointer :: cullCell

      logical, pointer :: on_a_sphere
      integer, pointer :: nCells

      integer :: iCell

      iErr = 0

      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)

      if ( on_a_sphere ) call mpas_log_write( 'WARNING: Can only mark north boundaries of planar meshes. ' &
                                           // 'Skipping marking of cells...')

      call mpas_pool_get_array(meshPool, 'xCell', xCell)
      call mpas_pool_get_array(meshPool, 'cullCell', cullCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      if ( associated(cullCell) ) then
         do iCell = 1, nCells
            if ( xCell(iCell) < xMin + 0.8_RKIND * edgeMin ) then
               cullCell(iCell) = 1
            end if
         end do
      end if

   end subroutine ocn_mark_west_boundary!}}}

   !***********************************************************************
   !
   !  routine ocn_mark_maxlevelcell
   !
   !> \brief   MaxLevelCell cell marker
   !> \author  Doug Jacobsen
   !> \date    03/31/2015
   !> \details
   !>  This routine marks cells for removal that have maxLevelCell <= 0.
   !
   !-----------------------------------------------------------------------
   subroutine ocn_mark_maxlevelcell(meshPool, iErr)!{{{
      implicit none

      type (mpas_pool_type), intent(in) :: meshPool
      integer, intent(out) :: iErr

      integer, dimension(:), pointer :: cullCell, maxLevelCell

      logical, pointer :: on_a_sphere

      integer, pointer :: nCells

      integer :: iCell

      iErr = 0

      call mpas_pool_get_array(meshPool, 'cullCell', cullCell)
      call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      if ( associated(cullCell) ) then
         do iCell = 1, nCells
            if ( maxLevelCell(iCell) <= 0 ) then
               cullCell(iCell) = 1
            end if
         end do
      end if

   end subroutine ocn_mark_maxlevelcell!}}}

!***********************************************************************

end module ocn_init_cell_markers

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
