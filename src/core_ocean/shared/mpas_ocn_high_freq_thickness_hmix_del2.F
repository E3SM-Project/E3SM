! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_high_freq_thickness_hmix_del2
!
!> \brief MPAS ocean horizontal high_freq_thickness mixing driver
!> \author  Mark Petersen
!> \date    July 2013
!> \details
!>  This module contains the main driver routine for computing 
!>  horizontal mixing tendencies for high frequency thickness mixing
!
!-----------------------------------------------------------------------

module ocn_high_freq_thickness_hmix_del2

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

   public :: ocn_high_freq_thickness_hmix_del2_tend, &
             ocn_high_freq_thickness_hmix_del2_init


   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical, pointer :: config_use_highFreqThick_del2
   real (kind=RKIND), pointer :: config_highFreqThick_del2

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_high_freq_thickness_hmix_del2_tend
!
!> \brief   Computes Laplacian tendency term for horizontal highFreqThickness mixing
!> \author  Mark Petersen
!> \date    July 2013
!> \details 
!>  This routine computes the horizontal mixing tendency for 
!>  high frequency thickness
!>  based on current state using a Laplacian parameterization.
!
!-----------------------------------------------------------------------

   subroutine ocn_high_freq_thickness_hmix_del2_tend(meshPool, highFreqThickness, tend_highFreqThickness, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         meshPool          !< Input: mesh information

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         highFreqThickness !< Input: high frequency thickness

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(inout) :: &
         tend_highFreqThickness          !< Input/Output: high freq thickness tendency

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

      integer :: iCell, iEdge, cell1, cell2, i, k
      integer, pointer :: nCells, nVertLevels
      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell, hhf_turb_flux, flux, r_tmp

      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2

      err = 0

      if(.not.config_use_highFreqThick_del2) return

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)

      call mpas_pool_get_array(meshPool, 'maxLevelEdgeTop', maxLevelEdgeTop)
      call mpas_pool_get_array(meshPool, 'cellsOnEdge', cellsOnEdge)
      call mpas_pool_get_array(meshPool, 'edgeMask', edgeMask)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'dvEdge', dvEdge)
      call mpas_pool_get_array(meshPool, 'dcEdge', dcEdge)
      call mpas_pool_get_array(meshPool, 'meshScalingDel2', meshScalingDel2)

      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgesOnCell', edgesOnCell)
      call mpas_pool_get_array(meshPool, 'edgeSignOnCell', edgeSignOnCell)

      do iCell = 1, nCells
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOncell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            r_tmp = meshScalingDel2(iEdge) * config_highFreqThick_del2 * dvEdge(iEdge) / dcEdge(iEdge)
           
            do k = 1, maxLevelEdgeTop(iEdge)
               ! \nabla h^{hf} on edge
               hhf_turb_flux = highFreqThickness(k,cell2) - highFreqThickness(k,cell1)

               ! div(\kappa_{hf} \nabla h^{hf}) at cell center
               flux = hhf_turb_flux * edgeMask(k, iEdge) * r_tmp

               tend_highFreqThickness(k, iCell) = tend_highFreqThickness(k, iCell) - edgeSignOnCell(i, iCell) * flux * invAreaCell
            end do

         end do
      end do

   end subroutine ocn_high_freq_thickness_hmix_del2_tend!}}}

!***********************************************************************
!
!  routine ocn_high_freq_thickness_hmix_del2_init
!
!> \brief   Initializes horizontal highFreqThickness mixing
!> \author  Mark Petersen
!> \date    July 2013
!> \details 
!>  This routine initializes the module for horizontal mixing of 
!>  high frequency thickness
!>  
!
!-----------------------------------------------------------------------

   subroutine ocn_high_freq_thickness_hmix_del2_init(err)!{{{

      integer, intent(out) :: err !< Output: error flag

      err = 0

      call mpas_pool_get_config(ocnConfigs, 'config_use_highFreqThick_del2', config_use_highFreqThick_del2)
      call mpas_pool_get_config(ocnConfigs, 'config_highFreqThick_del2', config_highFreqThick_del2)

   end subroutine ocn_high_freq_thickness_hmix_del2_init!}}}

end module ocn_high_freq_thickness_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
