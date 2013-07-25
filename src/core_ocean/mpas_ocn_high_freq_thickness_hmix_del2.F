! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
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

   public :: ocn_high_freq_thickness_hmix_del2_tend

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

   subroutine ocn_high_freq_thickness_hmix_del2_tend(grid, highFreqThickness, tend_highFreqThickness, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         grid          !< Input: grid information

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

      integer :: iCell, iEdge, nCells, nVertLevels, cell1, cell2, i, k
      integer, dimension(:), pointer :: maxLevelEdgeTop, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnEdge, edgeMask, edgesOnCell, edgeSignOnCell

      real (kind=RKIND) :: invAreaCell, hhf_turb_flux, flux, r_tmp

      real (kind=RKIND), dimension(:), pointer :: areaCell, dvEdge, dcEdge
      real (kind=RKIND), dimension(:), pointer :: meshScalingDel2

      err = 0

      nCells = grid % nCells
      nVertLevels = grid % nVertLevels

      maxLevelEdgeTop => grid % maxLevelEdgeTop % array
      cellsOnEdge => grid % cellsOnEdge % array
      edgeMask => grid % edgeMask % array
      areaCell => grid % areaCell % array
      dvEdge => grid % dvEdge % array
      dcEdge => grid % dcEdge % array
      meshScalingDel2 => grid % meshScalingDel2 % array

      nEdgesOnCell => grid % nEdgesOnCell % array
      edgesOnCell => grid % edgesOnCell % array
      edgeSignOnCell => grid % edgeSignOnCell % array

      do iCell = 1, nCells
         invAreaCell = 1.0 / areaCell(iCell)
         do i = 1, nEdgesOncell(iCell)
            iEdge = edgesOnCell(i, iCell)
            cell1 = cellsOnEdge(1,iEdge)
            cell2 = cellsOnEdge(2,iEdge)

            r_tmp = meshScalingDel2(iEdge) * config_hhf_del2 * dvEdge(iEdge) / dcEdge(iEdge)
           
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

end module ocn_high_freq_thickness_hmix_del2

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
