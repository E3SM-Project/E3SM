! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS) (LA-CC-13-047)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_thick_ale
!
!> \brief MPAS ocean ALE thickness driver
!> \author Mark Petersen
!> \date   August 2013
!> \details
!>  This module contains the routines for computing
!>  diagnostic variables, and other quantities such as vertTransportVelocityTop.
!
!-----------------------------------------------------------------------

module ocn_thick_ale

   use mpas_grid_types
   use mpas_configure
   use mpas_constants
   use mpas_timer

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

   public :: ocn_ALE_thickness, &
             ocn_thick_ale_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   integer :: use_freq_filtered_thickness

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_ALE_thickness
!
!> \brief   Computes desired ALE thickness at new time
!> \author  Mark Petersen
!> \date    August 2013
!> \details 
!>  This routine computes the desired Arbitrary Lagrangian-Eulerian (ALE)
!>  thickness at the new time. It uses the ALE formulation, and includes 
!>  contributions from SSH variations (z-star), high-frequency divergence
!>  (z-tilde), and imposes a minimum layer thickness.
!
!-----------------------------------------------------------------------
   subroutine ocn_ALE_thickness(mesh, verticalMesh, oldSSH, div_hu_btr, newHighFreqThickness, dt, ALE_thickness, err)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mesh_type), intent(in) :: &
         mesh           !< Input: horizonal mesh information

      type (verticalMesh_type), intent(in) :: &
         verticalMesh   !< Input: vertical mesh information

      real (kind=RKIND), dimension(:), intent(in) :: &
         oldSSH,   &!< Input: sea surface height at old time
         div_hu_btr !< Input: thickness-weighted barotropic divergence

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         newHighFreqThickness   !< Input: high frequency thickness.  Alters ALE thickness.

      real (kind=RKIND), intent(in) :: &
         dt     !< Input: time step

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:), intent(out) :: &
         ALE_thickness     !< Output: desired thickness at new time

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, k, i, nCells, nVertLevels, kMax
      integer, dimension(:), pointer :: maxLevelCell

      real (kind=RKIND) :: thicknessSum, newSSH, remainder, newThickness, thicknessWithRemainder
      real (kind=RKIND), dimension(:), pointer :: vertCoordMovementWeights
      real (kind=RKIND), dimension(:), allocatable :: &
         SSH_ALE_thickness,      & !> ALE thickness alteration due to SSH (z-star)
         prelim_ALE_Thickness,   & !> ALE thickness at new time
         min_ALE_thickness_down, & !> ALE thickness alteration due to min/max thickness
         min_ALE_thickness_up      !> ALE thickness alteration due to min/max thickness
      real (kind=RKIND), dimension(:,:), pointer :: &
         restingThickness   !>  Layer thickness when the ocean is at rest, i.e. without SSH or internal perturbations.

      err = 0

      maxLevelCell      => mesh % maxLevelCell % array
      vertCoordMovementWeights => mesh % vertCoordMovementWeights % array

      restingThickness => verticalMesh % restingThickness % array

      nCells      = mesh % nCells
      nVertLevels = mesh % nVertLevels

      allocate(SSH_ALE_thickness(nVertLevels), prelim_ALE_thickness(nVertLevels), &
         min_ALE_thickness_down(nVertLevels), min_ALE_thickness_up(nVertLevels))

      !
      ! ALE thickness alteration due to SSH (z-star)
      !
      do iCell=1,nCells
         kMax = maxLevelCell(iCell)

         newSSH = oldSSH(iCell) - dt*div_hu_btr(iCell)
         thicknessSum = 1e-14  
         do k = 1, kMax
            SSH_ALE_Thickness(k) = newSSH * vertCoordMovementWeights(k) * restingThickness(k, iCell)
            thicknessSum = thicknessSum + vertCoordMovementWeights(k) * restingThickness(k, iCell)
         end do
         SSH_ALE_Thickness = SSH_ALE_Thickness / thicknessSum

         ! Note that restingThickness is nonzero, and remaining terms are perturbations about zero.
         ALE_Thickness(1:kMax, iCell) = & 
             restingThickness(1:kMax,iCell) &
           + SSH_ALE_Thickness(1:kMax) &
           + use_freq_filtered_thickness * newHighFreqThickness(1:kMax,iCell)

      enddo

      !
      ! ALE thickness alteration due to minimum and maximum thickness
      !
      if (config_use_min_max_thickness) then

         do iCell=1,nCells
            kMax = maxLevelCell(iCell)

            ! go down the column:
            prelim_ALE_Thickness(1:kMax) = ALE_Thickness(1:kMax, iCell)
            remainder = 0.0
            do k = 1, kMax
               newThickness = max(min(prelim_ALE_Thickness(k)+remainder,config_max_thickness_factor*restingThickness(k,iCell)),config_min_thickness)
               min_ALE_thickness_down(k) = newThickness - prelim_ALE_Thickness(k) 
               remainder = remainder - min_ALE_thickness_down(k)
            end do

            ! go back up the column:
            min_ALE_thickness_up(kMax) = 0.0
            prelim_ALE_Thickness(1:kMax) = prelim_ALE_Thickness(1:kMax) + min_ALE_thickness_down(1:kMax)
            do k = kMax-1,1,-1
               newThickness = max(min(prelim_ALE_Thickness(k)+remainder,config_max_thickness_factor*restingThickness(k,iCell)),config_min_thickness)
               min_ALE_thickness_up(k) = newThickness - prelim_ALE_Thickness(k) 
               remainder = remainder - min_ALE_thickness_up(k)
            end do
            min_ALE_thickness_up(1) = min_ALE_thickness_up(1) + remainder

            ALE_Thickness(1:kMax, iCell) = ALE_Thickness(1:kMax, iCell) + min_ALE_thickness_down(1:kMax) + min_ALE_thickness_up(1:kMax)

         enddo

      endif ! config_use_min_max_thickness

      deallocate(SSH_ALE_Thickness, prelim_ALE_thickness, min_ALE_thickness_down, min_ALE_thickness_up)

   end subroutine ocn_ALE_thickness!}}}


!***********************************************************************
!
!  routine ocn_thick_ale_init
!
!> \brief   Initializes flags used within diagnostics routines.
!> \author  Mark Petersen
!> \date    August 2013
!> \details 
!>  This routine initializes flags related to quantities computed within
!>  other diagnostics routines.
!
!-----------------------------------------------------------------------
   subroutine ocn_thick_ale_init(err)!{{{
      integer, intent(out) :: err !< Output: Error flag

      err = 0

      if (config_use_freq_filtered_thickness) then
         use_freq_filtered_thickness = 1
      else
         use_freq_filtered_thickness = 0
      endif

    end subroutine ocn_thick_ale_init!}}}

end module ocn_thick_ale

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
