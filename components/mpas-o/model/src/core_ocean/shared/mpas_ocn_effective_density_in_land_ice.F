! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_effective_density_in_land_ice
!
!> \brief MPAS ocean effective density in land ice
!> \author Xylar Asay-Davis
!> \date   10/03/2015
!> \details
!>  This module contains routines for computing the effective seawater
!>  density in land ice using Arhimedes' principle.
!
!-----------------------------------------------------------------------

module ocn_effective_density_in_land_ice

   use mpas_constants
   use mpas_kind_types
   use mpas_derived_types
   use mpas_pool_routines

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

   public :: ocn_effective_density_in_land_ice_update

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_effective_density_in_land_ice_update
!
!> \brief updates effective density in land ice
!> \author Xylar Asay-Davis
!> \date   10/03/2015
!> \details
!>  This routine updates the value of the effective seawater density
!>  displaced by land ice, based on Archimedes' principle.  The effective
!>  density is smoothed and extrapolated by averaging with nearest neighbors
!>  (cellsOnCell).
!
!-----------------------------------------------------------------------

   subroutine ocn_effective_density_in_land_ice_update(meshPool, forcingPool, statePool, scratchPool, ierr)!{{{

      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(in) :: meshPool !< Input: mesh information
      type (mpas_pool_type), intent(in) :: forcingPool !< Input: Forcing information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------
      type (mpas_pool_type), intent(inout) :: statePool !< Input/Output: state information
      type (mpas_pool_type), intent(inout) :: scratchPool !< Input/Output: scratch information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: ierr !< Output: Error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      character (len=StrKIND), pointer :: config_land_ice_flux_mode

      real (kind=RKIND), dimension(:), pointer :: landIcePressure, &
                                                  ssh, &
                                                  effectiveDensityCur, &
                                                  effectiveDensityNew, &
                                                  effectiveDensityScratch

      integer, dimension(:), pointer :: landIceMask

      type (field1DReal), pointer :: effectiveDensityField

      real (kind=RKIND) :: weightSum

      integer :: iCell, cell2, i
      integer, pointer :: nCells, nEdges

      integer, dimension(:,:), pointer :: cellsOnCell, cellMask

      integer, dimension(:), pointer :: nEdgesOnCell

      ierr = 0

      call mpas_pool_get_config(ocnConfigs, 'config_land_ice_flux_mode', config_land_ice_flux_mode)
      if ( (trim(config_land_ice_flux_mode) .ne. 'coupled') ) then
         return
      end if

      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'cellMask', cellMask)

      call mpas_pool_get_array(forcingPool, 'landIceMask', landIceMask)
      call mpas_pool_get_array(forcingPool, 'landIcePressure', landIcePressure)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 2)
      call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', effectiveDensityCur, 1)
      call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', effectiveDensityNew, 2)

      call mpas_pool_get_field(scratchPool, 'effectiveDensityScratch', effectiveDensityField)
      call mpas_allocate_scratch_field(effectiveDensityField, .true.)
      effectiveDensityScratch => effectiveDensityField % array

      !$omp do schedule(runtime)
      do iCell = 1, nCells
         ! TODO: should only apply to floating land ice, once wetting/drying is supported
         if(landIceMask(iCell) == 1) then
            ! land ice is present to update the effective density
            effectiveDensityScratch(iCell) = -landIcePressure(iCell)/(ssh(iCell)*gravity)
         else
            ! we copy the previous effective density
            effectiveDensityScratch(iCell) = effectiveDensityCur(iCell)
         end if
      end do
      !$omp end do

      !$omp do schedule(runtime) private(weightSum, i, cell2)
      do iCell = 1, nCells
         ! smooth/extrapolate by averaging with nearest neighbors
         weightSum = 1.0_RKIND
         effectiveDensityNew(iCell) = effectiveDensityScratch(iCell)
         do i = 1, nEdgesOnCell(iCell)
            cell2 = cellsOnCell(i,iCell)
            effectiveDensityNew(iCell) = effectiveDensityNew(iCell) &
               + cellMask(1,cell2)*effectiveDensityScratch(cell2)
            weightSum = weightSum + cellMask(1,cell2)
         end do
         effectiveDensityNew(iCell) = effectiveDensityNew(iCell)/weightSum
      end do
      !$omp end do
      call mpas_deallocate_scratch_field(effectiveDensityField, .true.)

   !--------------------------------------------------------------------

   end subroutine ocn_effective_density_in_land_ice_update !}}}

!***********************************************************************

end module ocn_effective_density_in_land_ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
