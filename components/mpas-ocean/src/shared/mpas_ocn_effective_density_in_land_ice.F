! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.io/license.html
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
   use ocn_config
   use ocn_mesh

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
!>  displaced by land ice, based on Archimedes' principle.  The
!>  effective density is smoothed and extrapolated by averaging with
!>  nearest neighbors (cellsOnCell).
!
!-----------------------------------------------------------------------

   subroutine ocn_effective_density_in_land_ice_update(forcingPool, &
                                                    statePool, ierr)!{{{

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: &
         forcingPool !< [in] Forcing information

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: &
         statePool !< [inout] state information

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      integer, intent(out) :: ierr !< [out] Error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------

      integer :: &
         i, iCell,  &! loop indices for neighbor loop and cell loop
         cell2       ! neighbor cell address

      real (kind=RKIND) :: &
         weightSum   ! local sum of weights for masked points

      real (kind=RKIND), dimension(:), pointer :: &
         ssh,                 &! sea surface height
         landIcePressure,     &! land ice pressure
         effectiveDensityCur, &! effective density at current time
         effectiveDensityNew   ! effective density at new time level
                                                  
      integer, dimension(:), pointer :: &
         landIceFloatingMask   ! mask for land ice presence on cell

      ! Scratch Arrays
      ! effectiveDensityScratch: effective seawater density in land
      ! ice before horizontal averaging units: kg m^{-3}
      real (kind=RKIND), dimension(:), allocatable :: &
         effectiveDensityScratch

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Set default return code and exit if not needed
      ierr = 0
      if ( (trim(config_land_ice_flux_mode) /= 'coupled') ) return

      ! Retrieve forcing and state variables
      call mpas_pool_get_array(forcingPool, 'landIceFloatingMask', &
                                             landIceFloatingMask)
      call mpas_pool_get_array(forcingPool, 'landIcePressure', &
                                             landIcePressure)
      call mpas_pool_get_array(statePool, 'ssh', ssh, 2)
      call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', &
                                           effectiveDensityCur, 1)
      call mpas_pool_get_array(statePool, 'effectiveDensityInLandIce', &
                                           effectiveDensityNew, 2)
      !$acc enter data copyin(landIceFloatingMask, landIcePressure, ssh, &
      !$acc                   effectiveDensityCur, effectiveDensityNew)

      allocate(effectiveDensityScratch(nCellsAll))
      !$acc enter data create(effectiveDensityScratch)

      ! Compute effective density in each cell
      ! TODO: should only apply to floating land ice,
      !       once wetting/drying is supported
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(effectiveDensityScratch, effectiveDensityCur, &
      !$acc            landIcePressure, landIceFloatingMask, ssh)
#else
      !$omp parallel
      !$omp do schedule(runtime)
#endif
      do iCell = 1, nCellsAll
         if (landIceFloatingMask(iCell) == 1) then
            ! land ice is present to update the effective density
            effectiveDensityScratch(iCell) = -landIcePressure(iCell)/ &
                                              (ssh(iCell)*gravity)
         else
            ! we copy the previous effective density
            effectiveDensityScratch(iCell) = effectiveDensityCur(iCell)
         end if
      end do
#ifndef MPAS_OPENACC
      !$omp end do
#endif

      ! smooth/extrapolate density by averaging with nearest neighbors
#ifdef MPAS_OPENACC
      !$acc parallel loop &
      !$acc    present(effectiveDensityNew, effectiveDensityScratch, &
      !$acc            nEdgesOnCell, cellsOnCell, cellMask) &
      !$acc    private(weightSum)
#else
      !$omp do schedule(runtime) private(weightSum, i, cell2)
#endif
      do iCell = 1, nCellsAll
         weightSum = 1.0_RKIND
         effectiveDensityNew(iCell) = effectiveDensityScratch(iCell)
         do i = 1, nEdgesOnCell(iCell)
            cell2 = cellsOnCell(i,iCell)
            effectiveDensityNew(iCell) = effectiveDensityNew(iCell) &
                                       + cellMask(1,cell2)* &
                                         effectiveDensityScratch(cell2)
            weightSum = weightSum + cellMask(1,cell2)
         end do
         effectiveDensityNew(iCell) = effectiveDensityNew(iCell)/ &
                                      weightSum
      end do
#ifndef MPAS_OPENACC
      !$omp end do
      !$omp end parallel
#endif

      !$acc exit data delete(effectiveDensityScratch)
      !$acc exit data delete(landIceFloatingMask, landIcePressure, ssh, &
      !$acc                  effectiveDensityCur, effectiveDensityNew)
      deallocate(effectiveDensityScratch)

   !--------------------------------------------------------------------

   end subroutine ocn_effective_density_in_land_ice_update !}}}

!***********************************************************************

end module ocn_effective_density_in_land_ice

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
