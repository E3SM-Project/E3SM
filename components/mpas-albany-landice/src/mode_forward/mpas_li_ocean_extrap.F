! Copyright (c) 2013-2018,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  li_ocean_extrap
!
!> \MPAS land-ice ocean-data extrapolation driver
!> \author Holly Han
!> \date   January 2022
!> \details
!>  This module contains the routines for extrapolating
!>  ocean data (e.g., temperature, salinity, thermal forcing)
!>  into ice draft
!
!-----------------------------------------------------------------------

module li_ocean_extrap

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_dmpar
   use mpas_log
   use li_mask
   use li_setup

   implicit none
   private

   !--------------------------------------------------------------------
   ! Public parameters
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! Public member functions
   !--------------------------------------------------------------------

   public :: li_ocean_extrap_solve

   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine li_ocean_extrap_solve
!
!> \brief   Initializes ocean extrapolation scheme
!> \author  Holly Han
!> \date    12 Jan 2023
!> \details
!>  This routine performs horizontal and vertical extrapolation
!>  of ocean data (e.g., temperature, salinity, thermal forcing)
!
!-----------------------------------------------------------------------

   subroutine li_ocean_extrap_solve(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      logical, pointer :: config_ocean_data_extrapolation
      real, pointer :: config_sea_level
      type (block_type), pointer :: block
	  type (mpas_pool_type), pointer :: scratchPool, geometryPool, meshPool, extrapOceanDataPool
	  real (kind=RKIND), dimension(:,:), pointer :: TFoceanOld, TFoceanNew
	  real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing
	  real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_zBndsOcean
	  integer, dimension(:), pointer :: origOceanMaskHoriz
	  type (field1DInteger), pointer :: origOceanMaskHorizField
      !integer, dimension(:,:) pointer :: validOceanMaskField,  availOceanMaskField
      integer, dimension(:,:), pointer ::  validOceanMask, availOceanMask !masks to pass to flood-fill routine
	  real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography
	  real (kind=RKIND) :: layerTop
	  integer, dimension(:), allocatable :: seedOceanMaskHorizOld,seedOceanMaskHoriz, growOceanMaskHoriz
      integer, pointer :: nCells, nCellsSolve, nISMIP6OceanLayers, nCellsExtra
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
	  integer :: iCell, jCell, iLayer, iNeighbor, iter


      ! No init is needed.
      err = 0

      call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrapolation', config_ocean_data_extrapolation)

      if ( config_ocean_data_extrapolation ) then
         ! call the extrapolation scheme
         call mpas_log_write('ocean data will be extrapolated into the MALI ice draft')

			! initialize the ocean data and mask fields
			call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
		    call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
			call mpas_pool_get_subpool(block % structs, 'extrapOceanData', extrapOceanDataPool)
			call mpas_pool_get_array(meshPool, 'nCells', nCells)
			call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
			call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
			call mpas_pool_get_array(geometryPool, 'nCellsSolve', nCellsSolve)
			call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
			call mpas_pool_get_array(geometryPool, 'thickness', thickness)
			call mpas_pool_get_array(geometryPool, 'bedTropography', bedTopography)
			call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_3dThermalForcing', TFoceanOld)
			call mpas_pool_get_array(extrapOceanDataPool, 'origOceanMaskHoriz', origOceanMaskHoriz)
			call mpas_pool_get_araay(extrapOceanDataPool, 'validOceanMask', validOceanMask)
			call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)

			validOceanMask(:,:) = 0
			call mpas_pool_get_array(extrapOceanDataPool, 'availOceanMask', availOceanMask)
			availOceanMask(:,:) = 0

            ! -------- Switch between horizontal and vertical extrapolation
            ! update availOceanMask
            ! create a 2D mask based on the open ocean + floating ice +  grounded ice mask and n-extra cells into the grounded ice
            seedOceanMaskHoriz(:) = 0
            seedOceanMaskHorizOld(:) = 0
            growOceanMaskHoriz(:) = 0

            ! define seedOceanMaskHoriz and growOceanMaskHoriz for horizontal floodfill
            do iCell = 1, nCellsSolve
               if ( .not. li_mask_is_grounded_ice(cellMask(iCell)) .and. bedTopography(iCell) < config_sea_level ) then
                  seedOceanMaskHoriz(iCell) = 1
               endif
               if ( li_mask_is_grounded_ice(cellMask(iCell)) .and. bedTopography(iCell) < config_sea_level ) then
                  growOceanMaskHoriz(iCell) = 1
               endif
            enddo

            ! Start horizontal floodfill
            seedOceanMaskHorizOld(:) = seedOceanMaskHoriz(:)

            ! go through the loop to get nCells extra into grounded ice
            do iter = 1, nCellsExtra
               do iCell = 1, nCellsSolve
                  if ( growOceanMaskHoriz(iCell) == 1 ) then
                     do iNeighbor = 1, nEdgesOnCell(iCell)
                        jCell = cellsOnCell(iNeighbor, iCell)
                        if ( seedOceanMaskHorizOld(iCell) == 0 .and. seedOceanMaskHoriz(jCell) == 1 ) then
                           seedOceanMaskHoriz(iCell) = 1
                        endif
                     enddo
                  endif
               enddo
               seedOceanMaskHorizOld(:) = seedOceanMaskHoriz(:)
               ! Update halos
               call mpas_timer_start("halo updates")
               call mpas_dmpar_field_halo_exch(domain, 'growOceanMaskHoriz')
               call mpas_dmpar_field_halo_exch(domain, 'seedOceanMaskHoriz')
               call mpas_timer_stop("halo updates")
            enddo
            deallocate(seedOceanMaskHorizOld)

            ! make it a 3D mask based on the topography (loop through nISMIP6OceanLayers)
            call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
            call mpas_pool_get_array(extrapOceanDataPool, 'ismip6shelfMelt_zBndsOcean', ismip6shelfMelt_zBndsOcean)

            availOceanMask(:,:) = 0
            validOceanMask(:,:) = 0
            do iCell = 1, nCellsSolve
               do iLayer = 1, nISMIP6OceanLayers
                  layerTop = ismip6shelfMelt_zBndsOcean(1, iLayer)
                  if ( (bedTopography(iCell) <  layerTop) .and. (seedOceanMaskHoriz(iCell) == 1) ) then
                     availOceanMask(iLayer,iCell) = 1
                  endif
                  if ( (bedTopography(iCell) <  layerTop) .and. (origOceanMaskHoriz(iCell) == 1) ) then
                     validOceanMask(iLayer,iCell) = 1
                  endif
               enddo
            enddo



            ! initialize validOceanMask for the current timestep
            ! define validOceanMask to indicate where ocean data exists at the current timestep (original valid ocean data masked out by the current ocean bathymetry => find where ice thickness is zero and where bedrock is below the sea level)
            ! initialize the ocean data field based on the valid ocean mask
            ! update halo for availOceanMask, validOceanMask and ocean data
            ! -------- Start iteration until validOceanMask is unchanged
                ! perform horizontal extrapolation
                ! perform vertical extrapolation
                ! check if validOceanMask is unchanged

			!do iCell = 1, nCells
				!if ( (bedTopography(iCell) < config_sea_level .and. thickness(iCell) == 0) ) then
					!availOceanMask(iCell) = 1
					!validOceanMask(iCell) =  TFoceanOld(iCell)

			! define availOceanMask to indicate where the ocean data can be extrapolated onto (where ocean => ice thickness is zero and bedrock is below the sea level)

			! call flood fill routine 

			! update the ocean data field with the updated seedMask
			!TFoceanNew(:) = validOceanMaskField(:)

         !! TFoceanOld(:) = TFoceanNew(:)
         ! call vertical extrapolation 
         ! switch between horizontal and vertical extrap ??				

      else
         ! do nothing
         call mpas_log_write('ocean data will NOT be extrapolated into the MALI ice draft')
      endif
   !--------------------------------------------------------------------

   end subroutine li_ocean_extrap_solve

end module li_ocean_extrap

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||