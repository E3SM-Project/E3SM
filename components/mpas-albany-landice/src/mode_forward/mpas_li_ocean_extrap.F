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
      real (kind=RKIND) :: layerTop
      real (kind=RKIND), dimension(:,:), pointer :: TFocean, TFoceanOld
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_zOcean
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing, ismip6shelfMelt_zBndsOcean
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography
      integer, dimension(:), pointer :: origOceanMaskHoriz
      integer, dimension(:,:), pointer :: validOceanMask, availOceanMask !masks to pass to flood-fill routine
      integer, dimension(:), pointer :: seedOceanMaskHoriz, growOceanMaskHoriz, seedOceanMaskHorizInit
      integer, dimension(:), allocatable :: seedOceanMaskHorizOld
      integer, pointer :: nCells, nCellsSolve, nISMIP6OceanLayers, nCellsExtra
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      integer :: iCell, jCell, iLayer, iNeighbor, iter, err_tmp
      integer :: newMaskCountLocalAccum_horiz

      ! No init is needed.
      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrapolation', config_ocean_data_extrapolation)

      if ( config_ocean_data_extrapolation ) then
         ! call the extrapolation scheme
         call mpas_log_write('ocean data will be extrapolated into the MALI ice draft')
         block => domain % blocklist

        ! initialize the ocean data and mask fields
         call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
         call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrap_ncells_extra', nCellsExtra)
         call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
         call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
         call mpas_pool_get_subpool(block % structs, 'extrapOceanData', extrapOceanDataPool)
         call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
         call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
         call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
         call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
         call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
         call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
         call mpas_pool_get_array(geometryPool, 'thickness', thickness)
         call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_3dThermalForcing', ismip6shelfMelt_3dThermalForcing)
         call mpas_pool_get_array(geometryPool, 'origOceanMaskHoriz', origOceanMaskHoriz)
         call mpas_pool_get_array(extrapOceanDataPool, 'validOceanMask', validOceanMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'availOceanMask', availOceanMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'seedOceanMaskHoriz', seedOceanMaskHoriz)
         call mpas_pool_get_array(extrapOceanDataPool, 'growOceanMaskHoriz', growOceanMaskHoriz)
         call mpas_pool_get_array(extrapOceanDataPool, 'seedOceanMaskHorizInit', seedOceanMaskHorizInit)
         call mpas_pool_get_array(extrapOceanDataPool, 'TFoceanOld', TFoceanOld)
         call mpas_pool_get_array(extrapOceanDataPool, 'TFocean', TFocean)

         ! create a 2D mask based on the open ocean + floating ice +  grounded ice mask and n-extra cells into the grounded ice
         allocate(seedOceanMaskHorizOld(nCells+1))
         seedOceanMaskHorizOld(:) = 0
         seedOceanMaskHoriz(:) = 0
         growOceanMaskHoriz(:) = 0
         ! define seedOceanMaskHoriz and growOceanMaskHoriz for horizontal floodfill
         do iCell = 1, nCellsSolve
            if ( .not. li_mask_is_grounded_ice(cellMask(iCell)) .and. bedTopography(iCell) < config_sea_level ) then
               seedOceanMaskHoriz(iCell) = 1
            endif
            if ( bedTopography(iCell) < config_sea_level ) then
               growOceanMaskHoriz(iCell) = 1
            endif
         enddo

         ! Start horizontal floodfill
         seedOceanMaskHorizInit(:) = seedOceanMaskHoriz(:) !HH: benchmark/debug purpose
         seedOceanMaskHorizOld(:) = seedOceanMaskHoriz(:)
         ! go through the loop to get nCells extra into grounded ice
         do iter = 1, nCellsExtra
            do iCell = 1, nCellsSolve
               if ( growOceanMaskHoriz(iCell) == 1 .and. seedOceanMaskHorizOld(iCell) == 0 ) then
                  do iNeighbor = 1, nEdgesOnCell(iCell)
                     jCell = cellsOnCell(iNeighbor, iCell)
                     if ( seedOceanMaskHorizOld(jCell) == 1 ) then
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
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zOcean', ismip6shelfMelt_zOcean)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zBndsOcean', ismip6shelfMelt_zBndsOcean)
         availOceanMask(:,:) = 0
         validOceanMask(:,:) = 0
         do iCell = 1, nCellsSolve
            do iLayer = 1, nISMIP6OceanLayers
               layerTop = ismip6shelfMelt_zBndsOcean(1, iLayer)
               if ( (seedOceanMaskHoriz(iCell) == 1) .and. (bedTopography(iCell) <  layerTop) ) then
                  availOceanMask(iLayer,iCell) = 1
               endif
               if ( (origOceanMaskHoriz(iCell) == 1) .and. (bedTopography(iCell) <  layerTop) ) then
                  validOceanMask(iLayer,iCell) = 1
               endif
            enddo
         enddo
         call mpas_log_write('==HH: updating halos for the avail/valid ocean masks')
         ! Update halos
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'availOceanMask')
         call mpas_dmpar_field_halo_exch(domain, 'validOceanMask')
         call mpas_timer_stop("halo updates")

         ! call the horizontal extrapolation routine
         call mpas_log_write('==HH: calling the horizontal extrapolation routine')
         call horizontal_extrapolation(domain, availOceanMask, validOceanMask, newMaskCountLocalAccum_horiz, TFocean, err_tmp)
         err = ior(err, err_tmp)
         call mpas_log_write('==HH:Horizontal extrapolation done successfully!')


         ! initialize the local loop and count for validOceanMask
         localLoopCount = 0
         newMaskCountLocal = 1
         newMaskCountLocalAccum = 0
         !third variable to keep track of the number of valid ocean mask
         TFocean(:,:) = TFoceanOld(:,:)
         call mpas_log_write('==HH: starting the wrapper loop through newMaskCountLocal')
         do while ( newMaskCountLocal > 0 )
            localLoopCount = localLoopCount + 1
            newMaskCountLocal = 0
            ! Algorithm 4: perform vertical extrapolation of the ocean data
            validOceanMaskOld(:,:) = validOceanMask(:,:)
            do while ( newMaskCountLocal > 0 )
               do iCell = 1, nCellsSolve
                  do iLayer = 2, nISMIP6OceanLayers
                     if ( (availOceanMask(iLayer,iCell) == 1) .and. (validOceanMaskOld(iLayer,iCell) == 0) ) then
                     ! ========= placeholder for vertical averaging ===========
                     endif
                     ! update validOceanMask and new mask count
                     if ( validOceanMaskOld(iLayer-1,iCell) == 1 ) then
                        validOceanMask(iLayer,iCell) = 1
                     endif
                  enddo
                  validOceanMaskOld(:,:) = validOceanMask(:,:)
               enddo
            enddo
         enddo
      else
         ! do nothing
         call mpas_log_write('ocean data will NOT be extrapolated into the MALI ice draft')
      endif
   !--------------------------------------------------------------------

   end subroutine li_ocean_extrap_solve
!-----------------------------------------------------------------------


!***********************************************************************
!***********************************************************************
! Private subroutines:
!***********************************************************************
!***********************************************************************



!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine horizontal_extrapolation
!
!> \brief Extrapolate validOceanMask horizontally
!> \author Holly Kyeore Han
!> \date   November 2023
!> \details
!> This routine extrapolates takes the initialized availOceanMask
!> and validOceanMask and extrapolates validOceanMask in horizontal
!> direction until the local new mask count stops updating.
!> The output of the routine is an updated validOceanMask field and
!> newMaskCountLocal.

!-----------------------------------------------------------------------

   subroutine horizontal_extrapolation(domain, availOceanMask, validOceanMask, newMaskCountGlobal, TFocean, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      integer, dimension(:,:), pointer, intent(in) :: availOceanMask

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object
      integer, dimension(:,:), pointer, intent(inout) :: validOceanMask
      integer, intent(inout) :: err !< Output: error flag
      real (kind=RKIND), dimension(:,:), pointer, intent(inout) :: TFocean

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: newMaskCountGlobal

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: scratchPool, geometryPool, meshPool, extrapOceanDataPool
      real (kind=RKIND) :: layerTop, TFsum, areaSum
      integer, dimension(:,:), allocatable :: validOceanMaskOld
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_zOcean
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing
      real (kind=RKIND), dimension(:,:), pointer :: TFoceanOld
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography, areaCell
      integer, pointer :: nCells, nCellsSolve, nISMIP6OceanLayers, nCellsExtra
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      integer :: iCell, jCell, iLayer, iNeighbor, iter
      integer :: localLoopCount, newValidCount, newMaskCountLocalAccum

      err = 0

      ! initialize the ocean data and mask fields
      block => domain % blocklist
      call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrap_ncells_extra', nCellsExtra)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'extrapOceanData', extrapOceanDataPool)
      call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_3dThermalForcing', ismip6shelfMelt_3dThermalForcing)
      call mpas_pool_get_array(extrapOceanDataPool, 'TFoceanOld', TFoceanOld)


      ! get initial 3D valid data based on the original ISMIP6 field
      !TFoceanOld(:,:) = ismip6shelfMelt_3dThermalForcing(:,:) * validOceanMask(:,:)
      TFoceanOld(:,:) = 1.0 !HH: for now, have the ocean TF field as the field of ones to make it easy to verify the horizonal/vertical averaging.
      ! perform horizontal extrapolation until the validOceanMask is unchanged
      allocate(validOceanMaskOld(nISMIP6OceanLayers,nCells+1))
      validOceanMaskOld(:,:) = validOceanMask(:,:)

      ! initialize the local loop and count for validOceanMask
      localLoopCount = 0
      newMaskCountGlobal = 1
      do while ( newMaskCountGlobal > 0 )
         localLoopCount = localLoopCount + 1
         newMaskCountLocalAccum = 0
         do iCell = 1, nCellsSolve
            do iLayer = 1, nISMIP6OceanLayers
               if ( (availOceanMask(iLayer,iCell) == 1) .and. (validOceanMaskOld(iLayer,iCell) == 0) ) then
                  TFsum = 0.0
                  areaSum = 0.0
                  newValidCount = 0
                  do iNeighbor = 1, nEdgesOnCell(iCell)
                     jCell = cellsOnCell(iNeighbor, iCell)
                     if ( validOceanMaskOld(iLayer,jCell) == 1 ) then
                        validOceanMask(iLayer,iCell) = 1
                        newValidCount = 1
                        TFsum = TFsum + (TFoceanOld(iLayer,jCell) * areaCell(jCell))
                        areaSum = areaSum + areaCell(jCell)
                     endif
                  enddo
                  ! Accumulate cells added locally until we do the next global reduce
                  newMaskCountLocalAccum = newMaskCountLocalAccum + newValidCount
                  ! perform averaging of the thermal forcing field
                  if ( areaSum == 0.0 ) then
                     TFocean(iLayer,iCell) = 0.0
                  else
                     TFocean(iLayer,iCell) = TFsum / areaSum ! HH: we might want to do this regardless.. i.e. move this out of the if loop
                  endif
               endif
            enddo
         enddo
         validOceanMaskOld(:,:) = validOceanMask(:,:)
         TFoceanOld(:,:) = TFocean(:,:) !HH: is this what we want?

         ! update halo for validOceanMask and ocean data
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'validOceanMask')
         call mpas_dmpar_field_halo_exch(domain, 'TFocean')
         call mpas_timer_stop("halo updates")
         ! update count of cells added to mask globally
         call mpas_dmpar_sum_int(domain % dminfo, newMaskCountLocalAccum, newMaskCountGlobal)
      enddo

      deallocate(validOceanMaskOld)

   end subroutine horizontal_extrapolation
!-----------------------------------------------------------------------
end module li_ocean_extrap

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
