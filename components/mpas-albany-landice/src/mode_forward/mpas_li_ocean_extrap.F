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
!> \author Holly Han, modified by Alex Hager
!> \date   January 2022 (modified September 2024)
!> \details
!>  This module contains the routines for extrapolating
!>  ocean data (e.g., temperature and salinity, or thermal forcing)
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
   
   ! Coefficients in Jenkins 2011 shelf melt parameterizations
   real (kind=RKIND), parameter :: gamma1 = -5.7E-2
   real (kind=RKIND), parameter :: gamma2 = 8.8E-2
   real (kind=RKIND), parameter :: gamma3 = 7.61E-4

   real (kind=RKIND), parameter :: cp = 3974 !specific heat (J/kg), using value from Cowton 2015
   real (kind=RKIND), parameter :: L = 3.34E5 !latent heat of melting
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
      use li_calving, only: li_flood_fill

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
      real(kind=RKIND), pointer :: invalid_value_TF
      logical, pointer :: config_recalculate_thermal_forcing
      character (len=StrKIND), pointer :: config_TF_param_type
      real(kind=RKIND), pointer :: config_sea_level
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: scratchPool, geometryPool, meshPool, extrapOceanDataPool
      real (kind=RKIND) :: layerTop
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing, ismip6shelfMelt_zBndsOcean
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_zOcean
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography
      integer, dimension(:,:), pointer :: orig3dOceanMask
      integer, dimension(:,:), pointer :: validOceanMask, validOceanMaskOrig, availOceanMask !masks to pass to flood-fill routine
      integer, dimension(:), pointer :: seedOceanMaskHoriz, growOceanMaskHoriz
      real (kind=RKIND), dimension(:,:), pointer :: MPAS_3dOceanSalinity, MPAS_3dOceanTemperature
      integer, dimension(:), allocatable :: seedOceanMaskHorizOld
      integer, pointer :: nCells, nCellsSolve, nISMIP6OceanLayers, nCellsExtra
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      integer :: iCell, jCell, iLayer, iNeighbor, iter, err_tmp
      integer :: GlobalLoopCount, newMaskCountGlobal
      type (field1dInteger), pointer :: seedMaskField
      type (field1dInteger), pointer :: growMaskField
      integer, dimension(:), pointer :: connectedMarineMask, growMask !masks to pass to flood-fill routine
      real (kind=RKIND), dimension(:,:,:), allocatable :: oceanField
      
      ! No init is needed.
      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrapolation', config_ocean_data_extrapolation)

      if ( config_ocean_data_extrapolation ) then
              ! << TO DO >>: If using "retreat" athermal forcing paramterization, then need to translate T/S instead of extrap routine

         ! call the extrapolation scheme
         call mpas_log_write('ocean data will be extrapolated into the MALI ice draft')
         block => domain % blocklist

        ! initialize the ocean data and mask fields
         call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
         call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrap_ncells_extra', nCellsExtra)
         call mpas_pool_get_config(liConfigs, 'config_recalculate_thermal_forcing', config_recalculate_thermal_forcing)
         call mpas_pool_get_config(liConfigs, 'config_TF_param_type', config_TF_param_type)
         call mpas_pool_get_config(liConfigs, 'config_invalid_value_TF', invalid_value_TF)
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
         call mpas_pool_get_array(extrapOceanDataPool, 'orig3dOceanMask', orig3dOceanMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'validOceanMask', validOceanMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'availOceanMask', availOceanMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'seedOceanMaskHoriz', seedOceanMaskHoriz)
         call mpas_pool_get_array(extrapOceanDataPool, 'growOceanMaskHoriz', growOceanMaskHoriz)
         
         allocate(oceanField(nISMIP6OceanLayers,nCells+1,2))
         
         ! create a 2D mask based on the open ocean + floating ice + grounded ice mask to be used for defining the location of the n-extra cells into the grounded ice
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

         ! Start horizontal floodfill to define locations of n-extra cells
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
            ! Note: halo updates could be skipped for a number of iterations equal to the
            ! number of halo rows as a possible performance improvement in the future.
            ! (And the loop above would need to be changed to from nCellsSolve to nCells)
            call mpas_timer_start("halo updates")
            call mpas_dmpar_field_halo_exch(domain, 'growOceanMaskHoriz')
            call mpas_dmpar_field_halo_exch(domain, 'seedOceanMaskHoriz')
            call mpas_timer_stop("halo updates")
         enddo
         deallocate(seedOceanMaskHorizOld)

         ! Calculate mask of connected ocean
         call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
         call mpas_pool_get_field(scratchPool, 'seedMask', seedMaskField)
         call mpas_allocate_scratch_field(seedMaskField, single_block_in = .true.)
         connectedMarineMask => seedMaskField % array
         connectedMarineMask(:) = 0
         call mpas_pool_get_field(scratchPool, 'growMask', growMaskField)
         call mpas_allocate_scratch_field(growMaskField, single_block_in = .true.)
         growMask => growMaskField % array
         growMask(:) = 0

         do iCell = 1, nCells
            ! seedMask = open ocean cells in contact with the domain boundary
            if ((bedTopography(iCell) < config_sea_level) .and. (thickness(iCell) == 0.0_RKIND)) then
               do iNeighbor = 1, nEdgesOnCell(iCell)
                  if (cellsOnCell(iNeighbor, iCell) == nCells + 1) then
                     connectedMarineMask(iCell) = 1
                     exit  ! no need to keep checking neighbors
                  endif
               enddo
            endif
            ! growMask - all marine cells
            if (bedTopography(iCell) < config_sea_level) then
               growMask(iCell) = 1
            endif
         enddo
         ! now create mask of all marine locations connected to open ocean - to be used below to screen out lakes
         call li_flood_fill(connectedMarineMask, growMask, domain)

         ! make it a 3D mask based on the topography (loop through nISMIP6OceanLayers)
         call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zOcean', ismip6shelfMelt_zOcean)
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zBndsOcean', ismip6shelfMelt_zBndsOcean)
         ! check for valid data
         do iLayer = 1, nISMIP6OceanLayers
            if (ismip6shelfMelt_zOcean(iLayer) >= 0.0_RKIND) then
               call mpas_log_write("ismip6shelfMelt_zOcean has invalid value of $r in layer $i", MPAS_LOG_ERR, &
                       realArgs=(/ismip6shelfMelt_zOcean(iLayer)/), intArgs=(/iLayer/))
               call mpas_log_write("ismip6shelfMelt_zOcean must have negative values because they represent " // &
                                   "depths below sea level.", MPAS_LOG_ERR)
               err = ior(err, 1)
            endif
            if ((ismip6shelfMelt_zBndsOcean(1,iLayer) > 0.0_RKIND) .or. &
                    (ismip6shelfMelt_zBndsOcean(1,iLayer) < ismip6shelfMelt_zOcean(iLayer))) then
               call mpas_log_write("ismip6shelfMelt_zBndsOcean(1,:) has invalid value of $r in layer $i", MPAS_LOG_ERR, &
                       realArgs=(/ismip6shelfMelt_zBndsOcean(1,iLayer)/), intArgs=(/iLayer/))
               call mpas_log_write("ismip6shelfMelt_zBndsOcean(1,:) must be less than or equal to zero " // & 
                                   "because it represents the upper bound of an ocean layer", MPAS_LOG_ERR)
               err = ior(err, 1)
            endif
            if ((ismip6shelfMelt_zBndsOcean(2,iLayer) >= 0.0_RKIND) .or. &
                    (ismip6shelfMelt_zBndsOcean(2,iLayer) > ismip6shelfMelt_zOcean(iLayer))) then
               call mpas_log_write("ismip6shelfMelt_zBndsOcean(2,:) has invalid value of $r in layer $i", MPAS_LOG_ERR, &
                       realArgs=(/ismip6shelfMelt_zBndsOcean(2,iLayer)/), intArgs=(/iLayer/))
               call mpas_log_write("ismip6shelfMelt_zBndsOcean(2,:) must be less than zero " // & 
                                   "because it represents the lower bound of an ocean layer", MPAS_LOG_ERR)
               err = ior(err, 1)
            endif
         enddo
         availOceanMask(:,:) = 0
         do iCell = 1, nCells
            do iLayer = 1, nISMIP6OceanLayers
               layerTop = ismip6shelfMelt_zBndsOcean(1, iLayer)
               if ( (seedOceanMaskHoriz(iCell) == 1) .and. (connectedMarineMask(iCell) == 1)) then
                  if (bedTopography(iCell) < layerTop) then
                     availOceanMask(iLayer,iCell) = 1
                  else
                     ! If recalculating thermal forcing with ISMIP6retreat or AMretreat methods, then ignore
                     ! topography and allow for translation of water properties directly to glacier face. Otherwise,
                     ! cut off availOceanMask at first cell below seaflooor; this cell is necessary to keep in
                     ! order to ensure interpolation between cells directly above and below sea level is possible 
                    
                     if ((config_recalculate_thermal_forcing) .and. &
                            ( (config_TF_param_type == 'ISMIP6retreat') .or. (config_TF_param_type == 'AMretreat'))) then
                        availOceanMask(iLayer,iCell) = 1
                     else
                        availOceanMask(iLayer,iCell) = 1
                        exit
                     endif
                  endif
               endif
            enddo
         enddo

         ! Make a copy of original mask to use for extending the mask during extrapolation
         validOceanMask(:,:) = orig3dOceanMask(:,:)

         ! initialize the extrapolated fields
         if (config_recalculate_thermal_forcing) then
            
            call mpas_pool_get_array(extrapOceanDataPool, 'MPAS_3dOceanTemperature', MPAS_3dOceanTemperature)
            call mpas_pool_get_array(extrapOceanDataPool, 'MPAS_3dOceanSalinity', MPAS_3dOceanSalinity)

            !combine fields into one variable to expedite extrapolation
            oceanField(:,:,1) = MPAS_3dOceanTemperature(:,:) * validOceanMask(:,:)
            oceanField(:,:,2) = MPAS_3dOceanSalinity(:,:) * validOceanMask(:,:)
         else

            oceanField(:,:,1) = ismip6shelfMelt_3dThermalForcing * validOceanMask(:,:)
            oceanField(:,:,2) = invalid_value_TF
         endif
         
         ! initialize the invalid data locations with fill value
         do iCell = 1, nCellsSolve
            do iLayer = 1, nISMIP6OceanLayers
               if ((connectedMarineMask(iCell) == 0) .and. (bedTopography(iCell) < config_sea_level)) then
                  ! Don't assign invalid value to lakes/inland seas disconnected from global ocean
                  ! Let them retain the existing value:
                  ! This will take on the valid ocean data value where it exists or
                  ! zero where valid ocean data does not exist
               elseif (validOceanMask(iLayer,iCell) == 0) then
                  ! everywhere else where valid ocean data does not exist, insert invalid value outside of validOceanMask
                  oceanField(iLayer,iCell,:) = invalid_value_TF
               endif
            enddo
         enddo

         ! deallocate scratch fields used for flood fill
         call mpas_deallocate_scratch_field(seedMaskField, single_block_in=.true.)
         call mpas_deallocate_scratch_field(growMaskField, single_block_in=.true.)

         ! flood-fill the valid ocean mask and TF field through
         ! horizontal and vertial extrapolation
         ! get initial 3D valid data based on the original ISMIP6 field
         newMaskCountGlobal = 1
         GlobalLoopCount = 0
         do while (newMaskCountGlobal > 0)
            newMaskCountGlobal = 0
            GlobalLoopCount = GlobalLoopCount + 1
            if (GlobalLoopCount > 20) then
               ! There will only be an additional time through this loop for each sill behind a previous sill
               ! so it should not need to alternate very many times
               call mpas_log_write("Ocean extrapolation has alternated between horizontal and vertical " // &
                       "extrapolation more than $i times.  Aborting", MPAS_LOG_ERR, intArgs=(/GlobalLoopCount/))
               err = ior(err, 1)
            endif
            ! call the horizontal extrapolation routine
            call mpas_timer_start("horizontal scheme")
            call horizontal_extrapolation(domain, availOceanMask, validOceanMask, orig3dOceanMask, oceanField, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop("horizontal scheme")
            ! call the vertical extrapolation routine
            call mpas_timer_start("vertical scheme")
            call vertical_extrapolation(domain, availOceanMask, validOceanMask, newMaskCountGlobal, oceanField, err_tmp)
            err = ior(err, err_tmp)
            call mpas_timer_stop("vertical scheme")

            if (err > 0) then
               call mpas_log_write("Ocean extrapolation main iteration loop has encountered an error", MPAS_LOG_ERR)
               return
            endif
         enddo

         ! Reassign extrapolated TF back to primary TF field
         if (config_recalculate_thermal_forcing) then
            MPAS_3dOceanTemperature(:,:) = oceanField(:,:,1)
            MPAS_3dOceanSalinity(:,:) = oceanField(:,:,2)
         else 
            ismip6shelfMelt_3dThermalForcing(:,:) = oceanField(:,:,1)
         endif

         deallocate(oceanField)

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
!> \author Holly Kyeore Han, modified by Alex Hager
!> \date   November 2023 (modified September 2024)
!> \details
!> This routine extrapolates takes the initialized availOceanMask
!> and validOceanMask and extrapolates validOceanMask in horizontal
!> direction until the local new mask count stops updating.
!> The output of the routine is an updated validOceanMask field and
!> newMaskCountLocal.

!-----------------------------------------------------------------------

   subroutine horizontal_extrapolation(domain, availOceanMask, validOceanMask, validOceanMaskOrig, oceanField, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------
      integer, dimension(:,:), pointer, intent(in) :: availOceanMask, validOceanMaskOrig

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object
      integer, dimension(:,:), pointer, intent(inout) :: validOceanMask
      integer, intent(inout) :: err !< Output: error flag
      real (kind=RKIND), dimension(:,:,:), allocatable, intent(inout) :: oceanField

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: scratchPool, geometryPool, meshPool, extrapOceanDataPool
      real (kind=RKIND) :: layerTop, fieldSum, areaSum, weightCellLocal
      real (kind=RKIND), pointer :: weightCell
      integer, dimension(:,:), allocatable :: validOceanMaskOld
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography, areaCell
      integer, pointer :: nCells, nCellsSolve, nISMIP6OceanLayers, nCellsExtra
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:), pointer :: indexToCellID
      integer, dimension(:,:), pointer :: cellsOnCell
      integer :: iCell, jCell, iLayer, iNeighbor, iter, iField
      integer :: localLoopCount
      integer :: nValidNeighb, newValidCount, newMaskCountLocalAccum, newMaskCountGlobal
      integer :: newMaskCountTotal
      logical, pointer :: config_recalculate_thermal_forcing
      real (kind=RKIND), dimension(:,:,:), allocatable :: oceanFieldOld

      err = 0

      ! initialize the ocean data and mask fields
      block => domain % blocklist
      call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrap_ncells_extra', nCellsExtra)
      call mpas_pool_get_config(liConfigs,'config_weight_value_cell', weightCell)
      call mpas_pool_get_config(liConfigs, 'config_recalculate_thermal_forcing', config_recalculate_thermal_forcing)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'extrapOceanData', extrapOceanDataPool)
      call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)

      ! perform horizontal extrapolation until the validOceanMask is unchanged
      allocate(validOceanMaskOld(nISMIP6OceanLayers,nCells+1))
      allocate(oceanFieldOld(nISMIP6OceanLayers,nCells+1,2))
      validOceanMaskOld(:,:) = validOceanMask(:,:)
      oceanFieldOld(:,:,:) = oceanField(:,:,:)

      ! initialize the local loop and count for validOceanMask
      localLoopCount = 0
      newMaskCountTotal = 0
      newMaskCountGlobal = 1
      call mpas_log_write('Weight given to the cell with valid data from extrapolation: $r', realArgs=(/weightCell/))
      do while ( newMaskCountGlobal > 0 )
         localLoopCount = localLoopCount + 1
         newMaskCountLocalAccum = 0
         do iCell = 1, nCellsSolve
            do iLayer = 1, nISMIP6OceanLayers
               do iField = 1, 2  
                  if ( (availOceanMask(iLayer,iCell) == 1) .and. (validOceanMaskOrig(iLayer,iCell) == 0) ) then
                     fieldSum = 0.0
                     areaSum = 0.0
                     newValidCount = 0
                     nValidNeighb = 0
                     do iNeighbor = 1, nEdgesOnCell(iCell)
                        jCell = cellsOnCell(iNeighbor, iCell)
                        if ( validOceanMaskOld(iLayer,jCell) == 1 ) then
                           if ( oceanFieldOld(iLayer,jCell,iField) > 1.0e6_RKIND) then
                              ! raise error if an invalid ocean data value is used
                              call mpas_log_write("ocean data value used for extrapolation is invalid", &
                                   MPAS_LOG_ERR)
                              err = ior(err,1)
                           else
                              fieldSum = fieldSum + (oceanFieldOld(iLayer,jCell,iField) * areaCell(jCell))
                           endif
                            areaSum = areaSum + areaCell(jCell)
                            nValidNeighb = nValidNeighb + 1
                        endif
                     enddo
                     if ( validOceanMaskOld(iLayer,iCell) == 0 .and. nValidNeighb > 0 ) then
                        ! if current cell is not valid, set its weight to zero
                        weightCellLocal = 0.0_RKIND
                        validOceanMask(iLayer,iCell) = 1
                        newValidCount = 1
                     else
                        weightCellLocal = weightCell
                     endif
                     ! perform area-weighted averaging of ocean fields
                     if ( nValidNeighb == 0 ) then
                        oceanField(iLayer,iCell,iField) = oceanFieldOld(iLayer,iCell,iField)
                     else
                        oceanField(iLayer,iCell,iField) = ( weightCellLocal * oceanFieldOld(iLayer,iCell,iField) * areaCell(iCell) + &
                                                 & ((1.0_RKIND - weightCellLocal) * (fieldSum / nValidNeighb)) ) / &
                                                 & ( weightCellLocal * areaCell(iCell) + &
                                                 & (1.0_RKIND - weightCellLocal) * (areaSum / nValidNeighb) )
                     endif
                     ! Accumulate cells added locally until we do the next global reduce
                     newMaskCountLocalAccum = newMaskCountLocalAccum + newValidCount
                  endif
                  if (.not. config_recalculate_thermal_forcing) then
                     ! Don't repeat loop if only extrapolating thermal forcing
                     exit 
                  endif
               enddo
            enddo
         enddo

         ! update halo for validOceanMask
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'validOceanMask')
         call mpas_timer_stop("halo updates")

         validOceanMaskOld(:,:) = validOceanMask(:,:)
         oceanFieldOld(:,:,:) = oceanField(:,:,:)

         ! update count of cells added to mask globally
         call mpas_dmpar_sum_int(domain % dminfo, newMaskCountLocalAccum, newMaskCountGlobal)
         newMaskCountTotal = newMaskCountTotal + newMaskCountGlobal
         !call mpas_log_write('Horizontal extrap: Added total $i new cells to validOceanMask', intArgs=(/newMaskCountGlobal/))
      enddo
      call mpas_log_write('Horizontal extrapolation done after $i iterations.  Added total of $i cells across all processors', &
              intArgs=(/localLoopCount, newMaskCountTotal/))
      deallocate(validOceanMaskOld)
      deallocate(oceanFieldOld)

   end subroutine horizontal_extrapolation
!-----------------------------------------------------------------------


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine vertical_extrapolation
!
!> \brief Extrapolate validOceanMask vertically
!> \author Holly Kyeore Han, modified by Alex Hager
!> \date   November 2023 (modified September 2024)
!> \details
!> This routine extrapolates the horizontally extrapolated
!> validOceanMask through the vertical layers of the ocean.
!> The vertical extrapolation is completed once local new mask count
!> stops updating. The output of the routine is an updated
!> validOceanMask, ocean forcing fields (temperature and salinity, or 
!> thermal forcing) and newMaskCountLocal.

!-----------------------------------------------------------------------

   subroutine vertical_extrapolation(domain, availOceanMask, validOceanMask, newMaskCountGlobal, oceanField, err)

      use li_constants, only: oceanFreezingTempDepthDependence

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
      real (kind=RKIND), dimension(:,:,:), allocatable, intent(inout) :: oceanField

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------
      integer, intent(out) :: newMaskCountGlobal

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: scratchPool, geometryPool, meshPool, extrapOceanDataPool
      real (kind=RKIND) :: layerTop, fieldSum, areaSum
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_zOcean
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography, areaCell
      integer, pointer :: nCells, nCellsSolve, nISMIP6OceanLayers
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:), pointer :: indexToCellID
      integer, dimension(:,:), pointer :: cellsOnCell
      integer :: iCell, jCell, iLayer, iNeighbor, iter, iField
      integer :: localLoopCount, newMaskCountLocalAccum
      logical, pointer :: config_recalculate_thermal_forcing

      err = 0

      ! initialize the ocean data and mask fields
      block => domain % blocklist
      call mpas_pool_get_config(liConfigs, 'config_recalculate_thermal_forcing', config_recalculate_thermal_forcing)
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'extrapOceanData', extrapOceanDataPool)
      call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
      call mpas_pool_get_dimension(meshPool, 'nCellsSolve', nCellsSolve)
      call mpas_pool_get_array(meshPool, 'cellsOnCell', cellsOnCell)
      call mpas_pool_get_array(meshPool, 'areaCell', areaCell)
      call mpas_pool_get_array(meshPool, 'nEdgesOnCell', nEdgesOnCell)
      call mpas_pool_get_array(meshPool, 'indexToCellID', indexToCellID)
      call mpas_pool_get_array(geometryPool, 'cellMask', cellMask)
      call mpas_pool_get_array(geometryPool, 'thickness', thickness)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zOcean', ismip6shelfMelt_zOcean)

      ! initialize the local loop and count for validOceanMask
      newMaskCountGlobal = 0
      newMaskCountLocalAccum = 0
      do iCell = 1, nCellsSolve
         do iLayer = 2, nISMIP6OceanLayers
            do iField = 1, 2
               if ( (availOceanMask(iLayer,iCell) == 1) .and. (validOceanMask(iLayer,iCell) == 0) ) then
                  if ( validOceanMask(iLayer-1,iCell) == 1 ) then
                     if (oceanField(iLayer-1,iCell,iField) > 1.0e6_RKIND) then
                        ! raise error if an invalid ocean data value is used
                        call mpas_log_write("ocean data value used for extrapolation is invalid", &
                               MPAS_LOG_ERR)
                        err = ior(err,1)
                     else
                        if (config_recalculate_thermal_forcing) then
                           !if extrapolating temp/salt, assign value equal to overlying cell
                           oceanField(iLayer,iCell,iField) = oceanField(iLayer-1,iCell,iField)
                        else
                           !if extrapolating thermal forcing, need to account for depth-dependent difference in freezing temp
                           oceanField(iLayer,iCell,iField) = oceanField(iLayer-1,iCell,iField) - &
                              (ismip6shelfMelt_zOcean(iLayer-1) - ismip6shelfMelt_zOcean(iLayer)) * &
                               oceanFreezingTempDepthDependence
                        endif
                     endif
                     validOceanMask(iLayer,iCell) = 1
                     newMaskCountLocalAccum = newMaskCountLocalAccum + 1
                  endif
               endif
               
               if (.not. config_recalculate_thermal_forcing) then
                  exit
               endif
            
            enddo
         enddo
      enddo

      ! update halo for validOceanMask and ocean data
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'validOceanMask')
      call mpas_timer_stop("halo updates")
      ! update count of cells added to mask globally
      call mpas_dmpar_sum_int(domain % dminfo, newMaskCountLocalAccum, newMaskCountGlobal)
      call mpas_log_write('Vertical extrap: Added total $i new cells to the validOceanMask', intArgs=(/newMaskCountGlobal/))

   end subroutine vertical_extrapolation
!-----------------------------------------------------------------------


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine recalculate_thermal_forcing
!
!> \brief recalculate thermal forcing from extrapolated temp/salt
!> \author Alex Hager
!> \date   September 2024
!> \details
!> Use extrapolated temperature and salt to recalculate thermal forcing 
!> throughout the extrapolated region. User specifies which thermal forcing 
!> parameterization from Hager et al. (2023) to use with
!> config_TF_param_type. Output thermal forcing fields are 2D.
!-----------------------------------------------------------------------

   subroutine recalculate_thermal_forcing(domain, err)

      !-----------------------------------------------------------------
      ! input variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! input/output variables
      !-----------------------------------------------------------------
      type (domain_type), intent(inout) :: domain  !< Input/Output: domain object
      integer, intent(inout) :: err !< Output: error flag

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: geometryPool, meshPool, extrapOceanDataPool
      character (len=StrKIND), pointer :: config_TF_param_type
      real (kind=RKIND), pointer :: config_TF_iceberg_depth
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing
      real (kind=RKIND), dimension(:), pointer :: ismip6_2dThermalForcing
      real (kind=RKIND), dimension(:,:), pointer :: MPAS_3dOceanTemperature
      real (kind=RKIND), dimension(:,:), pointer :: MPAS_3dOceanSalinity
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_zOcean
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_zBndsOcean
      real (kind=RKIND), dimension(:), pointer :: icebergFjordMask
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), pointer :: dzAccumulated, dz, total, gadeSlope
      integer, pointer :: nISMIP6OceanLayers, nCells, iCell, iLayer, counter
      integer, dimension(:), allocatable :: icebergCellVertID
      integer, dimension(:), allocatable :: oceanCellVertID
      real (kind=RKIND), dimension(:), pointer :: freezingTemperature
      allocate(icebergCellVertID(nISMIP6OceanLayers))
      allocate(oceanCellVertID(nISMIP6OceanLayers))
      allocate(freezingTemperature(nISMIP6OceanLayers))

      block => domain % blocklist
      do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_config(liConfigs, 'config_TF_param_type', config_TF_param_type)
      call mpas_pool_get_config(liConfigs, 'config_TF_iceberg_depth', config_TF_iceberg_depth)
      call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_3dThermalForcing', ismip6shelfMelt_3dThermalForcing)
      call mpas_pool_get_array(geometryPool, 'ismip6_2dThermalForcing', ismip6_2dThermalForcing)
      call mpas_pool_get_array(geometryPool, 'MPAS_3dOceanTemperature', MPAS_3dOceanTemperature)
      call mpas_pool_get_array(geometryPool, 'MPAS_3dOceanSalinity', MPAS_3dOceanSalinity)
      call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zOcean', ismip6shelfMelt_zOcean)
      call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zBndsOcean', ismip6shelfMelt_zBndsOcean)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(geometryPool, 'icebergFjordMask', icebergFjordMask)
      call mpas_pool_get_array(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
      call mpas_pool_get_array(meshPool, 'nCells', nCells)

      do iCell = 1, nCells
        !Main difference between these parameterization is how we collapse into 2d, which happens in mpas_li_ice_shelt_melt.F
        !Difference between '...melt' and '...retreat' parameterizations happens in main body of mpas_li_ocean_extrap.F
        if ( (config_TF_param_type == 'AMretreat') .or. (config_TF_param_type == 'ISMIP6retreat') &
                .or. (config_TF_param_type == 'AMmelt') .or. (config_TF_param_type == 'ISMIP6melt') &
                .or. ((config_TF_param_type == 'AMfit') .and. (icebergFjordMask(iCell)) == 1) ) then

            ! From Jenkins 2011 and Slater 2020
            ismip6shelfMelt_3dThermalForcing(:,iCell) = MPAS_3dOceanTemperature(:,iCell) - (gamma1 * MPAS_3dOceanSalinity(:,iCell) + gamma2 + gamma3 * ismip6shelfMelt_zOcean)

        elseif ( (config_TF_param_type == 'AMberg') .or. ( (config_TF_param_type == 'AMfit') .and. (icebergFjordMask(iCell) == 1)) ) then
            ! Alter 3d profile following the iceberg Parameterization from Hager 2023 for fjord where icebergs are prevalent
            ! find cell just above config_TF_iceberg_depth
            icebergCellVertID(:) = 0
            do iLayer = 1, nISMIP6OceanLayers
               if ( (ismip6shelfMelt_zOcean(iLayer) <= config_TF_iceberg_depth) & 
                       .and. (bedTopography(iCell) >= ismip6shelfMelt_zOcean(iLayer)) )then
                  icebergCellVertID(iLayer) = iLayer
               endif
            enddo

            ! calculate freezing temperature. 
            ! << NOTE >>: Again we are using the simplified equation from Jenkins 2011 to do this. May eventually want to 
            ! switch to method in gsw toolbox to do this to be consistent with Hager 2023
            freezingTemperature = (gamma1 * MPAS_3dOceanSalinity(:,iCell) + gamma2 + gamma3 * ismip6shelfMelt_zOcean)
         
            ! Calculate Gade slope (slope of submarine melt mixing line) using cell just above iceberg depth as reference point. Assumes ice at pressure melting temperature
            gadeSlope = (1/MPAS_3dOceanSalinity(maxval(icebergCellVertID),iCell)) * (L/cp - (freezingTemperature(maxval(icebergCellVertID)) & 
                    - MPAS_3dOceanTemperature(maxval(icebergCellVertID),iCell)))

            ! Now that we have the slope of the melt line, alter ocean temperatures along submarine melt line
            do iLayer = 1, nISMIP6OceanLayers
               if (icebergCellVertID(iLayer) /= 0) then
                  MPAS_3dOceanTemperature(iLayer, iCell) = gadeSlope * (MPAS_3dOceanSalinity(iLayer, iCell) - MPAS_3dOceanSalinity(maxval(icebergCellVertID),iCell)) &
                          + MPAS_3dOceanTemperature(iLayer, iCell)
                  
                  ! Ensure no artificially adjusted temperatures fall below the freezing temperature
                  if (MPAS_3dOceanTemperature(iLayer, iCell) < freezingTemperature(iLayer)) then
                     MPAS_3dOceanTemperature(iLayer, iCell) = freezingTemperature(iLayer)
                  endif

                  ismip6shelfMelt_3dThermalForcing(iLayer, iCell) = MPAS_3dOceanTemperature(iLayer, iCell) - freezingTemperature(iLayer)
               endif
            enddo
        endif
      enddo
      
      deallocate(icebergCellVertID)
      deallocate(freezingTemperature)

      block => block % next
      enddo
   end subroutine recalculate_thermal_forcing
end module li_ocean_extrap


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
