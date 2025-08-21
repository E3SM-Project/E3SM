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
   use li_constants

   implicit none
   private

   !--------------------------------------------------------------------
   ! Public parameters
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   ! Public member functions
   !--------------------------------------------------------------------

   public :: li_ocean_extrap_solve
   public :: li_is_tf_invalid
   
   !--------------------------------------------------------------------
   ! Private module variables
   !--------------------------------------------------------------------
   
   ! Coefficients in Jenkins 2011 shelf melt parameterizations
   real (kind=RKIND), parameter :: gamma1 = -5.7E-2
   real (kind=RKIND), parameter :: gamma2 = 8.8E-2
   real (kind=RKIND), parameter :: gamma3 = 7.61E-4

   real (kind=RKIND), parameter :: default_salinity = 32.0
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
      character(len=StrKIND), pointer :: config_TF_vertical_calc
      logical, pointer :: config_calculate_thermal_forcing
      logical, pointer :: config_extrap_sill_blocking
      real(kind=RKIND), pointer :: config_sea_level
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: scratchPool, geometryPool, meshPool, extrapOceanDataPool
      real (kind=RKIND) :: layerCenter
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing, ismip6shelfMelt_zBndsOcean
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_zOcean
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography
      integer, dimension(:,:), pointer :: orig3dOceanMask
      integer, dimension(:,:), pointer :: validOceanMask, validOceanMaskOrig, availOceanMask !masks to pass to flood-fill routine
      integer, dimension(:), pointer :: seedOceanMaskHoriz, growOceanMaskHoriz
      real (kind=RKIND), dimension(:,:), pointer :: oceanSalinity, oceanTemperature
      integer, dimension(:), allocatable :: seedOceanMaskHorizOld
      integer, pointer :: nCells, nCellsSolve, nISMIP6OceanLayers, nCellsExtra
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:,:), pointer :: cellsOnCell
      integer :: iCell, jCell, iLayer, iNeighbor, iter, err_tmp
      integer :: GlobalLoopCount, newMaskCountGlobal, oceanOnLand
      type (field1dInteger), pointer :: seedMaskField
      type (field1dInteger), pointer :: growMaskField
      integer, dimension(:), pointer :: seedMask, growMask, connectedMarineMask, altConnectedMarineMask
      real (kind=RKIND), dimension(:,:,:), pointer :: oceanField
      integer, dimension(:), allocatable :: alteredCells
      ! No init is needed.
      err = 0
      err_tmp = 0

      call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrapolation', config_ocean_data_extrapolation)

      if ( config_ocean_data_extrapolation ) then

         ! call the extrapolation scheme
         call mpas_log_write('Beginning ocean data extrapolation')
         block => domain % blocklist

        ! initialize the ocean data and mask fields
         call mpas_pool_get_config(liConfigs, 'config_sea_level', config_sea_level)
         call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrap_ncells_extra', nCellsExtra)
         call mpas_pool_get_config(liConfigs, 'config_calculate_thermal_forcing', config_calculate_thermal_forcing)
         call mpas_pool_get_config(liConfigs, 'config_TF_vertical_calc', config_TF_vertical_calc) 
         call mpas_pool_get_config(liConfigs, 'config_extrap_sill_blocking', config_extrap_sill_blocking)
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
         call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zOcean', ismip6shelfMelt_zOcean)
         call mpas_pool_get_array(extrapOceanDataPool, 'orig3dOceanMask', orig3dOceanMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'validOceanMask', validOceanMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'availOceanMask', availOceanMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'seedOceanMaskHoriz', seedOceanMaskHoriz)
         call mpas_pool_get_array(extrapOceanDataPool, 'growOceanMaskHoriz', growOceanMaskHoriz)
         call mpas_pool_get_array(extrapOceanDataPool, 'oceanField', oceanField)
         call mpas_pool_get_array(extrapOceanDataPool, 'connectedMarineMask', connectedMarineMask)
         call mpas_pool_get_array(extrapOceanDataPool, 'altConnectedMarineMask', altConnectedMarineMask)

         if (config_calculate_thermal_forcing) then
            call mpas_pool_get_array(extrapOceanDataPool, 'oceanTemperature', oceanTemperature)
            call mpas_pool_get_array(extrapOceanDataPool, 'oceanSalinity', oceanSalinity)
         endif

         ! check parameterization compatibility. TF_vertical_calc = '200-500'
         ! assumes no knowledge of fjord geometry, so is not compatible with
         ! sill blocking 
         if ( (config_TF_vertical_calc == '200-500') .and. (config_extrap_sill_blocking) ) then
            call mpas_log_write("config_TF_vertical_calc set to '200-500' is not compatible with config_extrap_sill_blocking", MPAS_LOG_ERR)         
            err = ior(err, 1)
         endif

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
            ! Update halos
            ! Note: halo updates could be skipped for a number of iterations equal to the
            ! number of halo rows as a possible performance improvement in the future.
            ! (And the loop above would need to be changed to from nCellsSolve to nCells)
            call mpas_timer_start("halo updates")
            call mpas_dmpar_field_halo_exch(domain, 'growOceanMaskHoriz')
            call mpas_dmpar_field_halo_exch(domain, 'seedOceanMaskHoriz')
            call mpas_timer_stop("halo updates")
            seedOceanMaskHorizOld(:) = seedOceanMaskHoriz(:)
         enddo
         deallocate(seedOceanMaskHorizOld)

         ! Calculate mask of bedTopography below sea level that is connected to
         ! the ocean
         call mpas_pool_get_subpool(domain % blocklist % structs, 'scratch', scratchPool)
         call mpas_pool_get_field(scratchPool, 'seedMask', seedMaskField)
         call mpas_allocate_scratch_field(seedMaskField, single_block_in = .true.)
         seedMask => seedMaskField % array
         seedMask(:) = 0 
         
         call mpas_pool_get_field(scratchPool, 'growMask', growMaskField)
         call mpas_allocate_scratch_field(growMaskField, single_block_in = .true.)
         growMask => growMaskField % array
         growMask(:) = 0
         
         do iCell = 1, nCells
            ! seedMask = open ocean cells in contact with the domain boundary
            if ((bedTopography(iCell) < config_sea_level) .and. (thickness(iCell) == 0.0_RKIND)) then
               do iNeighbor = 1, nEdgesOnCell(iCell)
                  if (cellsOnCell(iNeighbor, iCell) == nCells + 1) then
                     seedMask(iCell) = 1
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
         call li_flood_fill(seedMask, growMask, domain)
         
         connectedMarineMask(:) = seedMask(:)
         
         ! halo update
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'connectedMarineMask')
         call mpas_timer_stop("halo updates")

         ! Calculate mask of ice-free cells with bedTopography below sea level
         ! and that are connected the ocean (again through only ice-free cells).
         ! Used in combination with connectedMarineMask to remove designation of 
         ! inland seas as connected ocean cells
         seedMask(:) = 0
         growMask(:) = 0

         do iCell = 1, nCells
            ! seedMask = open ocean cells in contact with the domain boundary
            if ((bedTopography(iCell) < config_sea_level) .and. (thickness(iCell) == 0.0_RKIND)) then
               do iNeighbor = 1, nEdgesOnCell(iCell)
                  if (cellsOnCell(iNeighbor, iCell) == nCells + 1) then
                     seedMask(iCell) = 1
                     exit  ! no need to keep checking neighbors
                  endif
               enddo
            endif
            ! altGrowMask - all marine cells without ice
            if ( (bedTopography(iCell) < config_sea_level) .and. (thickness(iCell) == 0.0_RKIND) ) then
               growMask(iCell) = 1
            endif
         enddo
        
         ! now create mask of all marine locations connected to open ocean
         ! through ice-free cells - to be used below to screen out lakes
         call li_flood_fill(seedMask, growMask, domain)
         
         altConnectedMarineMask(:) = seedMask(:)
         
         ! halo update
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'altConnectedMarineMask')
         call mpas_timer_stop("halo updates")
         
         ! alter orig3dOceanMask if necessary to avoid errors
         oceanOnLand = 0
         do iCell = 1, nCells
            ! Ensure MPAS ocean cells have not been interpolated onto non-ocean cells on the MALI mesh
            do iLayer = 1, nISMIP6OceanLayers
                if ( (bedTopography(iCell) > 0.0_RKIND) .and. (orig3dOceanMask(iLayer,iCell) == 1) ) then
                   orig3dOceanMask(iLayer,iCell) = 0
                   ismip6shelfMelt_3dThermalForcing(iLayer,iCell) = 0.0_RKIND
                   connectedMarineMask(iCell) = 0
                 
                   if (config_calculate_thermal_forcing) then
                      oceanTemperature(iLayer,iCell) = 0.0_RKIND
                      oceanSalinity(iLayer,iCell) = 0.0_RKIND
                   endif
                   oceanOnLand = 1
               endif
            enddo
         enddo
         
         if (oceanOnLand == 1) then
            call mpas_log_write("Non-ocean cells removed from orig3dOceanMask", MPAS_LOG_WARN)
         endif

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
            if ((config_calculate_thermal_forcing) .and. (.not. config_extrap_sill_blocking) &
                    .and. (bedTopography(iCell) < 0.0_RKIND)) then
                        
                  ! If recalculating thermal forcing without sill
                  ! blocking, ignore bathymetry and extrapolate water
                  ! properties to all cells below sea level. 
                  availOceanMask(:,iCell) = 1
            else
               do iLayer = 1, nISMIP6OceanLayers
                  layerCenter = ismip6shelfMelt_zOcean(iLayer)
                  
                  if ( ( (seedOceanMaskHoriz(iCell) == 1) .and. (connectedMarineMask(iCell) == 1) ) &
                          .and. ( (bedTopography(iCell) < layerCenter) &
                          .or. ((iLayer == 1) .and. (bedTopography(iCell) < 0.0_RKIND)) ) ) then
                     
                     ! Only assign surface cells and other ocean cells with cell centers above bed topography
                     ! to availOceanMask.  
                     availOceanMask(iLayer,iCell) = 1
                  endif  
               enddo
            endif
         enddo

         ! Make a copy of original mask to use for extending the mask during extrapolation
         validOceanMask(:,:) = orig3dOceanMask(:,:)

         ! initialize the extrapolated fields
         if (config_calculate_thermal_forcing) then
            
            !combine fields into one variable to expedite extrapolation
            oceanField(1,:,:) = oceanTemperature(:,:) * validOceanMask(:,:)
            oceanField(2,:,:) = oceanSalinity(:,:) * validOceanMask(:,:)
         else

            oceanField(1,:,:) = ismip6shelfMelt_3dThermalForcing * validOceanMask(:,:)
            oceanField(2,:,:) = invalid_value_TF
         endif
         
         ! initialize the invalid data locations with fill value
         do iCell = 1, nCells
            do iLayer = 1, nISMIP6OceanLayers
               if ((connectedMarineMask(iCell) == 0) .and. (bedTopography(iCell) < config_sea_level)) then
                  ! Don't assign invalid value to lakes/inland seas disconnected from global ocean
                  ! Let them retain the existing value:
                  ! This will take on the valid ocean data value where it exists or
                  ! zero where valid ocean data does not exist
               elseif (validOceanMask(iLayer,iCell) == 0) then
                  ! everywhere else where valid ocean data does not exist, insert invalid value outside of validOceanMask
                  oceanField(:,iLayer,iCell) = invalid_value_TF
               endif
            enddo
         enddo

         ! halo update
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'oceanField')
         call mpas_timer_stop("halo updates")

         ! deallocate scratch fields used for flood fill
         call mpas_deallocate_scratch_field(seedMaskField, single_block_in=.true.)

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
            ! if using the 200-500 TF parameterization, perform an initial vertical
            ! extrapolation to bottom z layer before extrapolating horizontally.
            if (config_TF_vertical_calc == '200-500') then
               call mpas_timer_start("vertical scheme")
               call vertical_extrapolation(domain, availOceanMask, validOceanMask, newMaskCountGlobal, oceanField, err_tmp)
               err = ior(err, err_tmp)
               call mpas_timer_stop("vertical scheme")
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
         if (config_calculate_thermal_forcing) then
            oceanTemperature(:,:) = oceanField(1,:,:)
            oceanSalinity(:,:) = oceanField(2,:,:)

            call calculate_thermal_forcing(domain, err)

         else 
            ismip6shelfMelt_3dThermalForcing(:,:) = oceanField(1,:,:)
         endif
            
         ! As a last step, identify all inland seas and lakes and set TF to zero in these
         ! locations. Otherwise, these cells are given an invalid_TF value that
         ! could result in runaway melting. Need a scratch Registry variable to
         ! iterate through and perform halo updates on, so borrow growMask,
         ! which is no longer used at this point in the routine

         allocate(alteredCells(nCells+1))
         alteredCells(:) = 0
         growMask(:) = 0

         do iCell = 1, nCells
            if ( (altConnectedMarineMask(iCell) /= connectedMarineMask(iCell)) .and. (thickness(iCell) == 0.0_RKIND) ) then
               do iLayer = 1, nISMIP6OceanLayers 
                  if ((iLayer == 1) .or. (ismip6shelfMelt_zOcean(iLayer) >= bedTopography(iCell))) then
                  
                     if (config_calculate_thermal_forcing) then
                        ! for ice-dammed inland seas, assume default salinity, and back-calculate sea temperature to ensure zero thermal forcing
                        oceanSalinity(iLayer,iCell) = default_salinity
                        oceanTemperature(iLayer,iCell) = (gamma1 * oceanSalinity(iLayer,iCell) + gamma2 + gamma3 * ismip6shelfMelt_zOcean(iLayer))
                     endif
                     ismip6shelfMelt_3dThermalForcing(iLayer,iCell) = 0.0_RKIND
                  endif
               enddo
               growMask(iCell) = 1
            endif
         enddo

         call mpas_timer_start("halo_updates")
         call mpas_dmpar_field_halo_exch(domain,'growMask')
         call mpas_timer_stop("halo_updates")

         do iter = 1, nCellsExtra
            ! Iterate through neighboring cells and alter neighbors within
            ! config_ocean_data_extrap_ncells_extra that have ice and are
            ! below sea level. GL melt is calculated from TF in first cell inside
            ! GL, and we need a buffer zone in case the glacier retreats to invalid
            ! TF values in the interior.

            alteredCells = growMask
            growMask(:) = 0

            do iCell = 1, nCells
               if (alteredCells(iCell) == 1) then
                  do iNeighbor = 1, nEdgesOnCell(iCell)
                     jCell = cellsOnCell(iNeighbor,iCell)
                     do iLayer = 1, nISMIP6OceanLayers
                        if ( ((thickness(jCell) > 0) .and. (bedTopography(jCell) < 0)) &
                              .and. ( (ismip6shelfMelt_zOcean(iLayer) >= bedTopography(jCell)) &
                              .or. (iLayer == 1))) then 
                           if (config_calculate_thermal_forcing) then
                                   ! Change temperature and/or salinity if
                                   ! outside of physical range
                              
                              if (li_is_tf_invalid(oceanSalinity(iLayer,jCell))) then
                                 oceanSalinity(iLayer,jCell) = default_salinity
                                 growMask(jCell) = 1
                              endif
                              if (li_is_tf_invalid(oceanTemperature(iLayer,jCell))) then
                                 oceanTemperature(iLayer,jCell) = (gamma1 * oceanSalinity(iLayer,jCell) + gamma2 + gamma3 * ismip6shelfMelt_zOcean(iLayer))
                                 growMask(jCell) = 1
                              endif
                           endif
                           if (li_is_tf_invalid(ismip6shelfMelt_3dThermalForcing(iLayer,jCell))) then
                              ismip6shelfMelt_3dThermalForcing(iLayer,jCell) = 0.0_RKIND
                              growMask(jCell) = 1
                           endif
                        endif
                     enddo
                  enddo
               endif
            enddo
            call mpas_timer_start("halo_updates")
            call mpas_dmpar_field_halo_exch(domain,'oceanSalinity')
            call mpas_dmpar_field_halo_exch(domain,'oceanTemperature')
            call mpas_dmpar_field_halo_exch(domain,'ismip6shelfMelt_3dThermalForcing')
            call mpas_dmpar_field_halo_exch(domain,'growMask')
            call mpas_timer_stop("halo_updates")
         enddo
         
         deallocate(alteredCells)
         call mpas_deallocate_scratch_field(growMaskField, single_block_in=.true.)
         
         call mpas_log_write('Ocean data extrapolation complete')
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
      real (kind=RKIND), dimension(:,:,:), pointer, intent(inout) :: oceanField

      !-----------------------------------------------------------------
      ! output variables
      !-----------------------------------------------------------------

      !-----------------------------------------------------------------
      ! local variables
      !-----------------------------------------------------------------
      type (block_type), pointer :: block
      type (mpas_pool_type), pointer :: scratchPool, geometryPool, meshPool, extrapOceanDataPool
      real (kind=RKIND) :: layerTop, areaSum, weightCellLocal
      real (kind=RKIND), pointer :: weightCell
      integer, dimension(:,:), allocatable :: validOceanMaskOld
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing
      real (kind=RKIND), dimension(:), pointer :: thickness, bedTopography, areaCell
      integer, pointer :: nCells, nCellsSolve, nISMIP6OceanLayers, nCellsExtra
      integer, dimension(:), pointer :: cellMask, nEdgesOnCell
      integer, dimension(:), pointer :: indexToCellID
      integer, dimension(:,:), pointer :: cellsOnCell
      integer :: iCell, jCell, iLayer, iNeighbor, iter, iField
      integer :: localLoopCount, nFields, newMaskCountTotal
      integer :: nValidNeighb, newValidCount, newMaskCountLocalAccum, newMaskCountGlobal
      logical, pointer :: config_calculate_thermal_forcing
      real (kind=RKIND), dimension(:,:,:), allocatable :: oceanFieldOld
      real (kind=RKIND), dimension(:), allocatable :: fieldSum

      err = 0

      ! initialize the ocean data and mask fields
      block => domain % blocklist
      call mpas_pool_get_config(liConfigs, 'config_ocean_data_extrap_ncells_extra', nCellsExtra)
      call mpas_pool_get_config(liConfigs,'config_weight_value_cell', weightCell)
      call mpas_pool_get_config(liConfigs, 'config_calculate_thermal_forcing', config_calculate_thermal_forcing)
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
      allocate(oceanFieldOld(2,nISMIP6OceanLayers,nCells+1))
      validOceanMaskOld(:,:) = validOceanMask(:,:)
      oceanFieldOld(:,:,:) = oceanField(:,:,:)

      if (.not. config_calculate_thermal_forcing) then
         nFields = 1
      else
         nFields = 2
      endif

      allocate(fieldSum(nFields))

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
               if ( (availOceanMask(iLayer,iCell) == 1) .and. (validOceanMaskOrig(iLayer,iCell) == 0) ) then
                  fieldSum(:) = 0.0
                  areaSum = 0.0
                  nValidNeighb = 0
                  newValidCount = 0
                  do iNeighbor = 1, nEdgesOnCell(iCell)
                     jCell = cellsOnCell(iNeighbor, iCell)
                     if ( validOceanMaskOld(iLayer,jCell) == 1 ) then
                         
                        do iField = 1, nFields  
                           if (li_is_tf_invalid(oceanFieldOld(iField,iLayer,jCell))) then
                              ! raise error if an invalid ocean data value is used
                              call mpas_log_write("ocean data value used for extrapolation is invalid" // &
                                   "in horizontal extrapolations: cell id=$i, iLayer=$i, iField=$i, value=$r", &
                                   MPAS_LOG_ERR, intArgs=(/indexToCellID(iCell), iLayer-1, iField/), realArgs=(/oceanField(iField,iLayer-1,iCell)/))
                              err = ior(err,1)
                           else
                              fieldSum(iField) = fieldSum(iField) + (oceanFieldOld(iField,iLayer,jCell) * areaCell(jCell))
                           endif
                        enddo

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

                  do iField = 1, nFields
                     ! perform area-weighted averaging of ocean fields
                     if ( nValidNeighb == 0 ) then
                        oceanField(iField,iLayer,iCell) = oceanFieldOld(iField,iLayer,iCell)
                     else
                        oceanField(iField,iLayer,iCell) = ( weightCellLocal * oceanFieldOld(iField,iLayer,iCell) * areaCell(iCell) + &
                           & ((1.0_RKIND - weightCellLocal) * (fieldSum(iField) / nValidNeighb)) ) / &
                                            & ( weightCellLocal * areaCell(iCell) + &
                                            & (1.0_RKIND - weightCellLocal) * (areaSum / nValidNeighb) )
                     endif
                  enddo

                  ! Accumulate cells added locally until we do the next global reduce
                  newMaskCountLocalAccum = newMaskCountLocalAccum + newValidCount
               endif
            enddo
         enddo

         ! update halo for validOceanMask
         call mpas_timer_start("halo updates")
         call mpas_dmpar_field_halo_exch(domain, 'validOceanMask')
         call mpas_dmpar_field_halo_exch(domain,'oceanField')
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
      deallocate(fieldSum)

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
      real (kind=RKIND), dimension(:,:,:), pointer, intent(inout) :: oceanField

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
      integer :: nFields, localLoopCount, newMaskCountLocalAccum
      logical, pointer :: config_calculate_thermal_forcing

      err = 0

      ! initialize the ocean data and mask fields
      block => domain % blocklist
      call mpas_pool_get_config(liConfigs, 'config_calculate_thermal_forcing', config_calculate_thermal_forcing)
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

      if (.not. config_calculate_thermal_forcing) then
         nFields = 1
      else
         nFields = 2
      endif

      ! initialize the local loop and count for validOceanMask
      newMaskCountGlobal = 0
      newMaskCountLocalAccum = 0
      do iCell = 1, nCellsSolve
         do iLayer = 2, nISMIP6OceanLayers
         ! Extrapolate water properties downward and stop 1 layer below bed topography. Don't update validOceanMask if cell is below bedTopography.
         ! This ensures partial cells with cell centers below bed topography are accounted for. 

            if ( (validOceanMask(iLayer,iCell) == 0) .and. (validOceanMask(iLayer-1,iCell) == 1) &
                    .and. (availOceanMask(iLayer-1,iCell) == 1) .and. &
                    (ismip6shelfMelt_zOcean(iLayer) >= bedTopography(iCell)) ) then
                  do iField = 1, nFields
                     if (li_is_tf_invalid(oceanField(iField,iLayer-1,iCell))) then
                        ! raise error if an invalid ocean data value is used
                        call mpas_log_write("ocean data value used for extrapolation is invalid" // &
                           "in vertical extrapolations: cell id=$i, iLayer=$i, TF=$r", &
                           MPAS_LOG_ERR, intArgs=(/indexToCellID(iCell), iLayer-1/), realArgs=(/oceanField(iField,iLayer-1,iCell)/))
                        err = ior(err,1)
                     else
                        if (config_calculate_thermal_forcing) then
                           !if extrapolating temp/salt, assign value equal to overlying cell
                           oceanField(iField,iLayer,iCell) = oceanField(iField,iLayer-1,iCell)
                        else
                           !if extrapolating thermal forcing, need to account for depth-dependent difference in freezing temp
                           oceanField(iField,iLayer,iCell) = oceanField(iField,iLayer-1,iCell) - &
                              (ismip6shelfMelt_zOcean(iLayer-1) - ismip6shelfMelt_zOcean(iLayer)) * &
                              oceanFreezingTempDepthDependence
                           exit
                        endif
                     endif
                  enddo
                  ! update validOceanMask if cell is defined in availOceanMask
                  if (availOceanMask(iLayer,iCell) == 1) then
                     validOceanMask(iLayer,iCell) = 1
                     newMaskCountLocalAccum = newMaskCountLocalAccum + 1
                  endif
            endif
         enddo
      enddo

      ! update halo for validOceanMask and ocean data
      call mpas_timer_start("halo updates")
      call mpas_dmpar_field_halo_exch(domain, 'validOceanMask')
      call mpas_dmpar_field_halo_exch(domain, 'oceanField')
      call mpas_timer_stop("halo updates")
      ! update count of cells added to mask globally
      call mpas_dmpar_sum_int(domain % dminfo, newMaskCountLocalAccum, newMaskCountGlobal)
      call mpas_log_write('Vertical extrap: Added total $i new cells to the validOceanMask', intArgs=(/newMaskCountGlobal/))

   end subroutine vertical_extrapolation
!-----------------------------------------------------------------------


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  !  routine calculate_thermal_forcing
!
!> \brief calculates thermal forcing from temp/salt
!> \author Alex Hager
!> \date   September 2024
!> \details
!> Use temperature and salinity to calculate thermal forcing. 
!> If config_TF_iceberg_melt is true, then ocean temperature 
!> in the upper water column will first be artificially cooled 
!> to imitate iceberg melt-driven cooling.
!-----------------------------------------------------------------------

   subroutine calculate_thermal_forcing(domain, err)

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
      character (len=StrKIND), pointer :: config_TF_iceberg_melt
      real (kind=RKIND), pointer :: config_TF_iceberg_depth
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_3dThermalForcing
      real (kind=RKIND), dimension(:,:), pointer :: oceanTemperature
      real (kind=RKIND), dimension(:,:), pointer :: oceanSalinity
      real (kind=RKIND), dimension(:), pointer :: ismip6shelfMelt_zOcean
      real (kind=RKIND), dimension(:,:), pointer :: ismip6shelfMelt_zBndsOcean
      integer, dimension(:), pointer :: icebergFjordMask
      real (kind=RKIND), dimension(:), pointer :: bedTopography
      real (kind=RKIND), pointer :: dzAccumulated, dz, total
      integer, pointer :: nISMIP6OceanLayers, nCells
      integer :: iCell, iLayer, counter, icebergCellVertID
      real (kind=RKIND), dimension(:), allocatable :: freezingTemperature
      real (kind=RKIND) :: inverseSalt, phaseChange, tempDiff, saltDiff, newTemp, gadeSlope
      real (kind=RKIND), pointer :: invalid_value_TF

      block => domain % blocklist
      do while (associated(block))
      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
      call mpas_pool_get_subpool(block % structs, 'geometry', geometryPool)
      call mpas_pool_get_subpool(block % structs, 'extrapOceanData', extrapOceanDataPool)
      call mpas_pool_get_config(liConfigs, 'config_TF_iceberg_melt', config_TF_iceberg_melt)
      call mpas_pool_get_config(liConfigs, 'config_TF_iceberg_depth', config_TF_iceberg_depth)
      call mpas_pool_get_config(liConfigs, 'config_invalid_value_TF', invalid_value_TF)
      call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_3dThermalForcing', ismip6shelfMelt_3dThermalForcing)
      call mpas_pool_get_array(extrapOceanDataPool, 'oceanTemperature', oceanTemperature)
      call mpas_pool_get_array(extrapOceanDataPool, 'oceanSalinity', oceanSalinity)
      call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zOcean', ismip6shelfMelt_zOcean)
      call mpas_pool_get_array(geometryPool, 'ismip6shelfMelt_zBndsOcean', ismip6shelfMelt_zBndsOcean)
      call mpas_pool_get_array(geometryPool, 'bedTopography', bedTopography)
      call mpas_pool_get_array(extrapOceanDataPool, 'icebergFjordMask', icebergFjordMask)
      call mpas_pool_get_dimension(meshPool, 'nISMIP6OceanLayers', nISMIP6OceanLayers)
      call mpas_pool_get_dimension(meshPool, 'nCells', nCells)

      allocate(freezingTemperature(nISMIP6OceanLayers))
     
      if ( (.not. config_TF_iceberg_melt == 'off') .and. (.not. config_TF_iceberg_melt == 'on') &
              .and. (.not. config_TF_iceberg_melt == 'mask') ) then
         call mpas_log_write("Unknown config_TF_iceberg_melt specified. Options are: 'on', 'off', or 'mask'", MPAS_LOG_ERR)
         err = ior(err, 1)
      endif
      
      ! Main difference between TF parameterization is how we collapse into 2d, which happens in mpas_li_ice_shelt_melt.F
      do iCell = 1, nCells
         if (config_TF_iceberg_melt == 'off') then

            ! From Jenkins 2011 and Slater 2020
            do iLayer = 1, nISMIP6OceanLayers
               if (li_is_tf_invalid(oceanTemperature(iLayer,iCell)) .or. li_is_tf_invalid(oceanSalinity(iLayer,iCell))) then
                  ismip6shelfMelt_3dThermalForcing(iLayer,iCell) = invalid_value_TF
               else 
                  ismip6shelfMelt_3dThermalForcing(iLayer,iCell) = oceanTemperature(iLayer,iCell) - (gamma1 * oceanSalinity(iLayer,iCell) + gamma2 + gamma3 * ismip6shelfMelt_zOcean(iLayer))
               endif
            enddo
      
         elseif ( (config_TF_iceberg_melt == 'on') .or. ( (config_TF_iceberg_melt == 'mask') .and. (icebergFjordMask(iCell) == 1)) ) then
      
            ! Alter 3d profile following the iceberg Parameterization from Hager 2023 for fjord where icebergs are prevalent
            ! find cell just above config_TF_iceberg_depth
            icebergCellVertID = 0
            freezingTemperature(:) = 0.0_RKIND

            do iLayer = 1, nISMIP6OceanLayers
               if ( ( (ismip6shelfMelt_zOcean(iLayer) >= config_TF_iceberg_depth) & 
                       .and. (ismip6shelfMelt_zOcean(iLayer) >= bedTopography(iCell)) ) & 
                       .or. ( (iLayer == 1) .and. (bedTopography(iCell) >= config_TF_iceberg_depth) &
                       .and. (bedTopography(iCell) < 0.0_RKIND) ) ) then
                  icebergCellVertID = iLayer

               endif
               ! calculate freezing temperature. 
               ! << NOTE >>: Again we are using the simplified equation from Jenkins 2011 to do this. May eventually want to 
               ! switch to method in gsw toolbox to do this to be consistent with Hager 2023
               freezingTemperature(iLayer) = (gamma1 * oceanSalinity(iLayer,iCell) + gamma2 + gamma3 * ismip6shelfMelt_zOcean(iLayer))
            enddo

            !Start alteration just above iceberg depth or bedTopography and work upwards
            do iLayer = icebergCellVertID, 1, -1 
               if (.not. (li_is_tf_invalid(oceanTemperature(iLayer,iCell)) .or. li_is_tf_invalid(oceanSalinity(iLayer,iCell)))) then
               
                  ! Calculate Gade slope (slope of submarine melt mixing line) at each z-level. Assumes ice at pressure melting temperature
                  if (oceanSalinity(iLayer,iCell) == 0) then
                     gadeSlope = 0.0_RKIND
                  else
                     inverseSalt = (1.0_RKIND / oceanSalinity(iLayer,iCell)) 
                     phaseChange = latent_heat_ice/cp_ice 
                     tempDiff = (freezingTemperature(iLayer) - oceanTemperature(iLayer,iCell))

                     gadeSlope = inverseSalt * (phaseChange - tempDiff)
                  endif

                  ! Now that we have the slope of the melt line, alter ocean temperatures along submarine melt line.
                  if (ismip6shelfMelt_zOcean(iLayer + 1) < bedTopography(iCell)) then
                     ! Not enough information to alter cells immediately above bed topography
                     saltDiff = 0.0_RKIND
                     newTemp = oceanTemperature(iLayer, iCell)
                  else
                     saltDiff = oceanSalinity(iLayer,iCell) - oceanSalinity(iLayer + 1,iCell)     
                     newTemp = (gadeSlope * saltDiff) + oceanTemperature(iLayer + 1, iCell)
                  endif
               
                  ! Ensure iceberg melting only cools ocean. Uncommon but possible to warm ocean if already a sharp thermocline in the ocean.
                  ! In reality, melt-driven upwelling can warm ocean but we are not accounting for this process.
                  oceanTemperature(iLayer,iCell) = min(newTemp,oceanTemperature(iLayer,iCell))
           
                  ! Ensure no artificially adjusted temperatures fall below the freezing temperature
                  oceanTemperature(iLayer, iCell) = max(freezingTemperature(iLayer), oceanTemperature(iLayer,iCell))
               endif
            enddo
            
            do iLayer = 1, nISMIP6OceanLayers
               if (li_is_tf_invalid(oceanTemperature(iLayer,iCell)) .or. li_is_tf_invalid(oceanSalinity(iLayer,iCell))) then
                  ismip6shelfMelt_3dThermalForcing(iLayer,iCell) = invalid_value_TF
               else 
                  ismip6shelfMelt_3dThermalForcing(iLayer, iCell) = oceanTemperature(iLayer, iCell) - freezingTemperature(iLayer)
               endif
            enddo
         endif
      enddo
      
      ! Ensure no TF where no T/S
      do iLayer = 1, nISMIP6OceanLayers
         if ( (oceanTemperature(iLayer,iCell) == 0.0_RKIND) .and. (oceanSalinity(iLayer,iCell) == 0.0_RKIND) ) then
            ismip6shelfMelt_3dThermalForcing(iLayer,iCell) = 0.0_RKIND
         endif
      enddo
      
      deallocate(freezingTemperature)

      block => block % next
      enddo
   end subroutine calculate_thermal_forcing

   !-----------------------------------------------------------------------

   !|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
   ! 
   !  function li_is_tf_invalid

   !> \brief Determines if a cell is invalid if it is within a threshold away from
   !  invalid_value_TF
   !> \author Matt Hoffman
   !> \date Jul 10 2025
   !-----------------------------------------------------------------------

   logical function li_is_tf_invalid(value)
      
      real(kind=RKIND), intent(in) :: value     
      real(kind=RKIND), pointer :: invalid_value
      real(kind=RKIND) :: eps
           
      call mpas_pool_get_config(liConfigs, 'config_invalid_value_TF', invalid_value)
           
      eps = 10.0_RKIND * epsilon(1.0_RKIND) * max(abs(value), abs(invalid_value))
      li_is_tf_invalid = (abs(value - invalid_value) < eps)
      return

   end function li_is_tf_invalid

end module li_ocean_extrap


!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
