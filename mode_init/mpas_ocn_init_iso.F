! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_init_iso
!
!> \brief MPAS ocean initialize case -- Idealized Southern Ocean (ISO)
!> \author Juan A. Saenz, based on idealized_acc and others
!> \date   12/08/2014
!> \details
!>  This module contains the routines for initializing the
!>  the idealized Southern Ocean (ISO) test case
!
!-----------------------------------------------------------------------

module ocn_init_iso

   use mpas_kind_types
   use mpas_io_units
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants

   use ocn_constants
   use ocn_init_vertical_grids
   use ocn_init_cell_markers

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

   public :: ocn_init_setup_iso, &
             ocn_init_validate_iso

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_init_setup_iso
!
!> \brief   Setup for ISO test case
!> \author  Juan A. Saenz
!> \date    02/26/2014
!> \details
!>  This routine sets up the initial conditions for the
!>  Idealized Southern Ocean configuration.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_setup_iso(domain, iErr)!{{{

   !--------------------------------------------------------------------

      type (domain_type), intent(inout) :: domain
      integer, intent(out) :: iErr

      ! local work variables
      type (block_type), pointer :: block_ptr

      type (mpas_pool_type), pointer :: meshPool, verticalMeshPool, statePool, forcingPool, tracersPool
      type (mpas_pool_type), pointer :: tracersSurfaceRestoringFieldsPool, tracersInteriorRestoringFieldsPool

      integer :: iCell, k, idx

      real (kind=RKIND) :: distance, xDistance, yDistance, zMid, sphereRadius
      real (kind=RKIND) :: currentLon, currentLat
      real (kind=RKIND) :: location, amplitude
      real (kind=RKIND) :: Tbottom, Tmin, TminGlobal
      real (kind=RKIND) :: depth, contSlopeWidthRad, widthWindASFRad, windStress
      real (kind=RKIND) :: widthQSouth, widthQMiddle, widthQNorth, heatFluxZonal, heatFlux1, heatFlux2
      real (kind=RKIND) :: temperature
      real (kind=RKIND), dimension(30) :: featureDepth

      real (kind=RKIND), dimension(:), pointer :: interfaceLocations


      ! Define config variable pointers
      character (len=StrKIND), pointer :: config_init_configuration, config_vertical_grid

      integer, pointer :: config_iso_vert_levels
      real (kind=RKIND), pointer :: config_iso_main_channel_depth
      real (kind=RKIND), pointer :: config_iso_north_wall_lat
      real (kind=RKIND), pointer :: config_iso_south_wall_lat
      logical, pointer :: config_iso_ridge_flag
      real (kind=RKIND), pointer :: config_iso_ridge_center_lon
      real (kind=RKIND), pointer :: config_iso_ridge_height
      real (kind=RKIND), pointer :: config_iso_ridge_width
      logical, pointer :: config_iso_plateau_flag
      real (kind=RKIND), pointer :: config_iso_plateau_center_lon
      real (kind=RKIND), pointer :: config_iso_plateau_center_lat
      real (kind=RKIND), pointer :: config_iso_plateau_height
      real (kind=RKIND), pointer :: config_iso_plateau_radius
      real (kind=RKIND), pointer :: config_iso_plateau_slope_width
      logical, pointer :: config_iso_shelf_flag
      real (kind=RKIND), pointer :: config_iso_shelf_depth
      real (kind=RKIND), pointer :: config_iso_shelf_width
      logical, pointer :: config_iso_cont_slope_flag
      real (kind=RKIND), pointer :: config_iso_max_cont_slope
      logical, pointer :: config_iso_embayment_flag
      real (kind=RKIND), pointer :: config_iso_embayment_radius
      real (kind=RKIND), pointer :: config_iso_embayment_depth
      real (kind=RKIND), pointer :: config_iso_embayment_center_lon
      real (kind=RKIND), pointer :: config_iso_embayment_center_lat
      logical, pointer :: config_iso_depression_flag
      real (kind=RKIND), pointer :: config_iso_depression_width
      real (kind=RKIND), pointer :: config_iso_depression_depth
      real (kind=RKIND), pointer :: config_iso_depression_center_lon
      real (kind=RKIND), pointer :: config_iso_depression_south_lat
      real (kind=RKIND), pointer :: config_iso_depression_north_lat
      real (kind=RKIND), pointer :: config_iso_salinity
      real (kind=RKIND), pointer :: config_iso_wind_stress_max
      real (kind=RKIND), pointer :: config_iso_asf_wind
      real (kind=RKIND), pointer :: config_iso_acc_wind
      real (kind=RKIND), pointer :: config_iso_wind_trans
      real (kind=RKIND), pointer :: config_iso_heat_flux_south
      real (kind=RKIND), pointer :: config_iso_heat_flux_middle
      real (kind=RKIND), pointer :: config_iso_heat_flux_north
      real (kind=RKIND), pointer :: config_iso_heat_flux_lat_ss
      real (kind=RKIND), pointer :: config_iso_heat_flux_lat_sm
      real (kind=RKIND), pointer :: config_iso_heat_flux_lat_mn
      real (kind=RKIND), pointer :: config_iso_surface_temperature_piston_velocity
      real (kind=RKIND), pointer :: config_iso_initial_temp_t1
      real (kind=RKIND), pointer :: config_iso_initial_temp_t2
      real (kind=RKIND), pointer :: config_iso_initial_temp_h0
      real (kind=RKIND), pointer :: config_iso_initial_temp_h1
      real (kind=RKIND), pointer :: config_iso_initial_temp_mt
      real (kind=RKIND), pointer :: config_iso_initial_temp_latS
      real (kind=RKIND), pointer :: config_iso_initial_temp_latN
      real (kind=RKIND), pointer :: config_iso_region1_center_lon
      real (kind=RKIND), pointer :: config_iso_region1_center_lat
      real (kind=RKIND), pointer :: config_iso_region2_center_lon
      real (kind=RKIND), pointer :: config_iso_region2_center_lat
      real (kind=RKIND), pointer :: config_iso_region3_center_lon
      real (kind=RKIND), pointer :: config_iso_region3_center_lat
      real (kind=RKIND), pointer :: config_iso_region4_center_lon
      real (kind=RKIND), pointer :: config_iso_region4_center_lat
      logical, pointer :: config_iso_heat_flux_region1_flag
      real (kind=RKIND), pointer :: config_iso_heat_flux_region1
      real (kind=RKIND), pointer :: config_iso_heat_flux_region1_radius
      logical, pointer :: config_iso_heat_flux_region2_flag
      real (kind=RKIND), pointer :: config_iso_heat_flux_region2
      real (kind=RKIND), pointer :: config_iso_heat_flux_region2_radius
      real (kind=RKIND), pointer :: config_iso_temperature_sponge_t1
      real (kind=RKIND), pointer :: config_iso_temperature_sponge_h1
      real (kind=RKIND), pointer :: config_iso_temperature_sponge_l1
      real (kind=RKIND), pointer :: config_iso_temperature_sponge_tau1

      logical, pointer :: config_iso_temperature_restore_region1_flag
      real (kind=RKIND), pointer :: config_iso_temperature_restore_t1
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcx1
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcy1
      logical, pointer :: config_iso_temperature_restore_region2_flag
      real (kind=RKIND), pointer :: config_iso_temperature_restore_t2
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcx2
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcy2
      logical, pointer :: config_iso_temperature_restore_region3_flag
      real (kind=RKIND), pointer :: config_iso_temperature_restore_t3
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcx3
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcy3
      logical, pointer :: config_iso_temperature_restore_region4_flag
      real (kind=RKIND), pointer :: config_iso_temperature_restore_t4
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcx4
      real (kind=RKIND), pointer :: config_iso_temperature_restore_lcy4

      ! Define dimension pointers
      integer, pointer :: nVertLevels, nCells, nVertLevelsP1
      integer, pointer :: index_temperature, index_salinity, index_tracer1

      ! Define variable pointers
      logical, pointer :: on_a_sphere
      integer, dimension(:), pointer :: maxLevelCell
      real (kind=RKIND), dimension(:), pointer :: refBottomDepth, bottomCell, refZMid
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness, restingThickness
      real (kind=RKIND), pointer :: sphere_radius
      real (kind=RKIND), dimension(:), pointer :: lonCell, latCell, bottomDepth
      real (kind=RKIND), dimension(:,:,:), pointer :: activeTracers, debugTracers
      real (kind=RKIND), dimension(:), pointer :: sensibleHeatFlux
      real (kind=RKIND), dimension(:), pointer :: windStressZonal, windStressMeridional
      real (kind=RKIND), dimension(:, :), pointer ::    activeTracersPistonVelocity, activeTracersSurfaceRestoringValue
      real (kind=RKIND), dimension(:, :, :), pointer :: activeTracersInteriorRestoringValue, activeTracersInteriorRestoringRate

      ! Define variables for the config_iso_ variables
      real (kind=RKIND) :: mainChannelDepth, northWallLat, southWallLat
      logical :: ridgeFlag
      real (kind=RKIND) :: ridgeCenterLon, ridgeHeight, ridgeWidth
      logical :: plateauFlag
      real (kind=RKIND) :: plateauCenterLon, plateauCenterLat
      real (kind=RKIND) :: plateauHeight, plateauRadius, plateauSlopeWidth
      logical :: shelfFlag
      real (kind=RKIND) :: shelfDepth, shelfWidth
      logical :: contSlopeFlag
      real (kind=RKIND) :: maxContSlope
      logical :: embaymentFlag
      real (kind=RKIND) :: embaymentRadius, embaymentDepth, embaymentCenterLon, embaymentCenterLat
      logical :: depressionFlag
      real (kind=RKIND) :: depressionWidth, depressionDepth
      real (kind=RKIND) :: depressionCenterLon, depressionSouthLat, depressionNorthLat
      real (kind=RKIND) :: salinity0
      real (kind=RKIND) :: windStressMax, windASF, windACC, latWindTrans
      real (kind=RKIND) :: QSouth, QNorth, QMiddle, transSS, transSM, transMN
      real (kind=RKIND) :: tempPistonVel, tempT1, tempT2, temph1, tempmT, temph0, tempLatS, tempLatN
      real (kind=RKIND) :: regionCenterLat1, regionCenterLon1, regionCenterLat2, regionCenterLon2
      real (kind=RKIND) :: regionCenterLat3, regionCenterLon3, regionCenterLat4, regionCenterLon4
      logical :: heatRegionFlag1
      real (kind=RKIND) :: heatRegion1flux, heatRegion1Radius
      logical :: heatRegionFlag2
      real (kind=RKIND) :: heatRegion2flux, heatRegion2Radius
      real (kind=RKIND) :: tempSpongeT1, tempSpongeh1, tempSpongeWeightL1, tempSpongeTau1
      logical :: tempRestoreFlag1
      real (kind=RKIND) :: tempRestoreT1, tempRestoreLcx1, tempRestoreLcy1
      logical :: tempRestoreFlag2
      real (kind=RKIND) :: tempRestoreT2, tempRestoreLcx2, tempRestoreLcy2
      logical :: tempRestoreFlag3
      real (kind=RKIND) :: tempRestoreT3, tempRestoreLcx3, tempRestoreLcy3
      logical :: tempRestoreFlag4
      real (kind=RKIND) :: tempRestoreT4, tempRestoreLcx4, tempRestoreLcy4

      iErr = 0

      call mpas_pool_get_config(domain % configs, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('iso')) return

      ! get config variables
      call mpas_pool_get_config(domain % configs, 'config_vertical_grid', config_vertical_grid)
      call mpas_pool_get_config(domain % configs, 'config_iso_vert_levels', config_iso_vert_levels)
      call mpas_pool_get_config(domain % configs, 'config_iso_main_channel_depth', config_iso_main_channel_depth)
      call mpas_pool_get_config(domain % configs, 'config_iso_north_wall_lat', config_iso_north_wall_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_south_wall_lat', config_iso_south_wall_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_ridge_flag', config_iso_ridge_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_ridge_center_lon', config_iso_ridge_center_lon)
      call mpas_pool_get_config(domain % configs, 'config_iso_ridge_height', config_iso_ridge_height)
      call mpas_pool_get_config(domain % configs, 'config_iso_ridge_width', config_iso_ridge_width)
      call mpas_pool_get_config(domain % configs, 'config_iso_plateau_flag', config_iso_plateau_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_plateau_center_lon', config_iso_plateau_center_lon)
      call mpas_pool_get_config(domain % configs, 'config_iso_plateau_center_lat', config_iso_plateau_center_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_plateau_height', config_iso_plateau_height)
      call mpas_pool_get_config(domain % configs, 'config_iso_plateau_radius', config_iso_plateau_radius)
      call mpas_pool_get_config(domain % configs, 'config_iso_plateau_slope_width', config_iso_plateau_slope_width)
      call mpas_pool_get_config(domain % configs, 'config_iso_shelf_flag', config_iso_shelf_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_shelf_depth', config_iso_shelf_depth)
      call mpas_pool_get_config(domain % configs, 'config_iso_shelf_width', config_iso_shelf_width)
      call mpas_pool_get_config(domain % configs, 'config_iso_cont_slope_flag', config_iso_cont_slope_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_max_cont_slope', config_iso_max_cont_slope)
      call mpas_pool_get_config(domain % configs, 'config_iso_embayment_flag', config_iso_embayment_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_embayment_radius', config_iso_embayment_radius)
      call mpas_pool_get_config(domain % configs, 'config_iso_embayment_depth', config_iso_embayment_depth)
      call mpas_pool_get_config(domain % configs, 'config_iso_embayment_center_lon', config_iso_embayment_center_lon)
      call mpas_pool_get_config(domain % configs, 'config_iso_embayment_center_lat', config_iso_embayment_center_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_depression_flag', config_iso_depression_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_depression_width', config_iso_depression_width)
      call mpas_pool_get_config(domain % configs, 'config_iso_depression_depth', config_iso_depression_depth)
      call mpas_pool_get_config(domain % configs, 'config_iso_depression_center_lon', config_iso_depression_center_lon)
      call mpas_pool_get_config(domain % configs, 'config_iso_depression_south_lat', config_iso_depression_south_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_depression_north_lat', config_iso_depression_north_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_salinity', config_iso_salinity)
      call mpas_pool_get_config(domain % configs, 'config_iso_wind_stress_max', config_iso_wind_stress_max)
      call mpas_pool_get_config(domain % configs, 'config_iso_asf_wind', config_iso_asf_wind)
      call mpas_pool_get_config(domain % configs, 'config_iso_acc_wind', config_iso_acc_wind)
      call mpas_pool_get_config(domain % configs, 'config_iso_wind_trans', config_iso_wind_trans)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_south', config_iso_heat_flux_south)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_middle', config_iso_heat_flux_middle)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_north', config_iso_heat_flux_north)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_lat_ss', config_iso_heat_flux_lat_ss)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_lat_sm', config_iso_heat_flux_lat_sm)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_lat_mn', config_iso_heat_flux_lat_mn)
      call mpas_pool_get_config(domain % configs, 'config_iso_surface_temperature_piston_velocity', &
                                config_iso_surface_temperature_piston_velocity)
      call mpas_pool_get_config(domain % configs, 'config_iso_initial_temp_t1', config_iso_initial_temp_t1)
      call mpas_pool_get_config(domain % configs, 'config_iso_initial_temp_t2', config_iso_initial_temp_t2)
      call mpas_pool_get_config(domain % configs, 'config_iso_initial_temp_h0', config_iso_initial_temp_h0)
      call mpas_pool_get_config(domain % configs, 'config_iso_initial_temp_h1', config_iso_initial_temp_h1)
      call mpas_pool_get_config(domain % configs, 'config_iso_initial_temp_mt', config_iso_initial_temp_mt)
      call mpas_pool_get_config(domain % configs, 'config_iso_initial_temp_latS', config_iso_initial_temp_latS)
      call mpas_pool_get_config(domain % configs, 'config_iso_initial_temp_latN', config_iso_initial_temp_latN)
      call mpas_pool_get_config(domain % configs, 'config_iso_region1_center_lon', config_iso_region1_center_lon)
      call mpas_pool_get_config(domain % configs, 'config_iso_region1_center_lat', config_iso_region1_center_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_region2_center_lon', config_iso_region2_center_lon)
      call mpas_pool_get_config(domain % configs, 'config_iso_region2_center_lat', config_iso_region2_center_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_region3_center_lon', config_iso_region3_center_lon)
      call mpas_pool_get_config(domain % configs, 'config_iso_region3_center_lat', config_iso_region3_center_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_region4_center_lon', config_iso_region4_center_lon)
      call mpas_pool_get_config(domain % configs, 'config_iso_region4_center_lat', config_iso_region4_center_lat)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_region1_flag', config_iso_heat_flux_region1_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_region1', config_iso_heat_flux_region1)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_region1_radius', config_iso_heat_flux_region1_radius)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_region2_flag', config_iso_heat_flux_region2_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_region2', config_iso_heat_flux_region2)
      call mpas_pool_get_config(domain % configs, 'config_iso_heat_flux_region2_radius', config_iso_heat_flux_region2_radius)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_sponge_t1', config_iso_temperature_sponge_t1)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_sponge_h1', config_iso_temperature_sponge_h1)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_sponge_l1', config_iso_temperature_sponge_l1)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_sponge_tau1', config_iso_temperature_sponge_tau1)

      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_region1_flag', &
                                config_iso_temperature_restore_region1_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_t1', config_iso_temperature_restore_t1)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_lcx1', config_iso_temperature_restore_lcx1)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_lcy1', config_iso_temperature_restore_lcy1)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_region2_flag', &
                                config_iso_temperature_restore_region2_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_t2', config_iso_temperature_restore_t2)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_lcx2', config_iso_temperature_restore_lcx2)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_lcy2', config_iso_temperature_restore_lcy2)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_region3_flag', &
                                config_iso_temperature_restore_region3_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_t3', config_iso_temperature_restore_t3)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_lcx3', config_iso_temperature_restore_lcx3)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_lcy3', config_iso_temperature_restore_lcy3)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_region4_flag', &
                                config_iso_temperature_restore_region4_flag)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_t4', config_iso_temperature_restore_t4)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_lcx4', config_iso_temperature_restore_lcx4)
      call mpas_pool_get_config(domain % configs, 'config_iso_temperature_restore_lcy4', config_iso_temperature_restore_lcy4)

      call mpas_pool_get_subpool(domain % blocklist % structs, 'mesh', meshPool)
      call mpas_pool_get_dimension(meshPool, 'nVertLevels', nVertLevels)
      call mpas_pool_get_dimension(meshPool, 'nVertLevelsP1', nVertLevelsP1)
      call mpas_pool_get_config(meshPool, 'on_a_sphere', on_a_sphere)
      call mpas_pool_get_config(meshPool, 'sphere_radius', sphere_radius)
      sphereRadius = sphere_radius


      if(.not. on_a_sphere) then
        call mpas_log_write( 'ISO test case can only be defined on a spherical mesh.', MPAS_LOG_CRIT)
        iErr = 1
        return
      else
       call mpas_log_write( 'ISO test case using spherical radius of size: $f', realArgs=(/ sphereRadius /) )
      end if

      ! Define interface locations
      allocate( interfaceLocations( nVertLevelsP1 ) )
      call ocn_generate_vertical_grid( config_vertical_grid, interfaceLocations, domain % configs )

      ! assign config variables
      nVertLevels  = config_iso_vert_levels
      mainChannelDepth  = config_iso_main_channel_depth
      northWallLat = config_iso_north_wall_lat * pii/180.0_RKIND
      southWallLat = config_iso_south_wall_lat * pii/180.0_RKIND
      ridgeFlag = config_iso_ridge_flag
      ridgeCenterLon = config_iso_ridge_center_lon * pii/180.0_RKIND
      ridgeHeight  = config_iso_ridge_height
      ridgeWidth   = config_iso_ridge_width
      plateauFlag = config_iso_plateau_flag
      plateauCenterLon = config_iso_plateau_center_lon * pii/180.0_RKIND
      plateauCenterLat = config_iso_plateau_center_lat * pii/180.0_RKIND
      plateauHeight = config_iso_plateau_height
      plateauRadius = config_iso_plateau_radius
      plateauSlopeWidth = config_iso_plateau_slope_width
      shelfFlag = config_iso_shelf_flag
      shelfDepth   = config_iso_shelf_depth
      shelfWidth   = config_iso_shelf_width
      contSlopeFlag = config_iso_cont_slope_flag
      maxContSlope = config_iso_max_cont_slope
      embaymentFlag = config_iso_embayment_flag
      embaymentRadius = config_iso_embayment_radius
      embaymentDepth  = config_iso_embayment_depth
      embaymentCenterLon = config_iso_embayment_center_lon * pii/180.0_RKIND
      embaymentCenterLat = config_iso_embayment_center_lat * pii/180.0_RKIND
      depressionFlag = config_iso_depression_flag
      depressionWidth = config_iso_depression_width
      depressionDepth  = config_iso_depression_depth
      depressionCenterLon = config_iso_depression_center_lon * pii/180.0_RKIND
      depressionSouthLat = config_iso_depression_south_lat * pii/180.0_RKIND
      depressionNorthLat = config_iso_depression_north_lat * pii/180.0_RKIND
      salinity0 = config_iso_salinity
      windStressMax = config_iso_wind_stress_max
      windASF = config_iso_asf_wind
      windACC = config_iso_acc_wind
      latWindTrans = config_iso_wind_trans * pii/180.0_RKIND
      QSouth = config_iso_heat_flux_south
      QMiddle = config_iso_heat_flux_middle
      QNorth = config_iso_heat_flux_north
      transSS = config_iso_heat_flux_lat_ss * pii/180.0_RKIND
      transSM = config_iso_heat_flux_lat_sm * pii/180.0_RKIND
      transMN = config_iso_heat_flux_lat_mn * pii/180.0_RKIND
      tempPistonVel = config_iso_surface_temperature_piston_velocity
      tempT1 = config_iso_initial_temp_t1
      tempT2 = config_iso_initial_temp_t2
      temph0 = config_iso_initial_temp_h0
      temph1 = config_iso_initial_temp_h1
      tempmT = config_iso_initial_temp_mt
      tempLatS = config_iso_initial_temp_latS * pii/180.0_RKIND
      tempLatN = config_iso_initial_temp_latN * pii/180.0_RKIND
      regionCenterLon1 = config_iso_region1_center_lon * pii/180.0_RKIND
      regionCenterLat1 = config_iso_region1_center_lat * pii/180.0_RKIND
      regionCenterLon2 = config_iso_region2_center_lon * pii/180.0_RKIND
      regionCenterLat2 = config_iso_region2_center_lat * pii/180.0_RKIND
      regionCenterLon3 = config_iso_region3_center_lon * pii/180.0_RKIND
      regionCenterLat3 = config_iso_region3_center_lat * pii/180.0_RKIND
      regionCenterLon4 = config_iso_region4_center_lon * pii/180.0_RKIND
      regionCenterLat4 = config_iso_region4_center_lat * pii/180.0_RKIND
      heatRegionFlag1 = config_iso_heat_flux_region1_flag
      heatRegion1flux = config_iso_heat_flux_region1
      heatRegion1Radius = config_iso_heat_flux_region1_radius
      heatRegionFlag2 = config_iso_heat_flux_region2_flag
      heatRegion2flux = config_iso_heat_flux_region2
      heatRegion2Radius = config_iso_heat_flux_region2_radius
      tempSpongeT1 = config_iso_temperature_sponge_t1
      tempSpongeh1 = config_iso_temperature_sponge_h1
      tempSpongeWeightL1 = config_iso_temperature_sponge_l1
      tempSpongeTau1 = config_iso_temperature_sponge_tau1

      tempRestoreFlag1 = config_iso_temperature_restore_region1_flag
      tempRestoreT1 = config_iso_temperature_restore_t1
      tempRestoreLcx1 = config_iso_temperature_restore_lcx1
      tempRestoreLcy1 = config_iso_temperature_restore_lcy1
      tempRestoreFlag2 = config_iso_temperature_restore_region2_flag
      tempRestoreT2 = config_iso_temperature_restore_t2
      tempRestoreLcx2 = config_iso_temperature_restore_lcx2
      tempRestoreLcy2 = config_iso_temperature_restore_lcy2
      tempRestoreFlag3 = config_iso_temperature_restore_region3_flag
      tempRestoreT3 = config_iso_temperature_restore_t3
      tempRestoreLcx3 = config_iso_temperature_restore_lcx3
      tempRestoreLcy3 = config_iso_temperature_restore_lcy3
      tempRestoreFlag4 = config_iso_temperature_restore_region4_flag
      tempRestoreT4 = config_iso_temperature_restore_t4
      tempRestoreLcx4 = config_iso_temperature_restore_lcx4
      tempRestoreLcy4 = config_iso_temperature_restore_lcy4


      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Setup the vertical grid
      !!!!!!!!!!!!!!!!!!!!!!!!!

      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_dimension(meshPool, 'nCells', nCells)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_array(verticalMeshPool, 'restingThickness', restingThickness)
        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)

        ! Set layerThickness and restingThickness
        do k = 1, nVertLevels
           layerThickness(k, :) = config_iso_main_channel_depth * ( interfaceLocations(k+1) - interfaceLocations(k) )
           restingThickness(k, :) = layerThickness(k, :)
        end do

         ! Set refBottomDepth
         do k = 1, nVertLevels
            refBottomDepth(k) = config_iso_main_channel_depth * interfaceLocations(k+1)
            refZMid(k) = -config_iso_main_channel_depth * (interfaceLocations(k)+interfaceLocations(k+1))/2.0_RKIND
         end do

        block_ptr => block_ptr % next

      end do

      !!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set Topography
      !!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'setting up topography'
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_array(meshPool, 'maxLevelCell', maxLevelCell)
        call mpas_pool_get_array(meshPool, 'bottomDepth', bottomDepth)
        call mpas_pool_get_array(meshPool, 'refBottomDepth', refBottomDepth)

        ! calculate the width of the continental slope,
        ! based on the specified max value of the slope of the continental slope, maxContSlope
        contSlopeWidthRad = &
           pii * 0.5_RKIND * (-shelfDepth + mainChannelDepth) / maxContSlope / sphereRadius

        do iCell = 1, nCells
           currentLon = lonCell(iCell)
           currentLat = latCell(iCell)

           bottomDepth(iCell) = 0.0_RKIND

           !!!!!!!!!!!!!!!!!!!!!!!!!
           ! Main channel
           if (currentLat <= northWallLat .and. currentLat >= southWallLat) then
              bottomDepth(iCell) = mainChannelDepth
           endif

           !!!!!!!!!!!!!!!!!!!!!!!!!
           ! set up fill-in features
           featureDepth = 1.0E6_RKIND

           ! feature 1: Add Ridge
           if (ridgeFlag) then
             distance = (currentLon - ridgeCenterLon) &
               *sphereRadius*cos(currentLat)
             if ( abs(distance) <= 0.6_RKIND * ridgeWidth ) then
                featureDepth(1) = mainChannelDepth - &
                  ridgeHeight * exp(-2.0_RKIND*(distance / ridgeWidth / 0.4_RKIND)**2)
             endif
           endif

           ! feature 2: Add Plateau
           if (plateauFlag) then
             distance = sqrt( &
               ( (currentLon - plateauCenterLon) * sphereRadius*cos(currentLat) )**2 &
               + ( (currentLat - plateauCenterLat) * sphereRadius )**2 &
               )
             if (abs(distance) <= plateauRadius) then
               featureDepth(2) = mainChannelDepth - plateauHeight
             else if (abs(distance) > plateauRadius .and. abs(distance) < plateauSlopeWidth) then
               featureDepth(2) = mainChannelDepth - plateauHeight * &
                    exp( -2 * ( (abs(distance)-plateauRadius) / plateauSlopeWidth / 0.4_RKIND ) **2 )
             endif
           endif

           ! feature 3: Add continental slope, or continental shelf break
           if (contSlopeFlag) then
             zMid = 0.5_RKIND*(mainChannelDepth+shelfDepth)
             amplitude = 0.5_RKIND*(-shelfDepth+mainChannelDepth)
             if (currentLat <= southWallLat + contSlopeWidthRad&
                 .and. currentLat > southWallLat) then
               featureDepth(3) = zMid - amplitude * sin( 0.5_RKIND*pii + pii/contSlopeWidthRad &
                 *(currentLat-southWallLat) )
             endif
           endif

           ! choose the shallowest
           bottomDepth(iCell) = min(minval(featureDepth), bottomDepth(iCell))


           !!!!!!!!!!!!!!!!!!!!!!!!!
           ! Set up dig-out features
           featureDepth = 0.0_RKIND

           ! feature 1: Continental shelf
           if (shelfFlag) then
             if (currentLat <= southWallLat .and. currentLat >= southWallLat-shelfWidth/sphereRadius) then
               featureDepth(1) = shelfDepth
             endif
           endif

           ! feature 2: Embayment
           if (embaymentFlag) then
             distance = sqrt( &
               ( (currentLon - embaymentCenterLon) * sphereRadius*cos(currentLat) )**2 &
               + ( (currentLat - embaymentCenterLat) * sphereRadius )**2 &
               )
             if(distance <= embaymentRadius .and. currentLat < embaymentCenterLat) then
               featureDepth(2) = embaymentDepth
             endif
           endif

           ! feature 3: depression
           if (depressionFlag) then
             distance = (currentLon - depressionCenterLon) * sphereRadius*cos(currentLat)
             if( abs(distance) <= 0.5_RKIND*depressionWidth &
               .and. currentLat >= depressionSouthLat .and. currentLat <= depressionNorthLat ) &
               then
               featureDepth(3) = depressionDepth
             endif
           endif

           ! choose the deepest one
           bottomDepth(iCell) = max(maxval(featureDepth), bottomDepth(iCell))


           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Set maxLevelCell to -1 for cells to be culled
           if (bottomDepth(iCell) > 0.0_RKIND) then
             maxLevelCell(iCell) = 1
           else
             maxLevelCell(iCell) = -1
           endif


           !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           ! Determine maxLevelCell based on bottomDepth and refBottomDepth
           ! Also set botomDepth based on refBottomDepth, since
           ! above bottomDepth was set with continuous analytical functions,
           ! and needs to be discrete
           if (maxLevelCell(iCell) > 0) then
             maxLevelCell(iCell) = nVertLevels
             if (nVertLevels .gt. 1) then
                do k = 1, nVertLevels
                  if (bottomDepth(iCell) < refBottomDepth(k) ) then
                      maxLevelCell(iCell) = k-1
                      bottomDepth(iCell) = refBottomDepth(k-1)
                      exit
                  end if
                end do
             end if
           end if

        enddo ! Looping through with iCell

        block_ptr => block_ptr % next
      enddo ! done setting topography

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! mark cells for culling
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
       block_ptr => domain % blocklist
       do while (associated(block_ptr))
          call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
          call ocn_mark_maxlevelcell(meshPool, iErr)
          block_ptr => block_ptr % next
       end do


      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      ! Set forcing boundary conditions and initial conditions
      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      write(*,*) 'setting up forcing and boundary conditions'
      block_ptr => domain % blocklist
      do while(associated(block_ptr))
        call mpas_pool_get_subpool(block_ptr % structs, 'mesh', meshPool)
        call mpas_pool_get_array(meshPool, 'lonCell', lonCell)
        call mpas_pool_get_array(meshPool, 'latCell', latCell)
        call mpas_pool_get_subpool(block_ptr % structs, 'state', statePool)
        call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
        call mpas_pool_get_dimension(tracersPool, 'index_temperature', index_temperature)
        call mpas_pool_get_dimension(tracersPool, 'index_salinity', index_salinity)
        call mpas_pool_get_dimension(tracersPool, 'index_tracer1', index_tracer1)
        call mpas_pool_get_array(tracersPool, 'activeTracers', activeTracers, 1)
        call mpas_pool_get_array(tracersPool, 'debugTracers', debugTracers, 1)
        call mpas_pool_get_subpool(block_ptr % structs, 'verticalMesh', verticalMeshPool)
        call mpas_pool_get_array(verticalMeshPool, 'refZMid', refZMid)
        call mpas_pool_get_subpool(block_ptr % structs, 'forcing', forcingPool)
        call mpas_pool_get_array(forcingPool, 'sensibleHeatFlux', sensibleHeatFlux)
        call mpas_pool_get_array(forcingPool, 'windStressZonal', windStressZonal)
        call mpas_pool_get_array(forcingPool, 'windStressMeridional', windStressMeridional, 1)
        call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceRestoringFields', tracersSurfaceRestoringFieldsPool)
        call mpas_pool_get_subpool(forcingPool, 'tracersInteriorRestoringFields', tracersInteriorRestoringFieldsPool)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersPistonVelocity', activeTracersPistonVelocity, 1)
        call mpas_pool_get_array(tracersSurfaceRestoringFieldsPool, 'activeTracersSurfaceRestoringValue', &
                                 activeTracersSurfaceRestoringValue, 1)
        call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringRate', &
                                 activeTracersInteriorRestoringRate, 1)
        call mpas_pool_get_array(tracersInteriorRestoringFieldsPool, 'activeTracersInteriorRestoringValue', &
                                 activeTracersInteriorRestoringValue, 1)

        activeTracersInteriorRestoringRate(:,:,:) = 0.0_RKIND
        activeTracersInteriorRestoringValue(:,:,:) = 0.0_RKIND
        activeTracersPistonVelocity(:,:)  = 0.0_RKIND
        activeTracersSurfaceRestoringValue(:,:) = 0.0_RKIND

        do iCell = 1, nCells
           currentLon = lonCell(iCell)
           currentLat = latCell(iCell)

           ! Set initial temperature
           idx = index_temperature
           do k = 1, nVertLevels
              zMid = refZMid(k)
              !temperature = tempT1 + tempT2*tanh(zMid/temph1) + tempmT * zMid
              temperature = (tempT1 + tempT2*tanh((zMid+temph0)/temph1) + tempmT * zMid) &
                * (-tempLatS+currentLat)*( 1.0_RKIND/(-tempLatS+tempLatN) )
              activeTracers(idx, k, iCell) = temperature
           enddo

           ! Set initial salinity
           idx = index_salinity
           activeTracers(idx, :, iCell) = salinity0

           ! Set up debugging tracers
           idx = index_tracer1
           if ( associated(debugTracers) ) then
              debugTracers(idx, :, iCell) = 1.0_RKIND
           end if

           ! Heat fluxes
           heatFluxZonal = 0.0_RKIND
           heatFlux1 = 0.0_RKIND
           heatFlux2 = 0.0_RKIND

           ! Setup zonally constant surface heat fluxes
           widthQSouth  = transSM - transSS
           widthQMiddle = transMN - transSM
           widthQNorth  = northWallLat - transMN
           if (currentLat > transSS .and. currentLat < transSM) then
               heatFluxZonal = QSouth*sin(pii*(currentLat-transSM)/widthQSouth)**2
           elseif (currentLat > transSM .and. currentLat < transMN) then
               heatFluxZonal = QMiddle*sin(pii*(currentLat-transMN)/widthQMiddle)**2
           elseif (currentLat > transMN .and. currentLat < northWallLat) then
               heatFluxZonal = QNorth*sin(pii*(currentLat-northWallLat)/widthQNorth)**2
           endif

           ! Setup heat flux over localized region 1
           if (heatRegionFlag1) then
              distance = sqrt( &
                ( (currentLon - regionCenterLon1) * sphereRadius*cos(currentLat) )**2 &
                + ( (currentLat - regionCenterLat1) * sphereRadius )**2 &
                )
              if (abs(distance) <= heatRegion1Radius) then
                heatFlux1 = heatRegion1flux * exp(-2.0_RKIND*(distance / 2.0_RKIND / heatRegion1Radius / 0.4_RKIND)**2)
              endif
           endif
           ! Setup heat flux over localized region 2
           if (heatRegionFlag2) then
              distance = sqrt( &
                ( (currentLon - regionCenterLon2) * sphereRadius*cos(currentLat) )**2 &
                + ( (currentLat - regionCenterLat2) * sphereRadius )**2 &
                )
              if (abs(distance) <= heatRegion2Radius) then
                heatFlux2 = heatRegion2flux * exp(-2.0_RKIND*(distance / 2.0_RKIND / heatRegion2Radius / 0.4_RKIND)**2)
              endif
           endif

           if (currentLat < transSM) then
             sensibleHeatFlux(iCell) = min(heatFluxZonal, heatFlux1, heatFlux2)
           else
             sensibleHeatFlux(iCell) = heatFluxZonal
           endif

           ! Set interior restoring
           do k = 1, nVertLevels
              zMid = refZMid(k)

              !Temperature
              !Interior restoring along northern wall
              distance = sphereRadius * ( currentLat - northWallLat)
              if(abs(distance) <= 3.0_RKIND*tempSpongeWeightL1) then
                 idx = index_temperature
                 temperature = tempSpongeT1 * exp(zMid/tempSpongeh1)
                 activeTracersInteriorRestoringValue(idx, k, iCell) = temperature

                 idx = index_temperature
                 activeTracersInteriorRestoringRate(idx, k, iCell) = exp(-abs(distance)/tempSpongeWeightL1) * ( 1.0_RKIND &
                                                                   / (tempSpongeTau1*86400.0_RKIND))
              endif

              ! Interior restoring at localized region 1
              if (tempRestoreFlag1) then
                 xDistance = (currentLon - regionCenterLon1) * sphereRadius*cos(currentLat)
                 yDistance = (currentLat - regionCenterLat1) * sphereRadius
                 if (abs(yDistance) <= tempRestoreLcy1 .and. abs(xDistance) <= tempRestoreLcx1) then
                    idx = index_temperature
                    activeTracersInteriorRestoringValue(idx, k, iCell) = TempRestoreT1

                    idx = index_temperature
                    activeTracersInteriorRestoringRate(idx, k, iCell) = ( 1.0_RKIND / (tempSpongeTau1*86400.0_RKIND)) * &
                       exp(-(2.0_RKIND*xDistance/tempRestoreLcx1)**2 - (2.0_RKIND*yDistance/tempRestoreLcy1)**2 )
                 endif
              endif

              ! Interior restoring at localized region 2
              if (tempRestoreFlag2) then
                 xDistance = (currentLon - regionCenterLon2) * sphereRadius*cos(currentLat)
                 yDistance = (currentLat - regionCenterLat2) * sphereRadius
                 if (abs(yDistance) <= tempRestoreLcy2 .and. abs(xDistance) <= tempRestoreLcx2) then
                    idx = index_temperature
                    activeTracersInteriorRestoringValue(idx, k, iCell) = TempRestoreT2

                    idx = index_temperature
                    activeTracersInteriorRestoringRate(idx, k, iCell) = ( 1.0_RKIND / (tempSpongeTau1*86400.0_RKIND)) * &
                       exp(-(2.0_RKIND*xDistance/tempRestoreLcx2)**2 - (2.0_RKIND*yDistance/tempRestoreLcy2)**2 )
                 endif
              endif

              ! Interior restoring at localized region 3
              if (tempRestoreFlag3) then
                 xDistance = (currentLon - regionCenterLon3) * sphereRadius*cos(currentLat)
                 yDistance = (currentLat - regionCenterLat3) * sphereRadius
                 if (abs(yDistance) <= tempRestoreLcy3 .and. abs(xDistance) <= tempRestoreLcx3) then
                    idx = index_temperature
                    activeTracersInteriorRestoringValue(idx, k, iCell) = TempRestoreT3

                    idx = index_temperature
                    activeTracersInteriorRestoringRate(idx, k, iCell) = ( 1.0_RKIND / (tempSpongeTau1*86400.0_RKIND)) * &
                       exp(-(2.0_RKIND*xDistance/tempRestoreLcx3)**2 - (2.0_RKIND*yDistance/tempRestoreLcy3)**2 )
                 endif
              endif

              ! Interior restoring at localized region 4
              if (tempRestoreFlag4) then
                 xDistance = (currentLon - regionCenterLon4) * sphereRadius*cos(currentLat)
                 yDistance = (currentLat - regionCenterLat4) * sphereRadius
                 if (abs(yDistance) <= tempRestoreLcy4 .and. abs(xDistance) <= tempRestoreLcx4) then
                    idx = index_temperature
                    activeTracersInteriorRestoringValue(idx, k, iCell) = TempRestoreT4

                    idx = index_temperature
                    activeTracersInteriorRestoringRate(idx, k, iCell) = ( 1.0_RKIND / (tempSpongeTau1*86400.0_RKIND)) * &
                       exp(-(2.0_RKIND*xDistance/tempRestoreLcx4)**2 - (2.0_RKIND*yDistance/tempRestoreLcy4)**2 )
                endif
               endif

              ! Salinity
              idx = index_salinity
              activeTracersInteriorRestoringValue(idx, k, iCell) = salinity0
              idx = index_salinity
              activeTracersInteriorRestoringRate(idx, k, iCell) = 0.0_RKIND

           enddo ! k = 1, nVertLevels, interior restoring loop


        end do ! iCell = 1, nCells

        ! fill activeTracersSurfaceRestoringValue surface restoring values
        ! fill activeTracersPistonVelocity with surface restoring rate
        activeTracersSurfaceRestoringValue(index_temperature,:) = activeTracers(index_temperature, 1, :)
        activeTracersPistonVelocity(index_temperature,:) = tempPistonVel
        activeTracersSurfaceRestoringValue(index_salinity,:) = activeTracers(index_salinity, 1, :)
        activeTracersPistonVelocity(index_salinity,:) = 0.0_RKIND

        ! Set wind stress
        widthWindASFRad = 1.1_RKIND*contSlopeWidthRad
        do iCell = 1, nCells
           currentLon = lonCell(iCell)
           currentLat = latCell(iCell)
           windStress = 0.0_RKIND

            ! Set wind stress over the ACC, or main channel
            if (currentLat > latWindTrans) then
              windStress = windACC * &
                sin( pii * (currentLat - latWindTrans) &
                    / (northWallLat-latWindTrans) )**2
            ! Set the wind over the continental slope front, over continental slope region
            else if (currentLat > latWindTrans - widthWindASFRad .and. currentLat < latWindTrans) then
              windStress = windASF * sin( pii * (latWindTrans-currentLat) / widthWindASFRad )**2
            endif
            windStressZonal(iCell) = windStress
            windStressMeridional(iCell) = 0.0_RKIND
        end do


        block_ptr => block_ptr % next
      end do

      write(*,*) 'exiting ocn_init_setup_iso'

   !--------------------------------------------------------------------

   end subroutine ocn_init_setup_iso!}}}

!***********************************************************************
!
!  routine ocn_init_validate_iso
!
!> \brief   Validation for ISO test case
!> \author  Juan A. Saenz
!> \date    02/26/2014
!> \details
!>  This routine validates the configuration options for the ISO test case.
!
!-----------------------------------------------------------------------

   subroutine ocn_init_validate_iso(configPool, packagePool, iocontext, iErr)!{{{

   !--------------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: configPool, packagePool
      type (mpas_io_context_type), intent(inout) :: iocontext
      integer, intent(out) :: iErr

      character (len=StrKIND), pointer :: config_init_configuration
      integer, pointer :: config_vert_levels, config_iso_vert_levels

      iErr = 0

      call mpas_pool_get_config(configPool, 'config_init_configuration', config_init_configuration)
      if(config_init_configuration .ne. trim('iso')) return

      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_iso_vert_levels', config_iso_vert_levels)

      if(config_vert_levels <= 0 .and. config_iso_vert_levels > 0) then
         config_vert_levels = config_iso_vert_levels
      else if (config_vert_levels <= 0) then
         call mpas_log_write( 'Validation failed for ISO. Not given a usable value for vertical levels.', MPAS_LOG_CRIT)
         iErr = 1
      end if

   !--------------------------------------------------------------------

   end subroutine ocn_init_validate_iso!}}}

end module ocn_init_iso





!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
