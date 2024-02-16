   function ocn_setup_derived_dimensions(readDimensions, dimensionPool, configPool) result(iErr)

      use mpas_derived_types
      use mpas_pool_routines
      use mpas_io_units
      use mpas_log, only : mpas_log_write

      implicit none

      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from
      type (mpas_pool_type), intent(inout) :: configPool !< Input: Pool containing namelist options with configs
      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into

      integer :: iErr, errLevel

      integer, pointer :: nCells
      integer, pointer :: nCellsP1
      integer, pointer :: nEdges
      integer, pointer :: nEdgesP1
      integer, pointer :: maxEdges
      integer, pointer :: maxEdges2
      integer, pointer :: nVertices
      integer, pointer :: nVerticesP1
      integer, pointer :: bnds
      integer, pointer :: TWO
      integer, pointer :: R3
      integer, pointer :: FOUR
      integer, pointer :: FIVE
      integer, pointer :: SIX
      integer, pointer :: FIFTEEN
      integer, pointer :: TWENTYONE
      integer, pointer :: vertexDegree
      integer, pointer :: nVertLevels
      integer, pointer :: nVertLevelsP1
      integer, pointer :: nPoints
      integer, pointer :: nPointGroups
      integer, pointer :: maxPointsInGroup
      integer, pointer :: nRegions
      integer, pointer :: nRegionGroups
      integer, pointer :: maxRegionsInGroup
      integer, pointer :: nTransects
      integer, pointer :: nTransectGroups
      integer, pointer :: maxTransectsInGroup
      integer, pointer :: maxEdgesInTransect
      integer, pointer :: nForcingGroupsMax
      integer, pointer :: maxTidalConstituents
      integer, pointer :: maxTidalConstituentsX2
      integer, pointer :: nStokesDriftWavenumbers
      integer, pointer :: nSfcAreaWeightedAvgFields
      integer, pointer :: nOceanRegions
      integer, pointer :: nTemperatureBins
      integer, pointer :: nSalinityBins
      integer, pointer :: nTemperatureBinsP1
      integer, pointer :: nSalinityBinsP1
      integer, pointer :: nOceanRegionsTmpCensus
      integer, pointer :: nLayerVolWeightedAvgFields
      integer, pointer :: nOceanRegionsTmp
      integer, pointer :: nZonalMeanBins
      integer, pointer :: nZonalMeanBinsP1
      integer, pointer :: nMerHeatTransBins
      integer, pointer :: nMerHeatTransBinsP1
      integer, pointer :: nParticles
      integer, pointer :: nBuoyancySurfaces
      integer, pointer :: nBuoyancyLayers
      integer, pointer :: nBuoyancyLayersP1
      integer, pointer :: nMocStreamfunctionBins
      integer, pointer :: nMocStreamfunctionBinsP1
      integer, pointer :: nDepthTracerIC
      integer, pointer :: nDepthEcosysIC
      integer, pointer :: config_vert_levels
      integer, pointer :: config_AM_zonalMean_num_bins
      integer, pointer :: config_AM_meridionalHeatTransport_num_bins
      integer, pointer :: config_AM_eliassenPalm_nBuoyancyLayers
      integer, pointer :: config_AM_mocStreamfunction_num_bins
      integer, pointer :: config_global_ocean_tracer_vert_levels
      integer, pointer :: config_global_ocean_ecosys_vert_levels

      iErr = 0
      errLevel = mpas_pool_get_error_level()
      call mpas_pool_set_error_level(MPAS_POOL_SILENT)

      nullify(config_vert_levels)
      call mpas_pool_get_config(configPool, 'config_vert_levels', config_vert_levels)
      nullify(config_AM_zonalMean_num_bins)
      call mpas_pool_get_config(configPool, 'config_AM_zonalMean_num_bins', config_AM_zonalMean_num_bins)
      nullify(config_AM_meridionalHeatTransport_num_bins)
      call mpas_pool_get_config(configPool, 'config_AM_meridionalHeatTransport_num_bins', &
config_AM_meridionalHeatTransport_num_bins)
      nullify(config_AM_eliassenPalm_nBuoyancyLayers)
      call mpas_pool_get_config(configPool, 'config_AM_eliassenPalm_nBuoyancyLayers', config_AM_eliassenPalm_nBuoyancyLayers)
      nullify(config_AM_mocStreamfunction_num_bins)
      call mpas_pool_get_config(configPool, 'config_AM_mocStreamfunction_num_bins', config_AM_mocStreamfunction_num_bins)
      nullify(config_global_ocean_tracer_vert_levels)
      call mpas_pool_get_config(configPool, 'config_global_ocean_tracer_vert_levels', config_global_ocean_tracer_vert_levels)
      nullify(config_global_ocean_ecosys_vert_levels)
      call mpas_pool_get_config(configPool, 'config_global_ocean_ecosys_vert_levels', config_global_ocean_ecosys_vert_levels)

      nullify(nCells)
      call mpas_pool_get_dimension(dimensionPool, 'nCells', nCells)
      nullify(nCellsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nCellsP1', nCellsP1)
      nullify(nEdges)
      call mpas_pool_get_dimension(dimensionPool, 'nEdges', nEdges)
      nullify(nEdgesP1)
      call mpas_pool_get_dimension(dimensionPool, 'nEdgesP1', nEdgesP1)
      nullify(maxEdges)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdges', maxEdges)
      nullify(maxEdges2)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdges2', maxEdges2)
      nullify(nVertices)
      call mpas_pool_get_dimension(dimensionPool, 'nVertices', nVertices)
      nullify(nVerticesP1)
      call mpas_pool_get_dimension(dimensionPool, 'nVerticesP1', nVerticesP1)
      nullify(bnds)
      call mpas_pool_get_dimension(dimensionPool, 'bnds', bnds)
      nullify(TWO)
      call mpas_pool_get_dimension(dimensionPool, 'TWO', TWO)
      nullify(R3)
      call mpas_pool_get_dimension(dimensionPool, 'R3', R3)
      nullify(FOUR)
      call mpas_pool_get_dimension(dimensionPool, 'FOUR', FOUR)
      nullify(FIVE)
      call mpas_pool_get_dimension(dimensionPool, 'FIVE', FIVE)
      nullify(SIX)
      call mpas_pool_get_dimension(dimensionPool, 'SIX', SIX)
      nullify(FIFTEEN)
      call mpas_pool_get_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
      nullify(TWENTYONE)
      call mpas_pool_get_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
      nullify(vertexDegree)
      call mpas_pool_get_dimension(dimensionPool, 'vertexDegree', vertexDegree)
      nullify(nVertLevels)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevels', nVertLevels)
      nullify(nVertLevelsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
      nullify(nPoints)
      call mpas_pool_get_dimension(dimensionPool, 'nPoints', nPoints)
      nullify(nPointGroups)
      call mpas_pool_get_dimension(dimensionPool, 'nPointGroups', nPointGroups)
      nullify(maxPointsInGroup)
      call mpas_pool_get_dimension(dimensionPool, 'maxPointsInGroup', maxPointsInGroup)
      nullify(nRegions)
      call mpas_pool_get_dimension(dimensionPool, 'nRegions', nRegions)
      nullify(nRegionGroups)
      call mpas_pool_get_dimension(dimensionPool, 'nRegionGroups', nRegionGroups)
      nullify(maxRegionsInGroup)
      call mpas_pool_get_dimension(dimensionPool, 'maxRegionsInGroup', maxRegionsInGroup)
      nullify(nTransects)
      call mpas_pool_get_dimension(dimensionPool, 'nTransects', nTransects)
      nullify(nTransectGroups)
      call mpas_pool_get_dimension(dimensionPool, 'nTransectGroups', nTransectGroups)
      nullify(maxTransectsInGroup)
      call mpas_pool_get_dimension(dimensionPool, 'maxTransectsInGroup', maxTransectsInGroup)
      nullify(maxEdgesInTransect)
      call mpas_pool_get_dimension(dimensionPool, 'maxEdgesInTransect', maxEdgesInTransect)
      nullify(nForcingGroupsMax)
      call mpas_pool_get_dimension(dimensionPool, 'nForcingGroupsMax', nForcingGroupsMax)
      nullify(maxTidalConstituents)
      call mpas_pool_get_dimension(dimensionPool, 'maxTidalConstituents', maxTidalConstituents)
      nullify(maxTidalConstituentsX2)
      call mpas_pool_get_dimension(dimensionPool, 'maxTidalConstituentsX2', maxTidalConstituentsX2)
      nullify(nStokesDriftWavenumbers)
      call mpas_pool_get_dimension(dimensionPool, 'nStokesDriftWavenumbers', nStokesDriftWavenumbers)
      nullify(nSfcAreaWeightedAvgFields)
      call mpas_pool_get_dimension(dimensionPool, 'nSfcAreaWeightedAvgFields', nSfcAreaWeightedAvgFields)
      nullify(nOceanRegions)
      call mpas_pool_get_dimension(dimensionPool, 'nOceanRegions', nOceanRegions)
      nullify(nTemperatureBins)
      call mpas_pool_get_dimension(dimensionPool, 'nTemperatureBins', nTemperatureBins)
      nullify(nSalinityBins)
      call mpas_pool_get_dimension(dimensionPool, 'nSalinityBins', nSalinityBins)
      nullify(nTemperatureBinsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nTemperatureBinsP1', nTemperatureBinsP1)
      nullify(nSalinityBinsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nSalinityBinsP1', nSalinityBinsP1)
      nullify(nOceanRegionsTmpCensus)
      call mpas_pool_get_dimension(dimensionPool, 'nOceanRegionsTmpCensus', nOceanRegionsTmpCensus)
      nullify(nLayerVolWeightedAvgFields)
      call mpas_pool_get_dimension(dimensionPool, 'nLayerVolWeightedAvgFields', nLayerVolWeightedAvgFields)
      nullify(nOceanRegionsTmp)
      call mpas_pool_get_dimension(dimensionPool, 'nOceanRegionsTmp', nOceanRegionsTmp)
      nullify(nZonalMeanBins)
      call mpas_pool_get_dimension(dimensionPool, 'nZonalMeanBins', nZonalMeanBins)
      nullify(nZonalMeanBinsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nZonalMeanBinsP1', nZonalMeanBinsP1)
      nullify(nMerHeatTransBins)
      call mpas_pool_get_dimension(dimensionPool, 'nMerHeatTransBins', nMerHeatTransBins)
      nullify(nMerHeatTransBinsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nMerHeatTransBinsP1', nMerHeatTransBinsP1)
      nullify(nParticles)
      call mpas_pool_get_dimension(dimensionPool, 'nParticles', nParticles)
      nullify(nBuoyancySurfaces)
      call mpas_pool_get_dimension(dimensionPool, 'nBuoyancySurfaces', nBuoyancySurfaces)
      nullify(nBuoyancyLayers)
      call mpas_pool_get_dimension(dimensionPool, 'nBuoyancyLayers', nBuoyancyLayers)
      nullify(nBuoyancyLayersP1)
      call mpas_pool_get_dimension(dimensionPool, 'nBuoyancyLayersP1', nBuoyancyLayersP1)
      nullify(nMocStreamfunctionBins)
      call mpas_pool_get_dimension(dimensionPool, 'nMocStreamfunctionBins', nMocStreamfunctionBins)
      nullify(nMocStreamfunctionBinsP1)
      call mpas_pool_get_dimension(dimensionPool, 'nMocStreamfunctionBinsP1', nMocStreamfunctionBinsP1)
      nullify(nDepthTracerIC)
      call mpas_pool_get_dimension(dimensionPool, 'nDepthTracerIC', nDepthTracerIC)
      nullify(nDepthEcosysIC)
      call mpas_pool_get_dimension(dimensionPool, 'nDepthEcosysIC', nDepthEcosysIC)

call mpas_log_write('Assigning remaining dimensions from definitions in Registry.xml ...')
      if ( .not. associated(nCells) ) then
         allocate(nCells)
         nCells = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nCells', nCells)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'nCellsP1', nCellsP1)
      if ( .not. associated(nCellsP1) ) then
         allocate(nCellsP1)
         nCellsP1 = nCells+1
call mpas_log_write('       nCellsP1 = $i', intArgs=(/nCells+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nCellsP1', nCellsP1)
          else if ( nCellsP1 == MPAS_MISSING_DIM ) then
         nCellsP1 = nCells+1
          end if

      if ( .not. associated(nEdges) ) then
         allocate(nEdges)
         nEdges = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nEdges', nEdges)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'nEdgesP1', nEdgesP1)
      if ( .not. associated(nEdgesP1) ) then
         allocate(nEdgesP1)
         nEdgesP1 = nEdges+1
call mpas_log_write('       nEdgesP1 = $i', intArgs=(/nEdges+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nEdgesP1', nEdgesP1)
          else if ( nEdgesP1 == MPAS_MISSING_DIM ) then
         nEdgesP1 = nEdges+1
          end if

      if ( .not. associated(maxEdges) ) then
         allocate(maxEdges)
         maxEdges = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'maxEdges', maxEdges)
      end if

      if ( .not. associated(maxEdges2) ) then
         allocate(maxEdges2)
         maxEdges2 = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'maxEdges2', maxEdges2)
      end if

      if ( .not. associated(nVertices) ) then
         allocate(nVertices)
         nVertices = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nVertices', nVertices)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'nVerticesP1', nVerticesP1)
      if ( .not. associated(nVerticesP1) ) then
         allocate(nVerticesP1)
         nVerticesP1 = nVertices+1
call mpas_log_write('       nVerticesP1 = $i', intArgs=(/nVertices+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nVerticesP1', nVerticesP1)
          else if ( nVerticesP1 == MPAS_MISSING_DIM ) then
         nVerticesP1 = nVertices+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'bnds', bnds)
      if ( .not. associated(bnds) ) then
         allocate(bnds)
         bnds = 2
call mpas_log_write('       bnds = $i', intArgs=(/2/))
         call mpas_pool_add_dimension(dimensionPool, 'bnds', bnds)
          else if ( bnds == MPAS_MISSING_DIM ) then
         bnds = 2
          end if

      call mpas_pool_get_dimension(dimensionPool, 'TWO', TWO)
      if ( .not. associated(TWO) ) then
         allocate(TWO)
         TWO = 2
call mpas_log_write('       TWO = $i', intArgs=(/2/))
         call mpas_pool_add_dimension(dimensionPool, 'TWO', TWO)
          else if ( TWO == MPAS_MISSING_DIM ) then
         TWO = 2
          end if

      call mpas_pool_get_dimension(dimensionPool, 'R3', R3)
      if ( .not. associated(R3) ) then
         allocate(R3)
         R3 = 3
call mpas_log_write('       R3 = $i', intArgs=(/3/))
         call mpas_pool_add_dimension(dimensionPool, 'R3', R3)
          else if ( R3 == MPAS_MISSING_DIM ) then
         R3 = 3
          end if

      call mpas_pool_get_dimension(dimensionPool, 'FOUR', FOUR)
      if ( .not. associated(FOUR) ) then
         allocate(FOUR)
         FOUR = 4
call mpas_log_write('       FOUR = $i', intArgs=(/4/))
         call mpas_pool_add_dimension(dimensionPool, 'FOUR', FOUR)
          else if ( FOUR == MPAS_MISSING_DIM ) then
         FOUR = 4
          end if

      call mpas_pool_get_dimension(dimensionPool, 'FIVE', FIVE)
      if ( .not. associated(FIVE) ) then
         allocate(FIVE)
         FIVE = 5
call mpas_log_write('       FIVE = $i', intArgs=(/5/))
         call mpas_pool_add_dimension(dimensionPool, 'FIVE', FIVE)
          else if ( FIVE == MPAS_MISSING_DIM ) then
         FIVE = 5
          end if

      call mpas_pool_get_dimension(dimensionPool, 'SIX', SIX)
      if ( .not. associated(SIX) ) then
         allocate(SIX)
         SIX = 6
call mpas_log_write('       SIX = $i', intArgs=(/6/))
         call mpas_pool_add_dimension(dimensionPool, 'SIX', SIX)
          else if ( SIX == MPAS_MISSING_DIM ) then
         SIX = 6
          end if

      call mpas_pool_get_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
      if ( .not. associated(FIFTEEN) ) then
         allocate(FIFTEEN)
         FIFTEEN = 15
call mpas_log_write('       FIFTEEN = $i', intArgs=(/15/))
         call mpas_pool_add_dimension(dimensionPool, 'FIFTEEN', FIFTEEN)
          else if ( FIFTEEN == MPAS_MISSING_DIM ) then
         FIFTEEN = 15
          end if

      call mpas_pool_get_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
      if ( .not. associated(TWENTYONE) ) then
         allocate(TWENTYONE)
         TWENTYONE = 21
call mpas_log_write('       TWENTYONE = $i', intArgs=(/21/))
         call mpas_pool_add_dimension(dimensionPool, 'TWENTYONE', TWENTYONE)
          else if ( TWENTYONE == MPAS_MISSING_DIM ) then
         TWENTYONE = 21
          end if

      if ( .not. associated(vertexDegree) ) then
         allocate(vertexDegree)
         vertexDegree = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'vertexDegree', vertexDegree)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'nVertLevels', nVertLevels)
      if ( .not. associated(nVertLevels) ) then
         allocate(nVertLevels)
         nVertLevels = config_vert_levels
call mpas_log_write('       nVertLevels = $i (config_vert_levels)', intArgs=(/config_vert_levels/))
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevels', nVertLevels)
          else if ( nVertLevels == MPAS_MISSING_DIM ) then
         nVertLevels = config_vert_levels
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
      if ( .not. associated(nVertLevelsP1) ) then
         allocate(nVertLevelsP1)
         nVertLevelsP1 = nVertLevels+1
call mpas_log_write('       nVertLevelsP1 = $i', intArgs=(/nVertLevels+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nVertLevelsP1', nVertLevelsP1)
          else if ( nVertLevelsP1 == MPAS_MISSING_DIM ) then
         nVertLevelsP1 = nVertLevels+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nPoints', nPoints)
      if ( .not. associated(nPoints) ) then
         allocate(nPoints)
         nPoints = 1
call mpas_log_write('       nPoints = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'nPoints', nPoints)
          else if ( nPoints == MPAS_MISSING_DIM ) then
         nPoints = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nPointGroups', nPointGroups)
      if ( .not. associated(nPointGroups) ) then
         allocate(nPointGroups)
         nPointGroups = 1
call mpas_log_write('       nPointGroups = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'nPointGroups', nPointGroups)
          else if ( nPointGroups == MPAS_MISSING_DIM ) then
         nPointGroups = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'maxPointsInGroup', maxPointsInGroup)
      if ( .not. associated(maxPointsInGroup) ) then
         allocate(maxPointsInGroup)
         maxPointsInGroup = 1
call mpas_log_write('       maxPointsInGroup = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'maxPointsInGroup', maxPointsInGroup)
          else if ( maxPointsInGroup == MPAS_MISSING_DIM ) then
         maxPointsInGroup = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nRegions', nRegions)
      if ( .not. associated(nRegions) ) then
         allocate(nRegions)
         nRegions = 1
call mpas_log_write('       nRegions = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'nRegions', nRegions)
          else if ( nRegions == MPAS_MISSING_DIM ) then
         nRegions = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nRegionGroups', nRegionGroups)
      if ( .not. associated(nRegionGroups) ) then
         allocate(nRegionGroups)
         nRegionGroups = 1
call mpas_log_write('       nRegionGroups = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'nRegionGroups', nRegionGroups)
          else if ( nRegionGroups == MPAS_MISSING_DIM ) then
         nRegionGroups = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'maxRegionsInGroup', maxRegionsInGroup)
      if ( .not. associated(maxRegionsInGroup) ) then
         allocate(maxRegionsInGroup)
         maxRegionsInGroup = 1
call mpas_log_write('       maxRegionsInGroup = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'maxRegionsInGroup', maxRegionsInGroup)
          else if ( maxRegionsInGroup == MPAS_MISSING_DIM ) then
         maxRegionsInGroup = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nTransects', nTransects)
      if ( .not. associated(nTransects) ) then
         allocate(nTransects)
         nTransects = 1
call mpas_log_write('       nTransects = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'nTransects', nTransects)
          else if ( nTransects == MPAS_MISSING_DIM ) then
         nTransects = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nTransectGroups', nTransectGroups)
      if ( .not. associated(nTransectGroups) ) then
         allocate(nTransectGroups)
         nTransectGroups = 1
call mpas_log_write('       nTransectGroups = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'nTransectGroups', nTransectGroups)
          else if ( nTransectGroups == MPAS_MISSING_DIM ) then
         nTransectGroups = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'maxTransectsInGroup', maxTransectsInGroup)
      if ( .not. associated(maxTransectsInGroup) ) then
         allocate(maxTransectsInGroup)
         maxTransectsInGroup = 1
call mpas_log_write('       maxTransectsInGroup = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'maxTransectsInGroup', maxTransectsInGroup)
          else if ( maxTransectsInGroup == MPAS_MISSING_DIM ) then
         maxTransectsInGroup = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'maxEdgesInTransect', maxEdgesInTransect)
      if ( .not. associated(maxEdgesInTransect) ) then
         allocate(maxEdgesInTransect)
         maxEdgesInTransect = 1
call mpas_log_write('       maxEdgesInTransect = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'maxEdgesInTransect', maxEdgesInTransect)
          else if ( maxEdgesInTransect == MPAS_MISSING_DIM ) then
         maxEdgesInTransect = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nForcingGroupsMax', nForcingGroupsMax)
      if ( .not. associated(nForcingGroupsMax) ) then
         allocate(nForcingGroupsMax)
         nForcingGroupsMax = 4
call mpas_log_write('       nForcingGroupsMax = $i', intArgs=(/4/))
         call mpas_pool_add_dimension(dimensionPool, 'nForcingGroupsMax', nForcingGroupsMax)
          else if ( nForcingGroupsMax == MPAS_MISSING_DIM ) then
         nForcingGroupsMax = 4
          end if

      call mpas_pool_get_dimension(dimensionPool, 'maxTidalConstituents', maxTidalConstituents)
      if ( .not. associated(maxTidalConstituents) ) then
         allocate(maxTidalConstituents)
         maxTidalConstituents = 8
call mpas_log_write('       maxTidalConstituents = $i', intArgs=(/8/))
         call mpas_pool_add_dimension(dimensionPool, 'maxTidalConstituents', maxTidalConstituents)
          else if ( maxTidalConstituents == MPAS_MISSING_DIM ) then
         maxTidalConstituents = 8
          end if

      call mpas_pool_get_dimension(dimensionPool, 'maxTidalConstituentsX2', maxTidalConstituentsX2)
      if ( .not. associated(maxTidalConstituentsX2) ) then
         allocate(maxTidalConstituentsX2)
         maxTidalConstituentsX2 = maxTidalConstituents+maxTidalConstituents
call mpas_log_write('       maxTidalConstituentsX2 = $i', intArgs=(/maxTidalConstituents+maxTidalConstituents/))
         call mpas_pool_add_dimension(dimensionPool, 'maxTidalConstituentsX2', maxTidalConstituentsX2)
          else if ( maxTidalConstituentsX2 == MPAS_MISSING_DIM ) then
         maxTidalConstituentsX2 = maxTidalConstituents+maxTidalConstituents
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nStokesDriftWavenumbers', nStokesDriftWavenumbers)
      if ( .not. associated(nStokesDriftWavenumbers) ) then
         allocate(nStokesDriftWavenumbers)
         nStokesDriftWavenumbers = 6
call mpas_log_write('       nStokesDriftWavenumbers = $i', intArgs=(/6/))
         call mpas_pool_add_dimension(dimensionPool, 'nStokesDriftWavenumbers', nStokesDriftWavenumbers)
          else if ( nStokesDriftWavenumbers == MPAS_MISSING_DIM ) then
         nStokesDriftWavenumbers = 6
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nSfcAreaWeightedAvgFields', nSfcAreaWeightedAvgFields)
      if ( .not. associated(nSfcAreaWeightedAvgFields) ) then
         allocate(nSfcAreaWeightedAvgFields)
         nSfcAreaWeightedAvgFields = 36
call mpas_log_write('       nSfcAreaWeightedAvgFields = $i', intArgs=(/36/))
         call mpas_pool_add_dimension(dimensionPool, 'nSfcAreaWeightedAvgFields', nSfcAreaWeightedAvgFields)
          else if ( nSfcAreaWeightedAvgFields == MPAS_MISSING_DIM ) then
         nSfcAreaWeightedAvgFields = 36
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nOceanRegions', nOceanRegions)
      if ( .not. associated(nOceanRegions) ) then
         allocate(nOceanRegions)
         nOceanRegions = 7
call mpas_log_write('       nOceanRegions = $i', intArgs=(/7/))
         call mpas_pool_add_dimension(dimensionPool, 'nOceanRegions', nOceanRegions)
          else if ( nOceanRegions == MPAS_MISSING_DIM ) then
         nOceanRegions = 7
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nTemperatureBins', nTemperatureBins)
      if ( .not. associated(nTemperatureBins) ) then
         allocate(nTemperatureBins)
         nTemperatureBins = 32
call mpas_log_write('       nTemperatureBins = $i', intArgs=(/32/))
         call mpas_pool_add_dimension(dimensionPool, 'nTemperatureBins', nTemperatureBins)
          else if ( nTemperatureBins == MPAS_MISSING_DIM ) then
         nTemperatureBins = 32
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nSalinityBins', nSalinityBins)
      if ( .not. associated(nSalinityBins) ) then
         allocate(nSalinityBins)
         nSalinityBins = 20
call mpas_log_write('       nSalinityBins = $i', intArgs=(/20/))
         call mpas_pool_add_dimension(dimensionPool, 'nSalinityBins', nSalinityBins)
          else if ( nSalinityBins == MPAS_MISSING_DIM ) then
         nSalinityBins = 20
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nTemperatureBinsP1', nTemperatureBinsP1)
      if ( .not. associated(nTemperatureBinsP1) ) then
         allocate(nTemperatureBinsP1)
         nTemperatureBinsP1 = nTemperatureBins+1
call mpas_log_write('       nTemperatureBinsP1 = $i', intArgs=(/nTemperatureBins+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nTemperatureBinsP1', nTemperatureBinsP1)
          else if ( nTemperatureBinsP1 == MPAS_MISSING_DIM ) then
         nTemperatureBinsP1 = nTemperatureBins+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nSalinityBinsP1', nSalinityBinsP1)
      if ( .not. associated(nSalinityBinsP1) ) then
         allocate(nSalinityBinsP1)
         nSalinityBinsP1 = nSalinityBins+1
call mpas_log_write('       nSalinityBinsP1 = $i', intArgs=(/nSalinityBins+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nSalinityBinsP1', nSalinityBinsP1)
          else if ( nSalinityBinsP1 == MPAS_MISSING_DIM ) then
         nSalinityBinsP1 = nSalinityBins+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nOceanRegionsTmpCensus', nOceanRegionsTmpCensus)
      if ( .not. associated(nOceanRegionsTmpCensus) ) then
         allocate(nOceanRegionsTmpCensus)
         nOceanRegionsTmpCensus = 7
call mpas_log_write('       nOceanRegionsTmpCensus = $i', intArgs=(/7/))
         call mpas_pool_add_dimension(dimensionPool, 'nOceanRegionsTmpCensus', nOceanRegionsTmpCensus)
          else if ( nOceanRegionsTmpCensus == MPAS_MISSING_DIM ) then
         nOceanRegionsTmpCensus = 7
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nLayerVolWeightedAvgFields', nLayerVolWeightedAvgFields)
      if ( .not. associated(nLayerVolWeightedAvgFields) ) then
         allocate(nLayerVolWeightedAvgFields)
         nLayerVolWeightedAvgFields = 36
call mpas_log_write('       nLayerVolWeightedAvgFields = $i', intArgs=(/36/))
         call mpas_pool_add_dimension(dimensionPool, 'nLayerVolWeightedAvgFields', nLayerVolWeightedAvgFields)
          else if ( nLayerVolWeightedAvgFields == MPAS_MISSING_DIM ) then
         nLayerVolWeightedAvgFields = 36
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nOceanRegionsTmp', nOceanRegionsTmp)
      if ( .not. associated(nOceanRegionsTmp) ) then
         allocate(nOceanRegionsTmp)
         nOceanRegionsTmp = 7
call mpas_log_write('       nOceanRegionsTmp = $i', intArgs=(/7/))
         call mpas_pool_add_dimension(dimensionPool, 'nOceanRegionsTmp', nOceanRegionsTmp)
          else if ( nOceanRegionsTmp == MPAS_MISSING_DIM ) then
         nOceanRegionsTmp = 7
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nZonalMeanBins', nZonalMeanBins)
      if ( .not. associated(nZonalMeanBins) ) then
         allocate(nZonalMeanBins)
         nZonalMeanBins = config_AM_zonalMean_num_bins
call mpas_log_write('       nZonalMeanBins = $i (config_AM_zonalMean_num_bins)', intArgs=(/config_AM_zonalMean_num_bins/))
         call mpas_pool_add_dimension(dimensionPool, 'nZonalMeanBins', nZonalMeanBins)
          else if ( nZonalMeanBins == MPAS_MISSING_DIM ) then
         nZonalMeanBins = config_AM_zonalMean_num_bins
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nZonalMeanBinsP1', nZonalMeanBinsP1)
      if ( .not. associated(nZonalMeanBinsP1) ) then
         allocate(nZonalMeanBinsP1)
         nZonalMeanBinsP1 = nZonalMeanBins+1
call mpas_log_write('       nZonalMeanBinsP1 = $i', intArgs=(/nZonalMeanBins+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nZonalMeanBinsP1', nZonalMeanBinsP1)
          else if ( nZonalMeanBinsP1 == MPAS_MISSING_DIM ) then
         nZonalMeanBinsP1 = nZonalMeanBins+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nMerHeatTransBins', nMerHeatTransBins)
      if ( .not. associated(nMerHeatTransBins) ) then
         allocate(nMerHeatTransBins)
         nMerHeatTransBins = config_AM_meridionalHeatTransport_num_bins
call mpas_log_write('       nMerHeatTransBins = $i (config_AM_meridionalHeatTransport_num_bins)', &
intArgs=(/config_AM_meridionalHeatTransport_num_bins/))
         call mpas_pool_add_dimension(dimensionPool, 'nMerHeatTransBins', nMerHeatTransBins)
          else if ( nMerHeatTransBins == MPAS_MISSING_DIM ) then
         nMerHeatTransBins = config_AM_meridionalHeatTransport_num_bins
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nMerHeatTransBinsP1', nMerHeatTransBinsP1)
      if ( .not. associated(nMerHeatTransBinsP1) ) then
         allocate(nMerHeatTransBinsP1)
         nMerHeatTransBinsP1 = nMerHeatTransBins+1
call mpas_log_write('       nMerHeatTransBinsP1 = $i', intArgs=(/nMerHeatTransBins+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nMerHeatTransBinsP1', nMerHeatTransBinsP1)
          else if ( nMerHeatTransBinsP1 == MPAS_MISSING_DIM ) then
         nMerHeatTransBinsP1 = nMerHeatTransBins+1
          end if

      if ( .not. associated(nParticles) ) then
         allocate(nParticles)
         nParticles = MPAS_MISSING_DIM
         call mpas_pool_add_dimension(dimensionPool, 'nParticles', nParticles)
      end if

      call mpas_pool_get_dimension(dimensionPool, 'nBuoyancySurfaces', nBuoyancySurfaces)
      if ( .not. associated(nBuoyancySurfaces) ) then
         allocate(nBuoyancySurfaces)
         nBuoyancySurfaces = 1
call mpas_log_write('       nBuoyancySurfaces = $i', intArgs=(/1/))
         call mpas_pool_add_dimension(dimensionPool, 'nBuoyancySurfaces', nBuoyancySurfaces)
          else if ( nBuoyancySurfaces == MPAS_MISSING_DIM ) then
         nBuoyancySurfaces = 1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nBuoyancyLayers', nBuoyancyLayers)
      if ( .not. associated(nBuoyancyLayers) ) then
         allocate(nBuoyancyLayers)
         nBuoyancyLayers = config_AM_eliassenPalm_nBuoyancyLayers
call mpas_log_write('       nBuoyancyLayers = $i (config_AM_eliassenPalm_nBuoyancyLayers)', &
intArgs=(/config_AM_eliassenPalm_nBuoyancyLayers/))
         call mpas_pool_add_dimension(dimensionPool, 'nBuoyancyLayers', nBuoyancyLayers)
          else if ( nBuoyancyLayers == MPAS_MISSING_DIM ) then
         nBuoyancyLayers = config_AM_eliassenPalm_nBuoyancyLayers
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nBuoyancyLayersP1', nBuoyancyLayersP1)
      if ( .not. associated(nBuoyancyLayersP1) ) then
         allocate(nBuoyancyLayersP1)
         nBuoyancyLayersP1 = nBuoyancyLayers+1
call mpas_log_write('       nBuoyancyLayersP1 = $i', intArgs=(/nBuoyancyLayers+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nBuoyancyLayersP1', nBuoyancyLayersP1)
          else if ( nBuoyancyLayersP1 == MPAS_MISSING_DIM ) then
         nBuoyancyLayersP1 = nBuoyancyLayers+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nMocStreamfunctionBins', nMocStreamfunctionBins)
      if ( .not. associated(nMocStreamfunctionBins) ) then
         allocate(nMocStreamfunctionBins)
         nMocStreamfunctionBins = config_AM_mocStreamfunction_num_bins
call mpas_log_write('       nMocStreamfunctionBins = $i (config_AM_mocStreamfunction_num_bins)', &
intArgs=(/config_AM_mocStreamfunction_num_bins/))
         call mpas_pool_add_dimension(dimensionPool, 'nMocStreamfunctionBins', nMocStreamfunctionBins)
          else if ( nMocStreamfunctionBins == MPAS_MISSING_DIM ) then
         nMocStreamfunctionBins = config_AM_mocStreamfunction_num_bins
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nMocStreamfunctionBinsP1', nMocStreamfunctionBinsP1)
      if ( .not. associated(nMocStreamfunctionBinsP1) ) then
         allocate(nMocStreamfunctionBinsP1)
         nMocStreamfunctionBinsP1 = nMocStreamfunctionBins+1
call mpas_log_write('       nMocStreamfunctionBinsP1 = $i', intArgs=(/nMocStreamfunctionBins+1/))
         call mpas_pool_add_dimension(dimensionPool, 'nMocStreamfunctionBinsP1', nMocStreamfunctionBinsP1)
          else if ( nMocStreamfunctionBinsP1 == MPAS_MISSING_DIM ) then
         nMocStreamfunctionBinsP1 = nMocStreamfunctionBins+1
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nDepthTracerIC', nDepthTracerIC)
      if ( .not. associated(nDepthTracerIC) ) then
         allocate(nDepthTracerIC)
         nDepthTracerIC = config_global_ocean_tracer_vert_levels
call mpas_log_write('       nDepthTracerIC = $i (config_global_ocean_tracer_vert_levels)', &
intArgs=(/config_global_ocean_tracer_vert_levels/))
         call mpas_pool_add_dimension(dimensionPool, 'nDepthTracerIC', nDepthTracerIC)
          else if ( nDepthTracerIC == MPAS_MISSING_DIM ) then
         nDepthTracerIC = config_global_ocean_tracer_vert_levels
          end if

      call mpas_pool_get_dimension(dimensionPool, 'nDepthEcosysIC', nDepthEcosysIC)
      if ( .not. associated(nDepthEcosysIC) ) then
         allocate(nDepthEcosysIC)
         nDepthEcosysIC = config_global_ocean_ecosys_vert_levels
call mpas_log_write('       nDepthEcosysIC = $i (config_global_ocean_ecosys_vert_levels)', &
intArgs=(/config_global_ocean_ecosys_vert_levels/))
         call mpas_pool_add_dimension(dimensionPool, 'nDepthEcosysIC', nDepthEcosysIC)
          else if ( nDepthEcosysIC == MPAS_MISSING_DIM ) then
         nDepthEcosysIC = config_global_ocean_ecosys_vert_levels
          end if

      call mpas_log_write(' ')
      call mpas_log_write(' ----- done assigning dimensions from Registry.xml -----')
      call mpas_log_write(' ')
      call mpas_log_write(' ')
      call mpas_pool_set_error_level(errLevel)

   end function ocn_setup_derived_dimensions


   function ocn_setup_decomposed_dimensions(block, manager, readDimensions, dimensionPool, totalBlocks) result(iErr)

      use mpas_derived_types
      use mpas_decomp
      use mpas_pool_routines
      use mpas_io_units
      use mpas_abort, only : mpas_dmpar_global_abort
      use mpas_log, only : mpas_log_write

      implicit none

      type (block_type), intent(inout) :: block !< Input: Pointer to block
      type (mpas_streamManager_type), intent(inout) :: manager !< Input: Stream manager
      type (mpas_pool_type), intent(inout) :: readDimensions !< Input: Pool to pull read dimensions from
      type (mpas_pool_type), intent(inout) :: dimensionPool !< Input/Output: Pool to add dimensions into
      integer, intent(in) :: totalBlocks !< Input: Number of blocks

      integer :: iErr
      type (field1DInteger), pointer :: ownedIndices
      procedure (mpas_decomp_function), pointer :: decompFunc

      integer, pointer :: nParticles

      iErr = 0
      call mpas_log_write('Processing decomposed dimensions ...')

      call mpas_pool_get_dimension(readDimensions, 'nParticles', nParticles)
      if ( .not. associated(nParticles)) then
         call mpas_log_write('Dimension ''nParticles'' was not defined, and cannot be decomposed.', MPAS_LOG_WARN)
      else
         call mpas_decomp_get_method(block % domain % decompositions, 'uniform', decompFunc, iErr)
         if ( iErr /= MPAS_DECOMP_NOERR ) then
            call mpas_dmpar_global_abort('ERROR: Decomposition method ''uniform'' used by dimension ''nParticles'' does not '&
//'exist.')
         end if

         allocate(ownedIndices)
         ownedIndices % hasTimeDimension = .false.
         ownedIndices % isActive = .true.
         ownedIndices % isVarArray = .false.
         ownedIndices % isDecomposed = .false.
         ownedIndices % isPersistent = .true.
         ownedIndices % defaultValue = 0
         ownedIndices % fieldName = 'nParticlesOwnedIndices'
         ownedIndices % dimNames(1) = 'nParticles'
         iErr = decompFunc(block, manager, nParticles, totalBlocks, ownedIndices % array)
         ownedIndices % dimSizes(1) = size(ownedIndices % array, dim=1)
         call mpas_pool_add_field(block % allFields, 'nParticlesOwnedIndices', ownedIndices)
         call mpas_pool_get_dimension(block % dimensions, 'nParticles', nParticles)
         nParticles = size(ownedIndices % array, dim=1)
         call mpas_log_write('       nParticles => $i indices owned by block $i', intArgs=(/nParticles, block % blockID/))
      end if

      call mpas_log_write(' ')
      call mpas_log_write(' ----- done processing decomposed dimensions -----')
      call mpas_log_write(' ')
      call mpas_log_write(' ')

   end function ocn_setup_decomposed_dimensions
