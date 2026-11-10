//===-- analysis/analysisGroups/MOC.cpp -------------------------*- C++ -*-===//
//
// Implementation of MOC constructor and helper methods. Reads configuration
// for binning parameters, regions, transect boundary conditions, and output
// specifications. Builds operator chains for each region to compute the MOC
// streamfunction, stores metadata for stream creation, and invokes base class
// method to create IOStreams.
//
//===----------------------------------------------------------------------===//

#include "analysisGroups/MOC.h"
#include "AnalysisOperator.h"
#include <iostream>
#include <set>

namespace OMEGA {

//------------------------------------------------------------------------------
// Constructs a MOC analysis group by reading configuration, building operator
// chains for each region (with required transect BC for non-Global regions),
// and creating IOStreams.
// Latitude-binned MOC streamfunction is computed by:
//   1. Converting VerticalPseudoVelocity to geometric coords
//   (PseudoToGeometric)
//   2. Computing vertical flux (BinaryMultiply with AreaCell)
//   3. Optionally masking to a region (ExtractRegion)
//   4. Accumulating into latitude bins (BinnedAccumulator)
//   5. Integrating horizontally south-to-north (PrefixSum, seeded by transect
//      BC for non-Global regions)
//   6. Converting to Sverdrups (ScalarMultiply × 1e-6)
// The Transects list is index-paired with non-Global entries in Regions.
// Every non-Global region must have a corresponding transect.
MOC::MOC(const std::string &GroupName, Config &AnalysisGroupOptions,
         Analysis *AnalysisManager) {

   Error Err;

   // Read binning parameters with defaults
   I4 NumBins = 180; // Default: 1-degree bins
   Err        = AnalysisGroupOptions.get("NumBins", NumBins);
   if (Err.isFail()) {
      LOG_INFO("MOC: NumBins not specified, using default: {}", NumBins);
   }

   Real MinLat = -90.0; // Default: full latitude range
   Err         = AnalysisGroupOptions.get("MinLat", MinLat);
   if (Err.isFail()) {
      LOG_INFO("MOC: MinLat not specified, using default: {}", MinLat);
   }

   Real MaxLat = 90.0;
   Err         = AnalysisGroupOptions.get("MaxLat", MaxLat);
   if (Err.isFail()) {
      LOG_INFO("MOC: MaxLat not specified, using default: {}", MaxLat);
   }

   // Read optional region list for regional MOC
   std::vector<std::string> RegionList;
   Err = AnalysisGroupOptions.get("Regions", RegionList);
   if (Err.isFail()) {
      // Default: compute global MOC only
      RegionList = {"Global"};
      LOG_INFO("MOC: Regions not specified, computing Global MOC only");
   }

   // Read transect list (index-paired with non-Global regions)
   // Every non-Global region requires a transect for its BC.
   std::vector<std::string> TransectList;
   Err = AnalysisGroupOptions.get("Transects", TransectList);
   if (Err.isFail()) {
      TransectList = {};
   }

   // Count non-Global regions and validate transect pairing
   size_t NumNonGlobal = 0;
   for (const auto &R : RegionList) {
      if (R != "Global")
         ++NumNonGlobal;
   }
   if (NumNonGlobal > 0 && TransectList.size() != NumNonGlobal) {
      ABORT_ERROR("MOC: {} non-Global region(s) specified but {} transect(s) "
                  "provided. Each non-Global region requires a transect.",
                  NumNonGlobal, TransectList.size());
   }

   parseTemporalPeriods(AnalysisGroupOptions);

   // Build operator chains for each region
   std::vector<std::string> ChainStems;
   size_t TransectIdx = 0;
   for (const auto &RegionName : RegionList) {

      // Determine mask field name (empty = global, no mask)
      std::string RegionMaskName;
      if (RegionName != "Global") {
         RegionMaskName = "RegionMask" + RegionName;
      }

      // Build internal transect BC chain for non-Global regions
      std::string BCFieldName;
      if (RegionName != "Global") {
         BCFieldName =
             buildTransectBCChain(TransectList[TransectIdx++], AnalysisManager);
      }

      // Build the MOC operator chain for this region (pass degrees; radians
      // conversion is done inside buildMOCChain)
      std::string ChainStem =
          buildMOCChain(RegionName, RegionMaskName, NumBins, MinLat, MaxLat,
                        BCFieldName, AnalysisManager);
      ChainStems.push_back(ChainStem);

      // Build the per-region depth coordinate chain (area-weighted horizontal
      // mean of GeomZInterface, masked to the region for non-Global regions).
      std::string DepthFieldName =
          buildDepthCoordChain(RegionName, RegionMaskName, AnalysisManager);
      DepthFieldNames.push_back(DepthFieldName);
   }

   // Build temporal chains with MOC-specific operator configurations
   buildTemporalChains(ChainStems, AnalysisGroupOptions, AnalysisManager,
                       ChainConfigs);

   // Stamp a readable long_name on every streamfunction output field. The
   // operators derive the units, 1e+06 m3 s-1 (Sverdrups) after the final
   // ScalarMultiply(1e-6); "Sv" would be read by udunits as the sievert.
   for (size_t i = 0; i < ChainStems.size(); ++i) {
      const auto &StemStr    = ChainStems[i];
      const auto &RegionName = RegionList[i];
      const std::string LongName =
          RegionName + " Meridional Overturning Streamfunction";

      std::vector<std::string> SFFieldNames = {StemStr};
      for (const auto &ReductionPeriod : ReductionPeriodList) {
         SFFieldNames.push_back(StemStr + "_TimeMean" + ReductionPeriod);
      }
      for (const auto &SFName : SFFieldNames) {
         if (Field::exists(SFName)) {
            auto SFField = Field::get(SFName);
            SFField->updateMetadata("long_name", LongName);
         }
      }
   }

   // Register each per-region depth coordinate into OpChainInfos.
   // The depth field is always output as an instantaneous output
   for (const auto &DepthFieldName : DepthFieldNames) {
      if (DepthFieldName.empty())
         continue;
      for (const auto &SnapshotPeriod : SnapshotPeriodList) {
         OpChainInfos.push_back(
             OpChainInfo{DepthFieldName, SnapshotPeriod, false});
      }
      for (const auto &ReductionPeriod : ReductionPeriodList) {
         OpChainInfos.push_back(
             OpChainInfo{DepthFieldName, ReductionPeriod, true});
      }
   }

   // Create IOStreams organized by output frequency and associate operators
   createAnalysisGroupStreams(GroupName, AnalysisGroupOptions, AnalysisManager);

   // Add the static bin boundary coordinate field to every output stream so
   // that the latitude axis values are available in the output file.
   if (!BinBoundaryFieldName.empty()) {
      std::set<std::string> StreamNamesUsed;
      for (const auto &Info : OpChainInfos) {
         std::string SName = GroupName + "_" + Info.FreqStr +
                             (Info.IsTimeReduction ? "TimeStats" : "Instants");
         StreamNamesUsed.insert(SName);
      }
      for (const auto &SName : StreamNamesUsed) {
         IOStream::get(SName)->addField(BinBoundaryFieldName);
      }

      // Add each region's depth coordinate field to every output stream so the
      // vertical (depth) axis values are available in the output file.
      for (const auto &SName : StreamNamesUsed) {
         for (const auto &DepthFieldName : DepthFieldNames) {
            if (!DepthFieldName.empty()) {
               IOStream::get(SName)->addField(DepthFieldName);
            }
         }
      }
   }

} // end MOC constructor

//------------------------------------------------------------------------------
// Builds the complete operator chain string for computing MOC streamfunction
// for a single region. The chain consists of:
//   Initialization (once):
//     0. LatCell_CoordinateBinning: Assigns cells to latitude bins (static)
//   Runtime (each timestep):
//     1. VerticalPseudoVelocity_PseudoToGeometric: Convert to geometric coords
//     2. VerticalVelocity_BinaryMultiply(AreaCell): Compute vertical flux
//     3. ExtractRegion(RegionName): Apply regional mask (regional MOC only)
//     4. VerticalFlux_BinnedAccumulator(LatCell_BinIndex): Accumulate into bins
//     5. BinnedFlux_PrefixSum: Horizontal integration (south to north)
//     6. MOC_ScalarMultiply: Convert to Sverdrups (Sv)
//
// This method constructs the complete chain string that will later be parsed
// by buildTemporalChains(), which appends temporal operators and calls
// parseChainAndBuildOps().
//
// Returns the complete chain string (stem) for temporal operator appending
std::string MOC::buildMOCChain(const std::string &RegionName,
                               const std::string &RegionMaskName, I4 NumBins,
                               Real MinLat, Real MaxLat,
                               const std::string &BCFieldName,
                               Analysis *AnalysisManager) {

   // -------------------------------------------------------------------------
   // Build the complete operator chain string
   // -------------------------------------------------------------------------
   // The complete chain encodes all transformation steps:
   //
   // For Global MOC:
   // VerticalPseudoVelocity_PseudoToGeometric_BinaryMultiply(AreaCell)_
   // BinnedAccumulator(LatCell_BinIndex)_PrefixSum_ScalarMultiply(1e-6)
   //
   // For Regional MOC:
   // VerticalPseudoVelocity_PseudoToGeometric_BinaryMultiply(AreaCell)_
   // ExtractRegion(RegionMask)_BinnedAccumulator(LatCell_BinIndex)_
   // PrefixSum(BC=BCName)_ScalarMultiply(1e-6)
   //
   // This will be parsed by Analysis::parseChainAndBuildOps() which:
   // 1. Splits by underscores to extract field names and operator types
   // 2. Creates operator instances with appropriate configurations
   // 3. Chains them together in sequence

   // Start with base vertical velocity field and conversion to geometric coords
   std::string ChainStr = "VerticalPseudoVelocity_PseudoToGeometric";

   // Multiply by AreaCell to get flux
   ChainStr += "_BinaryMultiply(AreaCell)";

   // Insert ExtractRegionOp for regional MOC after flux computation
   if (!RegionMaskName.empty()) {
      ChainStr += "_ExtractRegion(" + RegionMaskName + ")";
   }

   // Accumulate into latitude bins (uses LatCell_BinIndex from initialization)
   ChainStr += "_BinnedAccumulator(LatCell_BinIndex)";

   // Horizontal integration (south to north across latitude bins)
   // Non-Global regions seed with transect BC field
   if (!BCFieldName.empty()) {
      ChainStr += "_PrefixSum(BC=" + BCFieldName + ")";
   } else {
      ChainStr += "_PrefixSum";
   }

   // Convert to Sverdrups
   ChainStr += "_ScalarMultiply(1.0e-6)";

   // -------------------------------------------------------------------------
   // Store configuration for this chain
   // -------------------------------------------------------------------------
   // Each operator in the chain needs configuration. Build a merged config
   // containing parameters for all operators using makeOpConfig/opParam
   // helpers.

   // Build config with all operator parameters
   Config ChainConfig = makeOpConfig(
       opParam("NumBins", NumBins), // BinnedAccumulator: number of bins
       opParam(
           "Dimension",
           0), // PrefixSum: integrate along bin/latitude dimension (horizontal)
       opParam("Reverse",
               false), // PrefixSum: forward integration (south-to-north)
       opParam("IOName", std::string("MOC_streamfunction_") + RegionName));

   // Append config to vector (indexed in same order as RegionList/ChainStems)
   ChainConfigs.push_back(ChainConfig);

   // -------------------------------------------------------------------------
   // Build initialization chains (shared across all regions)
   // -------------------------------------------------------------------------
   // These are built immediately since they don't need temporal operators

   // Only build LatCell binning once (first region)
   static bool BinningBuilt = false;
   if (!BinningBuilt) {
      std::string BinningChainName = "LatCell_CoordinateBinning";

      // Convert latitude bounds to radians for the CoordinateBinningOp
      // (LatCell is stored in radians in the mesh)
      const Real DegToRad  = M_PI / 180.0;
      const Real MinLatRad = MinLat * DegToRad;
      const Real MaxLatRad = MaxLat * DegToRad;

      // Build binning config using makeOpConfig/opParam helpers
      Config BinningConfig = makeOpConfig(
          opParam("NumBins", NumBins),  // Number of latitude bins
          opParam("MinBin", MinLatRad), // Minimum latitude (radians)
          opParam("MaxBin", MaxLatRad)  // Maximum latitude (radians)
      );

      AnalysisManager->parseChainAndBuildOps(BinningChainName, BinningConfig);

      // -----------------------------------------------------------------------
      // Create the static latitude bin boundary coordinate field.
      // This is a 1D array of size NumBins+1 holding the southern edge of each
      // bin (plus the northern edge of the last bin) in degrees. Equivalent to
      // MPAS's binBoundaryMocStreamfunction. MinLat/MaxLat are in degrees as
      // read from config (the radians conversion is done separately for ops).
      // -----------------------------------------------------------------------
      BinBoundaryFieldName = "MOCLatBinBoundaries";

      // Match the margin CoordinateBinningOp applies before binning (see
      // CoordinateBinningOp.h: MinBin/MaxBin widened by (MaxBin-MinBin)*1.0e-6)
      // so the written boundaries are the edges the operator actually bins on.
      const Real Margin    = (MaxLat - MinLat) * 1.0e-6;
      const Real AdjMinLat = MinLat - Margin;
      const Real AdjMaxLat = MaxLat + Margin;
      const Real BinWidthDeg =
          (AdjMaxLat - AdjMinLat) / static_cast<Real>(NumBins);
      const I4 NBounds = NumBins + 1;

      // Create the dimension for NumBins+1 bin boundaries
      auto BinBoundDim = Dimension::create("NMocLatBinBoundaries", NBounds);

      std::vector<std::string> BinBoundDimNames = {"NMocLatBinBoundaries"};
      auto BinBoundField =
          Field::create(BinBoundaryFieldName,
                        "Southern latitude boundary of each MOC bin (plus "
                        "northern edge of last bin)", // Description
                        "degrees_north",              // Units
                        "latitude",                   // Standard name
                        AdjMinLat,                    // ValidMin
                        AdjMaxLat,                    // ValidMax
                        1,                            // Rank
                        BinBoundDimNames,             // Dimension names
                        false                         // Not time dependent
          );

      // Allocate and attach the boundary array
      Array1DReal BinBoundData(BinBoundaryFieldName + "_data", NBounds);
      BinBoundField->attachData<Array1DReal>(BinBoundData);

      // Compute and populate the real bin boundary values
      auto BinBoundHost = Kokkos::create_mirror_view(BinBoundData);
      for (I4 I = 0; I <= NumBins; ++I) {
         BinBoundHost(I) = AdjMinLat + I * BinWidthDeg;
      }
      Kokkos::deep_copy(BinBoundData, BinBoundHost);

      BinningBuilt = true;
   }

   // Return the complete chain string as the stem for temporal operators
   return ChainStr;

} // end buildMOCChain

//------------------------------------------------------------------------------
// Builds an internal transect BC chain. The chain computes transport through
// the transect and integrates vertically (bottom to top), producing a 1D field
// (NVertLevels) in m³/s. This field is used as the BC seed for the paired
// regional MOC PrefixSum. The chain consists of:
//   Runtime (each timestep):
//     1. PseudoThickness_PseudoToGeometric: Convert to geometric layer
//     thickness
//     2. LayerThickness_BinaryMultiply(NormalVelocity): thickness × velocity
//     3. EdgeTransport_BinaryMultiply(DvEdge): transport × edge width
//     4. TransportField_TransectAccumulator(TransectName): Accumulate across
//     transect (output registered as "TransectTransport(TransectName)" via
//     OutputName)
//     5. TransectTransport(TransectName)_PrefixSum: Vertical integration
//     (bottom to top)
//
// Registers ops immediately via parseChainAndBuildOps. The TransectAccumulator
// output is overridden to "TransectTransportTransectName" via the OutputName
// config key, keeping the BC field name short and the PrefixSum output name
// "TransectTransportTransectName_PrefixSum". Returns that short alias.
std::string MOC::buildTransectBCChain(const std::string &TransectName,
                                      Analysis *AnalysisManager) {

   // -------------------------------------------------------------------------
   // Build the complete operator chain string
   // -------------------------------------------------------------------------
   // The complete chain encodes all transformation steps:
   //
   // MeanPseudoThickEdge_PseudoToGeometric_BinaryMultiply(NormalVelocity)_
   // BinaryMultiply(DvEdge)_TransectAccumulator(TransectName)_PrefixSum
   //
   // This will be parsed by Analysis::parseChainAndBuildOps() which:
   // 1. Splits by underscores to extract field names and operator types
   // 2. Creates operator instances with appropriate configurations
   // 3. Chains them together in sequence
   std::string ChainStr = "MeanPseudoThickEdge_PseudoToGeometric";
   ChainStr += "_BinaryMultiply(NormalVelocity)";
   ChainStr += "_BinaryMultiply(DvEdge)";
   ChainStr += "_TransectAccumulator(" + TransectName + ")";
   ChainStr += "_PrefixSum";

   // Short alias for the final output of the TransectBC chain.
   std::string BCAlias = "TransectTransport" + TransectName;

   Config BCConfig = makeOpConfig(
       opParam("TransectName", TransectName),
       opParam("OutputName",
               BCAlias),          // override TransectAccumulator output name
       opParam("Dimension", 0),   // scan along vertical (dim 0 of 1D result)
       opParam("Reverse", true)); // bottom-to-top integration

   // Register immediately — not added to ChainStems, not written to output
   AnalysisManager->parseChainAndBuildOps(ChainStr, BCConfig);

   // Return the final output field name: BCAlias
   return BCAlias;

} // end buildTransectBCChain

//------------------------------------------------------------------------------
// Builds the per-region depth coordinate chain. Computes an area-weighted
// horizontal mean of the geometric interface height (GeomZInterface) to produce
// a representative interface depth for the region. For non-Global regions an
// ExtractRegion step masks GeomZInterface to the region first, so the HorzMean
// denominator (and hence the reported depth) reflects only cells in the region.
//
// The chain strings are:
//   Global:   GeomZInterface_HorzMean
//   Regional: GeomZInterface_<RegionMaskName>_HorzMean
//
// Returns the deterministic output field name so it can be added to the MOC
// output streams as the depth axis coordinate.
std::string MOC::buildDepthCoordChain(const std::string &RegionName,
                                      const std::string &RegionMaskName,
                                      Analysis *AnalysisManager) {

   // Base geometric interface height field from VertCoord
   const std::string GeomZFieldName = "GeomZInterface";

   // Build the chain string, inserting an ExtractRegion step for non-Global
   // regions. The chain parser reconstructs the running prefix at each node,
   // so the HorzMean output name is the full upstream prefix + "_HorzMean".
   std::string ChainStr = GeomZFieldName;
   if (!RegionMaskName.empty()) {
      ChainStr += "_ExtractRegion(" + RegionMaskName + ")";
   }
   ChainStr += "_HorzMean";

   // The HorzMean output field name is the upstream prefix + "_HorzMean".
   // ExtractRegion names its output as Input_MaskName, so the resulting depth
   // field name is deterministic and can be returned for stream attachment.
   std::string DepthFieldName = GeomZFieldName;
   if (!RegionMaskName.empty()) {
      DepthFieldName += "_" + RegionMaskName;
   }
   DepthFieldName += "_HorzMean";

   // Empty config
   Config Cfg;

   // Register the chain immediately (parseChainAndBuildOps handles duplicate
   // intermediate fields shared with other chains).
   AnalysisManager->parseChainAndBuildOps(ChainStr, Cfg);

   // Attach a descriptive long_name for output (once per unique field). The
   // HorzMean operator derives the units from GeomZInterface.
   if (Field::exists(DepthFieldName)) {
      auto DepthField = Field::get(DepthFieldName);
      DepthField->updateMetadata(
          "long_name",
          std::string("Area-weighted representative interface depth (") +
              RegionName + ")");
   }

   return DepthFieldName;

} // end buildDepthCoordChain

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
