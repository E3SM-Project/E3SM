//===-- analysis/analysisGroups/GlobalStats.cpp - ----------------*- C++-*-===//
//
// Implementation of GlobalStats constructor. Reads configuration for field
// list, spatial statistics, and temporal specifications (reduction periods
// and/or instantaneous output intervals). Builds operator chains for all
// field/statistic combinations, stores metadata for stream creation, and
// invokes base class method to create IOStreams organized by output frequency
// and type.
//
//===----------------------------------------------------------------------===//

#include "analysisGroups/GlobalStats.h"
#include <iostream>

namespace OMEGA {

//------------------------------------------------------------------------------
// Constructs a GlobalStats analysis group by reading configuration, building
// operator chains for all field/statistic combinations, and creating IOStreams.
// For each field and statistic, creates chains with optional temporal averaging
// (ReductionPeriod) and/or instantaneous output (SnapshotPeriod). Chains are
// grouped by their output characteristics and associated with appropriate
// IOStreams.
GlobalStats::GlobalStats(const std::string &GroupName,
                         Config &AnalysisGroupOptions,
                         Analysis *AnalysisManager) {

   Error Err1;
   Error Err2;

   // Read required field list from configuration
   std::vector<std::string> VarList;
   Err1 = AnalysisGroupOptions.get("Fields", VarList);
   CHECK_ERROR_ABORT(Err1, "GlobalStats: Fields list not found in Config");

   // Read required spatial statistics list from configuration
   // Each statistic name (e.g., "Mean", "Max") is prefixed with "Spatial"
   std::vector<std::string> OpList;
   Err1 = AnalysisGroupOptions.get("SpatialStats", OpList);
   CHECK_ERROR_ABORT(Err1,
                     "GlobalStats: SpatialStats list not found in Config");

   // Read optional temporal reduction periods (e.g., "1Day", "1Month")
   // If present, creates time-averaged output
   std::vector<std::string> ReductionPeriodList;
   Err1 = AnalysisGroupOptions.get("ReductionPeriod", ReductionPeriodList);

   // Read optional instantaneous output intervals (e.g., "6Hour")
   // If present, creates instantaneous snapshot output
   std::vector<std::string> SnapshotPeriodList;
   Err2 = AnalysisGroupOptions.get("SnapshotPeriod", SnapshotPeriodList);

   // Validate that at least one temporal specification is provided
   if (Err1.isFail() and Err2.isFail()) {
      ABORT_ERROR("GlobalStats: Error reading both ReductionPeriod and "
                  "SnapshotPeriod from Config, at least one must be present");
   }

   // Build operator chains for all field/statistic combinations
   for (const auto &VarName : VarList) {
      for (const auto &OpName : OpList) {

         // Construct operator type name (e.g., "SpatialMean")
         std::string OperatorType = "Spatial" + OpName;
         std::string ChainStr;

         std::string NewOpChainName = VarName + "_Spatial" + OpName;

         // Create temporal reduction chains: Field -> SpatialOp -> TimeMean
         // These produce time-averaged statistics over specified periods
         for (const auto &ReductionPeriod : ReductionPeriodList) {
            // Build chain string (e.g., "Temperature_SpatialMean_TimeMean1day")
            ChainStr =
                VarName + "_" + OperatorType + "_TimeMean" + ReductionPeriod;

            // Store metadata for later stream creation and association
            OpChainInfos.push_back(
                OpChainInfo{ChainStr, ReductionPeriod, true});

            // Parse chain string and instantiate operators
            AnalysisManager->parseChainAndBuildOps(ChainStr);
         }

         // Create instantaneous output chains: Field -> SpatialOp (no temporal
         // reduction) These produce instantaneous snapshots at specified
         // frequencies
         if (!SnapshotPeriodList.empty()) {
            // Build chain string without temporal operator
            ChainStr = VarName + "_" + OperatorType;

            // Parse chain string and instantiate operators (if not already
            // built)
            AnalysisManager->parseChainAndBuildOps(ChainStr);
         }

         // Store metadata for each instantaneous output interval
         for (const auto &SnapshotPeriod : SnapshotPeriodList) {
            // Same chain string, but different frequency and no temporal
            // reduction
            OpChainInfos.push_back(
                OpChainInfo{ChainStr, SnapshotPeriod, false});
         }
      }
   }

   // Create IOStreams organized by output frequency and type, and associate
   // operators with the appropriate streams based on OpChainInfos metadata
   createAnalysisGroupStreams(GroupName, AnalysisGroupOptions, AnalysisManager);

} // end GlobalStats constructor

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
