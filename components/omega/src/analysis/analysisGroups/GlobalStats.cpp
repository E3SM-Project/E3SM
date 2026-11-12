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

   // Build stems: Field × SpatialStat combinations (no temporal operators yet)
   std::vector<std::string> ChainStems;
   for (const auto &VarName : VarList) {
      for (const auto &OpName : OpList) {

         // Construct operator type name (e.g., "SpatialMean")
         std::string OperatorType = "Spatial" + OpName;

         // Build stem chain string (e.g., "Temperature_SpatialMean")
         std::string StemStr = VarName + "_" + OperatorType;
         ChainStems.push_back(StemStr);
      }
   }

   // Add temporal operator names to chain, build operators and populate
   // OpChainInfos
   buildTemporalChains(ChainStems, AnalysisGroupOptions, AnalysisManager);

   // Create IOStreams organized by output frequency and type, and associate
   // operators with the appropriate streams based on OpChainInfos metadata
   createAnalysisGroupStreams(GroupName, AnalysisGroupOptions, AnalysisManager);

} // end GlobalStats constructor

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
