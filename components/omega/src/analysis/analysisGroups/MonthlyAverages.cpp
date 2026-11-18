//===-- analysis/analysisGroups/MonthlyAverages.cpp - ------------*- C++-*-===//
//
// Implementation of the MonthlyAverages constructor. Reads the required field
// list from config, applies a default temporal reduction period of "1Month"
// when none is supplied (an optional ReductionPeriod in config overrides the
// default), builds one TimeMean operator chain per field, and creates the
// time-averaged output stream via the AnalysisGroup base class machinery.
//
//===----------------------------------------------------------------------===//

#include "analysisGroups/MonthlyAverages.h"
#include <vector>

namespace OMEGA {

//------------------------------------------------------------------------------
// Constructs a MonthlyAverages analysis group. Each requested field is averaged
// over a temporal window. There is no spatial reduction, so the operator chain
// stems are simply the field names, to which the base class appends a TimeMean
// operator. The averaging period defaults to "1Month" but may be overridden by
// supplying a ReductionPeriod option in config.
MonthlyAverages::MonthlyAverages(const std::string &GroupName,
                                 Config &AnalysisGroupOptions,
                                 Analysis *AnalysisManager) {

   Error Err1;

   // Read required field list from configuration
   std::vector<std::string> VarList;
   Err1 = AnalysisGroupOptions.get("Fields", VarList);
   CHECK_ERROR_ABORT(Err1, "MonthlyAverages: Fields list not found in Config");

   // Apply a default temporal reduction period of one month when the user has
   // not supplied a ReductionPeriod. If ReductionPeriod is present in config,
   // it is used as-is, allowing the averaging period to be overridden.
   std::vector<std::string> UserPeriods;
   Error PeriodErr = AnalysisGroupOptions.get("ReductionPeriod", UserPeriods);
   if (PeriodErr.isFail() || UserPeriods.empty()) {
      std::vector<std::string> DefaultPeriod{"1Month"};
      AnalysisGroupOptions.remove("ReductionPeriod");
      AnalysisGroupOptions.add("ReductionPeriod", DefaultPeriod);
   }

   // This group produces temporal reductions, not instantaneous snapshots, so
   // remove any SnapshotPeriod option.
   std::vector<std::string> SnapshotPeriods;
   Error SnapshotErr =
       AnalysisGroupOptions.get("SnapshotPeriod", SnapshotPeriods);
   if (!SnapshotErr.isFail() && !SnapshotPeriods.empty()) {
      LOG_WARN("MonthlyAverages AnalysisGroup ignores SnapshotPeriod; "
               "no instantaneous snapshots will be produced");
      AnalysisGroupOptions.remove("SnapshotPeriod");
   }

   // No spatial operator: the chain stems are simply the field names. The base
   // class appends the TimeMean operator to each stem.
   std::vector<std::string> ChainStems = VarList;

   // Add temporal operator names to chain, build operators and populate
   // OpChainInfos
   buildTemporalChains(ChainStems, AnalysisGroupOptions, AnalysisManager);

   // Create the IOStream for the monthly time-averaged output and associate
   // operators with it based on OpChainInfos metadata
   createAnalysisGroupStreams(GroupName, AnalysisGroupOptions, AnalysisManager);

} // end MonthlyAverages constructor

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
