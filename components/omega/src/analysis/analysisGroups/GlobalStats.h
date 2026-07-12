#ifndef OMEGA_GLOBALSTATS_H
#define OMEGA_GLOBALSTATS_H

//===-- analysis/analysisGroups/GlobalStats.h - GlobalStats ----*- C++ -*-===//
//
/// \file
/// \brief Defines the GlobalStats analysis group for computing global
/// statistics
///
/// GlobalStats is a bundled AnalysisGroup that automates creation of spatial
/// reduction operators (mean, min, max, etc.) and optional temporal reduction
/// for a set of fields. Users specify:
/// - A list of field names to analyze
/// - A list of spatial statistics to compute (e.g., "Mean", "Max", "Min")
/// - Optional temporal reduction periods (e.g., "1Day", "1Month") for
///   statistics over a time window
/// - Optional discrete sampling intervals (e.g., "1Hour") for instantaneous
///   snapshots without temporal reduction
///
/// The constructor builds operator chains for all combinations of fields and
/// statistics, creates separate IOStreams for each unique period and type
/// (temporal reduction vs. instantaneous output) combination, and associates
/// operators with the appropriate streams.
///
/// Example usage in config:
/// \ConfigInput
/// GlobalStats:
///   Fields: [Temperature, Salinity]
///   SpatialStats: [Mean, Max, Min]
///   ReductionPeriod: [1Day, 1Month]
///   SnapshotPeriod: [6Hour]
/// \EndConfigInput
///
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AnalysisGroup.h"
#include "Config.h"
#include "operators/Ops.h"
#include <string>

namespace OMEGA {

/// GlobalStats is a bundled AnalysisGroup that automates the creation of
/// spatial reduction operators for global statistics (mean, min, max, sum,
/// variance, etc.) across a set of fields. It supports both temporal reduction
/// (time-averaged statistics over specified periods) and discrete sampling
/// (instantaneous snapshots at specified intervals).
///
/// The constructor reads configuration options, builds operator chains for
/// all field/statistic combinations, and creates IOStreams organized by
/// output interval and type (temporal reduction vs. instantaneous output).
class GlobalStats : public AnalysisGroup {
 public:
   /// Constructs a GlobalStats analysis group. Reads field list, spatial
   /// statistics list, and temporal specifications (reduction periods and/or
   /// instantaneous output intervals) from config. Builds operator chains for
   /// all field/statistic combinations, creates IOStreams for each unique
   /// (period, type) pair, and associates operators with streams.
   GlobalStats(const std::string &GroupName, ///< [in] name of this group
               Config &AnalysisGroupOptions, ///< [in] group configuration
               Analysis *AnalysisManager     ///< [in] analysis manager
   );

   /// Default destructor
   ~GlobalStats() = default;

}; // end class GlobalStats

} // end namespace OMEGA

#endif
