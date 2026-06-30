#ifndef OMEGA_GLOBALSTATS_H
#define OMEGA_GLOBALSTATS_H

//===-- analysis/analysisGroups/GlobalStats.h - GlobalStats ----*- C++ -*-===//
//
/// \file
/// \brief Defines the GlobalStats analysis group for computing global
/// statistics
///
/// GlobalStats is a bundled AnalysisGroup that automates creation of spatial
/// reduction operators (mean, min, max, sum, etc.) and optional temporal
/// averaging for a set of fields. Users specify:
/// - A list of field names to analyze
/// - A list of spatial statistics to compute (e.g., "Mean", "Max", "Min")
/// - Optional temporal reduction periods (e.g., "1day", "1month") for
///   time-averaged statistics
/// - Optional discrete sampling frequencies (e.g., "1hour") for instantaneous
///   snapshots without temporal averaging
///
/// The constructor builds operator chains for all combinations of fields and
/// statistics, creates separate IOStreams for each unique (period, type)
/// combination, and associates operators with the appropriate streams.
///
/// Example usage in config:
/// \code
/// GlobalStats:
///   Fields: [Temperature, Salinity]
///   SpatialStats: [Mean, Max, Min]
///   ReductionPeriod: [1day, 1month]
///   SampleFreq: [6hour]
/// \endcode
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
/// (instantaneous snapshots at specified frequencies).
///
/// The constructor reads configuration options, builds operator chains for
/// all field/statistic combinations, and creates IOStreams organized by
/// output frequency and type (temporal reduction vs. discrete sampling).
class GlobalStats : public AnalysisGroup {
 public:
   /// Constructs a GlobalStats analysis group. Reads field list, spatial
   /// statistics list, and temporal specifications (reduction periods and/or
   /// sample frequencies) from config. Builds operator chains for all
   /// field/statistic combinations, creates IOStreams for each unique
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
