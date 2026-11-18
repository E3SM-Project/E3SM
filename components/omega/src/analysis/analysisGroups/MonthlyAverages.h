#ifndef OMEGA_MONTHLYAVERAGES_H
#define OMEGA_MONTHLYAVERAGES_H

//===-- analysis/analysisGroups/MonthlyAverages.h ----------------*- C++-*-===//
//
/// \file
/// \brief Defines the MonthlyAverages analysis group for monthly time averages
///
/// MonthlyAverages is a bundled AnalysisGroup that automates creation of
/// monthly time-averaged output for a set of requested global fields. Unlike
/// GlobalStats, it applies no spatial reduction; each requested field is fed
/// directly into a temporal-reduction (TimeMean) operator. The averaging
/// period defaults to one month, but may be overridden via an optional
/// ReductionPeriod config option.
///
/// Users specify:
/// - A list of field names to average (required)
/// - An optional ReductionPeriod list to override the default "1Month" period
///
/// Any SnapshotPeriod specified in config for this group is ignored, since the
/// group produces temporal reductions rather than instantaneous snapshots. The
/// resulting output stream is named "MonthlyAverages_<Period>TimeStats" (e.g.,
/// "MonthlyAverages_1MonthTimeStats"), consistent with the convention used by
/// other bundled analysis groups.
///
/// Example usage in config:
/// \ConfigInput
/// MonthlyAverages:
///   Enable: true
///   Fields: [Temperature, Salinity, NormalVelocity, PseudoThickness]
///   ReductionPeriod: [1Month]   # optional; defaults to [1Month]
/// \EndConfigInput
///
//===----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AnalysisGroup.h"
#include "Config.h"
#include "operators/Ops.h"
#include <string>

namespace OMEGA {

/// MonthlyAverages is a bundled AnalysisGroup that automates the creation of
/// monthly time-averaged output for a set of requested fields. No spatial
/// reduction is applied; each field is passed directly to a TimeMean operator.
/// The averaging period defaults to one month, but may be overridden by an
/// optional ReductionPeriod config option.
///
/// The constructor reads the required field list from config, builds one
/// operator chain per field (field -> TimeMean<Period>), and creates a single
/// IOStream for the time-averaged output.
class MonthlyAverages : public AnalysisGroup {
 public:
   /// Constructs a MonthlyAverages analysis group. Reads the required field
   /// list from config, applies a default reduction period of "1Month" when
   /// none is supplied (an optional ReductionPeriod overrides the default,
   /// while any SnapshotPeriod is ignored), builds one TimeMean operator chain
   /// per field, and creates the time-averaged output IOStream.
   MonthlyAverages(const std::string &GroupName, ///< [in] name of this group
                   Config &AnalysisGroupOptions, ///< [in] group configuration
                   Analysis *AnalysisManager     ///< [in] analysis manager
   );

   /// Default destructor
   ~MonthlyAverages() = default;

}; // end class MonthlyAverages

} // end namespace OMEGA

#endif
