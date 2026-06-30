#ifndef OMEGA_GROUPS_H
#define OMEGA_GROUPS_H

//===-- analysis/analysisGroups/Groups.h - AnalysisGroup includes -*- C++
//-*-===//
//
/// \file
/// \brief Convenience header that includes all bundled AnalysisGroup classes
///
/// This header provides a single include point for all bundled AnalysisGroup
/// implementations. Bundled groups automate common analysis workflows by
/// creating pre-configured operator chains and organizing output streams.
/// Users can include this header to access all available bundled groups without
/// needing to know individual header locations.
///
/// Currently available bundled groups:
/// - GlobalStats: Spatial statistics (mean, min, max, etc.) with optional
///   temporal averaging
///
/// New bundled groups should be added to this file to maintain the convenience
/// of a single include point for analysis group functionality.
///
//===----------------------------------------------------------------------===//

#include "analysisGroups/GlobalStats.h"

#endif
