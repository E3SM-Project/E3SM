#ifndef OMEGA_OPS_H
#define OMEGA_OPS_H

//===-- analysis/operators/Ops.h - Operator includes -----------*- C++ -*-===//
//
/// \file
/// \brief Convenience header that includes all analysis operator classes
///
/// This header provides a single include point for all concrete analysis
/// operator implementations. Each operator performs a specific transformation
/// on input fields (spatial reductions, temporal averaging, etc.). Users can
/// include this header to access all available operators without needing to
/// know individual header locations.
///
/// See Analysis.md in docs/devGuide for a description of currently available
/// operators.
///
/// All operators are templated on Kokkos array type to support multiple scalar
/// types, ranks, and memory locations. The AnalysisOpFactory handles type-safe
/// dispatch based on Field metadata.
///
/// New operators should be added to this file to maintain the convenience of
/// a single include point for operator functionality.
///
//===----------------------------------------------------------------------===//

#include "operators/BinaryMultiplyOp.h"
#include "operators/BinnedAccumulatorOp.h"
#include "operators/CoordinateBinningOp.h"
#include "operators/ExtractRegionOp.h"
#include "operators/HorzMeanOp.h"
#include "operators/PrefixSumOp.h"
#include "operators/PseudoToGeometricOp.h"
#include "operators/ScalarMultiplyOp.h"
#include "operators/SpatialMaxOp.h"
#include "operators/SpatialMeanOp.h"
#include "operators/SpatialMinOp.h"
#include "operators/SpatialStdDevOp.h"
#include "operators/TimeMeanOp.h"
#include "operators/TransectAccumulatorOp.h"

#endif
