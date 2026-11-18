//===-- analysis/operators/Ops.cpp - Operator registration ------*- C++ -*-===//
//
// Implementation of the Analysis method that registers all base analysis
// operators with the AnalysisOpFactory. This registration occurs during
// Analysis initialization, making all operators available for runtime
// instantiation via the factory. Each registerAllArrayVariants() call
// registers dozens of type-specific variants for a single operator template,
// covering all combinations of scalar type, rank, and memory location.
//
//===----------------------------------------------------------------------===//

#include "operators/Ops.h"
#include "Analysis.h"
#include "AnalysisOpFactory.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Registers all base analysis operators with the AnalysisOpFactory. Called
// during Analysis initialization to populate the factory registry. Each
// registration call expands a single operator template over all supported
// array type combinations (scalar type, rank, memory location), enabling
// type-safe dispatch at operator creation time.
//
// This method is defined here in operators/Ops.cpp rather than in Analysis.cpp
// to maintain separation of concerns: Analysis.cpp handles core Analysis
// initialization and lifecycle management, while this file centralizes all
// operator-specific registration logic. Co-locating the registration function
// with the operator includes (Ops.h) keeps operator management centralized in
// the operators directory and allows adding new operators without modifying
// Analysis.cpp.
void Analysis::registerAllBaseAnalysisOperators() {

   AnalysisOpFactory::registerAllArrayVariants<BinnedAccumulatorOp>(
       "BinnedAccumulator");
   AnalysisOpFactory::registerAllArrayVariants<BinaryMultiplyOp>(
       "BinaryMultiply");
   AnalysisOpFactory::registerAllArrayVariants<ExtractRegionOp>(
       "ExtractRegion");
   AnalysisOpFactory::register2DRealVariants<HorzMeanOp>("HorzMean");
   AnalysisOpFactory::registerAllArrayVariants<PrefixSumOp>("PrefixSum");
   AnalysisOpFactory::registerAllArrayVariants<PseudoToGeometricOp>(
       "PseudoToGeometric");
   AnalysisOpFactory::registerAllArrayVariants<ScalarMultiplyOp>(
       "ScalarMultiply");
   AnalysisOpFactory::register1DVariants<CoordinateBinningOp>(
       "CoordinateBinning");
   AnalysisOpFactory::registerAllArrayVariants<SpatialMaxOp>("SpatialMax");
   AnalysisOpFactory::registerAllArrayVariants<SpatialMinOp>("SpatialMin");
   AnalysisOpFactory::registerAllArrayVariants<SpatialMeanOp>("SpatialMean");
   AnalysisOpFactory::registerAllArrayVariants<SpatialStdDevOp>(
       "SpatialStdDev");
   AnalysisOpFactory::registerAllArrayVariants<TimeMeanOp>("TimeMean");
   AnalysisOpFactory::register2DRealVariants<TransectAccumulatorOp>(
       "TransectAccumulator");

} // end registerAllBaseAnalysisOperators

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
