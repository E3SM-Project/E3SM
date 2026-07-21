//===-- analysis/AnalysisOpFactory.cpp - Factory implementation --*- C++-*-===//
//
// Implementation of AnalysisOpFactory static methods for operator registration
// and creation. The factory maintains a runtime registry (Meyer's singleton)
// that maps fully-qualified operator type keys to constructor functions.
// Registration occurs during program initialization via static initializers.
// Creation performs type-safe dispatch by inspecting Field metadata to select
// the appropriate templated operator specialization.
//
//===----------------------------------------------------------------------===//

#include "AnalysisOpFactory.h"
#include <iostream>

namespace OMEGA {

//------------------------------------------------------------------------------
// Registers an operator variant by adding it to the factory registry. Checks
// for duplicate registration and aborts if the key already exists, which
// indicates a programming error (likely duplicate registration in static
// initializers).
void AnalysisOpFactory::registerOperator(const std::string &Label,
                                         CreatorFunc Creator) {

   auto &Reg = registry();

   // Check for duplicate registration (programming error)
   if (Reg.find(Label) != Reg.end()) {
      ABORT_ERROR("AnalysisOpFactory: Operator type {} is already registered",
                  Label);
   }

   // Add constructor function to registry
   Reg[Label] = Creator;
} // end registerOperator

//------------------------------------------------------------------------------
// Creates an operator instance by inspecting the primary upstream Field's
// metadata to construct a fully-qualified type key, looking up the constructor
// in the registry, and invoking it. The primary Field (first in UpstreamNames)
// determines the array type template parameter for the operator.
std::unique_ptr<AnalysisOperator>
AnalysisOpFactory::createOp(const std::string &OpType,
                            const std::vector<std::string> &UpstreamNames,
                            Config Options) {

   // Validate that all upstream Fields exist before attempting creation
   for (const auto &FieldName : UpstreamNames) {
      auto FieldPtr = Field::get(FieldName);
      if (!FieldPtr) {
         ABORT_ERROR("Field '{}' not found for operator creation", FieldName);
      }
   }

   // Extract metadata from primary upstream Field (determines array type)
   auto FieldPtr       = Field::get(UpstreamNames[0]);
   ArrayDataType DType = FieldPtr->getType();
   int Rank            = FieldPtr->getNumDims();
   ArrayMemLoc MemLoc  = FieldPtr->getMemoryLocation();

   // Map metadata to array type name (e.g., "Array2DR8_Device")
   std::string arrayTypeName = getArrayTypeName(DType, Rank, MemLoc);

   // Build fully-qualified operator key (e.g., "SpatialMean_Array2DR8_Device")
   std::string FullOpType = OpType + "_" + arrayTypeName;

   // Look up constructor in registry
   auto &Reg = registry();
   auto it   = Reg.find(FullOpType);

   if (it == Reg.end()) {
      ABORT_ERROR("Operator type {} not found", FullOpType);
   }

   // Invoke the registered constructor function
   return it->second(UpstreamNames, Options);
} // end createOp

//------------------------------------------------------------------------------
// Checks whether an operator variant with the given key exists in the registry.
// Returns true if found, false otherwise.
bool AnalysisOpFactory::hasOperator(const std::string &Type) {
   auto &Reg = registry();
   return Reg.find(Type) != Reg.end();
} // end hasOperator

//------------------------------------------------------------------------------
// Constructs an array type name string from Field metadata by expanding the
// OMEGA_ANALYSIS_ARRAY_TYPES macro and checking each combination. Returns a
// string like "Array2DR8_Device" that is used as part of the operator variant
// key. Aborts if the metadata combination is not supported.
std::string AnalysisOpFactory::getArrayTypeName(ArrayDataType DType, int Rank,
                                                ArrayMemLoc MemLoc) {
// Define a macro to check if metadata matches this array type
#define TRY_ARRAY_TYPE(dt, r, ml, ArrayT)          \
   if (DType == dt && Rank == r && MemLoc == ml) { \
      return std::string(#ArrayT) + "_" + #ml;     \
   }

   // Expand macro over all supported array type combinations
   OMEGA_ANALYSIS_ARRAY_TYPES(TRY_ARRAY_TYPE)

#undef TRY_ARRAY_TYPE

   // If we reach here, the array type is not supported
   ABORT_ERROR(
       "Unsupported array type/Rank/location: DType={}, Rank={}, MemLoc={}",
       static_cast<int>(DType), Rank, static_cast<int>(MemLoc));

   return {};
} // end getArrayTypeName

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
