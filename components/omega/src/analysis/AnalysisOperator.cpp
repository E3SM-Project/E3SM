//===-- analysis/AnalysisOperator.cpp - AnalysisOperator impl ---*- C++ -*-===//
//
// Implementation of AnalysisOperator base class methods. Provides default
// implementations for constructors, getters, cache validation, and resource
// management. Derived classes override the pure virtual compute() method and
// may override initialize() to perform operator-specific setup.
//
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Default constructor - initializes cache tracking variables to indicate
// no computation has occurred yet.
AnalysisOperator::AnalysisOperator() {
   // Initialize cache tracking variables
   FieldComputed = false;
   LastComputed  = TimeInstant();
} // end default constructor

//------------------------------------------------------------------------------
// Constructor with operator type name - stores the type and initializes
// cache tracking variables. Derived classes call this to set OperatorTypeName.
AnalysisOperator::AnalysisOperator(const std::string &OperatorType) {

   // Store operator type (e.g., "SpatialMean", "TimeMean")
   OperatorTypeName = OperatorType;

   // Initialize cache tracking variables
   FieldComputed = false;
   LastComputed  = TimeInstant();
} // end constructor

//------------------------------------------------------------------------------
// Base class initialization - stores pointers to mesh, vertical coordinate,
// and MPI communicator. These are used by derived classes during compute().
// Derived classes may override this to perform additional setup, but should
// call this base implementation to ensure pointers are stored.
void AnalysisOperator::initialize(const MachEnv *Env, const HorzMesh *InMesh,
                                  const VertCoord *InVCoord, Config Options) {

   // Store pointers needed during compute()
   Mesh   = InMesh;
   VCoord = InVCoord;
   Comm   = Env->getComm();

} // end initialize

//------------------------------------------------------------------------------
// Destructor - removes this operator's output Fields from the Field registry.
// Checks if each output Field exists before attempting destruction to handle
// cases where Fields may have been removed elsewhere.
AnalysisOperator::~AnalysisOperator() {
   // Clean up output Fields registered by this operator
   for (const auto &OutputName : OutputNames) {
      if (Field::exists(OutputName)) {
         Field::destroy(OutputName);
      }
   }
} // end destructor

//------------------------------------------------------------------------------
// Returns the operator type name (e.g., "SpatialMax", "TimeMean")
const std::string AnalysisOperator::getOperatorType() {
   return OperatorTypeName;
} // end getOperatorType

//------------------------------------------------------------------------------
// Returns the unique instance name for this operator (e.g.,
// "Temperature_SpatialMean")
const std::string AnalysisOperator::getName() {
   return InstanceName;
} // end getName

//------------------------------------------------------------------------------
// Returns the list of input field names required by this operator
const std::vector<std::string> AnalysisOperator::getInputFieldNames() {
   return InputNames;
} // end getInputFieldNames

//------------------------------------------------------------------------------
// Returns the list of output field names produced by this operator
const std::vector<std::string> AnalysisOperator::getOutputFieldNames() {
   return OutputNames;
} // end getOutputFieldNames

//------------------------------------------------------------------------------
// Checks whether the operator's output is valid for the given timestamp.
// Returns true if the operator has been computed and the LastComputed
// timestamp matches the current timestamp (cache hit), false otherwise.
// This prevents redundant computation when multiple downstream operators
// share this intermediate result.
bool AnalysisOperator::isCacheValid(const TimeInstant &TimeStamp) {
   bool IsValid = false;

   // Cache is valid if we've computed and timestamp matches
   if (FieldComputed && LastComputed == TimeStamp) {
      IsValid = true;
   }

   return IsValid;
} // end isCacheValid

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
