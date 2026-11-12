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
// Parses a frequency string like "1day" or "6hour" into numeric and unit
// components. Returns a vector with [0] = numeric part, [1] = units part.
// Appends 's' to units if not already present for consistency with
// TimeInterval.
std::vector<std::string> parseFreqStr(const std::string &FreqStr) {

   std::string DigitStr;
   std::string UnitsStr;

   // Find where digits end and units begin
   size_t Pos = FreqStr.find_first_not_of("0123456789");
   if (Pos != std::string::npos) {
      DigitStr = FreqStr.substr(0, Pos);
      UnitsStr = FreqStr.substr(Pos);
   }

   // Validate that we found both components
   if (DigitStr.empty() || UnitsStr.empty()) {
      ABORT_ERROR("Analysis: Invalid frequency string found in Config: {}",
                  FreqStr);
   }

   // Ensure units end with 's' for consistency with TimeInterval
   if (UnitsStr.back() != 's') {
      UnitsStr += 's';
   }

   return {DigitStr, UnitsStr};

} // end parseFreqStr

//------------------------------------------------------------------------------
// Default constructor
AnalysisOperator::AnalysisOperator() = default;

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

   OMEGA_REQUIRE(Env, "Null MachEnv pointer in AnalysisOperator::initialize");
   OMEGA_REQUIRE(InMesh,
                 "Null HorzMesh pointer in AnalysisOperator::initialize");
   OMEGA_REQUIRE(InVCoord,
                 "Null VertCoord pointer in AnalysisOperator::initialize");

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
// Returns a string attribute of a Field, or an empty string if the Field has
// no such attribute. Analysis output is derived from fields that may lack
// units or a standard name, so absence is the empty string, not an error.
std::string AnalysisOperator::getFieldAttribute(const std::string &FieldName,
                                                const std::string &AttName) {

   auto FieldPtr = Field::get(FieldName);
   OMEGA_REQUIRE(FieldPtr, "AnalysisOperator: field {} not found", FieldName);

   std::string Value;
   if (FieldPtr->hasMetadata(AttName)) {
      Error Err = FieldPtr->getMetadata(AttName, Value);
      CHECK_ERROR_ABORT(Err,
                        "AnalysisOperator: attribute {} of field {} is "
                        "not a string",
                        AttName, FieldName);
   }
   return Value;
} // end getFieldAttribute

//------------------------------------------------------------------------------
// Collects the units, standard name and cell methods an output inherits from
// its input. A new cell method is appended to the input's cell_methods so the
// attribute lists the reductions in the order they were applied, as CF
// requires (e.g. a time mean of a spatial mean is "area: mean time: mean").
AnalysisOperator::InheritedMetadata
AnalysisOperator::inheritMetadata(const std::string &InputName,
                                  const std::string &CellMethod) {

   InheritedMetadata Meta;
   Meta.Units       = getFieldAttribute(InputName, "units");
   Meta.StdName     = getFieldAttribute(InputName, "standard_name");
   Meta.CellMethods = getFieldAttribute(InputName, "cell_methods");

   if (!CellMethod.empty()) {
      if (!Meta.CellMethods.empty())
         Meta.CellMethods += " ";
      Meta.CellMethods += CellMethod;
   }
   return Meta;
} // end inheritMetadata

//------------------------------------------------------------------------------
// Returns the cell method for a reduction over all owned mesh entities and
// layers of an input Field. The spatial operators treat a 1D field as
// horizontal only and the last dimension of any higher-rank field as
// vertical, so the method covers "area" alone or "area" and "depth". Both
// are names CF allows in cell_methods without a matching coordinate.
std::string AnalysisOperator::spatialCellMethod(const std::string &InputName,
                                                const std::string &Method) {

   auto FieldPtr = Field::get(InputName);
   OMEGA_REQUIRE(FieldPtr, "AnalysisOperator: field {} not found", InputName);

   if (FieldPtr->getNumDims() >= 2)
      return "area: depth: " + Method;
   return "area: " + Method;
} // end spatialCellMethod

//------------------------------------------------------------------------------
// Creates an output Field carrying the inherited CF metadata. The cell
// methods are stored as a std::string so the metadata can be retrieved and
// written as a string attribute.
std::shared_ptr<Field> AnalysisOperator::createOutputField(
    const std::string &OutputName, const std::string &Description,
    const InheritedMetadata &Meta, const std::any ValidMin,
    const std::any ValidMax, const int NumDims,
    const std::vector<std::string> &DimNames) {

   auto OutputField =
       Field::create(OutputName, Description, Meta.Units, Meta.StdName,
                     ValidMin, ValidMax, NumDims, DimNames);

   if (!Meta.CellMethods.empty())
      OutputField->addMetadata("cell_methods", Meta.CellMethods);

   return OutputField;
} // end createOutputField

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

//------------------------------------------------------------------------------
// If Options contains a non-empty "OutputName" key, overrides OutputNames[0]
// and InstanceName with that value. Derived class constructors call this after
// setting their default OutputNames to support caller-assigned short names.
void AnalysisOperator::applyOutputNameOverride(Config &Options) {
   std::string Override;
   if (OutputNames.empty()) {
      return;
   }
   Error Err = Options.get("OutputName", Override);
   if (Err.isSuccess() && !Override.empty()) {
      OutputNames[0] = Override;
      InstanceName   = Override;
   }
} // end applyOutputNameOverride

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
