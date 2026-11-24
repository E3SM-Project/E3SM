#ifndef OMEGA_SCALARMULTIPLYOP_H
#define OMEGA_SCALARMULTIPLYOP_H

//===-- analysis/operators/ScalarMultiplyOp.h -------------------*- C++ -*-===//
//
/// \file
/// \brief Defines the ScalarMultiplyOp operator for scalar multiplication
///
/// ScalarMultiplyOp multiplies all elements of an input field by a
/// user-specified scalar value. The operator is templated on the Kokkos array
/// type (ArrayT) of the input field, supporting 1D, 2D, and 3D+ fields. The
/// output has the same shape and dimensions as the input.
///
/// This operator is particularly useful for unit conversions (e.g., m³/s to
/// Sverdrups for ocean transport) and for scaling derived quantities.
///
/// Configuration:
/// - Scalar: The multiplicative factor (specified in parentheses in operator
/// chain)
/// - InPlace: If true, modifies input array directly (default: false)
///
/// Example usage in operator chain:
/// \code
///   Field_ScalarMultiply(1.0e-6)
/// \endcode
/// where ScalarMultiply(1.0e-6) multiplies the Field by 1.0e-6.
/// This is useful for unit conversions (e.g., m³/s to Sv).
/// When InPlace=true, no additional memory is allocated; the input array is
/// modified directly and the output Field references the same data.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"

#include <cmath>

namespace OMEGA {

/// ScalarMultiplyOp multiplies all elements of an input field by a scalar
/// value. The operator handles arrays of any rank (1D, 2D, 3D+) and preserves
/// the input dimensions in the output. The scalar multiplication is performed
/// element-wise using Kokkos parallel operations. Output type matches input
/// type unless input is integral, in which case output is Real.
template <typename ArrayT> class ScalarMultiplyOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type - same as input array type
   using OutputArrayT = ArrayT;

   /// Constructs a ScalarMultiplyOp operator. Creates output Field with the
   /// same dimensions as the input field, allocates output data array, and
   /// registers the output Field in the Field registry. The output Field name
   /// is constructed as InputName + "_ScalarMultiply". Reads the scalar
   /// multiplicative factor from the Options config (which is set by the parser
   /// from the parenthesized argument in the operator chain).
   ScalarMultiplyOp(const std::vector<std::string>
                        &UpstreamNames, ///< [in] input field names
                    Config Options      ///< [in] operator config
                    )
       : AnalysisOperator("ScalarMultiply") {

      // Store input field names
      InputNames = UpstreamNames;

      // Read required scalar parameter from configuration.
      // Read as string first to preserve the original representation for field
      // naming (std::to_string loses precision, e.g. "1.0e-6" -> "0.000001").
      // The units are scaled with the scalar in double precision, so that in
      // a single-precision build 1e-6 still gives units of 1e+06 m3 s-1.
      std::string ScalarStr;
      double ScalarR8;
      Error Err = Options.get("Scalar", ScalarStr);
      if (Err.isSuccess()) {
         ScalarR8 = std::stod(ScalarStr);
         Scalar   = static_cast<Real>(ScalarR8);
      } else {
         // Fallback: try reading directly as Real
         Err = Options.get("Scalar", Scalar);
         if (Err.isFail()) {
            ABORT_ERROR(
                "ScalarMultiplyOp: Required parameter 'Scalar' not found "
                "in configuration");
         }
         ScalarStr = std::to_string(Scalar);
         ScalarR8  = static_cast<double>(Scalar);
      }

      // Read optional InPlace parameter (default: false)
      Err = Options.get("InPlace", InPlace);
      if (Err.isFail()) {
         InPlace = false; // Default to allocating new array
      }

      // Retrieve input Field to get dimensions and metadata
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Get dimension info
      auto NDims = InputField->getNumDims();
      std::vector<std::string> DimNames;
      InputField->getDimNames(DimNames);

      // Construct output field name and set instance name.
      // Use the original scalar string to match the chain string
      // representation.
      std::string OutputFieldName =
          InputNames[0] + "_ScalarMultiply(" + ScalarStr + ")";
      OutputNames  = {OutputFieldName};
      InstanceName = OutputFieldName;

      // Get input metadata
      std::string InputDescr;
      ScalarT InputValidMin, InputValidMax;
      InputField->getMetadata("long_name", InputDescr);
      InputField->getMetadata("valid_min", InputValidMin);
      InputField->getMetadata("valid_max", InputValidMax);

      // Multiplying the values by the scalar scales the units by its
      // reciprocal, so a transport in m3 s-1 multiplied by 1e-6 is in
      // 1e+06 m3 s-1 (Sverdrups). A negative scalar also reverses the sign of
      // the quantity, which then no longer matches the input standard name. A
      // zero scalar gives zeros, which are valid in any units.
      auto Meta = inheritMetadata(InputNames[0]);
      if (ScalarR8 != 0)
         Meta.Units = CFUnits::scale(Meta.Units, 1.0 / std::abs(ScalarR8));
      if (ScalarR8 < 0)
         Meta.StdName = "";

      // Create output Field with same dimensions as input
      auto OutputField =
          createOutputField(OutputNames[0],
                            InputDescr + " multiplied by " +
                                std::to_string(Scalar), // Description
                            Meta,                       // CF metadata
                            InputValidMin,              // Min valid
                            InputValidMax,              // Max valid
                            NDims,                      // Rank
                            DimNames                    // Dimension names
          );

      // Store array size for parallel iteration
      ArraySize = static_cast<I4>(InputData.size());

      // Handle InPlace vs new allocation
      if (InPlace) {
         // InPlace: output Field references the same array as input
         // This modifies the input data directly
         OutputData = InputData;
      } else {
         // Allocate new output data array matching input layout
         OutputData = OutputArrayT(OutputNames[0] + "_out", InputData.layout());
      }

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

      // Propagate regional mask from input to output if present
      if (InputField->hasRegionalMask()) {
         OutputField->setRegionalMask(InputField->getRegionalMask());
      }

   } // end constructor

   /// Computes the scalar multiplication by retrieving input data and
   /// performing element-wise multiplication using Kokkos parallel_for. The
   /// operation uses flat indexing to handle arrays of any rank efficiently.
   /// Updates output data, timestamp, and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Create local scope reference to output array for kernel capture
      OMEGA_SCOPE(LocOutputData, OutputData);

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Perform element-wise multiplication using flat indexing
      // This works for arrays of any rank
      auto LocalScalar = Scalar;
      parallelFor(
          {ArraySize}, KOKKOS_LAMBDA(const int FlatIdx) {
             LocOutputData.data()[FlatIdx] =
                 InputData.data()[FlatIdx] * LocalScalar;
          });

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array holding the scaled field values
   OutputArrayT OutputData;

   /// The scalar multiplier read from configuration
   Real Scalar;

   /// Whether to modify input array in-place (true) or allocate new array
   /// (false)
   bool InPlace;

   /// Total size of the array for flat indexing
   I4 ArraySize;

}; // end class ScalarMultiplyOp

} // end namespace OMEGA

#endif
