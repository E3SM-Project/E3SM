#ifndef OMEGA_GLOBALMEANOP_H
#define OMEGA_GLOBALMEANOP_H

//===-- analysis/operators/SpatialMeanOp.h - SpatialMeanOp ------*- C++ -*-===//
//
/// \file
/// \brief Defines the SpatialMeanOp operator for computing spatial mean
///
/// SpatialMeanOp computes the weighted mean of a field across all owned mesh
/// entities (cells, edges, or vertices) and active layers, excluding halo
/// regions. A horizontal field is weighted by the area of each entity and a
/// field with a vertical dimension by mass (area times pseudo-thickness), so
/// the mean is independent of the mesh resolution; see SpatialWeights. The
/// operator computes the weighted sum of field values and the sum of weights,
/// then divides to get the mean. For 3D+ fields, the weight sum is multiplied
/// by the product of extra dimension sizes.
///
/// The operator is templated on the Kokkos array type (ArrayT) of the input
/// field, supporting 1D (horizontal only), 2D (horizontal + vertical), and 3D+
/// (extra dimensions + horizontal + vertical) fields. The output is a scalar
/// (1D array with single element) stored in a Field with dimension "Scalar".
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Reductions.h"
#include "SpatialWeights.h"

namespace OMEGA {

/// SpatialMeanOp computes the global spatial mean of a field across all owned
/// mesh entities and active vertical layers, weighted by area for a
/// horizontal field and by mass for a layered one. The operator handles 1D,
/// 2D, and 3D+ input fields, computes the weighted sum of values and the sum
/// of weights, and divides to get the mean. For 3D+ fields, accounts for
/// extra dimensions in the normalization. Output is a scalar Field.
template <typename ArrayT> class SpatialMeanOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Constructs a SpatialMeanOp operator. Creates output Field as scalar
   /// (1D array with single element), allocates output data array, and
   /// registers the output Field in the Field registry. The output Field
   /// name is constructed as InputName + "_SpatialMean".
   SpatialMeanOp(const std::vector<std::string>
                     &UpstreamNames, ///< [in] input field names
                 Config Options      ///< [in] operator config
                 )
       : AnalysisOperator("SpatialMean") {

      // Store input field names
      InputNames = UpstreamNames;

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_SpatialMean";
      OutputNames                 = {OutputFieldName};
      InstanceName                = OutputFieldName;

      // Allocate output data array (single scalar value, always Real type)
      OutputData = Array1DReal(OutputNames[0], 1);

      // Create scalar dimension for output Field
      I4 NDims = 1;
      std::vector<std::string> DimNames(NDims);
      DimNames[0]    = "Scalar";
      auto ScalarDim = Dimension::create(DimNames[0], 1);

      // Register output Field with metadata. A mean has the units and
      // standard name of the field it reduces; its cell_methods record the
      // reduction after any the input already carries.
      auto Meta = inheritMetadata(InputNames[0],
                                  spatialCellMethod(InputNames[0], "mean"));
      auto OutputField =
          createOutputField(OutputNames[0],
                            "Spatial mean of " + InputNames[0], // Description
                            Meta,                               // CF metadata
                            -std::numeric_limits<Real>::max(),  // Min valid
                            std::numeric_limits<Real>::max(),   // Max valid
                            NDims,                              // Rank
                            DimNames                            // Dim names
          );

      // Attach output data array to Field
      OutputField->template attachData<Array1DReal>(OutputData);

   } // end constructor

   /// Stores the mesh, vertical coordinate and communicator and sets up the
   /// weights for the input field
   void initialize(const MachEnv *Env, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {
      AnalysisOperator::initialize(Env, InMesh, InVCoord, Options);
      Weights.init(InputNames[0], Mesh, VCoord, Comm);
   }

   /// Computes the spatial mean by retrieving input data, updating the
   /// weights from the current state, computing the weighted sum of values
   /// over the owned entities and active layers and dividing by the sum of
   /// weights. For 3D+ fields, scales the weight sum by the product of extra
   /// dimension sizes. Updates output data, timestamp, and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Product of the extra dimension sizes over which the weights repeat
      Real Replication = 1;
      for (I4 I = 0; I < static_cast<I4>(InputData.rank) - 2; ++I)
         Replication *= InputData.extent(I);

      // Weights depend on the state for layered fields, so refresh them
      Weights.update();

      // Mean: weighted sum of values over the sum of weights
      Real ValSum = Weights.weightedSum(InputData);
      SpatialMean = ValSum / (Weights.weightSum() * Replication);

      // Write result to output array
      deepCopy(OutputData, SpatialMean);

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array holding the computed spatial mean (single scalar value)
   /// Always Real type regardless of input type
   Array1DReal OutputData;

   /// Temporary storage for the computed mean value before copying to
   /// OutputData
   Real SpatialMean;

   /// Area or mass weights of the input field's entities and layers
   SpatialWeights Weights;

}; // end class SpatialMeanOp

} // end namespace OMEGA

#endif
