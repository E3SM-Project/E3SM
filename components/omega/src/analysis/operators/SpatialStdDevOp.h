#ifndef OMEGA_STDDEV_H
#define OMEGA_STDDEV_H

//===-- analysis/operators/SpatialStdDevOp.h - SpatialStdDevOp --*- C++ -*-===//
//
/// \file
/// \brief Defines the SpatialStdDevOp operator for computing standard deviation
///
/// SpatialStdDevOp computes the weighted standard deviation of a field across
/// all owned mesh entities (cells, edges, or vertices) and active layers,
/// excluding halo regions, with the same area or mass weights as
/// SpatialMeanOp (see SpatialWeights). The operator requires the spatial mean
/// as input (computed by SpatialMeanOp), calculates the weighted sum of
/// squared deviations from the mean and the sum of weights, divides to get
/// the variance, and takes the square root to get the standard deviation.
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

/// SpatialStdDevOp computes the global spatial standard deviation of a field
/// across all owned mesh entities and active vertical layers, weighted by
/// area for a horizontal field and by mass for a layered one. The operator
/// requires the spatial mean as input, performs a weighted sum of squared
/// deviations, and takes the square root of the variance. Output is a scalar
/// Field.
template <typename ArrayT> class SpatialStdDevOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Constructs a SpatialStdDevOp operator. Declares two inputs: the field
   /// itself and its spatial mean (computed by SpatialMeanOp). Creates output
   /// Field as scalar (1D array with single element), allocates the output
   /// data array, and registers the output Field in the Field registry. The
   /// output Field name is constructed as InputName + "_SpatialStdDev".
   SpatialStdDevOp(const std::vector<std::string>
                       &UpstreamNames, ///< [in] input field names
                   Config Options      ///< [in] operator config
                   )
       : AnalysisOperator("SpatialStdDev") {

      // Declare two inputs: field and its spatial mean
      // Mean must be computed first (e.g., by SpatialMeanOp in the chain)
      InputNames = {UpstreamNames[0], UpstreamNames[0] + "_SpatialMean"};

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_SpatialStdDev";
      OutputNames                 = {OutputFieldName};
      InstanceName                = OutputFieldName;

      // Allocate output data array (single scalar value)
      OutputData = Array1DReal(OutputNames[0], 1);

      // Create scalar dimension for output Field
      I4 NDims = 1;
      std::vector<std::string> DimNames(NDims);
      DimNames[0]    = "Scalar";
      auto ScalarDim = Dimension::create(DimNames[0], 1);

      // Register output Field with metadata. A standard deviation has the
      // units and standard name of the field it reduces; its cell_methods
      // record the reduction after any the input already carries.
      auto Meta = inheritMetadata(
          InputNames[0],
          spatialCellMethod(InputNames[0], "standard_deviation"));
      auto OutputField = createOutputField(
          OutputNames[0],
          "Standard deviation of " + InputNames[0], // Description
          Meta,                                     // CF metadata
          static_cast<Real>(0),                     // Min valid value
          std::numeric_limits<Real>::max(),         // Max valid value
          NDims,                                    // Rank
          DimNames                                  // Dimension names
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

   /// Computes the spatial standard deviation by retrieving the input data
   /// and spatial mean, updating the weights from the current state,
   /// computing the weighted sum of squared deviations from the mean over
   /// the owned entities and active layers, dividing by the sum of weights
   /// to get the variance, and taking the square root. Updates output data,
   /// timestamp, and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Product of the extra dimension sizes over which the weights repeat
      Real Replication = 1;
      for (I4 I = 0; I < static_cast<I4>(InputData.rank) - 2; ++I)
         Replication *= InputData.extent(I);

      // Retrieve spatial mean value computed by upstream SpatialMeanOp
      auto MeanField = Field::get(InputNames[1]);
      auto MeanVal   = MeanField->template getDataArray<Array1DReal>();
      auto MeanHost  = createHostMirrorCopy(MeanVal);
      Real Mean      = MeanHost(0);

      // Weights depend on the state for layered fields, so refresh them
      Weights.update();

      // Variance: weighted sum of squared deviations over sum of weights
      Real DevSum   = Weights.weightedSquaredDeviation(InputData, Mean);
      Real Variance = DevSum / (Weights.weightSum() * Replication);

      // Compute standard deviation: square root of variance
      StdDev = std::sqrt(Variance);

      // Write result to output array
      deepCopy(OutputData, StdDev);

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array holding the computed standard deviation (single scalar
   /// value). Always Real type regardless of input type
   Array1DReal OutputData;

   /// Temporary storage for the computed standard deviation before copying to
   /// OutputData
   Real StdDev;

   /// Area or mass weights of the input field's entities and layers
   SpatialWeights Weights;

}; // end class SpatialStdDevOp

} // end namespace OMEGA

#endif
