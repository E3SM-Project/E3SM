#ifndef OMEGA_STDDEV_H
#define OMEGA_STDDEV_H

//===-- analysis/operators/SpatialStdDevOp.h - SpatialStdDevOp --*- C++ -*-===//
//
/// \file
/// \brief Defines the SpatialStdDevOp operator for computing standard deviation
///
/// SpatialStdDevOp computes the standard deviation of a field across all
/// owned mesh entities (cells, edges, or vertices), excluding halo regions.
/// The operator requires the spatial mean as input (computed by
/// SpatialMeanOp), computes squared differences from the mean in a work array,
/// calculates the masked sum of squared differences and sum of mask values,
/// divides to get variance, and takes the square root to get standard
/// deviation.
///
/// The operator is templated on the Kokkos array type (ArrayT) of the input
/// field, supporting 1D (horizontal only), 2D (horizontal + vertical), and 3D+
/// (extra dimensions + horizontal + vertical) fields. The output is a scalar
/// (1D array with single element) stored in a Field with dimension "Scalar".
///
/// For 1D inputs, the horizontal-only mask (k=0 layer of the full mask by
/// default) is used. For 2D+ inputs, the full 2D mask (horizontal × vertical)
/// is applied. Unlike simpler operators, SpatialStdDevOp allocates a work
/// array matching the input field layout to store squared differences before
/// reduction.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Reductions.h"

namespace OMEGA {

/// SpatialStdDevOp computes the global spatial standard deviation of a field
/// across all owned mesh entities and active vertical layers. The operator
/// requires spatial mean as input, computes squared differences in a work
/// array, performs masked sum reduction, and takes square root of the
/// variance. Output is a scalar Field.
template <typename ArrayT> class SpatialStdDevOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Constructs a SpatialStdDevOp operator. Declares two inputs: the field
   /// itself and its spatial mean (computed by SpatialMeanOp). Creates output
   /// Field as scalar (1D array with single element), allocates output data
   /// array and work array matching input layout, and registers the output
   /// Field in the Field registry. The output Field name is constructed as
   /// InputName + "_SpatialStdDev".
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

      // Register output Field with metadata
      auto OutputField =
          Field::create(OutputNames[0],
                        "Standard deviation of " + InputNames[0], // Description
                        "",                                       // Units
                        "",                                // Standard name
                        static_cast<Real>(0),              // Min valid value
                        std::numeric_limits<Real>::max(),  // Max valid value
                        -std::numeric_limits<Real>::max(), // Fill value
                        NDims,                             // Rank
                        DimNames                           // Dimension names
          );

      // Attach output data array to Field
      OutputField->template attachData<Array1DReal>(OutputData);

      // Allocate work array matching input field layout but always using Real
      // type Used to store squared differences: (x - mean)^2 Must be Real type
      // to preserve precision when computing variance
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Create Real-type array with same shape as input
      if constexpr (ArrayT::rank == 1) {
         WorkArray = Array1D_t<Real>(OutputNames[0] + "_work_array",
                                     InputData.extent(0));
      } else if constexpr (ArrayT::rank == 2) {
         WorkArray = Array2D_t<Real>(OutputNames[0] + "_work_array",
                                     InputData.extent(0), InputData.extent(1));
      } else if constexpr (ArrayT::rank == 3) {
         WorkArray = Array3D_t<Real>(OutputNames[0] + "_work_array",
                                     InputData.extent(0), InputData.extent(1),
                                     InputData.extent(2));
      }

   } // end constructor

   /// Computes the spatial standard deviation by retrieving input data and
   /// spatial mean, determining the appropriate mesh index space and vertical
   /// mask, constructing index ranges to exclude halo regions, filling work
   /// array with squared differences from mean, computing masked sum of squared
   /// differences and sum of mask values, dividing to get variance, and taking
   /// square root. Updates output data, timestamp, and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Create local scope reference to work array for kernel capture
      OMEGA_SCOPE(LocWorkArray, WorkArray);

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Get dimension names to determine array structure
      std::vector<std::string> InputDimNames;
      InputField->getDimNames(InputDimNames);
      I4 NDims = InputDimNames.size();

      // Determine mesh index space (cells/edges/vertices) from dimension name
      // For 1D: dimension is horizontal
      // For 2D+: second-to-last dimension is horizontal
      std::string IndexSpaceName = InputDimNames[std::max(0, NDims - 2)];

      // Get appropriate mask and owned entity count for this index space
      Array2DReal MaskArray;
      I4 NOwned      = 0;
      I4 NVertLayers = VCoord->NVertLayers;

      if (IndexSpaceName == "NCells") {
         MaskArray = VCoord->CellMask;
         NOwned    = Mesh->NCellsOwned;
      } else if (IndexSpaceName == "NEdges") {
         MaskArray = VCoord->EdgeMask;
         NOwned    = Mesh->NEdgesOwned;
      } else if (IndexSpaceName == "NVertices") {
         MaskArray = VCoord->VertexMask;
         NOwned    = Mesh->NVerticesOwned;
      } else {
         ABORT_ERROR("SpatialStdDevOp: Unknown index space {}", IndexSpaceName);
      }

      // Construct index range for input data to exclude halo cells and inactive
      // layers Format: [dim0_start, dim0_end, dim1_start, dim1_end, ...]
      std::vector<I4> IndxRange;

      if (NDims == 1) {
         // 1D array: horizontal dimension only
         IndxRange = {0, NOwned - 1};
      } else if (NDims == 2) {
         // 2D array: (horizontal, vertical)
         IndxRange = {0, NOwned - 1, 0, NVertLayers - 1};
      } else {
         // 3D+ array: (extra dims..., horizontal, vertical)
         IndxRange.resize(2 * NDims);

         // Extra dimensions: include full extent
         for (I4 I = 0; I < NDims - 2; ++I) {
            IndxRange[2 * I]     = 0;
            IndxRange[2 * I + 1] = InputData.extent(I) - 1;
         }

         // Horizontal dimension (second to last): exclude halo
         IndxRange[2 * (NDims - 2)]     = 0;
         IndxRange[2 * (NDims - 2) + 1] = NOwned - 1;

         // Vertical dimension (last): all layers
         IndxRange[2 * (NDims - 1)]     = 0;
         IndxRange[2 * (NDims - 1) + 1] = NVertLayers - 1;
      }

      // Index range for mask array (always 2D: horizontal × vertical)
      std::vector<I4> MaskIndxRange = {0, NOwned - 1, 0, NVertLayers - 1};

      // Retrieve spatial mean value computed by upstream SpatialMeanOp
      auto MeanField = Field::get(InputNames[1]);
      auto MeanVal   = MeanField->template getDataArray<Array1DReal>();

      // Fill work array with squared differences: (x - mean)^2
      // Mask will be applied later during reduction
      I4 NSize           = static_cast<I4>(InputData.size());
      const int Arr1Rank = InputData.rank;
      parallelFor(
          {NSize}, KOKKOS_LAMBDA(const int FlatIdx) {
             // Compute horizontal and vertical indices from flat index
             // (used for debugging/validation, though not needed for this
             // computation)
             int HorizIdx = 0;
             int VertIdx  = 0;

             if (Arr1Rank == 1) {
                HorizIdx = FlatIdx;
             } else if (Arr1Rank == 2) {
                HorizIdx = FlatIdx / InputData.extent(1);
                VertIdx  = FlatIdx % InputData.extent(1);
             } else {
                int IdxLastTwo = FlatIdx % (InputData.extent(Arr1Rank - 2) *
                                            InputData.extent(Arr1Rank - 1));
                HorizIdx       = IdxLastTwo / InputData.extent(Arr1Rank - 1);
                VertIdx        = IdxLastTwo % InputData.extent(Arr1Rank - 1);
             }

             // Compute squared difference and store in work array
             // Cast input to Real for computation, WorkArray is Real type
             auto Diff =
                 static_cast<Real>(InputData.data()[FlatIdx]) - MeanVal(0);
             LocWorkArray.data()[FlatIdx] = Diff * Diff;
          });

      // Compute masked sum of squared differences and sum of mask values
      // Cast WorkSum to Real immediately to ensure proper precision
      Real WorkSum;
      Real MaskSum;

      if (NDims == 1) {
         // For 1D arrays, use horizontal-only mask (k=0 column of 2D mask)
         // Copy to contiguous 1D array to avoid LayoutStride incompatibility
         if (Mask1D.size() == 0)
            Mask1D = Array1D_t<Real>("Mask1D", MaskArray.extent(0));

         auto LocalMaskArray = MaskArray;
         auto LocalMask1D    = Mask1D;
         parallelFor(
             {static_cast<I4>(MaskArray.extent(0))},
             KOKKOS_LAMBDA(int I) { LocalMask1D(I) = LocalMaskArray(I, 0); });

         WorkSum = static_cast<Real>(
             globalMaskedSum(WorkArray, Mask1D, Comm, &IndxRange));

         // Use cached mask sum if available, otherwise compute and cache it
         if (CachedMaskSum < 0) {
            CachedMaskSum = globalSum(Mask1D, Comm, &IndxRange);
         }
         MaskSum = CachedMaskSum;
      } else {
         // For 2D+ arrays, use full 2D mask
         WorkSum = static_cast<Real>(
             globalMaskedSum(WorkArray, MaskArray, Comm, &IndxRange));

         // Use cached mask sum if available, otherwise compute and cache it
         if (CachedMaskSum < 0) {
            CachedMaskSum = globalSum(MaskArray, Comm, &MaskIndxRange);

            // For 3D+ arrays, scale mask sum by product of extra dimension
            // sizes. This accounts for replication of the 2D mask across extra
            // dimensions
            if (NDims > 2) {
               I4 ExtraDimSize = 1;
               for (I4 I = 0; I < NDims - 2; ++I) {
                  ExtraDimSize *= InputData.extent(I);
               }
               CachedMaskSum *= ExtraDimSize;
            }
         }
         MaskSum = CachedMaskSum;
      }

      // Compute variance: sum of squared diffs / sum of mask
      Real Variance = WorkSum / MaskSum;

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

   /// Work array matching input field layout but always Real type to preserve
   /// precision Used to store squared differences (x - mean)^2 before masked
   /// reduction
   typename std::conditional<
       ArrayT::rank == 1, Array1D_t<Real>,
       typename std::conditional<ArrayT::rank == 2, Array2D_t<Real>,
                                 Array3D_t<Real>>::type>::type WorkArray;

   /// Temporary storage for the computed standard deviation before copying to
   /// OutputData
   Real StdDev;

   /// Contiguous 1D mask for horizontal-only operations (1D inputs).
   /// Stores k=0 column of the 2D mask. Allocated lazily on first compute
   /// to avoid LayoutStride subviews incompatible with reduction functions.
   Array1D_t<Real> Mask1D;

   /// Cached mask sum computed on first pass and reused for subsequent calls.
   /// The mask is constant in time, so this optimization avoids redundant
   /// global reduction operations. Initialized to -1.0 to indicate not yet
   /// computed.
   Real CachedMaskSum{static_cast<Real>(-1.0)};

}; // end class SpatialStdDevOp

} // end namespace OMEGA

#endif
