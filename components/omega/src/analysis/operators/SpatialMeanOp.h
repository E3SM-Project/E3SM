#ifndef OMEGA_GLOBALMEANOP_H
#define OMEGA_GLOBALMEANOP_H

//===-- analysis/operators/SpatialMeanOp.h - SpatialMeanOp ------*- C++ -*-===//
//
/// \file
/// \brief Defines the SpatialMeanOp operator for computing spatial mean
///
/// SpatialMeanOp computes the mean of a field across all owned mesh entities
/// (cells, edges, or vertices), excluding halo regions. The operator computes
/// the masked sum of field values and the sum of mask values, then divides to
/// get the mean. For 3D+ fields, the mask sum is multiplied by the product of
/// extra dimension sizes.
///
/// The operator is templated on the Kokkos array type (ArrayT) of the input
/// field, supporting 1D (horizontal only), 2D (horizontal + vertical), and 3D+
/// (extra dimensions + horizontal + vertical) fields. The output is a scalar
/// (1D array with single element) stored in a Field with dimension "Scalar".
///
/// For 1D inputs, the horizontal-only mask (k=0 column of the 2D mask) is
/// used. For 2D+ inputs, the full 2D mask (horizontal × vertical) is applied.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Reductions.h"

namespace OMEGA {

/// SpatialMeanOp computes the global spatial mean of a field across all owned
/// mesh entities and active vertical layers. The operator handles 1D, 2D, and
/// 3D+ input fields, computes masked sum of values and sum of mask values, and
/// divides to get the mean. For 3D+ fields, accounts for extra dimensions in
/// the normalization. Output is a scalar Field.
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

      // Allocate output data array (single scalar value)
      OutputData = Array1D_t<ScalarT>(OutputNames[0], 1);

      // Create scalar dimension for output Field
      I4 NDims = 1;
      std::vector<std::string> DimNames(NDims);
      DimNames[0]    = "Scalar";
      auto ScalarDim = Dimension::create(DimNames[0], 1);

      // Register output Field with metadata
      auto OutputField =
          Field::create(OutputNames[0],
                        "Spatial mean of " + InputNames[0],   // Description
                        "",                                   // Units
                        "",                                   // Standard name
                        -std::numeric_limits<ScalarT>::max(), // Min valid value
                        std::numeric_limits<ScalarT>::max(),  // Max valid value
                        -std::numeric_limits<ScalarT>::max(), // Fill value
                        NDims,                                // Rank
                        DimNames                              // Dimension names
          );

      // Attach output data array to Field
      OutputField->template attachData<Array1D_t<ScalarT>>(OutputData);

   } // end constructor

   /// Computes the spatial mean by retrieving input data, determining the
   /// appropriate mesh index space and vertical mask, constructing index ranges
   /// to exclude halo regions, computing the masked sum of values and sum of
   /// mask values, and dividing to get the mean. For 3D+ fields, scales mask
   /// sum by the product of extra dimension sizes. Updates output data,
   /// timestamp, and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

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
         ABORT_ERROR("SpatialMeanOp: Unknown index space {}", IndexSpaceName);
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

      // Compute masked sum of values and sum of mask values
      ScalarT ValSum;
      ScalarT MaskSum;

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

         ValSum = globalMaskedSum(InputData, Mask1D, Comm, &IndxRange);

         // Use cached mask sum if available, otherwise compute and cache it
         if (CachedMaskSum < 0) {
            CachedMaskSum = globalSum(Mask1D, Comm, &IndxRange);
         }
         MaskSum = CachedMaskSum;
      } else {
         // For 2D+ arrays, use full 2D mask
         ValSum = globalMaskedSum(InputData, MaskArray, Comm, &IndxRange);

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

      // Compute mean: sum of masked values / sum of mask values
      SpatialMean = ValSum / MaskSum;

      // Write result to output array
      deepCopy(OutputData, SpatialMean);

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array holding the computed spatial mean (single scalar value)
   Array1D_t<ScalarT> OutputData;

   /// Temporary storage for the computed mean value before copying to
   /// OutputData
   ScalarT SpatialMean;

   /// Contiguous 1D mask for horizontal-only operations (1D inputs).
   /// Stores k=0 column of the 2D mask. Allocated lazily on first compute
   /// to avoid LayoutStride subviews incompatible with reduction functions.
   Array1D_t<Real> Mask1D;

   /// Cached mask sum computed on first pass and reused for subsequent calls.
   /// The mask is constant in time, so this optimization avoids redundant
   /// global reduction operations. Initialized to -1 to indicate not yet
   /// computed.
   ScalarT CachedMaskSum{static_cast<ScalarT>(-1)};

}; // end class SpatialMeanOp

} // end namespace OMEGA

#endif
