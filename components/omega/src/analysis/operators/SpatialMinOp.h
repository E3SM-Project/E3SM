#ifndef OMEGA_GLOBALMINOP_H
#define OMEGA_GLOBALMINOP_H

//===-- analysis/operators/SpatialMinOp.h - SpatialMinOp --------*- C++ -*-===//
//
/// \file
/// \brief Defines the SpatialMinOp operator for computing spatial minimum
///
/// SpatialMinOp computes the global minimum value of a field across all owned
/// mesh entities (cells, edges, or vertices), excluding halo regions. The
/// operator uses MPI reduction to find the minimum across all MPI ranks,
/// respecting the vertical mask for active/inactive layers.
///
/// The operator is templated on the Kokkos array type (ArrayT) of the input
/// field, supporting 1D (horizontal only), 2D (horizontal + vertical), and 3D+
/// (extra dimensions + horizontal + vertical) fields. The output is a scalar
/// (1D array with single element) stored in a Field with dimension "Scalar".
///
/// For 1D inputs, the horizontal-only mask (k=0 column of the 2D mask) is used.
/// For 2D+ inputs, the full 2D mask (horizontal × vertical) is applied. Index
/// ranges are constructed to exclude halo cells and include only owned entities
/// and active vertical layers.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Reductions.h"

namespace OMEGA {

/// SpatialMinOp computes the global spatial minimum of a field across all owned
/// mesh entities and active vertical layers. The operator handles 1D, 2D, and
/// 3D+ input fields, constructs appropriate index ranges to exclude halos,
/// applies vertical masking, and performs MPI reduction to find the minimum
/// across all ranks. Output is a scalar Field.
template <typename ArrayT> class SpatialMinOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Constructs a SpatialMinOp operator. Creates output Field as scalar
   /// (1D array with single element), allocates output data array, and
   /// registers the output Field in the Field registry. The output Field
   /// name is constructed as InputName + "_SpatialMin".
   SpatialMinOp(const std::vector<std::string>
                    &UpstreamNames, ///< [in] input field names
                Config Options      ///< [in] operator config
                )
       : AnalysisOperator("SpatialMin") {

      // Store input field names
      InputNames = UpstreamNames;

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_SpatialMin";
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
                        "Spatial minimum of " + InputNames[0], // Description
                        "",                                    // Units
                        "",                                    // Standard name
                        -std::numeric_limits<ScalarT>::max(), // Min valid value
                        std::numeric_limits<ScalarT>::max(),  // Max valid value
                        NDims,                                // Rank
                        DimNames                              // Dimension names
          );

      // Attach output data array to Field
      OutputField->template attachData<Array1D_t<ScalarT>>(OutputData);

   } // end constructor

   /// Computes the spatial minimum by retrieving input data, determining the
   /// appropriate mesh index space (cells/edges/vertices) and vertical mask,
   /// constructing index ranges to exclude halo regions, and calling
   /// globalMaskedMin() with MPI reduction. For 1D inputs, uses horizontal-only
   /// mask; for 2D+ inputs, uses full 2D mask. Updates output data, timestamp,
   /// and computed flag.
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
         ABORT_ERROR("SpatialMinOp: Unknown index space {}", IndexSpaceName);
      }

      // Construct index range to exclude halo cells and inactive layers
      // Format: [dim0_start, dim0_end, dim1_start, dim1_end, ...]
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

      // Compute global masked minimum
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

         SpatialMin = globalMaskedMin(InputData, Mask1D, Comm, &IndxRange);
      } else {
         // For 2D+ arrays, use full 2D mask
         SpatialMin = globalMaskedMin(InputData, MaskArray, Comm, &IndxRange);
      }

      // Write result to output array
      deepCopy(OutputData, SpatialMin);

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array holding the computed spatial minimum (single scalar
   /// value)
   Array1D_t<ScalarT> OutputData;

   /// Temporary storage for the computed minimum value before copying to
   /// OutputData
   ScalarT SpatialMin;

   /// Contiguous 1D mask for horizontal-only operations (1D inputs).
   /// Stores k=0 column of the 2D mask. Allocated lazily on first compute
   /// to avoid LayoutStride subviews incompatible with reduction functions.
   Array1D_t<Real> Mask1D;

}; // end class SpatialMinOp

} // namespace OMEGA
#endif
