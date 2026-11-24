#ifndef OMEGA_HORZMEANOP_H
#define OMEGA_HORZMEANOP_H

//===-- analysis/operators/HorzMeanOp.h - HorzMeanOp -----------*- C++ -*-===//
//
/// \file
/// \brief Defines the HorzMeanOp operator for area-weighted horizontal means
///
/// HorzMeanOp reduces a 2D field (horizontal x vertical) over the horizontal
/// dimension, producing a 1D result with one value per vertical level or
/// interface. Unlike SpatialMeanOp (which collapses all spatial dimensions to
/// a single scalar), HorzMeanOp preserves the vertical dimension. This is used
/// to build a representative per-interface depth coordinate for MOC output.
///
/// For each vertical index k the operator computes the weighted mean over the
/// owned horizontal locations:
///
///   Output(k) = sum_i ( w_i * P_i * z_ik ) / sum_i ( w_i * P_i )
///
/// where w_i is an optional horizontal weight (default AreaCell), and P_i is
/// a presence mask that is 1 when the horizontal point is inside the region of
/// interest.
///
/// The mid-layer vs interface staggering is auto-detected in initialize() by
/// comparing the input field vertical extent to VertCoord NVertLayers, so the
/// operator handles both mid-layer- and interface-dimensioned inputs.
///
/// Configuration:
/// - WeightField: Name of a 1D horizontal weight field (optional). Defaults to
///                the mesh AreaCell array when omitted (area weighting).
///
/// Example usage in operator chains:
/// \code
///   GeomZInterface_HorzMean
///   GeomZInterface_ExtractRegion(Atlantic)_HorzMean
/// \endcode
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Reductions.h"

#include <limits>
#include <vector>

namespace OMEGA {

/// HorzMeanOp computes an (optionally area-weighted) horizontal mean of a 2D
/// field, preserving the vertical dimension. The output is a 1D Real Field with
/// one value per vertical level/interface of the input field.
template <typename ArrayT> class HorzMeanOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Constructs a HorzMeanOp operator. Reads the optional WeightField config
   /// parameter, builds the 1D output Field over the input's vertical
   /// dimension, allocates the output data array, and registers the output
   /// Field. The output Field name is constructed as InputName + "_HorzMean".
   HorzMeanOp(const std::vector<std::string>
                  &UpstreamNames, ///< [in] input field names
              Config Options      ///< [in] operator config
              )
       : AnalysisOperator("HorzMean") {

      Error Err;

      // Store input field names
      InputNames = UpstreamNames;

      // Read optional weight field name (default: AreaCell)
      Err = Options.get("WeightField", WeightFieldName);
      if (Err.isFail()) {
         WeightFieldName = "AreaCell";
      }

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_HorzMean";
      OutputNames                 = {OutputFieldName};
      InstanceName                = OutputFieldName;

      // Allow callers to override the output name
      applyOutputNameOverride(Options);

      // Determine the vertical dimension from the input Field (last dimension).
      // The output is 1D over that vertical/interface dimension.
      auto InputField = Field::get(InputNames[0]);
      std::vector<std::string> InputDimNames;
      InputField->getDimNames(InputDimNames);
      I4 InputRank = static_cast<I4>(InputDimNames.size());
      if (InputRank < 2) {
         ABORT_ERROR("HorzMeanOp: Input field '{}' must be at least 2D "
                     "(horizontal x vertical)",
                     InputNames[0]);
      }

      std::string VertDimName = InputDimNames[InputRank - 1];
      auto VertDim            = Dimension::get(VertDimName);
      if (VertDim == nullptr) {
         ABORT_ERROR("HorzMeanOp: Vertical dimension '{}' not found for input "
                     "field '{}'",
                     VertDimName, InputNames[0]);
      }
      NVertSize = VertDim->getLengthGlobal();

      // Allocate 1D output data array (always Real type)
      OutputData = Array1DReal(OutputNames[0], NVertSize);

      // Register output Field with metadata over the vertical dimension. A
      // horizontal mean has the units and standard name of the field it
      // averages; its cell_methods record the reduction after any the input
      // already carries.
      I4 NDims                          = 1;
      std::vector<std::string> DimNames = {VertDimName};
      auto Meta        = inheritMetadata(InputNames[0], "area: mean");
      auto OutputField = createOutputField(
          OutputNames[0], "Area-weighted horizontal mean of " + InputNames[0],
          Meta,                              // CF metadata
          -std::numeric_limits<Real>::max(), // Min valid value
          std::numeric_limits<Real>::max(),  // Max valid value
          NDims,                             // Rank
          DimNames                           // Dimension names
      );

      // Attach output data array to Field
      OutputField->template attachData<Array1DReal>(OutputData);

   } // end constructor

   /// Initializes the operator after all Fields exist. Determines the
   /// horizontal index space from the input's second-to-last dimension name,
   /// stores the owned horizontal count and the per-column MinLayer/MaxLayer
   /// active-layer bounds from VertCoord, detects whether the input is on
   /// mid-layers or interfaces, and resolves the horizontal weight array.
   void initialize(const MachEnv *InEnv, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {

      AnalysisOperator::initialize(InEnv, InMesh, InVCoord, Options);

      auto InputField = Field::get(InputNames[0]);
      std::vector<std::string> InputDimNames;
      InputField->getDimNames(InputDimNames);
      I4 InputRank = static_cast<I4>(InputDimNames.size());

      // Horizontal dimension is second-to-last for 2D+ fields
      std::string IndexSpaceName = InputDimNames[InputRank - 2];

      if (IndexSpaceName == "NCells") {
         MinLayer    = VCoord->MinLayerCell;
         MaxLayer    = VCoord->MaxLayerCell;
         NHorizOwned = Mesh->NCellsOwned;
      } else if (IndexSpaceName == "NEdges") {
         MinLayer    = VCoord->MinLayerEdgeBot;
         MaxLayer    = VCoord->MaxLayerEdgeTop;
         NHorizOwned = Mesh->NEdgesOwned;
      } else if (IndexSpaceName == "NVertices") {
         MinLayer    = VCoord->MinLayerVertexBot;
         MaxLayer    = VCoord->MaxLayerVertexTop;
         NHorizOwned = Mesh->NVerticesOwned;
      } else {
         ABORT_ERROR("HorzMeanOp: Unknown index space {}", IndexSpaceName);
      }

      // Detect vertical staggering of the input field. The MinLayer/MaxLayer
      // arrays are layer indices; a midlayer field has NVertLayers levels and
      // uses the active range [MinLayer, MaxLayer], while an interface field
      // has NVertLayers+1 levels and uses [MinLayer, MaxLayer+1] (interfaces
      // bracket the active layers, including the bottom interface).
      AtInterfaces = (NVertSize == VCoord->NVertLayersP1);

      // Resolve the horizontal weight array. Default is the mesh AreaCell; a
      // named 1D field may be supplied via the WeightField config parameter.
      if (WeightFieldName == "AreaCell") {
         WeightArray = Mesh->AreaCell;
      } else {
         auto WeightField = Field::get(WeightFieldName);
         if (WeightField == nullptr) {
            ABORT_ERROR("HorzMeanOp: Weight field '{}' not found in Field "
                        "registry",
                        WeightFieldName);
         }
         WeightArray = WeightField->template getDataArray<Array1DReal>();
      }

      if (static_cast<I4>(WeightArray.extent(0)) < NHorizOwned) {
         ABORT_ERROR("HorzMeanOp: Weight field '{}' size {} smaller than owned "
                     "horizontal count {}",
                     WeightFieldName, WeightArray.extent(0), NHorizOwned);
      }

      // Check if the primary input Field carries a 1D regional mask (attached
      // by an upstream ExtractRegionOp). If so, store it for use in compute().
      if (InputField->hasRegionalMask()) {
         RegionalMask = InputField->getRegionalMask();
         if (static_cast<I4>(RegionalMask.extent(0)) < NHorizOwned) {
            ABORT_ERROR(
                "HorzMeanOp: Regional mask on field '{}' has size {} which is "
                "smaller than the owned horizontal count {}",
                InputNames[0], RegionalMask.extent(0), NHorizOwned);
         }
         UseRegionalMask = true;
      }

   } // end initialize

   /// Computes the area-weighted horizontal mean per vertical level. For each
   /// vertical index the operator accumulates a weighted numerator and
   /// denominator over the owned horizontal locations, performs a single global
   /// reduction over the concatenated numerator/denominator vectors, and
   /// divides to obtain the mean (guarding against a zero denominator).
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      const I4 NVert  = NVertSize;
      const I4 NHoriz = NHorizOwned;

      // Lazily allocate the per-level local reduction result arrays. These are
      // kept as device arrays so that when a device-aware globalSum overload
      // becomes available the host sync below can be removed entirely.
      if (LocalNum.size() == 0) {
         LocalNum = Array1DReal("HorzMeanLocalNum", NVert);
         LocalDen = Array1DReal("HorzMeanLocalDen", NVert);
      }

      // Local device references captured by value into the kernel.
      OMEGA_SCOPE(LocInput, InputData);
      OMEGA_SCOPE(LocMinLayer, MinLayer);
      OMEGA_SCOPE(LocMaxLayer, MaxLayer);
      const bool LocAtInterfaces = AtInterfaces;
      OMEGA_SCOPE(LocWeight, WeightArray);
      OMEGA_SCOPE(LocNum, LocalNum);
      OMEGA_SCOPE(LocDen, LocalDen);
      OMEGA_SCOPE(LocRegMask, RegionalMask);
      const bool LocUseRegMask = UseRegionalMask;

      // Single on-device kernel: one team per vertical interface (outer over
      // the vertical dimension), with an inner reduction over the owned
      // horizontal locations accumulating the weighted numerator and the
      // weighted denominator simultaneously.
      parallelForOuter(
          "HorzMeanReduce", {NVert},
          KOKKOS_LAMBDA(int K, const TeamMember &Team) {
             Real SumNum, SumDen;
             parallelReduceInner(
                 Team, NHoriz,
                 INNER_LAMBDA(int IHoriz, Real &AccNum, Real &AccDen) {
                    // The MinLayer/MaxLayer arrays give the per-column active
                    // layer index range. For a midlayer input the active range
                    // is [MinLayer, MaxLayer]. For an interface input (one more
                    // level) the active range is extended by one level to
                    // [MinLayer, MaxLayer+1] so the bottom interface bracketing
                    // the deepest active layer is included.
                    const I4 KMin = LocMinLayer(IHoriz);
                    const I4 KMax =
                        LocMaxLayer(IHoriz) + (LocAtInterfaces ? 1 : 0);
                    const bool Active = (K >= KMin && K <= KMax);

                    const Real Value = LocInput(IHoriz, K);

                    // Presence: physically active cell AND inside the region.
                    // When a 1D regional mask is attached to the input Field
                    // (propagated from an upstream ExtractRegionOp) use it
                    // directly. Otherwise fall back to the legacy non-zero
                    // value check so that ExtractRegion-zeroing chains
                    // continue to work without a mask attachment.
                    const bool Present = LocUseRegMask
                                             ? (LocRegMask(IHoriz) != 0)
                                             : (Value != 0);
                    if (Active && Present) {
                       const Real Weight = LocWeight(IHoriz);
                       AccNum += Weight * Value;
                       AccDen += Weight;
                    }
                 },
                 SumNum, SumDen);

             // Inner reduction results are collective across the team, so every
             // thread has the final values; store once per team.
             LocNum(K) = SumNum;
             LocDen(K) = SumDen;
          });

      // TODO: replace this host sync + std::vector globalSum with a
      // device-array globalSum overload once it is supported.
      auto LocalNumHost = createHostMirrorCopy(LocalNum);
      auto LocalDenHost = createHostMirrorCopy(LocalDen);

      // Concatenated reduction buffer: [Num(0..NVert-1), Den(0..NVert-1)]
      std::vector<R8> LocalSums(2 * NVert, 0.0);
      for (I4 K = 0; K < NVert; ++K) {
         LocalSums[K]         = static_cast<R8>(LocalNumHost(K));
         LocalSums[NVert + K] = static_cast<R8>(LocalDenHost(K));
      }

      // Single global reduction over the concatenated numerator/denominator
      std::vector<R8> GlobalSums = globalSum(LocalSums, Comm);

      // Divide numerator by denominator (guarding zero denominator)
      auto OutputHost = Kokkos::create_mirror_view(OutputData);
      for (I4 K = 0; K < NVert; ++K) {
         const R8 Num = GlobalSums[K];
         const R8 Den = GlobalSums[NVert + K];
         OutputHost(K) =
             (Den != 0.0) ? static_cast<Real>(Num / Den) : static_cast<Real>(0);
      }
      Kokkos::deep_copy(OutputData, OutputHost);

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array holding the per-level horizontal mean (Real type)
   Array1DReal OutputData;

   /// Vertical (or interface) dimension size of the input/output field
   I4 NVertSize;

   /// Number of owned horizontal points (for MPI correctness)
   I4 NHorizOwned;

   /// Per-column active layer index bounds from VertCoord. The active vertical
   /// range for a column is [MinLayer, MaxLayer] for a midlayer input and
   /// [MinLayer, MaxLayer+1] for an interface input (see AtInterfaces).
   Array1DI4 MinLayer;
   Array1DI4 MaxLayer;

   /// True when the input field is dimensioned on interfaces (NVertLayersP1)
   /// rather than mid-layers (NVertLayers). Detected in initialize().
   bool AtInterfaces = false;

   /// Horizontal weight array (default AreaCell)
   Array1DReal WeightArray;

   /// Per-level local (owned-only) reduction results kept on the device.
   /// Lazily allocated on first compute.
   Array1DReal LocalNum;
   Array1DReal LocalDen;

   /// Name of the weight field (default "AreaCell")
   std::string WeightFieldName;

   /// 1D regional mask (horizontal dimension) propagated from the primary
   /// input Field by an upstream ExtractRegionOp. Mask value 0 excludes the
   /// horizontal point from both numerator and denominator; non-zero includes
   /// it. Empty (not set) when no mask was attached to the input Field.
   Array1DI4 RegionalMask;

   /// True when the primary input Field carried a regional mask at initialize()
   /// time, i.e. RegionalMask is valid and should be applied in compute().
   bool UseRegionalMask = false;

}; // end class HorzMeanOp

} // end namespace OMEGA

#endif
