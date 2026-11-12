#ifndef OMEGA_BINARYMULTIPLYOP_H
#define OMEGA_BINARYMULTIPLYOP_H

//===-- analysis/operators/BinaryMultiplyOp.h -------------------*- C++ -*-===//
//
/// \file
/// \brief Defines the BinaryMultiplyOp operator for element-wise field
/// multiplication
///
/// BinaryMultiplyOp performs element-wise multiplication of two fields with
/// matching dimensions. This operator is useful for computing products like:
/// - Flux calculations (velocity × area)
/// - Density-weighted quantities (field × density)
/// - Coordinate transformations (field1 × field2)
///
/// The operator supports two modes:
/// 1. Same-rank multiplication: Both fields have the same dimensions
///    Output(i,j,...) = Field1(i,j,...) × Field2(i,j,...)
///
/// 2. Vertical expansion: Field1 is 2D/3D and Field2 is 1D (horizontal only)
///    The 1D field value is replicated across all vertical layers
///    Output(i,j) = Field1(i,j) × Field2(i)
///
/// Example usage in operator chain:
/// \code
///   VerticalVelocity_BinaryMultiply(AreaCell)
/// \endcode
/// Computes the vertical flux by multiplying 2D vertical velocity by 1D cell
/// area, with the area value replicated across all vertical layers.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"

namespace OMEGA {

/// BinaryMultiplyOp performs element-wise multiplication of two fields.
/// Supports both same-rank multiplication and vertical expansion where 1D
/// fields are replicated across all vertical layers when multiplying with 2D/3D
/// fields. The output field has the same dimensions as the first input field.
template <typename ArrayT> class BinaryMultiplyOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type - same as input array types
   using OutputArrayT = ArrayT;

   /// Constructs a BinaryMultiplyOp operator. Retrieves both input field
   /// dimensions and metadata, validates that they match, creates output
   /// Field for the product, allocates output data array, and registers
   /// the output Field. The output Field name is constructed as
   /// Field1Name + "_" + Field2Name + "_Product".
   BinaryMultiplyOp(const std::vector<std::string>
                        &UpstreamNames, ///< [in] input field names
                    Config Options      ///< [in] operator config
                    )
       : AnalysisOperator("BinaryMultiply") {

      // Store input field names
      // UpstreamNames[0] = first field
      // UpstreamNames[1] = second field
      InputNames = UpstreamNames;

      if (InputNames.size() < 2) {
         ABORT_ERROR("BinaryMultiplyOp: Requires 2 input fields");
      }

      // Retrieve both input fields
      auto Field1 = Field::get(InputNames[0]);
      auto Field2 = Field::get(InputNames[1]);
      auto Data1  = Field1->template getDataArray<ArrayT>();

      // Get rank information for both fields
      const I4 Rank1 = ArrayT::rank;
      const I4 Rank2 = Field2->getNumDims();

      // Store whether Field2 is 1D (will be expanded vertically if Field1 is
      // 2D/3D)
      IsVerticalExpansion = (Rank1 > 1 && Rank2 == 1);

      // Validate dimension compatibility
      if (!IsVerticalExpansion && Rank1 != Rank2) {
         ABORT_ERROR("BinaryMultiplyOp: Input fields must have same rank or "
                     "second field must be 1D for vertical expansion, got "
                     "ranks {} and {}",
                     Rank1, Rank2);
      }

      // For same-rank case, validate all dimensions match
      if (!IsVerticalExpansion) {
         auto Data2 = Field2->template getDataArray<ArrayT>();
         for (I4 dim = 0; dim < Rank1; ++dim) {
            if (Data1.extent(dim) != Data2.extent(dim)) {
               ABORT_ERROR("BinaryMultiplyOp: Dimension {} mismatch: "
                           "Field1={}, Field2={}",
                           dim, Data1.extent(dim), Data2.extent(dim));
            }
         }
      } else {
         // For vertical expansion case, validate horizontal dimensions match
         auto Data2_1D = Field2->template getDataArray<Array1DReal>();
         if (Data1.extent(Rank1 - 2) != Data2_1D.extent(0)) {
            ABORT_ERROR("BinaryMultiplyOp: Horizontal dimension mismatch for "
                        "vertical expansion: "
                        "Field1 horizontal dim={}, Field2 size={}",
                        Data1.extent(Rank1 - 2), Data2_1D.extent(0));
         }
      }

      // Store dimensions based on rank
      if (Rank1 == 1) {
         // 1D: horizontal only (cells/edges/vertices), no vertical structure
         NHorizDim = Data1.extent(0);
         NVertSize = 1; // No vertical dimension
      } else if (Rank1 == 2) {
         // 2D: NHorizDim × NVertSize
         NHorizDim = Data1.extent(0);
         NVertSize = Data1.extent(1);
      } else if (Rank1 == 3) {
         // 3D: OuterDim × NHorizDim × NVertSize
         NHorizDim = Data1.extent(1);
         NVertSize = Data1.extent(2);
      } else {
         ABORT_ERROR("BinaryMultiplyOp: Unsupported rank {}. "
                     "Supports 1D, 2D, and 3D arrays",
                     Rank1);
      }

      // Get dimension info from first field
      auto NDims = Field1->getNumDims();
      std::vector<std::string> DimNames;
      Field1->getDimNames(DimNames);

      // Construct output field name and set instance name
      // Format: Input1_BinaryMultiply(Input2) to match chain syntax
      std::string OutputFieldName =
          InputNames[0] + "_BinaryMultiply(" + InputNames[1] + ")";
      OutputNames  = {OutputFieldName};
      InstanceName = OutputFieldName;

      // Get descriptions from both input fields
      std::string Descr1, Descr2;
      Field1->getMetadata("long_name", Descr1);
      Field2->getMetadata("long_name", Descr2);
      std::string OutputDescr = "Product of " + Descr1 + " and " + Descr2;

      // A product has the product of the input units (m s-1 times m2 is
      // m3 s-1), which are unknown if either input's are. It is a new
      // quantity, so it has no standard name or cell methods.
      InheritedMetadata Meta;
      Meta.Units = CFUnits::multiply(getFieldAttribute(InputNames[0], "units"),
                                     getFieldAttribute(InputNames[1], "units"));

      // Create output Field
      auto OutputField = createOutputField(
          OutputNames[0],
          OutputDescr,                          // Description
          Meta,                                 // CF metadata
          -std::numeric_limits<ScalarT>::max(), // Min valid value
          std::numeric_limits<ScalarT>::max(),  // Max valid value
          NDims,                                // Rank
          DimNames                              // Dimension names
      );

      // Allocate output data array matching input layout
      OutputData = OutputArrayT(OutputNames[0] + "_out", Data1.layout());

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

      // Propagate regional mask from first input to output if present
      // (Binary operations preserve the spatial structure of the first input)
      if (Field1->hasRegionalMask()) {
         OutputField->setRegionalMask(Field1->getRegionalMask());
      }

   } // end constructor

   /// Initializes the operator after all Fields exist. Determines the index
   /// space (cells, edges, or vertices) from the first input field's horizontal
   /// dimension name and stores the appropriate MinLayer/MaxLayer arrays
   /// from VertCoord for bounding the inner vertical loop. Also determines
   /// whether the horizontal dimension is mesh-distributed and stores the
   /// appropriate owned count for MPI-correct parallel loops.
   void initialize(const MachEnv *Env, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {

      AnalysisOperator::initialize(Env, InMesh, InVCoord, Options);

      constexpr I4 InputRank = ArrayT::rank;

      // For 1D arrays, check if it's a mesh dimension to determine owned count
      if constexpr (InputRank == 1) {
         auto Field1 = Field::get(InputNames[0]);
         std::vector<std::string> DimNames;
         Field1->getDimNames(DimNames);
         std::string IndexSpaceName = DimNames[0];

         IsMeshDimension =
             (IndexSpaceName == "NCells" || IndexSpaceName == "NEdges" ||
              IndexSpaceName == "NVertices");

         if (IsMeshDimension) {
            if (IndexSpaceName == "NCells") {
               NHorizOwned = Mesh->NCellsOwned;
            } else if (IndexSpaceName == "NEdges") {
               NHorizOwned = Mesh->NEdgesOwned;
            } else if (IndexSpaceName == "NVertices") {
               NHorizOwned = Mesh->NVerticesOwned;
            }
         } else {
            ABORT_ERROR("BinaryMultiplyOp: Unknown index space {}",
                        IndexSpaceName);
         }
      }

      if constexpr (InputRank > 1) {
         auto Field1 = Field::get(InputNames[0]);
         std::vector<std::string> DimNames;
         Field1->getDimNames(DimNames);
         // Horizontal dimension is 2nd-to-last for 2D/3D
         std::string IndexSpaceName = DimNames[InputRank - 2];

         IsMeshDimension =
             (IndexSpaceName == "NCells" || IndexSpaceName == "NEdges" ||
              IndexSpaceName == "NVertices");

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
            ABORT_ERROR("BinaryMultiplyOp: Unknown index space {}",
                        IndexSpaceName);
         }
      }

   } // end initialize

   /// Computes the element-wise multiplication by retrieving both input
   /// data arrays and performing parallel multiplication using hierarchical
   /// parallelism (outer loop over horizontal dimension, inner loop over
   /// vertical layers bounded by MinLayer/MaxLayer). Supports 1D (horizontal
   /// only), 2D (horizontal × vertical), and 3D (additional outer dimension ×
   /// horizontal × vertical) arrays. Also supports vertical expansion where a
   /// 1D field value is replicated across all vertical layers when multiplying
   /// with 2D/3D fields. Updates output data, timestamp, and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Get array rank for conditional logic
      constexpr I4 InputRank = ArrayT::rank;

      // Retrieve input Fields
      auto Field1 = Field::get(InputNames[0]);
      auto Field2 = Field::get(InputNames[1]);
      auto Data1  = Field1->template getDataArray<ArrayT>();

      auto Output = OutputData;

      if constexpr (InputRank == 1) {
         // 1D case: horizontal field only, no vertical structure
         auto Data2 = Field2->template getDataArray<ArrayT>();
         parallelFor(
             {NHorizOwned}, KOKKOS_LAMBDA(int IHoriz) {
                Output(IHoriz) =
                    static_cast<ScalarT>(static_cast<Real>(Data1(IHoriz)) *
                                         static_cast<Real>(Data2(IHoriz)));
             });

      } else if constexpr (InputRank == 2) {
         // 2D case: hierarchical parallelism over horizontal × vertical
         // Inner loop bounded by MinLayer/MaxLayer for partial columns

         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);

         if (IsVerticalExpansion) {
            // Vertical expansion: Data2 is 1D, replicate across vertical layers
            auto Data2_1D = Field2->template getDataArray<Array1DReal>();
            parallelForOuter(
                "BinaryMultiply2D_VertExpand", LaunchConfig({NHorizOwned}),
                KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                   const Real Data2Val = static_cast<Real>(Data2_1D(IHoriz));
                   const I4 KMin       = LocMinLayer(IHoriz);
                   const I4 KMax       = LocMaxLayer(IHoriz);
                   const I4 KRange     = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K        = KMin + KIdx;
                          Output(IHoriz, K) = static_cast<ScalarT>(
                              static_cast<Real>(Data1(IHoriz, K)) * Data2Val);
                       });
                });
         } else {
            // Same rank: element-wise multiplication
            auto Data2 = Field2->template getDataArray<ArrayT>();
            parallelForOuter(
                "BinaryMultiply2D", LaunchConfig({NHorizOwned}),
                KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz);
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K        = KMin + KIdx;
                          Output(IHoriz, K) = static_cast<ScalarT>(
                              static_cast<Real>(Data1(IHoriz, K)) *
                              static_cast<Real>(Data2(IHoriz, K)));
                       });
                });
         }

      } else if constexpr (InputRank == 3) {
         // 3D case: hierarchical parallelism over dim0 × horizontal × vertical
         // Inner loop bounded by MinLayer/MaxLayer for partial columns

         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         const I4 Dim0 = Data1.extent(0);

         if (IsVerticalExpansion) {
            // Vertical expansion: Data2 is 1D, replicate across vertical layers
            auto Data2_1D = Field2->template getDataArray<Array1DReal>();
            parallelForOuter(
                "BinaryMultiply3D_VertExpand",
                LaunchConfig({Dim0, NHorizOwned}),
                KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                   const Real Data2Val = static_cast<Real>(Data2_1D(IHoriz));
                   const I4 KMin       = LocMinLayer(IHoriz);
                   const I4 KMax       = LocMaxLayer(IHoriz);
                   const I4 KRange     = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K            = KMin + KIdx;
                          Output(I0, IHoriz, K) = static_cast<ScalarT>(
                              static_cast<Real>(Data1(I0, IHoriz, K)) *
                              Data2Val);
                       });
                });
         } else {
            // Same rank: element-wise multiplication
            auto Data2 = Field2->template getDataArray<ArrayT>();
            parallelForOuter(
                "BinaryMultiply3D", LaunchConfig({Dim0, NHorizOwned}),
                KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                   const I4 KMin   = LocMinLayer(IHoriz);
                   const I4 KMax   = LocMaxLayer(IHoriz);
                   const I4 KRange = vertRange(KMin, KMax);
                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K            = KMin + KIdx;
                          Output(I0, IHoriz, K) = static_cast<ScalarT>(
                              static_cast<Real>(Data1(I0, IHoriz, K)) *
                              static_cast<Real>(Data2(I0, IHoriz, K)));
                       });
                });
         }
      }

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array holding the product field values
   OutputArrayT OutputData;

   /// Number of points in horizontal dimension (cells, edges, or vertices)
   I4 NHorizDim;

   /// Number of owned points in horizontal dimension (for MPI correctness)
   I4 NHorizOwned;

   /// Vertical size of the array
   I4 NVertSize;

   /// Whether the horizontal dimension is a mesh-distributed dimension
   bool IsMeshDimension;

   /// Whether Field2 is 1D and should be replicated across vertical dimension
   bool IsVerticalExpansion;

   /// Min active layer index for each horizontal point (only for 2D/3D arrays)
   Array1DI4 MinLayer;

   /// Max active layer index for each horizontal point (only for 2D/3D arrays)
   Array1DI4 MaxLayer;

}; // end class BinaryMultiplyOp

} // end namespace OMEGA

#endif
