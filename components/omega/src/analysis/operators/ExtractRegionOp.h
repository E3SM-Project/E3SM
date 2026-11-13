#ifndef OMEGA_EXTRACTREGIONOP_H
#define OMEGA_EXTRACTREGIONOP_H

//===-- analysis/operators/ExtractRegionOp.h - ExtractRegionOp -*- C++ -*-===//
//
/// \file
/// \brief Defines the ExtractRegionOp operator for applying regional masks
///
/// ExtractRegionOp applies a named regional mask to a field by multiplying
/// the field values by the mask. The mask must already exist in the Field
/// registry as an Array1DI4 field. This operator enables regional statistics
/// by multiplying field values by mask values (where mask == 0 excludes points
/// by leaving output unchanged), which can then be processed by spatial
/// reduction operators that use the attached regional mask.
///
/// The operator multiplies each horizontal location by its mask value:
///   Output(i, k, ...) = Input(i, k, ...) * Mask(i)
///
/// where Mask(i) is typically 0 (exclude) or 1 (include).
///
/// Configuration:
/// - MaskName: Name of the mask field in the Field registry (required)
///
/// Example usage in operator chain:
/// \code
///   Temperature_ExtractRegion(Atlantic)_SpatialMean
/// \endcode
/// where `ExtractRegion(Atlantic)` applies the "Atlantic" mask to the
/// Temperature field, enabling regional spatial statistics.
///
/// The mask field should be Array1DI4 with dimension matching the horizontal
/// dimension of the input field (NCells, NEdges, or NVertices). Values are:
/// - 0: exclude this horizontal location (data not copied to output)
/// - 1: include this horizontal location (data multiplied by 1 and copied)
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"

namespace OMEGA {

/// ExtractRegionOp applies a named regional mask from the Field registry to
/// a field by multiplying field values with the mask. The operator copies
/// masked data to the output field and attaches the mask to enable downstream
/// spatial reduction operators to exclude masked-out points.
template <typename ArrayT> class ExtractRegionOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type - same as input array type
   using OutputArrayT = ArrayT;

   /// Constructs an ExtractRegionOp operator. Reads the mask field name from
   /// config, retrieves the mask from the Field registry, creates an output
   /// Field matching the input dimensions, and attaches the regional mask to
   /// the output Field. The output Field name is constructed as
   /// InputName + "_" + MaskName (without "Mask" suffix if present).
   ExtractRegionOp(const std::vector<std::string>
                       &UpstreamNames, ///< [in] input field names
                   Config Options      ///< [in] operator config
                   )
       : AnalysisOperator("ExtractRegion") {

      // Store input field names
      InputNames = UpstreamNames;

      // Read required MaskName parameter from configuration
      std::string MaskFieldName;
      Error Err = Options.get("MaskName", MaskFieldName);

      if (Err.isFail()) {
         ABORT_ERROR("ExtractRegion Op: Required parameter 'MaskName' not "
                     "found in configuration");
      }

      // Retrieve the mask field from registry
      auto MaskField = Field::get(MaskFieldName);
      if (MaskField == nullptr) {
         ABORT_ERROR("ExtractRegionOp: Mask field '{}' not found in Field "
                     "registry",
                     MaskFieldName);
      }

      // Verify mask is 1D integer array
      if (MaskField->getNumDims() != 1) {
         ABORT_ERROR("ExtractRegionOp: Mask field '{}' must be 1D, got {}D",
                     MaskFieldName, MaskField->getNumDims());
      }

      // Extract mask data
      RegionalMask = MaskField->template getDataArray<Array1DI4>();

      // Retrieve input Field to get dimensions and metadata
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Support 1D (horizontal), 2D (horizontal × vertical), and 3D arrays
      constexpr I4 InputRank = ArrayT::rank;

      // Get dimension names from input field
      std::vector<std::string> InputDimNames;
      InputField->getDimNames(InputDimNames);

      // Get dimension name from mask field
      std::vector<std::string> MaskDimNames;
      MaskField->getDimNames(MaskDimNames);

      // Determine the horizontal dimension name based on rank
      std::string InputHorizDimName;
      if constexpr (InputRank == 1) {
         // 1D: horizontal only
         InputHorizDimName = InputDimNames[0];
         NHorizDim         = InputData.extent(0);
         NVertSize         = 1; // No vertical dimension
      } else if constexpr (InputRank == 2) {
         // 2D: NHorizDim × NVertSize (horizontal is first dimension)
         InputHorizDimName = InputDimNames[0];
         NHorizDim         = InputData.extent(0);
         NVertSize         = InputData.extent(1);
      } else if constexpr (InputRank == 3) {
         // 3D: Dim0 × NHorizDim × NVertSize (horizontal is second dimension)
         InputHorizDimName = InputDimNames[1];
         Dim0              = InputData.extent(0);
         NHorizDim         = InputData.extent(1);
         NVertSize         = InputData.extent(2);
      } else {
         ABORT_ERROR("ExtractRegionOp: Unsupported rank {}. "
                     "Supports 1D, 2D, and 3D arrays",
                     InputRank);
      }

      // Verify mask dimension name matches input horizontal dimension name
      if (MaskDimNames[0] != InputHorizDimName) {
         ABORT_ERROR("ExtractRegionOp: Mask dimension name '{}' does not match "
                     "input horizontal dimension name '{}'",
                     MaskDimNames[0], InputHorizDimName);
      }

      // Construct output field name: InputName + "_" + RegionName
      // Strip "Mask" suffix from MaskFieldName if present for cleaner names
      std::string RegionName = MaskFieldName;
      if (RegionName.size() > 4 &&
          RegionName.substr(RegionName.size() - 4) == "Mask") {
         RegionName = RegionName.substr(0, RegionName.size() - 4);
      }

      std::string OutputFieldName =
          InputNames[0] + "_ExtractRegion(" + RegionName + ")";
      OutputNames  = {OutputFieldName};
      InstanceName = OutputFieldName;

      // Get dimension info from input
      auto NDims = InputField->getNumDims();
      std::vector<std::string> DimNames;
      InputField->getDimNames(DimNames);

      // Get input metadata. Masking a region leaves the units, standard name
      // and cell methods unchanged.
      std::string InputDescr;
      ScalarT InputValidMin, InputValidMax;

      InputField->getMetadata("long_name", InputDescr);
      InputField->getMetadata("valid_min", InputValidMin);
      InputField->getMetadata("valid_max", InputValidMax);
      auto Meta = inheritMetadata(InputNames[0]);

      // Create output Field with same dimensions as input
      auto OutputField = createOutputField(
          OutputNames[0],
          InputDescr + " (region: " + RegionName + ")", // Description
          Meta,                                         // CF metadata
          InputValidMin,                                // Min valid (unchanged)
          InputValidMax,                                // Max valid (unchanged)
          NDims,                                        // Rank
          DimNames                                      // Dimension names
      );

      // Allocate output data array matching input layout
      OutputData = OutputArrayT(OutputNames[0] + "_out", InputData.layout());

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

   } // end constructor

   /// Initializes the operator after all Fields exist. Determines the index
   /// space (cells, edges, or vertices) from the input field's horizontal
   /// dimension name and stores the appropriate MinLayer/MaxLayer arrays
   /// from VertCoord for bounding the inner vertical loop.
   void initialize(const MachEnv *Env, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {

      AnalysisOperator::initialize(Env, InMesh, InVCoord, Options);

      constexpr I4 InputRank = ArrayT::rank;

      // For 2D/3D arrays, get MinLayer/MaxLayer for partial columns
      if constexpr (InputRank > 1) {
         auto Field1 = Field::get(InputNames[0]);
         std::vector<std::string> DimNames;
         Field1->getDimNames(DimNames);
         // Horizontal dimension is 2nd-to-last for 2D/3D
         std::string IndexSpaceName = DimNames[InputRank - 2];

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
            ABORT_ERROR("ExtractRegionOp: Unknown index space {}",
                        IndexSpaceName);
         }
      } else {
         // For 1D arrays, determine owned count
         auto Field1 = Field::get(InputNames[0]);
         std::vector<std::string> DimNames;
         Field1->getDimNames(DimNames);
         std::string IndexSpaceName = DimNames[0];

         if (IndexSpaceName == "NCells") {
            NHorizOwned = Mesh->NCellsOwned;
         } else if (IndexSpaceName == "NEdges") {
            NHorizOwned = Mesh->NEdgesOwned;
         } else if (IndexSpaceName == "NVertices") {
            NHorizOwned = Mesh->NVerticesOwned;
         } else {
            ABORT_ERROR("ExtractRegionOp: Unknown index space {}",
                        IndexSpaceName);
         }
      }

      // Attach the regional mask to the output Field
      // If input already has a mask, compute intersection
      auto InputField = Field::get(InputNames[0]);
      if (InputField->hasRegionalMask()) {
         Array1DI4 InputMask = InputField->getRegionalMask();
         Array1DI4 IntersectionMask =
             Array1DI4("IntersectionMask", RegionalMask.extent(0));
         OMEGA_SCOPE(LocInputMask, InputMask);
         OMEGA_SCOPE(LocRegionalMask, RegionalMask);
         OMEGA_SCOPE(LocIntersectionMask, IntersectionMask);
         parallelFor(
             {static_cast<I4>(RegionalMask.extent(0))}, KOKKOS_LAMBDA(int I) {
                LocIntersectionMask(I) = LocInputMask(I) * LocRegionalMask(I);
             });
         auto OutputField = Field::get(OutputNames[0]);
         OutputField->setRegionalMask(IntersectionMask);
      } else {
         auto OutputField = Field::get(OutputNames[0]);
         OutputField->setRegionalMask(RegionalMask);
      }

   } // end initialize

   /// Computes the regional extraction by multiplying input data with the
   /// regional mask. Uses hierarchical parallelism with MinLayer/MaxLayer
   /// bounding for partial columns. The mask (1D horizontal) is applied
   /// across all vertical layers and extra dimensions.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Get array rank for conditional logic
      constexpr I4 InputRank = ArrayT::rank;

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      OMEGA_SCOPE(LocOutput, OutputData);
      OMEGA_SCOPE(LocMask, RegionalMask);

      if constexpr (InputRank == 1) {
         // 1D case: horizontal field only
         parallelFor(
             {NHorizOwned}, KOKKOS_LAMBDA(int iHoriz) {
                const Real MaskVal = static_cast<Real>(LocMask(iHoriz));
                if (MaskVal != 0) {
                   LocOutput(iHoriz) = static_cast<ScalarT>(
                       static_cast<Real>(InputData(iHoriz)) *
                       static_cast<Real>(MaskVal));
                }
             });

      } else if constexpr (InputRank == 2) {
         // 2D case: hierarchical parallelism over horizontal × vertical
         // Inner loop bounded by MinLayer/MaxLayer for partial columns
         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);

         parallelForOuter(
             "ExtractRegion2D", LaunchConfig({NHorizOwned}),
             KOKKOS_LAMBDA(int iHoriz, const TeamMember &Team) {
                const Real MaskVal = static_cast<Real>(LocMask(iHoriz));
                if (MaskVal != 0) {
                   const I4 KMin   = LocMinLayer(iHoriz);
                   const I4 KMax   = LocMaxLayer(iHoriz);
                   const I4 KRange = vertRange(KMin, KMax);

                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K           = KMin + KIdx;
                          LocOutput(iHoriz, K) = static_cast<ScalarT>(
                              static_cast<Real>(InputData(iHoriz, K)) *
                              MaskVal);
                       });
                }
             });

      } else if constexpr (InputRank == 3) {
         // 3D case: hierarchical parallelism over dim0 × horizontal
         // Inner loop bounded by MinLayer/MaxLayer for partial columns
         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         OMEGA_SCOPE(LocDim0, Dim0);

         parallelForOuter(
             "ExtractRegion3D", LaunchConfig({LocDim0, NHorizOwned}),
             KOKKOS_LAMBDA(int i0, int iHoriz, const TeamMember &Team) {
                const Real MaskVal = static_cast<Real>(LocMask(iHoriz));
                if (MaskVal != 0) {
                   const I4 KMin   = LocMinLayer(iHoriz);
                   const I4 KMax   = LocMaxLayer(iHoriz);
                   const I4 KRange = vertRange(KMin, KMax);

                   parallelForInner(
                       Team, KRange, INNER_LAMBDA(int KIdx) {
                          const I4 K               = KMin + KIdx;
                          LocOutput(i0, iHoriz, K) = static_cast<ScalarT>(
                              static_cast<Real>(InputData(i0, iHoriz, K)) *
                              MaskVal);
                       });
                }
             });
      }

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

 private:
   /// Output data array holding the masked field values
   OutputArrayT OutputData;

   /// The regional mask (1D horizontal)
   Array1DI4 RegionalMask;

   /// Number of points in horizontal dimension
   I4 NHorizDim;

   /// Number of owned points in horizontal dimension (for MPI correctness)
   I4 NHorizOwned;

   /// Vertical size of the field array
   I4 NVertSize;

   /// Size of first dimension for 3D arrays
   I4 Dim0;

   /// Min active layer index for each horizontal point (only for 2D/3D arrays)
   Array1DI4 MinLayer;

   /// Max active layer index for each horizontal point (only for 2D/3D arrays)
   Array1DI4 MaxLayer;

}; // end class ExtractRegionOp

} // end namespace OMEGA

#endif
