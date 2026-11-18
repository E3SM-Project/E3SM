#ifndef OMEGA_PREFIXSUMOP_H
#define OMEGA_PREFIXSUMOP_H

//===-- analysis/operators/PrefixSumOp.h - PrefixSumOp ----------*- C++ -*-===//
//
/// \file
/// \brief Defines the PrefixSumOp operator for cumulative summation
///
/// PrefixSumOp computes the cumulative sum (prefix sum, or scan) of an input
/// field along a specified dimension. The operator supports both forward
/// (inclusive scan from start to end) and reverse (inclusive scan from end to
/// start) directions.
///
/// Mathematical formulation:
/// - Forward scan: Output[i] = sum(Input[0] + Input[1] + ... + Input[i])
/// - Reverse scan: Output[i] = sum(Input[i] + Input[i+1] + ... + Input[n-1])
///
/// This is particularly useful for vertical integration in ocean models, where
/// the MOC (Meridional Overturning Circulation) stream function requires
/// cumulative integration of transport from bottom to top. The reverse scan
/// direction is typically used for this case, starting from the ocean bottom
/// (maximum depth index) and integrating upward to the surface.
///
/// The operator is templated on the Kokkos array type (ArrayT) of the input
/// field. Currently supports 2D arrays (e.g., latitude bins × depth levels).
/// The output has the same shape and dimensions as the input.
///
/// An optional boundary condition (BC) field may be provided as a second
/// upstream input. When present, the BC field is used to seed the accumulator
/// at index 0 of the scan dimension before processing the first input element.
/// This is used in regional MOC calculations to anchor the streamfunction to
/// the transport through a bounding transect at the southern boundary. The BC
/// field must be 1D with the same size as the non-scan dimension (e.g.,
/// NVertLevels for a scan over latitude bins). The BC field is expected to be
/// in the same units as the input field (m³/s, not Sverdrups).
///
/// Configuration:
/// - Dimension: The dimension along which to compute cumulative sum (required)
///              0 = first dimension, 1 = second dimension
///              Note: For 2D arrays with mesh dimensions (NCells, NEdges,
///              NVertices), scanning along dimension 0 is not allowed in
///              distributed (MPI) configurations. Horizontal scans are only
///              supported for non-distributed dimensions (e.g., bins).
/// - Reverse: If true, scan from end to start; if false, scan from start to
///            end (default: false)
///
/// Example usage in operator chain (with BC):
/// \code
///   BinnedFlux_PrefixSum(BC=TransectMOC_TransectName)
/// \endcode
/// where the BC field provides the southern boundary transport for a regional
/// MOC calculation.
///
/// Example usage in operator chain (without BC):
/// \code
///   BinnedTransport_PrefixSum
/// \endcode
/// where PrefixSum performs cumulative vertical integration for MOC, typically
/// with Dimension=1 and Reverse=true to integrate from ocean bottom to top.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"

namespace OMEGA {

/// PrefixSumOp computes the cumulative sum of an input field along a specified
/// dimension. Supports forward and reverse directions for flexibility. The
/// output array has the same shape as the input. Currently optimized for 2D
/// arrays common in MOC calculations.
template <typename ArrayT> class PrefixSumOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type - same as input array type
   using OutputArrayT = ArrayT;

   /// Constructs a PrefixSumOp operator. Creates output Field with the same
   /// dimensions as the input field, allocates output data array, and
   /// registers the output Field in the Field registry. The output Field name
   /// is constructed as InputName + "_PrefixSum". Reads the dimension and
   /// direction parameters from the Options config. If a second upstream name
   /// is provided it is treated as the boundary condition (BC) field name.
   PrefixSumOp(const std::vector<std::string>
                   &UpstreamNames, ///< [in] input field names
               Config Options      ///< [in] operator config
               )
       : AnalysisOperator("PrefixSum") {

      // Store input field names
      InputNames = UpstreamNames;

      // Optional second upstream: boundary condition field
      // Fetch the BC view
      if (InputNames.size() > 1) {
         BCFieldName  = InputNames[1];
         auto BCField = Field::get(BCFieldName);
         BCData       = BCField->template getDataArray<Array1DReal>();
      }

      // Read required dimension parameter from configuration
      Error Err = Options.get("Dimension", ScanDimension);
      if (Err.isFail()) {
         ABORT_ERROR("PrefixSumOp: Required parameter 'Dimension' not found in "
                     "configuration");
      }

      // Read optional reverse parameter (default: false = forward scan)
      Err = Options.get("Reverse", ReverseDirection);
      if (Err.isFail()) {
         ReverseDirection = false; // Default to forward scan
      }

      // Retrieve input Field to get dimensions and metadata
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Validate dimension parameter
      if (ScanDimension < 0 || ScanDimension >= static_cast<I4>(ArrayT::rank)) {
         ABORT_ERROR("PrefixSumOp: Dimension {} out of range for array rank {}",
                     ScanDimension, static_cast<int>(ArrayT::rank));
      }

      if constexpr (ArrayT::rank == 1) {
         if (ScanDimension != 0) {
            ABORT_ERROR(
                "PrefixSumOp: For 1D arrays, Dimension must be 0, got {}",
                ScanDimension);
         }
      }

      // Get dimension info
      auto NDims = InputField->getNumDims();
      std::vector<std::string> DimNames;
      InputField->getDimNames(DimNames);

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_PrefixSum";
      if (InputNames.size() > 1) {
         OutputFieldName += "(BC=" + BCFieldName + ")";
      }
      OutputNames  = {OutputFieldName};
      InstanceName = OutputFieldName;

      // Get input metadata. A cumulative sum has the units of the values
      // summed.
      std::string InputDescr;
      ScalarT InputValidMin, InputValidMax;
      InputField->getMetadata("long_name", InputDescr);
      InputField->getMetadata("valid_min", InputValidMin);
      InputField->getMetadata("valid_max", InputValidMax);
      auto Meta = inheritMetadata(InputNames[0]);

      // The BC field seeds the sum, so it must be in the same units as the
      // input. If either has unknown units, so does the sum.
      if (!BCFieldName.empty()) {
         CFUnits InputUnits, BCUnits;
         Error Err = CFUnits::parse(Meta.Units, InputUnits);
         CHECK_ERROR_ABORT(Err, "PrefixSumOp: invalid units of field {}",
                           InputNames[0]);
         Err = CFUnits::parse(getFieldAttribute(BCFieldName, "units"), BCUnits);
         CHECK_ERROR_ABORT(Err, "PrefixSumOp: invalid units of BC field {}",
                           BCFieldName);
         if (InputUnits.isUnknown() || BCUnits.isUnknown()) {
            Meta.Units = "";
         } else if (InputUnits != BCUnits) {
            ABORT_ERROR("PrefixSumOp: BC field {} has units '{}' but field {} "
                        "has units '{}'",
                        BCFieldName, BCUnits.str(), InputNames[0],
                        InputUnits.str());
         }
      }

      // Allow callers to override the output field name (e.g., for a short
      // alias) via the "OutputName" config key.
      applyOutputNameOverride(Options);

      // Create output Field with same dimensions as input
      auto OutputField =
          createOutputField(OutputNames[0],
                            "Cumulative sum of " + InputDescr, // Description
                            Meta,                              // CF metadata
                            InputValidMin,                     // Min valid
                            InputValidMax,                     // Max valid
                            NDims,                             // Rank
                            DimNames                           // Dim names
          );

      // Allocate output data array matching input layout
      OutputData = OutputArrayT(OutputNames[0] + "_out", InputData.layout());

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

   } // end constructor

   /// Initializes the operator after all Fields exist. For 2D arrays,
   /// examines the first dimension name to determine if it's a mesh dimension
   /// (NCells, NEdges, NVertices) or a non-mesh dimension (e.g., bins). For
   /// mesh dimensions, stores the appropriate MinLayer/MaxLayer arrays from
   /// VertCoord for bounding the inner vertical loop when scanning along the
   /// vertical dimension (dimension 1). For non-mesh dimensions (bins), no
   /// MinLayer/MaxLayer bounds are needed. The dimension name is used in
   /// compute() to determine whether to use owned count (mesh dimensions) or
   /// full extent (non-mesh dimensions) for MPI-correct parallel loops.
   void initialize(const MachEnv *Env, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {

      AnalysisOperator::initialize(Env, InMesh, InVCoord, Options);

      // Set Dim0Name for 1D arrays (no MinLayer/MaxLayer needed)
      if constexpr (ArrayT::rank == 1) {
         auto Field1 = Field::get(InputNames[0]);
         std::vector<std::string> DimNames;
         Field1->getDimNames(DimNames);
         Dim0Name        = DimNames[0];
         IsMeshDimension = false;
      }

      // Only need MinLayer/MaxLayer and dimension info for 2D arrays
      if constexpr (ArrayT::rank == 2) {
         auto Field1 = Field::get(InputNames[0]);
         std::vector<std::string> DimNames;
         Field1->getDimNames(DimNames);
         // First dimension can be mesh (NCells/NEdges/NVertices) or non-mesh
         // (bins)
         Dim0Name = DimNames[0];

         if (Dim0Name == "NCells") {
            MinLayer        = VCoord->MinLayerCell;
            MaxLayer        = VCoord->MaxLayerCell;
            IsMeshDimension = true;
         } else if (Dim0Name == "NEdges") {
            MinLayer        = VCoord->MinLayerEdgeBot;
            MaxLayer        = VCoord->MaxLayerEdgeTop;
            IsMeshDimension = true;
         } else if (Dim0Name == "NVertices") {
            MinLayer        = VCoord->MinLayerVertexBot;
            MaxLayer        = VCoord->MaxLayerVertexTop;
            IsMeshDimension = true;
         } else {
            // For non-mesh dimensions (e.g., bins), MinLayer/MaxLayer remain
            // uninitialized
            IsMeshDimension = false;
         }
      }

   } // end initialize

   /// Computes the cumulative sum along the specified dimension. For 2D
   /// arrays, uses hierarchical parallelism with parallelForOuter and
   /// parallelScanInner. Supports both forward (start to end) and reverse
   /// (end to start) directions. If a BC field is registered, seeds the
   /// accumulator at scan index 0 with the BC values before processing input.
   /// Updates output data, timestamp, and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Retrieve input Field and extract data array
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Compute cumulative sum based on array rank using constexpr to prevent
      // instantiation of incorrect branches at compile time
      if constexpr (ArrayT::rank == 1) {
         compute1D(InputData);
      } else if constexpr (ArrayT::rank == 2) {
         compute2D(InputData);
      } else {
         ABORT_ERROR(
             "PrefixSumOp: Only 1D and 2D arrays are supported. Rank: {}",
             static_cast<int>(ArrayT::rank));
      }

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

   /// Computes cumulative sum for 1D arrays using Kokkos::parallel_scan.
   /// Supports both forward (start to end) and reverse (end to start)
   /// directions over the full extent of the array.
   void compute1D(const ArrayT &InputData) {
      const I4 N = InputData.extent(0);
      OMEGA_SCOPE(LocalOutput, OutputData);

      if (ReverseDirection) {
         Kokkos::parallel_scan(
             "PrefixSum_1D_Reverse", N,
             KOKKOS_LAMBDA(int IRev, Real &Accum, bool IsFinal) {
                const int I = N - 1 - IRev;
                Accum += static_cast<Real>(InputData(I));
                if (IsFinal) {
                   LocalOutput(I) = static_cast<ScalarT>(Accum);
                }
             });
      } else {
         Kokkos::parallel_scan(
             "PrefixSum_1D_Forward", N,
             KOKKOS_LAMBDA(int I, Real &Accum, bool IsFinal) {
                Accum += static_cast<Real>(InputData(I));
                if (IsFinal) {
                   LocalOutput(I) = static_cast<ScalarT>(Accum);
                }
             });
      }
   } // end compute1D

   /// Computes cumulative sum for 2D arrays using hierarchical parallelism.
   /// Uses parallelForOuter for the non-scan dimension and parallelScanInner
   /// for the scan dimension. For vertical scans (dimension 1), bounds the
   /// loop with MinLayer/MaxLayer to handle partial columns. Horizontal scans
   /// (dimension 0) are only allowed for non-distributed dimensions like bins.
   /// When UseBC is true, the BCData array is added to the accumulator before
   /// processing the first element along the scan dimension (index 0), seeding
   /// the running sum with the boundary condition value.
   void compute2D(const ArrayT &InputData) {

      // Check if scanning along distributed dimension (not allowed)
      if (ScanDimension == 0 && (Dim0Name == "NCells" || Dim0Name == "NEdges" ||
                                 Dim0Name == "NVertices")) {
         ABORT_ERROR("PrefixSumOp: Cannot scan along distributed dimension {}. "
                     "Horizontal scans only supported for non-distributed "
                     "dimensions (e.g., bins).",
                     Dim0Name);
      }

      // Determine loop extent for dimension 0
      // For mesh dimensions, use owned count; for non-mesh dimensions, use full
      // extent
      I4 N0;
      if (IsMeshDimension) {
         if (Dim0Name == "NCells") {
            N0 = Mesh->NCellsOwned;
         } else if (Dim0Name == "NEdges") {
            N0 = Mesh->NEdgesOwned;
         } else if (Dim0Name == "NVertices") {
            N0 = Mesh->NVerticesOwned;
         } else {
            N0 = InputData.extent(0); // fallback
         }
      } else {
         N0 = InputData.extent(
             0); // bins or other non-distributed dims use full extent
      }
      const I4 N1 = InputData.extent(1);

      OMEGA_SCOPE(LocOutput, OutputData);

      if (ScanDimension == 0) {
         // Scan along first dimension (horizontal bins)
         // Outer loop over second dimension (vertical), inner scan over first
         OMEGA_SCOPE(LocBC, BCData);
         if (ReverseDirection) {
            // Reverse: scan from end to start
            // BC is seeded at the last index (N0-1), i.e. iRev==0 maps to
            // i=N0-1
            parallelForOuter(
                "PrefixSum_Dim0_Reverse", LaunchConfig({N1}),
                KOKKOS_LAMBDA(int J, const TeamMember &Team) {
                   parallelScanInner(
                       Team, N0,
                       INNER_LAMBDA(int IRev, Real &Accum, bool IsFinal) {
                          const int I = N0 - 1 - IRev;
                          if (IRev == 0 && LocBC.is_allocated() &&
                              J < LocBC.extent(0))
                             Accum += static_cast<Real>(LocBC(J));
                          Accum += static_cast<Real>(InputData(I, J));
                          if (IsFinal) {
                             LocOutput(I, J) = static_cast<ScalarT>(Accum);
                          }
                       });
                });
         } else {
            // Forward: scan from start to end
            // BC is seeded at index 0 before the first input element
            parallelForOuter(
                "PrefixSum_Dim0_Forward", LaunchConfig({N1}),
                KOKKOS_LAMBDA(int J, const TeamMember &Team) {
                   parallelScanInner(
                       Team, N0,
                       INNER_LAMBDA(int I, Real &Accum, bool IsFinal) {
                          if (I == 0 && LocBC.is_allocated() &&
                              J < LocBC.extent(0))
                             Accum += static_cast<Real>(LocBC(J));
                          Accum += static_cast<Real>(InputData(I, J));
                          if (IsFinal) {
                             LocOutput(I, J) = static_cast<ScalarT>(Accum);
                          }
                       });
                });
         }
      } else { // ScanDimension == 1
         // Scan along second dimension (vertical)
         // Outer loop over first dimension (bins/etc), inner scan over vertical

         if (IsMeshDimension) {
            // For mesh dimensions, use MinLayer/MaxLayer to bound vertical
            // loops for partial columns
            OMEGA_SCOPE(LocMinLayer, MinLayer);
            OMEGA_SCOPE(LocMaxLayer, MaxLayer);

            if (ReverseDirection) {
               // Reverse: scan from end to start
               // Iterate from KMax down to KMin, accumulating as we go
               // Result at K = sum from K to KMax
               parallelForOuter(
                   "PrefixSum_Dim1_Reverse", LaunchConfig({N0}),
                   KOKKOS_LAMBDA(int I, const TeamMember &Team) {
                      const I4 KMin   = LocMinLayer(I);
                      const I4 KMax   = LocMaxLayer(I);
                      const I4 KRange = vertRange(KMin, KMax);
                      parallelScanInner(
                          Team, KRange,
                          INNER_LAMBDA(int KIdx, Real &Accum, bool IsFinal) {
                             // Iterate from KMax down: KIdx 0→KRange-1 maps to
                             // K KMax→KMin
                             const I4 K = KMax - KIdx;
                             Accum += static_cast<Real>(InputData(I, K));
                             if (IsFinal) {
                                LocOutput(I, K) = static_cast<ScalarT>(Accum);
                             }
                          });
                   });
            } else {
               // Forward: scan from start to end (top to bottom for vertical)
               parallelForOuter(
                   "PrefixSum_Dim1_Forward", LaunchConfig({N0}),
                   KOKKOS_LAMBDA(int I, const TeamMember &Team) {
                      const I4 KMin   = LocMinLayer(I);
                      const I4 KMax   = LocMaxLayer(I);
                      const I4 KRange = vertRange(KMin, KMax);
                      parallelScanInner(
                          Team, KRange,
                          INNER_LAMBDA(int KIdx, Real &Accum, bool IsFinal) {
                             // Map KIdx to actual layer: 0→KMin, KRange-1→KMax
                             const I4 K = KMin + KIdx;
                             Accum += static_cast<Real>(InputData(I, K));
                             if (IsFinal) {
                                LocOutput(I, K) = static_cast<ScalarT>(Accum);
                             }
                          });
                   });
            }
         } else {
            // For non-mesh dimensions (e.g., bins), scan over full extent
            if (ReverseDirection) {
               // Reverse: scan from end to start
               parallelForOuter(
                   "PrefixSum_Dim1_Reverse_Full", LaunchConfig({N0}),
                   KOKKOS_LAMBDA(int I, const TeamMember &Team) {
                      parallelScanInner(
                          Team, N1,
                          INNER_LAMBDA(int KIdx, Real &Accum, bool IsFinal) {
                             // Iterate from N1-1 down to 0
                             const I4 K = N1 - 1 - KIdx;
                             Accum += static_cast<Real>(InputData(I, K));
                             if (IsFinal) {
                                LocOutput(I, K) = static_cast<ScalarT>(Accum);
                             }
                          });
                   });
            } else {
               // Forward: scan from start to end
               parallelForOuter(
                   "PrefixSum_Dim1_Forward_Full", LaunchConfig({N0}),
                   KOKKOS_LAMBDA(int I, const TeamMember &Team) {
                      parallelScanInner(
                          Team, N1,
                          INNER_LAMBDA(int KIdx, Real &Accum, bool IsFinal) {
                             // KIdx maps directly to K
                             const I4 K = KIdx;
                             Accum += static_cast<Real>(InputData(I, K));
                             if (IsFinal) {
                                LocOutput(I, K) = static_cast<ScalarT>(Accum);
                             }
                          });
                   });
            }
         }
      }

   } // end compute2D

 private:
   /// Output data array holding the cumulative sum
   OutputArrayT OutputData;

   /// The dimension along which to perform cumulative sum (0 or 1 for 2D)
   I4 ScanDimension;

   /// Whether to scan in reverse direction (true) or forward (false)
   bool ReverseDirection;

   /// Name of the optional boundary condition field (empty = no BC)
   std::string BCFieldName = "";

   /// Boundary condition view, fetched in constructor (unallocated = no BC)
   Array1DReal BCData;

   /// Min active layer index for each horizontal point (only for 2D arrays)
   Array1DI4 MinLayer;

   /// Max active layer index for each horizontal point (only for 2D arrays)
   Array1DI4 MaxLayer;

   /// Name of the first dimension (only for 2D arrays) to determine owned count
   std::string Dim0Name;

   /// Flag indicating whether the first dimension is a mesh dimension (NCells,
   /// NEdges, NVertices) If true, MinLayer/MaxLayer are used; if false, full
   /// extent is used
   bool IsMeshDimension = false;

}; // end class PrefixSumOp

} // end namespace OMEGA

#endif
