#ifndef OMEGA_TRANSECTACCUMULATOROP_H
#define OMEGA_TRANSECTACCUMULATOROP_H

//===-- analysis/operators/TransectAccumulatorOp.h --------------*- C++ -*-===//
//
/// \file
/// \brief Defines the TransectAccumulatorOp operator for transect accumulation
///
/// TransectAccumulatorOp accumulates edge-based transport values across a
/// defined transect. This operator is essential for transect-based MOC
/// (Meridional Overturning Circulation) calculations, where it accumulates
/// edge normal velocities multiplied by edge width and layer thickness across
/// a transect (a collection of edges). The operator performs local accumulation
/// followed by MPI reduction across all ranks to produce global transect
/// totals.
///
/// The operator takes three inputs:
/// 1. A transport field to accumulate (e.g., normalVelocity × dvEdge ×
/// layerThickness)
/// 2. A transect edge mask field (1 = on transect, 0 = not on transect)
/// 3. A transect edge sign field (+1 or -1 for direction)
///
/// The operator outputs a 1D array of accumulated transport values at each
/// vertical level for the specified transect.
///
/// Configuration:
/// - TransectName: Name of the transect for field naming (required)
///
/// Example usage in operator chain:
/// \code
///   TransportField_TransectAccumulator(TransectName)
/// \endcode
/// where TransectAccumulator(TransectName) sums transport values across
/// transect edges for MOC.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Reductions.h"
#include "mpi.h"

namespace OMEGA {

/// TransectAccumulatorOp accumulates edge-based transport field values across
/// a transect. Supports 2D input fields (NEdges × NVertLayers) and produces
/// 1D output (NVertLayers). Includes transect masking and directional signing,
/// and performs MPI reduction to compute global transect totals. This is a core
/// operator for transect-based MOC streamfunction calculations.
template <typename ArrayT>
class TransectAccumulatorOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type is always 1D (NVertLayers) for transect accumulation
   using OutputArrayT = Array1DReal;

   /// Constructs a TransectAccumulatorOp operator. Reads configuration
   /// (transect name), retrieves mask and sign fields, creates output Field
   /// for transect accumulation, allocates output and temporary arrays, and
   /// registers the output Field. The default output Field name is
   /// "InputField_TransectAccumulator(TransectName)".
   TransectAccumulatorOp(const std::vector<std::string>
                             &UpstreamNames, ///< [in] input field names
                         Config Options      ///< [in] operator config
                         )
       : AnalysisOperator("TransectAccumulator") {

      // Store input field names
      // UpstreamNames[0] = transport field to accumulate
      // UpstreamNames[1] = transect edge mask field
      // UpstreamNames[2] = transect edge sign field
      InputNames = UpstreamNames;

      if (InputNames.size() < 3) {
         ABORT_ERROR("TransectAccumulatorOp: Requires 3 inputs (transport "
                     "field, mask field, and sign field)");
      }

      // Read required TransectName parameter
      Error Err = Options.get("TransectName", TransectName);
      if (Err.isFail()) {
         ABORT_ERROR("TransectAccumulatorOp: Required parameter 'TransectName' "
                     "not found in configuration");
      }

      if (TransectName.empty()) {
         ABORT_ERROR("TransectAccumulatorOp: TransectName cannot be empty");
      }

      // Retrieve input transport Field to get dimensions
      auto TransportField = Field::get(InputNames[0]);
      auto TransportData  = TransportField->template getDataArray<ArrayT>();

      // Validate input rank (must be 2D: NEdges × NVertLayers)
      const I4 InputRank = ArrayT::rank;
      if (InputRank != 2) {
         ABORT_ERROR("TransectAccumulatorOp: Input must be 2D (NEdges × "
                     "NVertLayers), got rank {}",
                     InputRank);
      }

      // Get dimensions
      const I4 NEdges      = TransportData.extent(0);
      const I4 NVertLayers = TransportData.extent(1);

      // Store dimensions
      HorzSize = NEdges;
      VertSize = NVertLayers;

      // Get dimension info from transport field and validate horizontal
      // dimension
      std::vector<std::string> TransportDimNames;
      TransportField->getDimNames(TransportDimNames);

      if (TransportDimNames[0] != "NEdges") {
         ABORT_ERROR("TransectAccumulatorOp: Transport field horizontal "
                     "dimension must be 'NEdges', got '{}'",
                     TransportDimNames[0]);
      }

      // Validate mask field horizontal dimension
      auto MaskFieldCheck = Field::get(InputNames[1]);
      std::vector<std::string> MaskDimNames;
      MaskFieldCheck->getDimNames(MaskDimNames);
      if (MaskDimNames[0] != "NEdges") {
         ABORT_ERROR("TransectAccumulatorOp: Mask field horizontal dimension "
                     "must be 'NEdges', got '{}'",
                     MaskDimNames[0]);
      }

      // Validate sign field horizontal dimension
      auto SignFieldCheck = Field::get(InputNames[2]);
      std::vector<std::string> SignDimNames;
      SignFieldCheck->getDimNames(SignDimNames);
      if (SignDimNames[0] != "NEdges") {
         ABORT_ERROR("TransectAccumulatorOp: Sign field horizontal dimension "
                     "must be 'NEdges', got '{}'",
                     SignDimNames[0]);
      }

      // Construct output field name and set instance name.
      std::string OutputFieldName =
          InputNames[0] + "_TransectAccumulator(" + TransectName + ")";
      OutputNames  = {OutputFieldName};
      InstanceName = OutputFieldName;

      // Get input metadata. A signed sum of the transport through the
      // transect edges has the units of the transport.
      std::string TransportDescr;
      ScalarT TransportValidMin, TransportValidMax;
      TransportField->getMetadata("long_name", TransportDescr);
      TransportField->getMetadata("valid_min", TransportValidMin);
      TransportField->getMetadata("valid_max", TransportValidMax);
      auto Meta = inheritMetadata(InputNames[0]);

      // Create output Field dimensions (1D: NVertLayers only)
      std::vector<std::string> OutputDimNames = {TransportDimNames[1]};

      // Create output Field
      auto OutputField = createOutputField(
          OutputNames[0],
          "Transect " + TransectName + " accumulated transport", // Description
          Meta,                                                  // CF metadata
          TransportValidMin,                                     // Min valid
          TransportValidMax,                                     // Max valid
          1,                                                     // Rank (1D)
          OutputDimNames                                         // Dim names
      );

      // Allocate output and local accumulation arrays (both 1D)
      OutputData = OutputArrayT("TransectAccum_out", NVertLayers);
      LocalAccum = OutputArrayT("TransectAccum_local", NVertLayers);

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

   } // end constructor

   /// Initializes the operator after Fields are created. Stores mesh pointer
   /// and determines the number of owned edges to avoid double-counting
   /// halo edges in MPI reductions.
   void initialize(const MachEnv *Env, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {

      // Call base class initialization to store Mesh, VCoord, and Comm
      AnalysisOperator::initialize(Env, InMesh, InVCoord, Options);

      // Transect accumulation operates on edges
      NOwned   = Mesh->NEdgesOwned;
      MinLayer = VCoord->MinLayerEdgeBot;
      MaxLayer = VCoord->MaxLayerEdgeTop;

   } // end initialize

   /// Computes the transect accumulation by:
   /// 1. Zeroing local accumulation array
   /// 2. Iterating over all owned edges and accumulating transport values
   ///    multiplied by mask and sign into vertical levels
   /// 3. Performing MPI_Allreduce to sum local accumulations across all ranks
   /// 4. Storing global totals in output array
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Retrieve input fields
      auto TransportField = Field::get(InputNames[0]);
      auto TransportData  = TransportField->template getDataArray<ArrayT>();
      auto MaskField      = Field::get(InputNames[1]);
      auto MaskData       = MaskField->template getDataArray<Array1DI4>();
      auto SignField      = Field::get(InputNames[2]);
      auto SignData       = SignField->template getDataArray<Array1DI4>();

      // Compute transect accumulation
      computeTransectAccum(TransportData, MaskData, SignData);

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

   /// Computes transect accumulation for 2D edge-based input using hierarchical
   /// parallelism. Outer loop over owned edges, inner loop over active vertical
   /// layers bounded by MinLayerEdgeBot/MaxLayerEdgeTop.
   void computeTransectAccum(const ArrayT &TransportData,
                             const Array1DI4 &MaskData,
                             const Array1DI4 &SignData) {

      auto LocalAccumData = LocalAccum;
      auto LocalNOwned    = NOwned;
      auto LocalVertSize  = VertSize;

      // Zero local accumulation array
      parallelFor(
          {LocalVertSize}, KOKKOS_LAMBDA(int k) { LocalAccumData(k) = 0.0; });

      // Accumulate transport across transect edges using hierarchical
      // parallelism Outer loop: owned edges; inner loop: active vertical layers
      OMEGA_SCOPE(LocMinLayer, MinLayer);
      OMEGA_SCOPE(LocMaxLayer, MaxLayer);
      parallelForOuter(
          "TransectAccum2D", LaunchConfig({LocalNOwned}),
          KOKKOS_LAMBDA(int iEdge, const TeamMember &Team) {
             const I4 Mask = MaskData(iEdge);
             const I4 Sign = SignData(iEdge);

             if (Mask != 0) {
                const I4 KMin   = LocMinLayer(iEdge);
                const I4 KMax   = LocMaxLayer(iEdge);
                const I4 KRange = vertRange(KMin, KMax);
                parallelForInner(
                    Team, KRange, INNER_LAMBDA(int KIdx) {
                       const I4 K = KMin + KIdx;
                       Real Transport =
                           static_cast<Real>(TransportData(iEdge, K));
                       Real Contribution = Transport * static_cast<Real>(Mask) *
                                           static_cast<Real>(Sign);
                       Kokkos::atomic_add(&LocalAccumData(K), Contribution);
                    });
             }
          });

      // Perform MPI reduction to get global totals
      // Copy local accumulation from device to host
      auto LocalAccumHost = Kokkos::create_mirror_view(LocalAccum);
      deepCopy(LocalAccumHost, LocalAccum);

      // Create host mirror for output
      auto OutputDataHost = Kokkos::create_mirror_view(OutputData);

      // Sum local accumulations across all ranks element-wise. Use the
      // reproducible Omega globalSum (multifield) helper, which performs an
      // element-wise reduction across ranks using the stored communicator and
      // avoids a direct MPI_Allreduce call here.
      std::vector<R8> LocalLevels(VertSize);
      for (I4 K = 0; K < VertSize; ++K)
         LocalLevels[K] = static_cast<R8>(LocalAccumHost(K));

      std::vector<R8> GlobalLevels = globalSum(LocalLevels, Comm);

      for (I4 K = 0; K < VertSize; ++K)
         OutputDataHost(K) = GlobalLevels[K];

      // Copy global totals back to device
      deepCopy(OutputData, OutputDataHost);

   } // end computeTransectAccum

 private:
   /// Output data array holding global transect accumulation (1D: NVertLayers)
   OutputArrayT OutputData;

   /// Local accumulation array before MPI reduction (1D: NVertLayers)
   OutputArrayT LocalAccum;

   /// Name of the transect
   std::string TransectName;

   /// Horizontal dimension size (NEdges)
   I4 HorzSize;

   /// Number of owned edges (excludes halo edges for MPI)
   I4 NOwned;

   /// Min active layer index for each edge (MinLayerEdgeBot)
   Array1DI4 MinLayer;

   /// Max active layer index for each edge (MaxLayerEdgeTop)
   Array1DI4 MaxLayer;

   /// Vertical dimension size (NVertLayers)
   I4 VertSize;

}; // end class TransectAccumulatorOp

} // end namespace OMEGA

#endif
