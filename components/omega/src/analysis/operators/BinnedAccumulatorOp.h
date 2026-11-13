#ifndef OMEGA_BINNEDACCUMULATOROP_H
#define OMEGA_BINNEDACCUMULATOROP_H

//===-- analysis/operators/BinnedAccumulatorOp.h ----------------*- C++ -*-===//
//
/// \file
/// \brief Defines the BinnedAccumulatorOp operator for binned accumulation
///
/// BinnedAccumulatorOp accumulates field values into spatial bins. This
/// operator is essential for MOC (Meridional Overturning Circulation)
/// calculations, where it accumulates edge normal velocities into latitude
/// bins. The operator performs local accumulation followed by MPI reduction
/// across all ranks to produce global bin totals.
///
/// The operator takes two inputs:
/// 1. A value field to accumulate (1D or 2D)
/// 2. A bin index field (from CoordinateBinningOp) indicating which bin each
///    entity belongs to
///
/// If the value field has a Field-attached regional mask it is applied
/// automatically during accumulation (1.0 = include, 0.0 = exclude).
///
/// The output replaces the horizontal dimension with NumBins. For 2D input
/// (NCells x NVertLevels) the output is 2D (NumBins x NVertLevels); for 1D
/// input the output is 1D (NumBins).
///
/// Configuration:
/// - NumBins: Number of spatial bins (required)
///
/// Example usage in operator chain:
/// \code
///   NormalVelocity_BinnedAccumulator(LatCell_BinIndex)
/// \endcode
/// where LatCell_BinIndex is the bin index field produced by
/// CoordinateBinningOp.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Reductions.h"
#include "mpi.h"

namespace OMEGA {

/// BinnedAccumulatorOp accumulates field values into spatial bins. Supports
/// 1D and 2D input fields and produces output with the horizontal dimension
/// replaced by NumBins. If the input Field has a Field-attached regional mask
/// it is applied automatically. Uses hierarchical parallelism for 2D input.
/// Performs MPI reduction to compute global bin totals. This is the core
/// operator for MOC streamfunction calculations.
template <typename ArrayT> class BinnedAccumulatorOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type depends on input rank:
   /// For 2D input (e.g., NCells x NVertLevels), output is 2D (NumBins x
   /// NVertLevels) For 1D input (e.g., NCells), output is 1D (NumBins)
   using OutputArrayT =
       std::conditional_t<ArrayT::rank == 2, Array2DR8, Array1DR8>;

   /// Constructs a BinnedAccumulatorOp operator. Reads the required NumBins
   /// configuration parameter, creates the output Field for binned
   /// accumulation, allocates output and local accumulation arrays, and
   /// registers the output Field. The output Field name is constructed as
   /// InputName + "_BinnedAccum".
   BinnedAccumulatorOp(const std::vector<std::string>
                           &UpstreamNames, ///< [in] input field names
                       Config Options      ///< [in] operator config
                       )
       : AnalysisOperator("BinnedAccumulator") {

      // Store input field names:
      //   UpstreamNames[0] = value field to accumulate
      //   UpstreamNames[1] = bin index field (from CoordinateBinningOp)
      InputNames = UpstreamNames;

      if (InputNames.size() < 2) {
         ABORT_ERROR("BinnedAccumulatorOp: Requires at least 2 inputs (value "
                     "field and bin index field)");
      }

      // Read required NumBins parameter
      Error Err = Options.get("NumBins", NumBins);
      if (Err.isFail()) {
         ABORT_ERROR("BinnedAccumulatorOp: Required parameter 'NumBins' not "
                     "found in configuration");
      }

      if (NumBins <= 0) {
         ABORT_ERROR("BinnedAccumulatorOp: NumBins must be positive, got {}",
                     NumBins);
      }

      // Retrieve input value Field to get dimensions
      auto ValueField = Field::get(InputNames[0]);
      auto ValueData  = ValueField->template getDataArray<ArrayT>();

      // Validate input rank (1D or 2D)
      const I4 InputRank = ArrayT::rank;
      if (InputRank != 1 && InputRank != 2) {
         ABORT_ERROR("BinnedAccumulatorOp: Input must be 1D or 2D, got rank {}",
                     InputRank);
      }

      // Get dimensions
      const I4 Dim0 = ValueData.extent(0);
      I4 Dim1       = 1;
      if (InputRank == 2) {
         Dim1 = ValueData.extent(1);
      }

      // Store dimensions
      HorzSize = Dim0;
      VertSize = Dim1;

      // Get dimension info from value field
      auto NDims = ValueField->getNumDims();
      std::vector<std::string> ValueDimNames;
      ValueField->getDimNames(ValueDimNames);

      // Construct output field name and set instance name
      // Format: Input_BinnedAccumulator(BinIndexField) to match chain syntax
      std::string OutputFieldName =
          InputNames[0] + "_BinnedAccumulator(" + InputNames[1] + ")";
      OutputNames  = {OutputFieldName};
      InstanceName = OutputFieldName;

      // Get input metadata. A sum over the entities in each bin has the
      // units of the values summed.
      std::string ValueDescr;
      ScalarT ValueValidMin, ValueValidMax;
      ValueField->getMetadata("long_name", ValueDescr);
      ValueField->getMetadata("valid_min", ValueValidMin);
      ValueField->getMetadata("valid_max", ValueValidMax);
      auto Meta = inheritMetadata(InputNames[0]);

      // Create output Field dimensions
      std::string NumBinsDimName = "NumBins" + InputNames[1];
      auto NumBinsDim            = Dimension::create(NumBinsDimName, NumBins);

      // Replace horizontal dimension with bin dimension
      std::vector<std::string> OutputDimNames;
      if (InputRank == 1) {
         OutputDimNames = {NumBinsDimName};
      } else {
         OutputDimNames = {NumBinsDimName, ValueDimNames[1]};
      }

      // Create output Field
      auto OutputField = createOutputField(
          OutputNames[0],
          "Binned accumulation of " + ValueDescr, // Description
          Meta,                                   // CF metadata
          ValueValidMin,                          // Min valid
          ValueValidMax,                          // Max valid
          InputRank,                              // Rank
          OutputDimNames                          // Dimension names
      );

      // Allocate output and local accumulation arrays
      if (InputRank == 1) {
         OutputData = OutputArrayT("BinnedAccum_out", NumBins);
         LocalAccum = OutputArrayT("BinnedAccum_local", NumBins);
      } else {
         OutputData = OutputArrayT("BinnedAccum_out", NumBins, Dim1);
         LocalAccum = OutputArrayT("BinnedAccum_local", NumBins, Dim1);
      }

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

   } // end constructor

   /// Initializes the operator after Fields are created. Stores mesh pointer,
   /// determines the number of owned entities to avoid double-counting halo
   /// cells in MPI reductions, and stores MinLayer/MaxLayer arrays from
   /// VertCoord for bounding the inner vertical loop.
   void initialize(const MachEnv *Env, const HorzMesh *InMesh,
                   const VertCoord *InVCoord, Config Options) override {

      // Call base class initialization to store Mesh, VCoord, and Comm
      AnalysisOperator::initialize(Env, InMesh, InVCoord, Options);

      // Determine index space, set number of owned entities, and
      // MinLayer/MaxLayer
      auto ValueField = Field::get(InputNames[0]);
      std::vector<std::string> DimNames;
      ValueField->getDimNames(DimNames);
      std::string IndexSpaceName = DimNames[0];

      if (IndexSpaceName == "NCells") {
         NOwned   = Mesh->NCellsOwned;
         MinLayer = VCoord->MinLayerCell;
         MaxLayer = VCoord->MaxLayerCell;
      } else if (IndexSpaceName == "NEdges") {
         NOwned   = Mesh->NEdgesOwned;
         MinLayer = VCoord->MinLayerEdgeBot;
         MaxLayer = VCoord->MaxLayerEdgeTop;
      } else if (IndexSpaceName == "NVertices") {
         NOwned   = Mesh->NVerticesOwned;
         MinLayer = VCoord->MinLayerVertexBot;
         MaxLayer = VCoord->MaxLayerVertexTop;
      } else {
         ABORT_ERROR("BinnedAccumulatorOp: Unknown index space {}",
                     IndexSpaceName);
      }

   } // end initialize

   /// Computes the binned accumulation by:
   /// 1. Zeroing local accumulation array
   /// 2. Accumulating values into bins for owned entities only (excludes halo),
   ///    applying the Field-attached regional mask if present
   /// 3. Performing MPI_Allreduce to sum local accumulations across all ranks
   /// 4. Storing global totals in the output array
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Retrieve input fields
      auto ValueField    = Field::get(InputNames[0]);
      auto ValueData     = ValueField->template getDataArray<ArrayT>();
      auto BinIndexField = Field::get(InputNames[1]);
      auto BinIndexData  = BinIndexField->template getDataArray<Array1DI4>();

      // Apply Field-attached regional mask if present
      Array1DI4 MaskData;
      bool LocalUseMask = false;

      if (ValueField->hasRegionalMask()) {
         MaskData     = ValueField->getRegionalMask();
         LocalUseMask = true;
      }

      // Dispatch to rank-specific implementation using constexpr to avoid
      // instantiating incorrect branches at compile time
      if constexpr (ArrayT::rank == 1) {
         computeAccum1D(ValueData, BinIndexData, MaskData, LocalUseMask);
      } else if constexpr (ArrayT::rank == 2) {
         computeAccum2D(ValueData, BinIndexData, MaskData, LocalUseMask);
      } else {
         ABORT_ERROR("BinnedAccumulatorOp: Unsupported rank {}. "
                     "Supports 1D and 2D arrays only.",
                     static_cast<int>(ArrayT::rank));
      }

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

   /// Computes binned accumulation for 1D input. Zeroes the local array,
   /// then atomically accumulates owned-entity values into bins, and
   /// performs MPI_Allreduce to obtain global totals.
   void computeAccum1D(const ArrayT &ValueData, const Array1DI4 &BinIndexData,
                       const Array1DI4 &MaskData, bool LocalUseMask) {

      auto LocalAccumData = LocalAccum;
      auto LocalNumBins   = NumBins;
      auto LocalNOwned    = NOwned;

      // Zero local accumulation array
      parallelFor(
          {LocalNumBins},
          KOKKOS_LAMBDA(int IBin) { LocalAccumData(IBin) = 0.0; });

      // Accumulate values into bins (only owned entities, not halo)
      // Atomic adds handle multiple threads writing to the same bin
      parallelFor(
          {LocalNOwned}, KOKKOS_LAMBDA(int I) {
             I4 BinIdx = BinIndexData(I);

             if (BinIdx >= 0 && BinIdx < LocalNumBins) {
                Real Value = static_cast<Real>(ValueData(I));

                if (LocalUseMask) {
                   Value *= MaskData(I);
                }

                Kokkos::atomic_add(&LocalAccumData(BinIdx), Value);
             }
          });

      // Copy local accumulation from device to host for MPI reduction
      auto LocalAccumHost = createHostMirrorCopy(LocalAccum);

      auto OutputDataHost = createHostMirrorCopy(OutputData);

      // Sum local accumulations across all ranks element-wise. Use the
      // reproducible Omega globalSum (multifield) helper, which performs an
      // element-wise reduction across ranks using the same communicator and
      // avoids a direct MPI_Allreduce call here.
      std::vector<R8> LocalBins(NumBins);
      for (I4 IBin = 0; IBin < NumBins; ++IBin)
         LocalBins[IBin] = static_cast<R8>(LocalAccumHost(IBin));

      std::vector<R8> GlobalBins = globalSum(LocalBins, Comm);

      for (I4 IBin = 0; IBin < NumBins; ++IBin)
         OutputDataHost(IBin) = GlobalBins[IBin];

      // Copy global totals back to device
      deepCopy(OutputData, OutputDataHost);

   } // end computeAccum1D

   /// Computes binned accumulation for 2D input using hierarchical parallelism:
   /// outer loop over owned horizontal entities, inner loop over vertical
   /// levels. Zeroes the local array, accumulates with atomic adds, then
   /// performs MPI_Allreduce to obtain global totals.
   void computeAccum2D(const ArrayT &ValueData, const Array1DI4 &BinIndexData,
                       const Array1DI4 &MaskData, bool LocalUseMask) {

      auto LocalAccumData = LocalAccum;
      auto LocalNumBins   = NumBins;
      auto LocalNOwned    = NOwned;
      auto LocalVertSize  = VertSize;

      // Zero local accumulation array
      parallelFor(
          {LocalNumBins, LocalVertSize},
          KOKKOS_LAMBDA(int IBin, int K) { LocalAccumData(IBin, K) = 0.0; });

      // Accumulate values into bins (only owned entities, not halo)
      // Hierarchical parallelism: outer over horizontal, inner over vertical
      // Inner loop bounded by MinLayer/MaxLayer for partial columns
      // Atomic adds handle multiple threads writing to the same bin
      OMEGA_SCOPE(LocMinLayer, MinLayer);
      OMEGA_SCOPE(LocMaxLayer, MaxLayer);
      parallelForOuter(
          "BinnedAccum2D", LaunchConfig({LocalNOwned}),
          KOKKOS_LAMBDA(int I, const TeamMember &Team) {
             I4 BinIdx = BinIndexData(I);

             if (BinIdx >= 0 && BinIdx < LocalNumBins) {
                const Real MaskVal =
                    LocalUseMask ? static_cast<Real>(MaskData(I)) : 1.0;
                const I4 KMin   = LocMinLayer(I);
                const I4 KMax   = LocMaxLayer(I);
                const I4 KRange = vertRange(KMin, KMax);
                parallelForInner(
                    Team, KRange, INNER_LAMBDA(int KIdx) {
                       const I4 K = KMin + KIdx;
                       Real Value =
                           static_cast<Real>(ValueData(I, K)) * MaskVal;
                       Kokkos::atomic_add(&LocalAccumData(BinIdx, K), Value);
                    });
             }
          });

      // Copy local accumulation from device to host for MPI reduction
      auto LocalAccumHost = createHostMirrorCopy(LocalAccum);

      auto OutputDataHost = createHostMirrorCopy(OutputData);

      // Sum local accumulations across all ranks element-wise. Flatten the
      // (NumBins x VertSize) local array into a single vector, reduce with the
      // reproducible Omega globalSum (multifield) helper, then unflatten into
      // the output host mirror.
      I4 TotalSize = NumBins * VertSize;
      std::vector<R8> LocalFlat(TotalSize);
      for (I4 IBin = 0; IBin < NumBins; ++IBin)
         for (I4 K = 0; K < VertSize; ++K)
            LocalFlat[IBin * VertSize + K] =
                static_cast<R8>(LocalAccumHost(IBin, K));

      std::vector<R8> GlobalFlat = globalSum(LocalFlat, Comm);

      for (I4 IBin = 0; IBin < NumBins; ++IBin)
         for (I4 K = 0; K < VertSize; ++K)
            OutputDataHost(IBin, K) = GlobalFlat[IBin * VertSize + K];

      // Copy global totals back to device
      deepCopy(OutputData, OutputDataHost);

   } // end computeAccum2D

 private:
   /// Output data array holding global binned accumulation
   OutputArrayT OutputData;

   /// Local accumulation array (before MPI reduction)
   OutputArrayT LocalAccum;

   /// Number of spatial bins
   I4 NumBins;

   /// Horizontal dimension size (NCells, NEdges, or NVertices)
   I4 HorzSize;

   /// Number of owned entities (excludes halo cells for MPI reduction)
   I4 NOwned;

   /// Vertical dimension size (NVertLevels or 1 for 1D input)
   I4 VertSize;

   /// Min active layer index for each horizontal point (for 2D input)
   Array1DI4 MinLayer;

   /// Max active layer index for each horizontal point (for 2D input)
   Array1DI4 MaxLayer;

}; // end class BinnedAccumulatorOp

} // end namespace OMEGA

#endif
