#ifndef OMEGA_COORDINATEBINNINGOP_H
#define OMEGA_COORDINATEBINNINGOP_H

//===-- analysis/operators/CoordinateBinningOp.h ----------------*- C++ -*-===//
//
//
/// \file
/// \brief Defines the CoordinateBinningOp operator for coordinate binning
///
/// CoordinateBinningOp assigns mesh cells (or edges/vertices) to bins
/// based on a coordinate field (e.g., latitude). This operator is particularly
/// useful for meridional overturning circulation (MOC) calculations, where
/// cells are binned by latitude for subsequent accumulation and integration.
///
/// The operator takes a 1D coordinate field (e.g., latCell) and outputs a 1D
/// integer array containing the bin index for each cell. Bin boundaries are
/// either auto-computed from the mesh extent or specified via configuration.
/// The binning computation occurs once during initialization and the results
/// are cached, as bin assignments are constant over time.
///
/// Configuration:
/// - NumBins: Number of bins (required)
/// - MinBin: Minimum coordinate value for binning (optional, auto-computed if
///           not specified)
/// - MaxBin: Maximum coordinate value for binning (optional, auto-computed if
///           not specified)
///
/// Example usage in operator chain:
/// \code
///   LatCell_CoordinateBinning
/// \endcode
/// where CoordinateBinning assigns cells to latitude bins for MOC calculation.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "Reductions.h"

namespace OMEGA {

/// CoordinateBinningOp assigns mesh entities to bins based on a
/// coordinate field. Bin boundaries are either auto-computed from mesh extent
/// or user-specified. The binning is computed once during initialization and
/// cached. Output is a 1D integer array of bin indices.
template <typename ArrayT> class CoordinateBinningOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output is always I4 (integer bin indices)
   using OutputArrayT = Array1DI4;

   /// Constructs a CoordinateBinningOp operator. Reads binning configuration
   /// (number of bins, optional min/max bounds), creates output Field for bin
   /// indices, allocates output data array, and registers the output Field.
   /// The actual binning computation is deferred to initialize() when mesh
   /// information is available. The output Field name is constructed as
   /// InputName + "_BinIndex".
   CoordinateBinningOp(const std::vector<std::string>
                           &UpstreamNames, ///< [in] input field names
                       Config Options      ///< [in] operator config
                       )
       : AnalysisOperator("CoordinateBinning") {

      // Store input field names
      InputNames = UpstreamNames;

      // Read required NumBins parameter
      Error Err = Options.get("NumBins", NumBins);
      if (Err.isFail()) {
         ABORT_ERROR("CoordinateBinningOp: Required parameter 'NumBins' not "
                     "found in configuration");
      }

      if (NumBins <= 0) {
         ABORT_ERROR("CoordinateBinningOp: NumBins must be positive, got {}",
                     NumBins);
      }

      // Read optional MinBin parameter
      Err = Options.get("MinBin", MinBin);
      if (Err.isFail()) {
         AutoComputeMin = true;
      } else {
         AutoComputeMin = false;
      }

      // Read optional MaxBin parameter
      Err = Options.get("MaxBin", MaxBin);
      if (Err.isFail()) {
         AutoComputeMax = true;
      } else {
         AutoComputeMax = false;
      }

      // Retrieve input Field to get dimensions
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Validate input is 1D
      if (ArrayT::rank != 1) {
         ABORT_ERROR(
             "CoordinateBinningOp: Input must be 1D (coordinate array), "
             "got rank {}",
             static_cast<int>(ArrayT::rank));
      }

      // Get dimension info
      std::vector<std::string> InputDimNames;
      InputField->getDimNames(InputDimNames);

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_BinIndex";
      OutputNames                 = {OutputFieldName};
      InstanceName                = OutputFieldName;

      // Create output Field with same horizontal dimension as input
      // but with I4 type for bin indices
      auto OutputField =
          Field::create(OutputNames[0],
                        "Bin index for " + InputNames[0], // Description
                        "1",          // Units (dimensionless)
                        "",           // Standard name
                        0,            // Min valid value
                        NumBins - 1,  // Max valid value
                        1,            // Rank
                        InputDimNames // Dimension names
          );

      // Store array size
      ArraySize = InputData.extent(0);

      // Allocate output data array (bin indices for each cell)
      OutputData = OutputArrayT(OutputNames[0] + "_out", ArraySize);

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

      // Mark binning as not yet computed
      BinningComputed = false;

   } // end constructor

   /// Initializes the operator by computing bin boundaries and assigning cells
   /// to bins. If MinBin/MaxBin are not specified, computes them from the
   /// global min/max of the coordinate field. Bin assignments are cached and
   /// reused in subsequent compute() calls.
   void initialize(const MachEnv *Env,      ///< [in] machine environment
                   const HorzMesh *Mesh,    ///< [in] horizontal mesh
                   const VertCoord *VCoord, ///< [in] vertical coordinate
                   Config Options           ///< [in] operator-specific options
                   ) override {

      // Call base class initialize
      AnalysisOperator::initialize(Env, Mesh, VCoord, Options);

      // Retrieve input coordinate field
      auto InputField = Field::get(InputNames[0]);
      auto InputData  = InputField->template getDataArray<ArrayT>();

      // Auto-compute MinBin and/or MaxBin if not specified
      if (AutoComputeMin || AutoComputeMax) {
         // Determine index space from field dimension name
         std::vector<std::string> InputDimNames;
         InputField->getDimNames(InputDimNames);
         std::string IndexSpaceName = InputDimNames[0];

         I4 NOwned = 0;
         if (IndexSpaceName == "NCells") {
            NOwned = Mesh->NCellsOwned;
         } else if (IndexSpaceName == "NEdges") {
            NOwned = Mesh->NEdgesOwned;
         } else if (IndexSpaceName == "NVertices") {
            NOwned = Mesh->NVerticesOwned;
         } else {
            ABORT_ERROR("CoordinateBinningOp: Unknown index space {}",
                        IndexSpaceName);
         }

         // Compute global min/max over owned entities
         std::vector<I4> IndxRange = {0, NOwned - 1};

         if (AutoComputeMin) {
            MinBin =
                static_cast<Real>(globalMinVal(InputData, Comm, &IndxRange));
         }

         if (AutoComputeMax) {
            MaxBin =
                static_cast<Real>(globalMaxVal(InputData, Comm, &IndxRange));
         }
      }

      // Add small margins to avoid boundary issues
      const Real Margin = (MaxBin - MinBin) * 1.0e-6;
      MinBin -= Margin;
      MaxBin += Margin;

      // Compute bin width
      BinWidth = (MaxBin - MinBin) / NumBins;

      if (BinWidth <= 0) {
         ABORT_ERROR("CoordinateBinningOp: Invalid bin width {}. MinBin={}, "
                     "MaxBin={}, NumBins={}",
                     BinWidth, MinBin, MaxBin, NumBins);
      }

      // Compute bin assignments for all entities
      computeBinning(InputData);

      BinningComputed = true;

   } // end initialize

   /// Computes bin assignments. Since binning is constant in time, this
   /// simply validates that initialization has occurred. The actual binning
   /// is done once in initialize() and cached.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      if (!BinningComputed) {
         ABORT_ERROR("CoordinateBinningOp::compute called before initialize");
      }

      // Update cache validity markers
      // (binning is constant, so no recomputation needed)
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

   /// Computes bin index for each entity based on coordinate value
   void computeBinning(const ArrayT &CoordData) {

      auto LocalOutput   = OutputData;
      auto LocalMinBin   = MinBin;
      auto LocalBinWidth = BinWidth;
      auto LocalNumBins  = NumBins;

      parallelFor(
          {ArraySize}, KOKKOS_LAMBDA(int i) {
             Real CoordVal = static_cast<Real>(CoordData(i));

             // Compute bin index: floor((coord - minBin) / binWidth)
             I4 BinIdx =
                 static_cast<I4>((CoordVal - LocalMinBin) / LocalBinWidth);

             // Clamp to valid range [0, NumBins-1]
             if (BinIdx < 0) {
                BinIdx = 0;
             } else if (BinIdx >= LocalNumBins) {
                BinIdx = LocalNumBins - 1;
             }

             LocalOutput(i) = BinIdx;
          });

   } // end computeBinning

 private:
   /// Output data array holding bin indices for each entity
   OutputArrayT OutputData;

   /// Number of spatial bins
   I4 NumBins;

   /// Minimum coordinate value for binning
   Real MinBin;

   /// Maximum coordinate value for binning
   Real MaxBin;

   /// Bin width (derived from MinBin, MaxBin, NumBins)
   Real BinWidth;

   /// Whether to auto-compute MinBin from mesh
   bool AutoComputeMin;

   /// Whether to auto-compute MaxBin from mesh
   bool AutoComputeMax;

   /// Total size of the coordinate array
   I4 ArraySize;

   /// Whether binning has been computed (once in initialize())
   bool BinningComputed;

}; // end class CoordinateBinningOp

} // end namespace OMEGA

#endif
