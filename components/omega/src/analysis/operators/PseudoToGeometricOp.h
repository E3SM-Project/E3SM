#ifndef OMEGA_PSEUDOTOGEOMETRICOP_H
#define OMEGA_PSEUDOTOGEOMETRICOP_H

//===-- analysis/operators/PseudoToGeometricOp.h ----------------*- C++ -*-===//
//
/// \file
/// \brief Defines the PseudoToGeometricOp operator for coordinate conversion
///
/// PseudoToGeometricOp converts fields from pseudo-height coordinates to
/// geometric coordinates in Omega's non-Boussinesq formulation. This operator
/// is essential for diagnostics and analysis that require physical (geometric)
/// quantities rather than the prognostic pseudo-coordinate variables.
///
/// In Omega, the vertical coordinate is based on pseudo-height (proportional
/// to pressure), and prognostic variables include pseudo-thickness and
/// pseudo-velocity. The relationship between pseudo and geometric quantities
/// is:
///
///   geometric_quantity = (RhoSw * SpecVol) * pseudo_quantity
///
/// where RhoSw is the reference density and SpecVol is the in-situ specific
/// volume (SpecVol = 1/Rho). This is equivalent to (RhoSw/Rho) *
/// pseudo_quantity.
///
/// Vertical Grid Staggering:
/// The operator automatically handles vertical grid staggering. If the pseudo
/// field is at layer interfaces (NVertLayers+1) while specific volume is at
/// layer midpoints (NVertLayers), the operator interpolates specific volume to
/// interfaces using:
/// - Top interface (k=0): uses top layer specific volume
/// - Interior interfaces: averages adjacent layer specific volumes
/// - Bottom interface: uses bottom layer specific volume
///
/// The operator takes one input:
/// 1. A pseudo-coordinate field (e.g., VerticalPseudoVelocity, PseudoThickness)
///
/// The in-situ specific volume field (SpecVol) is fetched via the Field
/// registry
///
/// Example usage in operator chain:
/// \code
///   VerticalPseudoVelocity_PseudoToGeometric
/// \endcode
/// where PseudoToGeometric converts pseudo vertical velocity to geometric.
///
//===----------------------------------------------------------------------===//

#include "AnalysisOperator.h"
#include "GlobalConstants.h"

namespace OMEGA {

/// PseudoToGeometricOp converts pseudo-coordinate fields to geometric
/// coordinates using the specific volume scaling relationship from Omega's
/// non-Boussinesq formulation. Supports 1D (horizontal), 2D (horizontal ×
/// vertical), and 3D (structured horizontal × vertical) arrays. Works on any
/// horizontal index space (cells, edges, or vertices). The conversion formula
/// is: Output = (RhoSw * SpecVol) * Input.
template <typename ArrayT> class PseudoToGeometricOp : public AnalysisOperator {
 public:
   /// Scalar type extracted from the input array type
   using ScalarT = typename ArrayT::non_const_value_type;

   /// Output array type - same as input array type
   using OutputArrayT = ArrayT;

   /// Constructs a PseudoToGeometricOp operator. Retrieves input field
   /// dimensions and metadata, creates output Field for the
   /// geometric-coordinate quantity, allocates output data array, and registers
   /// the output Field. The output Field name is constructed as InputName +
   /// "_Geometric".
   PseudoToGeometricOp(const std::vector<std::string>
                           &UpstreamNames, ///< [in] input field names
                       Config Options      ///< [in] operator config
                       )
       : AnalysisOperator("PseudoToGeometric") {

      // Store input field names
      // UpstreamNames[0] = pseudo-coordinate field to convert
      InputNames = UpstreamNames;

      if (InputNames.size() < 1) {
         ABORT_ERROR("PseudoToGeometricOp: Requires 1 input (pseudo field)");
      }

      // Retrieve input pseudo-coordinate field to get dimensions
      auto PseudoField = Field::get(InputNames[0]);
      auto PseudoData  = PseudoField->template getDataArray<ArrayT>();

      // Get dimension info from pseudo field
      auto NDims = PseudoField->getNumDims();
      std::vector<std::string> DimNames;
      PseudoField->getDimNames(DimNames);

      // Construct output field name and set instance name
      std::string OutputFieldName = InputNames[0] + "_PseudoToGeometric";
      OutputNames                 = {OutputFieldName};
      InstanceName                = OutputFieldName;

      // Get input metadata. The conversion factor RhoSw * SpecVol is
      // dimensionless, so the geometric quantity has the units of the pseudo
      // quantity.
      std::string PseudoDescr;
      PseudoField->getMetadata("long_name", PseudoDescr);
      auto Meta = inheritMetadata(InputNames[0]);

      // Create output Field with updated metadata
      auto OutputField = createOutputField(
          OutputNames[0],
          "Geometric conversion of" + PseudoDescr, // Description
          Meta,                                    // CF metadata
          -std::numeric_limits<ScalarT>::max(),    // Min valid
          std::numeric_limits<ScalarT>::max(),     // Max valid
          NDims,                                   // Rank
          DimNames                                 // Dimension names
      );

      // Allocate output data array matching input layout
      OutputData = OutputArrayT(OutputNames[0] + "_out", PseudoData.layout());

      // Attach output data array to Field
      OutputField->template attachData<OutputArrayT>(OutputData);

   } // end constructor

   /// Initializes the operator after all Fields exist. Determines the index
   /// space (cells, edges, or vertices) from the pseudo field's horizontal
   /// dimension name and stores the appropriate MinLayer/MaxLayer arrays
   /// from VertCoord and connectivity arrays from Mesh for horizontal
   /// interpolation of SpecVol.
   void initialize(const MachEnv *Env,      ///< [in] machine environment
                   const HorzMesh *Mesh,    ///< [in] horizontal mesh
                   const VertCoord *VCoord, ///< [in] vertical coordinate
                   Config Options           ///< [in] operator config
                   ) override {

      // Call base class initialize to store Mesh, VCoord, etc.
      AnalysisOperator::initialize(Env, Mesh, VCoord, Options);

      // Get array rank for conditional logic
      constexpr I4 InputRank = ArrayT::rank;

      auto PseudoField = Field::get(InputNames[0]);
      auto PseudoData  = PseudoField->template getDataArray<ArrayT>();

      // Determine index space from horizontal dimension
      std::string IndexSpaceName;
      if constexpr (InputRank == 1) {
         // 1D: first (only) dimension
         std::vector<std::string> DimNames;
         PseudoField->getDimNames(DimNames);
         IndexSpaceName = DimNames[0];
         NVertSize      = 1;
      } else {
         // 2D/3D: horizontal is 2nd to last dimension
         std::vector<std::string> DimNames;
         PseudoField->getDimNames(DimNames);
         IndexSpaceName = DimNames[InputRank - 2];
         NVertSize      = PseudoData.extent(InputRank - 1);
      }

      // Set index space type and connectivity arrays
      if (IndexSpaceName == "NCells") {
         NHorizDim      = Mesh->NCellsOwned;
         IndexSpaceType = IndexSpace::Cell;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerCell;
            MaxLayer = VCoord->MaxLayerCell;
         }
      } else if (IndexSpaceName == "NEdges") {
         NHorizDim      = Mesh->NEdgesOwned;
         IndexSpaceType = IndexSpace::Edge;
         CellsOnEdge    = Mesh->CellsOnEdge;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerEdgeBot;
            MaxLayer = VCoord->MaxLayerEdgeTop;
         }
      } else if (IndexSpaceName == "NVertices") {
         NHorizDim      = Mesh->NVerticesOwned;
         IndexSpaceType = IndexSpace::Vertex;
         CellsOnVertex  = Mesh->CellsOnVertex;
         VertexDegree   = Mesh->VertexDegree;
         if constexpr (InputRank > 1) {
            MinLayer = VCoord->MinLayerVertexBot;
            MaxLayer = VCoord->MaxLayerVertexTop;
         }
      } else {
         ABORT_ERROR("PseudoToGeometricOp: Unknown index space {}",
                     IndexSpaceName);
      }

   } // end initialize

   /// Computes the conversion from pseudo to geometric coordinates by
   /// multiplying the pseudo field by (RhoSw * SpecVol) at each horizontal
   /// point and vertical layer. Uses hierarchical parallelism with outer loop
   /// over horizontal dimension and inner loop over vertical layers respecting
   /// MinLayer/MaxLayer bounds. For fields at layer interfaces (like vertical
   /// velocity), specific volume is interpolated from layer midpoints to
   /// interfaces. Updates output data, timestamp, and computed flag.
   void compute(const TimeInstant &TimeStamp ///< [in] current timestamp
                ) override {

      // Get array rank for conditional logic
      constexpr I4 InputRank = ArrayT::rank;

      // Retrieve input pseudo field
      auto PseudoField = Field::get(InputNames[0]);
      auto PseudoData  = PseudoField->template getDataArray<ArrayT>();

      // Retrieve SpecVol from Field registry
      auto SpecVolField = Field::get("SpecVol");
      auto SpecVolData  = SpecVolField->getDataArray<Array2DReal>();

      OMEGA_SCOPE(LocOutput, OutputData);

      // Create functor for horizontal interpolation
      HorizInterpFunctor getSpecVol{SpecVolData, IndexSpaceType, CellsOnEdge,
                                    CellsOnVertex, VertexDegree};

      if constexpr (InputRank == 1) {
         // 1D case: horizontal field only, no vertical structure
         parallelFor(
             {NHorizDim}, KOKKOS_LAMBDA(int IHoriz) {
                Real SpecVolHoriz = getSpecVol(IHoriz, 0);
                LocOutput(IHoriz) =
                    static_cast<ScalarT>((RhoSw * SpecVolHoriz) *
                                         static_cast<Real>(PseudoData(IHoriz)));
             });

      } else if constexpr (InputRank == 2) {
         // 2D case: hierarchical parallelism over horizontal × vertical

         // Check vertical staggering
         const I4 SpecVolNVertSize = SpecVolData.extent(1);
         const bool AtInterfaces   = (NVertSize == SpecVolNVertSize + 1);

         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         OMEGA_SCOPE(LocNVertSize, NVertSize);

         parallelForOuter(
             "PseudoToGeometric2D", LaunchConfig({NHorizDim}),
             KOKKOS_LAMBDA(int IHoriz, const TeamMember &Team) {
                const I4 KMin   = LocMinLayer(IHoriz);
                const I4 KMax   = LocMaxLayer(IHoriz);
                const I4 KRange = vertRange(KMin, KMax);

                parallelForInner(
                    Team, KRange, INNER_LAMBDA(int KIdx) {
                       const I4 K = KMin + KIdx;

                       // Horizontal interpolation first, then vertical if
                       // needed
                       Real SpecVolAtK;
                       if (AtInterfaces) {
                          if (K == 0) {
                             // Top interface: horizontally interpolate at top
                             // layer
                             SpecVolAtK = getSpecVol(IHoriz, 0);
                          } else if (K == LocNVertSize - 1) {
                             // Bottom interface: horizontally interpolate at
                             // bottom
                             SpecVolAtK =
                                 getSpecVol(IHoriz, SpecVolNVertSize - 1);
                          } else {
                             // Interior interface: horiz interp then vert
                             // average
                             Real SpecVolAbove = getSpecVol(IHoriz, K - 1);
                             Real SpecVolBelow = getSpecVol(IHoriz, K);
                             SpecVolAtK =
                                 0.5_Real * (SpecVolAbove + SpecVolBelow);
                          }
                       } else {
                          // Same vertical staggering: horizontal interpolation
                          // only
                          SpecVolAtK = getSpecVol(IHoriz, K);
                       }

                       LocOutput(IHoriz, K) = static_cast<ScalarT>(
                           (RhoSw * SpecVolAtK) *
                           static_cast<Real>(PseudoData(IHoriz, K)));
                    });
             });

      } else if constexpr (InputRank == 3) {
         // 3D case: hierarchical parallelism over tracers/dim0 × horizontal
         // Following the pattern from TimeStepper.cpp for tracer loops

         // Check vertical staggering
         const I4 SpecVolNVertSize = SpecVolData.extent(1);
         const bool AtInterfaces   = (NVertSize == SpecVolNVertSize + 1);

         OMEGA_SCOPE(LocMinLayer, MinLayer);
         OMEGA_SCOPE(LocMaxLayer, MaxLayer);
         OMEGA_SCOPE(LocNVertSize, NVertSize);

         const I4 Dim0 = PseudoData.extent(0);

         parallelForOuter(
             "PseudoToGeometric3D", LaunchConfig({Dim0, NHorizDim}),
             KOKKOS_LAMBDA(int I0, int IHoriz, const TeamMember &Team) {
                const I4 KMin   = LocMinLayer(IHoriz);
                const I4 KMax   = LocMaxLayer(IHoriz);
                const I4 KRange = vertRange(KMin, KMax);

                parallelForInner(
                    Team, KRange, INNER_LAMBDA(int KIdx) {
                       const I4 K = KMin + KIdx;

                       // Horizontal interpolation first, then vertical if
                       // needed
                       Real SpecVolAtK;
                       if (AtInterfaces) {
                          if (K == 0) {
                             // Top interface: horizontally interpolate at top
                             // layer
                             SpecVolAtK = getSpecVol(IHoriz, 0);
                          } else if (K == LocNVertSize - 1) {
                             // Bottom interface: horizontally interpolate at
                             // bottom
                             SpecVolAtK =
                                 getSpecVol(IHoriz, SpecVolNVertSize - 1);
                          } else {
                             // Interior interface: horiz interp then vert
                             // average
                             Real SpecVolAbove = getSpecVol(IHoriz, K - 1);
                             Real SpecVolBelow = getSpecVol(IHoriz, K);
                             SpecVolAtK =
                                 0.5_Real * (SpecVolAbove + SpecVolBelow);
                          }
                       } else {
                          // Same vertical staggering: horizontal interpolation
                          // only
                          SpecVolAtK = getSpecVol(IHoriz, K);
                       }

                       LocOutput(I0, IHoriz, K) = static_cast<ScalarT>(
                           (RhoSw * SpecVolAtK) *
                           static_cast<Real>(PseudoData(I0, IHoriz, K)));
                    });
             });
      }

      // Update cache validity markers
      LastComputed  = TimeStamp;
      FieldComputed = true;

   } // end compute

   /// Enum for horizontal index space type
   enum class IndexSpace { Cell, Edge, Vertex };

   /// This functor holds all necessary data to perform horizontal
   /// interpolation of SpecVol.
   struct HorizInterpFunctor {
      Array2DReal SpecVolData;
      IndexSpace IndexSpaceType;
      Array2DI4 CellsOnEdge;
      Array2DI4 CellsOnVertex;
      I4 VertexDegree;

      /// Horizontally interpolates SpecVol to the field's index space.
      /// For cell-based fields, returns SpecVol directly.
      /// For edge-based fields, averages from 2 adjacent cells.
      /// For vertex-based fields, averages from VertexDegree adjacent cells.
      KOKKOS_INLINE_FUNCTION
      Real operator()(const I4 IHoriz, ///< Horizontal index
                      const I4 K       ///< Vertical layer
      ) const {

         if (IndexSpaceType == IndexSpace::Cell) {
            // Cell-based: use directly
            return SpecVolData(IHoriz, K);

         } else if (IndexSpaceType == IndexSpace::Edge) {
            // Edge-based: average from 2 adjacent cells
            const I4 Cell1 = CellsOnEdge(IHoriz, 0);
            const I4 Cell2 = CellsOnEdge(IHoriz, 1);
            return 0.5_Real * (SpecVolData(Cell1, K) + SpecVolData(Cell2, K));

         } else { // IndexSpace::Vertex
            // Vertex-based: average from VertexDegree adjacent cells
            Real SpecVolSum = 0.0_Real;
            for (int j = 0; j < VertexDegree; ++j) {
               const I4 Cell = CellsOnVertex(IHoriz, j);
               SpecVolSum += SpecVolData(Cell, K);
            }
            return SpecVolSum / static_cast<Real>(VertexDegree);
         }
      }
   };

 private:
   /// Output data array holding the geometric-coordinate field
   OutputArrayT OutputData;

   /// Number of points in horizontal dimension (cells, edges, or vertices)
   I4 NHorizDim;

   /// Vertical size of the pseudo field array
   I4 NVertSize;

   /// Type of horizontal index space (cell, edge, or vertex)
   IndexSpace IndexSpaceType;

   /// Connectivity: cells adjacent to each edge (for edge-based fields)
   Array2DI4 CellsOnEdge;

   /// Connectivity: cells sharing each vertex (for vertex-based fields)
   Array2DI4 CellsOnVertex;

   /// Number of cells sharing each vertex (for vertex-based fields)
   I4 VertexDegree;

   /// Min active layer index for each horizontal point (only for 2D/3D arrays)
   Array1DI4 MinLayer;

   /// Max active layer index for each horizontal point (only for 2D/3D arrays)
   Array1DI4 MaxLayer;

}; // end class PseudoToGeometricOp

} // end namespace OMEGA

#endif
