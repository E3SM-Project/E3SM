//===-- analysis/SpatialWeights.cpp - weights for statistics ----*- C++ -*-===//
//
// Implementation of SpatialWeights. The area of each entity and its active
// mask are taken from the mesh and vertical coordinate at init(). The mass of
// each layer is the masked area times the current pseudo-thickness, which is
// read through the state's PseudoThickness Field at each update() so that it
// follows the time level the state is on.
//
//===----------------------------------------------------------------------===//

#include "SpatialWeights.h"
#include "Error.h"
#include "Field.h"
#include "Logging.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Reductions.h"

namespace OMEGA {

//------------------------------------------------------------------------------
// Default constructor
SpatialWeights::SpatialWeights()
    : Mesh(nullptr), VCoord(nullptr), Comm(MPI_COMM_NULL),
      Space(IndexSpace::Cell), Horizontal(true), Interface(false), NOwned(0),
      NVert(0), WeightSum(0), SumIsValid(false) {}

//------------------------------------------------------------------------------
// Sets up the weights for the named input Field
void SpatialWeights::init(const std::string &InName, const HorzMesh *InMesh,
                          const VertCoord *InVCoord, MPI_Comm InComm) {

   OMEGA_REQUIRE(InMesh, "SpatialWeights: null mesh for field {}", InName);
   OMEGA_REQUIRE(InVCoord,
                 "SpatialWeights: null vertical coordinate for "
                 "field {}",
                 InName);

   InputName = InName;
   Mesh      = InMesh;
   VCoord    = InVCoord;
   Comm      = InComm;

   auto InputField = Field::get(InputName);
   OMEGA_REQUIRE(InputField, "SpatialWeights: field {} not found", InputName);

   // The horizontal dimension is the only one of a rank-1 field and the
   // second-to-last one otherwise; the vertical dimension is the last
   std::vector<std::string> DimNames;
   InputField->getDimNames(DimNames);
   I4 NDims = DimNames.size();
   OMEGA_REQUIRE(NDims >= 1, "SpatialWeights: field {} has no dimensions",
                 InputName);

   Horizontal               = (NDims == 1);
   std::string HorizDimName = DimNames[std::max(0, NDims - 2)];
   I4 NEntities             = 0;
   Array1DReal EntityArea;
   Array2DReal EntityMask;

   if (HorizDimName == "NCells") {
      Space      = IndexSpace::Cell;
      NOwned     = Mesh->NCellsOwned;
      NEntities  = Mesh->NCellsSize;
      EntityArea = Mesh->AreaCell;
      EntityMask = VCoord->CellMask;
   } else if (HorizDimName == "NEdges") {
      Space      = IndexSpace::Edge;
      NOwned     = Mesh->NEdgesOwned;
      NEntities  = Mesh->NEdgesSize;
      EntityMask = VCoord->EdgeMask;
      // The area associated with an edge is the product of its two lengths
      EntityArea = Array1DReal("SpatialWeightsAreaEdge", NEntities);
      OMEGA_SCOPE(LocArea, EntityArea);
      OMEGA_SCOPE(LocDcEdge, Mesh->DcEdge);
      OMEGA_SCOPE(LocDvEdge, Mesh->DvEdge);
      parallelFor(
          {NEntities}, KOKKOS_LAMBDA(int IEdge) {
             LocArea(IEdge) = LocDcEdge(IEdge) * LocDvEdge(IEdge);
          });
   } else if (HorizDimName == "NVertices") {
      Space      = IndexSpace::Vertex;
      NOwned     = Mesh->NVerticesOwned;
      NEntities  = Mesh->NVerticesSize;
      EntityArea = Mesh->AreaTriangle;
      EntityMask = VCoord->VertexMask;
   } else {
      ABORT_ERROR("SpatialWeights: field {} has horizontal dimension {}, "
                  "expected NCells, NEdges or NVertices",
                  InputName, HorizDimName);
   }
   Area = EntityArea;
   Mask = EntityMask;

   I4 NVertLayers = VCoord->NVertLayers;

   Terms = Array1DReal("SpatialWeightsTerms_" + InputName, NEntities);

   if (Horizontal) {
      // Area weights for entities that are active at the surface
      NVert       = 1;
      AreaWeights = Array1DReal("SpatialWeightsArea_" + InputName, NEntities);
      Active =
          Array2DReal("SpatialWeightsActive_" + InputName, NEntities, NVert);
      OMEGA_SCOPE(LocWeights, AreaWeights);
      OMEGA_SCOPE(LocActive, Active);
      OMEGA_SCOPE(LocArea, Area);
      OMEGA_SCOPE(LocMask, Mask);
      parallelFor(
          {NEntities}, KOKKOS_LAMBDA(int I) {
             const bool IsActive = LocMask(I, 0) > 0;
             LocActive(I, 0)     = IsActive ? 1._Real : 0._Real;
             LocWeights(I)       = IsActive ? LocArea(I) : 0._Real;
          });
   } else {
      std::string VertDimName = DimNames[NDims - 1];
      if (VertDimName == "NVertLayers") {
         Interface = false;
         NVert     = NVertLayers;
      } else if (VertDimName == "NVertLayersP1") {
         Interface = true;
         NVert     = NVertLayers + 1;
      } else {
         ABORT_ERROR("SpatialWeights: field {} has vertical dimension {}, "
                     "expected NVertLayers or NVertLayersP1",
                     InputName, VertDimName);
      }
      MassWeights =
          Array2DReal("SpatialWeightsMass_" + InputName, NEntities, NVert);
      Active =
          Array2DReal("SpatialWeightsActive_" + InputName, NEntities, NVert);
      Layer = Array2DReal("SpatialWeightsLayer_" + InputName, NEntities,
                          NVertLayers);

      // A layer is active where the entity mask says so; an interface is
      // active if either adjacent layer is
      OMEGA_SCOPE(LocActive, Active);
      OMEGA_SCOPE(LocMask, Mask);
      const bool LocInterface = Interface;
      parallelFor(
          {NEntities, NVert}, KOKKOS_LAMBDA(int I, int K) {
             bool IsActive;
             if (LocInterface) {
                const bool Above = (K > 0) && (LocMask(I, K - 1) > 0);
                const bool Below = (K < NVertLayers) && (LocMask(I, K) > 0);
                IsActive         = Above || Below;
             } else {
                IsActive = LocMask(I, K) > 0;
             }
             LocActive(I, K) = IsActive ? 1._Real : 0._Real;
          });
   }

   SumIsValid = false;

} // end init

//------------------------------------------------------------------------------
// Computes the mass of each active layer of each entity: the entity area
// times the pseudo-thickness interpolated to the entity from the active cells
// around it. Inactive layers get zero without their values being read.
void SpatialWeights::computeLayerMass() {

   auto State = OceanState::getDefault();
   OMEGA_REQUIRE(State,
                 "SpatialWeights: no default ocean state to take the "
                 "pseudo-thickness from for field {}",
                 InputName);
   auto ThickField = Field::get(State->PseudoThicknessFldName);
   OMEGA_REQUIRE(ThickField,
                 "SpatialWeights: pseudo-thickness field {} not found",
                 State->PseudoThicknessFldName);
   Array2DReal ThickCell = ThickField->getDataArray<Array2DReal>();

   I4 NEntities   = Layer.extent(0);
   I4 NVertLayers = Layer.extent(1);

   OMEGA_SCOPE(LocLayer, Layer);
   OMEGA_SCOPE(LocArea, Area);
   OMEGA_SCOPE(LocMask, Mask);
   OMEGA_SCOPE(LocCellMask, VCoord->CellMask);

   if (Space == IndexSpace::Cell) {
      parallelFor(
          {NEntities, NVertLayers}, KOKKOS_LAMBDA(int ICell, int K) {
             LocLayer(ICell, K) = (LocMask(ICell, K) > 0)
                                      ? LocArea(ICell) * ThickCell(ICell, K)
                                      : 0._Real;
          });
   } else if (Space == IndexSpace::Edge) {
      OMEGA_SCOPE(LocCellsOnEdge, Mesh->CellsOnEdge);
      parallelFor(
          {NEntities, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int K) {
             Real Mass = 0._Real;
             if (LocMask(IEdge, K) > 0) {
                // The edge is active only where both cells are, so this is
                // the mean of two active thicknesses
                Real Thick = 0._Real;
                for (I4 J = 0; J < 2; ++J) {
                   const I4 JCell = LocCellsOnEdge(IEdge, J);
                   if (LocCellMask(JCell, K) > 0)
                      Thick += 0.5_Real * ThickCell(JCell, K);
                }
                Mass = LocArea(IEdge) * Thick;
             }
             LocLayer(IEdge, K) = Mass;
          });
   } else {
      OMEGA_SCOPE(LocCellsOnVertex, Mesh->CellsOnVertex);
      OMEGA_SCOPE(LocKiteAreas, Mesh->KiteAreasOnVertex);
      const I4 VertexDegree = Mesh->VertexDegree;
      parallelFor(
          {NEntities, NVertLayers}, KOKKOS_LAMBDA(int IVertex, int K) {
             Real Mass = 0._Real;
             if (LocMask(IVertex, K) > 0) {
                // The kite areas sum to the triangle area, so this is the
                // mass in the dual cell from its active kites
                for (I4 J = 0; J < VertexDegree; ++J) {
                   const I4 JCell = LocCellsOnVertex(IVertex, J);
                   if (LocCellMask(JCell, K) > 0)
                      Mass += LocKiteAreas(IVertex, J) * ThickCell(JCell, K);
                }
             }
             LocLayer(IVertex, K) = Mass;
          });
   }

} // end computeLayerMass

//------------------------------------------------------------------------------
// Recomputes the mass weights and the global weight sum
void SpatialWeights::update() {

   std::vector<I4> Range;

   if (Horizontal) {
      // Area weights are constant, so their sum is computed once
      if (!SumIsValid) {
         Range      = {0, NOwned - 1};
         WeightSum  = globalSum(AreaWeights, Comm, &Range);
         SumIsValid = true;
      }
      return;
   }

   computeLayerMass();

   I4 NEntities   = Layer.extent(0);
   I4 NVertLayers = Layer.extent(1);
   OMEGA_SCOPE(LocLayer, Layer);
   OMEGA_SCOPE(LocWeights, MassWeights);

   if (Interface) {
      // Each interface carries half the mass of each adjacent layer; the
      // surface and bottom interfaces have only one
      parallelFor(
          {NEntities, NVert}, KOKKOS_LAMBDA(int I, int K) {
             Real Above       = (K > 0) ? LocLayer(I, K - 1) : 0._Real;
             Real Below       = (K < NVertLayers) ? LocLayer(I, K) : 0._Real;
             LocWeights(I, K) = 0.5_Real * (Above + Below);
          });
   } else {
      deepCopy(MassWeights, Layer);
   }

   Range      = {0, NOwned - 1, 0, NVert - 1};
   WeightSum  = globalSum(MassWeights, Comm, &Range);
   SumIsValid = true;

} // end update

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
