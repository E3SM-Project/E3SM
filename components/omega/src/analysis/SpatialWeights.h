#ifndef OMEGA_SPATIALWEIGHTS_H
#define OMEGA_SPATIALWEIGHTS_H
//===-- analysis/SpatialWeights.h - weights for statistics ------*- C++ -*-===//
//
/// \file
/// \brief Defines the SpatialWeights class used by the spatial statistics
///
/// A global mean or standard deviation of a field must weight each mesh
/// entity by how much of the ocean it represents, or on a variable-resolution
/// mesh the statistic is biased toward wherever the cells are small.
/// SpatialWeights builds the weights for a reduction of a field over all of
/// its owned mesh entities and active layers:
///
/// - A horizontal field (rank 1) is weighted by the area of each entity:
///   AreaCell for cells, DcEdge times DvEdge for edges and AreaTriangle for
///   vertices, masked to active entities.
/// - A layered field (last dimension NVertLayers) is weighted by mass: the
///   entity area times the current PseudoThickness, which is the layer's
///   renormalized mass per unit area, masked to active layers. The
///   pseudo-thickness is the mean of the two cells on an edge and the
///   kite-area-weighted mean of the cells on a vertex, as in the auxiliary
///   variables.
/// - An interface field (last dimension NVertLayersP1) is weighted by half
///   the mass of each adjacent active layer, so the weights sum to the
///   column mass and the result is the mass-weighted mean of the field taken
///   as piecewise linear between interfaces.
///
/// Only active entries take part: the weights are built from the active
/// masks by selection rather than multiplication, and the reductions skip
/// inactive entries, so the fill values in inactive layers are never read.
/// Area weights are constant and computed once. Mass weights depend on the
/// state and are recomputed by update() before each reduction. The weighted
/// sums accumulate each owned entity's contribution over its layers and any
/// leading extra dimensions (e.g. tracers) in a fixed order and pass the
/// per-entity terms to globalSum, so they are as reproducible as Omega's
/// other global sums.
///
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Error.h"
#include "HorzMesh.h"
#include "OmegaKokkos.h"
#include "Reductions.h"
#include "VertCoord.h"
#include "mpi.h"

#include <string>
#include <vector>

namespace OMEGA {

/// Weights for a reduction of a field over its owned mesh entities and
/// active layers, by area for a horizontal field and by mass otherwise
class SpatialWeights {

 public:
   /// The mesh index space a field lives on
   enum class IndexSpace { Cell, Edge, Vertex };

   /// Default constructor; init() must be called before use
   SpatialWeights();

   /// Sets up the weights for the named input Field. The index space is
   /// determined from the horizontal dimension name (NCells, NEdges or
   /// NVertices), which is the only dimension of a rank-1 field and the
   /// second-to-last dimension otherwise, and the vertical dimension name
   /// (NVertLayers or NVertLayersP1) selects layer or interface mass
   /// weights. Aborts on any other dimension name. Computes the constant
   /// area weights.
   void init(const std::string &InputName, ///< [in] input field name
             const HorzMesh *Mesh,         ///< [in] horizontal mesh
             const VertCoord *VCoord,      ///< [in] vertical coordinate
             MPI_Comm Comm                 ///< [in] communicator for sums
   );

   /// Recomputes the mass weights from the current pseudo-thickness and the
   /// global weight sum. For a horizontal field the area weights are
   /// constant, so only the sum is computed, and only once.
   void update();

   /// Returns the global sum of the weights times the field over the owned
   /// entities, active layers and all leading extra dimensions of the input
   /// array. Call update() first.
   template <typename ArrayT>
   Real weightedSum(const ArrayT &Data ///< [in] input array
   ) {
      return moment(Data, static_cast<Real>(0), 1);
   }

   /// Returns the global sum of the weights times the squared deviation of
   /// the field from a given mean, over the same entries as weightedSum
   template <typename ArrayT>
   Real weightedSquaredDeviation(const ArrayT &Data, ///< [in] input array
                                 Real Mean           ///< [in] mean of field
   ) {
      return moment(Data, Mean, 2);
   }

   /// True for a horizontal (rank-1) input weighted by area alone
   bool isHorizontal() const { return Horizontal; }

   /// Global sum of the weights over owned entities and layers, valid after
   /// update(). A field with leading extra dimensions replicates the weights
   /// over them, so the caller multiplies by their product.
   Real weightSum() const { return WeightSum; }

   /// Area weights of a horizontal field (NEntitiesSize), masked
   Array1DReal AreaWeights;

   /// Mass weights of a layered or interface field (NEntitiesSize by
   /// NVertLayers or NVertLayers + 1), masked
   Array2DReal MassWeights;

   /// Computes the mass of each active layer of each entity from the
   /// current pseudo-thickness into Layer
   void computeLayerMass();

   /// Global sum over owned entities, active layers and leading extra
   /// dimensions of the weight times (Data - Center)^Power, for Power 1 or 2
   template <typename ArrayT>
   Real moment(const ArrayT &Data, Real Center, int Power) {

      checkExtents(Data);
      constexpr int Rank = ArrayT::rank;

      OMEGA_SCOPE(LocTerms, Terms);
      OMEGA_SCOPE(LocAreaWeights, AreaWeights);
      OMEGA_SCOPE(LocMassWeights, MassWeights);
      OMEGA_SCOPE(LocActive, Active);
      const I4 LocNVert = NVert;

      if constexpr (Rank == 1) {
         parallelFor(
             {NOwned}, KOKKOS_LAMBDA(int I) {
                Real Sum = 0;
                if (LocActive(I, 0) > 0) {
                   Real Dev = static_cast<Real>(Data(I)) - Center;
                   Sum = LocAreaWeights(I) * (Power == 2 ? Dev * Dev : Dev);
                }
                LocTerms(I) = Sum;
             });
      } else if constexpr (Rank == 2) {
         parallelFor(
             {NOwned}, KOKKOS_LAMBDA(int I) {
                Real Sum = 0;
                for (I4 K = 0; K < LocNVert; ++K) {
                   if (LocActive(I, K) > 0) {
                      Real Dev = static_cast<Real>(Data(I, K)) - Center;
                      Sum +=
                          LocMassWeights(I, K) * (Power == 2 ? Dev * Dev : Dev);
                   }
                }
                LocTerms(I) = Sum;
             });
      } else if constexpr (Rank == 3) {
         const I4 N0 = Data.extent(0);
         parallelFor(
             {NOwned}, KOKKOS_LAMBDA(int I) {
                Real Sum = 0;
                for (I4 L = 0; L < N0; ++L) {
                   for (I4 K = 0; K < LocNVert; ++K) {
                      if (LocActive(I, K) > 0) {
                         Real Dev = static_cast<Real>(Data(L, I, K)) - Center;
                         Sum += LocMassWeights(I, K) *
                                (Power == 2 ? Dev * Dev : Dev);
                      }
                   }
                }
                LocTerms(I) = Sum;
             });
      } else {
         const I4 N0 = Data.extent(0);
         const I4 N1 = Data.extent(1);
         parallelFor(
             {NOwned}, KOKKOS_LAMBDA(int I) {
                Real Sum = 0;
                for (I4 L0 = 0; L0 < N0; ++L0) {
                   for (I4 L1 = 0; L1 < N1; ++L1) {
                      for (I4 K = 0; K < LocNVert; ++K) {
                         if (LocActive(I, K) > 0) {
                            Real Dev =
                                static_cast<Real>(Data(L0, L1, I, K)) - Center;
                            Sum += LocMassWeights(I, K) *
                                   (Power == 2 ? Dev * Dev : Dev);
                         }
                      }
                   }
                }
                LocTerms(I) = Sum;
             });
      }

      std::vector<I4> Range = {0, NOwned - 1};
      return globalSum(Terms, Comm, &Range);
   }

 private:
   const HorzMesh *Mesh;    ///< horizontal mesh
   const VertCoord *VCoord; ///< vertical coordinate
   MPI_Comm Comm;           ///< communicator for global sums

   std::string InputName; ///< input field name, for messages
   IndexSpace Space;      ///< index space of the input field
   bool Horizontal;       ///< true for a rank-1 input
   bool Interface;        ///< true if the vertical dimension is interfaces
   I4 NOwned;             ///< number of owned entities
   I4 NVert;              ///< vertical extent of the weights

   Array1DReal Area;   ///< area of each entity (NEntitiesSize)
   Array2DReal Mask;   ///< active-layer mask of each entity
   Array2DReal Active; ///< mask of active weights (NEntitiesSize by NVert)
   Array2DReal Layer;  ///< work array: mass of each active layer
   Array1DReal Terms;  ///< work array: per-entity terms of a weighted sum

   Real WeightSum;  ///< global sum of the weights
   bool SumIsValid; ///< true once the sum of constant weights is known

   /// Checks that an input array has the extents the weights were built for
   template <typename ArrayT> void checkExtents(const ArrayT &Data) const {
      constexpr int Rank = ArrayT::rank;
      if (Horizontal) {
         OMEGA_REQUIRE(Rank == 1,
                       "SpatialWeights: field {} is horizontal but the "
                       "array has rank {}",
                       InputName, Rank);
      } else {
         OMEGA_REQUIRE(Rank >= 2,
                       "SpatialWeights: field {} is layered but the array "
                       "has rank {}",
                       InputName, Rank);
         OMEGA_REQUIRE(static_cast<I4>(Data.extent(Rank - 1)) == NVert,
                       "SpatialWeights: field {} has vertical extent {} but "
                       "the weights have {}",
                       InputName, Data.extent(Rank - 1), NVert);
      }
      OMEGA_REQUIRE(
          static_cast<I4>(Data.extent(Rank == 1 ? 0 : Rank - 2)) >= NOwned,
          "SpatialWeights: field {} has fewer entities than owned", InputName);
   }

}; // end class SpatialWeights

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_SPATIALWEIGHTS_H
