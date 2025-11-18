#ifndef OMEGA_HORZOPERATORS_H
#define OMEGA_HORZOPERATORS_H

#include "DataTypes.h"
#include "HorzMesh.h"
#include "HorzUtil.h"

namespace OMEGA {

class DivergenceOnCell {
 public:
   DivergenceOnCell(HorzMesh const *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &DivCell, int ICell,
                                   int KChunk,
                                   const Array2DReal &VecEdge) const {
      const int KStart       = KChunk * VecLength;
      const Real InvAreaCell = 1._Real / AreaCell(ICell);

      Real DivCellTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnCell(ICell); ++J) {
         const int JEdge = EdgesOnCell(ICell, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            DivCellTmp[KVec] -= DvEdge(JEdge) * EdgeSignOnCell(ICell, J) *
                                VecEdge(JEdge, K) * InvAreaCell;
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K       = KStart + KVec;
         DivCell(ICell, K) = DivCellTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array1DReal DvEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeSignOnCell;
};

class GradientOnEdge {
 public:
   GradientOnEdge(HorzMesh const *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &GradEdge, int IEdge,
                                   int KChunk,
                                   const Array2DReal &ScalarCell) const {
      const int KStart     = KChunk * VecLength;
      const Real InvDcEdge = 1._Real / DcEdge(IEdge);
      const auto JCell0    = CellsOnEdge(IEdge, 0);
      const auto JCell1    = CellsOnEdge(IEdge, 1);

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K = KStart + KVec;
         GradEdge(IEdge, K) =
             InvDcEdge * (ScalarCell(JCell1, K) - ScalarCell(JCell0, K));
      }
   }

 private:
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
};

class CurlOnVertex {
 public:
   CurlOnVertex(HorzMesh const *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &CurlVertex, int IVertex,
                                   int KChunk,
                                   const Array2DReal &VecEdge) const {
      const int KStart           = KChunk * VecLength;
      const Real InvAreaTriangle = 1._Real / AreaTriangle(IVertex);

      Real CurlVertexTmp[VecLength] = {0};

      for (int J = 0; J < VertexDegree; ++J) {
         const int JEdge = EdgesOnVertex(IVertex, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            CurlVertexTmp[KVec] += DcEdge(JEdge) *
                                   EdgeSignOnVertex(IVertex, J) *
                                   VecEdge(JEdge, K) * InvAreaTriangle;
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K            = KStart + KVec;
         CurlVertex(IVertex, K) = CurlVertexTmp[KVec];
      }
   }

 private:
   I4 VertexDegree;
   Array2DI4 EdgesOnVertex;
   Array1DReal DcEdge;
   Array1DReal AreaTriangle;
   Array2DReal EdgeSignOnVertex;
};

class TangentialReconOnEdge {
 public:
   TangentialReconOnEdge(HorzMesh const *Mesh);

   KOKKOS_FUNCTION void operator()(const Array2DReal &ReconEdge, int IEdge,
                                   int KChunk,
                                   const Array2DReal &VecEdge) const {
      const int KStart = KChunk * VecLength;

      Real ReconEdgeTmp[VecLength] = {0};

      for (int J = 0; J < NEdgesOnEdge(IEdge); ++J) {
         const int JEdge = EdgesOnEdge(IEdge, J);
         for (int KVec = 0; KVec < VecLength; ++KVec) {
            const int K = KStart + KVec;
            ReconEdgeTmp[KVec] += WeightsOnEdge(IEdge, J) * VecEdge(JEdge, K);
         }
      }

      for (int KVec = 0; KVec < VecLength; ++KVec) {
         const int K         = KStart + KVec;
         ReconEdge(IEdge, K) = ReconEdgeTmp[KVec];
      }
   }

 private:
   Array1DI4 NEdgesOnEdge;
   Array2DI4 EdgesOnEdge;
   Array2DReal WeightsOnEdge;
};

enum class InterpCellToEdgeOption { Anisotropic, Isotropic };

class InterpCellToEdge {
 public:
   InterpCellToEdge(const HorzMesh *Mesh);

   KOKKOS_FUNCTION Real operator()(int IEdge, const Array1DReal &ArrayCell,
                                   InterpCellToEdgeOption Option) const {
      switch (Option) {
      case InterpCellToEdgeOption::Anisotropic:
         return interpolateAnisotropic(IEdge, ArrayCell);
      case InterpCellToEdgeOption::Isotropic:
         return interpolateIsotropic(IEdge, ArrayCell);
      default:
         return interpolateIsotropic(IEdge, ArrayCell);
      }
   };

 private:
   KOKKOS_FUNCTION Real
   interpolateAnisotropic(int IEdge, const Array1DReal &ArrayCell) const {
      const int JCell0 = CellsOnEdge(IEdge, 0);
      const int JCell1 = CellsOnEdge(IEdge, 1);

      return 0.5_Real * (ArrayCell(JCell0) + ArrayCell(JCell1));
   };

   KOKKOS_FUNCTION Real
   interpolateIsotropic(int IEdge, const Array1DReal &ArrayCell) const {

      Real Accum     = 0;
      Real AreaAccum = 0;
      for (int J = 0; J < 2; ++J) {
         const int JVertex = VerticesOnEdge(IEdge, J);
         for (int L = 0; L < VertexDegree; ++L) {
            const Real KiteArea = KiteAreasOnVertex(JVertex, L);
            const int LCell     = CellsOnVertex(JVertex, L);

            Accum += ArrayCell(LCell) * KiteArea;
            AreaAccum += KiteArea;
         }
      }

      const Real InvAreaAccum = 1._Real / AreaAccum;
      return Accum * InvAreaAccum;
   };

   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;
   Array2DI4 CellsOnVertex;
   Array2DReal KiteAreasOnVertex;
   I4 VertexDegree;
};

class SecondDerivativeOnCell {
   /*
      \brief deriv two computation
      \author Doug Jacobsen, Bill Skamarock
      \date   03/09/12
      \details
       This routine precomputes the second derivative values for tracer
       advection. It computes cell coefficients for the polynomial fit
       as described in:
       Skamarock, W. C., & Gassmann, A. (2011).
         Conservative Transport Schemes for Spherical Geodesic Meshs:
         High-Order Flux Operators for ODE-Based Time Integration.
         Monthly Weather Review, 139(9), 2962-2975.
         doi:10.1175/MWR-D-10-05056.1
       This is performed during model initialization.
   */
 public:
   SecondDerivativeOnCell(HorzMesh const *Mesh);
   KOKKOS_FUNCTION ~SecondDerivativeOnCell() {}
   KOKKOS_FUNCTION void operator()(const Array3DReal &DerivTwo,
                                   int ICell) const {
      const int NEdges = NEdgesOnCell(ICell);
      if (MaxMaxEdges < NEdges)
         printf("Error: Number of edges on cell:%d exceeds maximum "
                "expected:%d for cell:%d",
                NEdges, MaxMaxEdges, ICell);

      // check to see if we are reaching outside the halo
      auto CellList = Kokkos::subview(CellListCell, ICell, Kokkos::ALL);
      bool doCell   = true;
      CellList[0]   = ICell;
      for (int J = 1; J <= NEdges; ++J)
         CellList[J] = CellsOnCell(ICell, J - 1);

      for (int I = 0; I <= NEdges; ++I)
         if (NCellsAll <= CellList[I])
            doCell = false;
      if (!doCell)
         return;

      auto XP      = Kokkos::subview(XPCell, ICell, Kokkos::ALL);
      auto YP      = Kokkos::subview(YPCell, ICell, Kokkos::ALL);
      auto Angle2D = Kokkos::subview(Angle2DCell, ICell, Kokkos::ALL);
      auto B       = Kokkos::subview(BCell, ICell, Kokkos::ALL, Kokkos::ALL);
      if (OnSphere) {
         const Array1DI4 edgesOnCell =
             Kokkos::subview(EdgesOnCell, ICell, Kokkos::ALL);
         DetermineSphericalPatchGeometry(
             NEdges, edgesOnCell, VerticesOnEdge, XCell, YCell, ZCell, XVertex,
             YVertex, ZVertex, CellList, XP, YP, Angle2D, ThetaAbs[ICell]);
      } else { // On an x-y plane
         const Array1DI4 edgesOnCell =
             Kokkos::subview(EdgesOnCell, ICell, Kokkos::ALL);
         DeterminePlanerPatchGeometry(ICell, NEdges, edgesOnCell, CellsOnEdge,
                                      AngleEdge, DcEdge, XP, YP, Angle2D);
      }
      LeastSquaresFit(XP, YP, NEdges, B);

      // fill second derivative stencil for rk advection
      for (int I = 0; I < NEdges; ++I) {
         const I4 IEdge   = EdgesOnCell(ICell, I);
         const I4 Ind     = (ICell == CellsOnEdge(IEdge, 0)) ? 0 : 1;
         const Real Theta = Angle2D[I];
         const Real x     = Kokkos::cos(Theta);
         const Real y     = Kokkos::sin(Theta);
         const Real xx    = x * x;
         const Real xy    = x * y;
         const Real yy    = y * y;
         for (int J = 0; J <= NEdges; ++J)
            DerivTwo(J, Ind, IEdge) = 2._Real * xx * B(3, J) +
                                      2._Real * xy * B(4, J) +
                                      2._Real * yy * B(5, J);
      }
   }

 private:
   // MaxMaxEdges is used to dimention arrays that include ICell and the
   // neighbor cells, so it is technically one more than MaxEdges.
   static const I4 MaxMaxEdges = 10;
   static constexpr R8 Pii     = 3.141592653589793_Real;

   const bool OnSphere;
   const I4 NCellsAll;
   const I4 MaxEdges;
   Array1DI4 NEdgesOnCell;
   Array2DI4 EdgesOnCell;
   Array2DI4 CellsOnCell;
   Array2DI4 CellsOnEdge;
   Array2DI4 VerticesOnEdge;

   Array1DReal XCell;
   Array1DReal YCell;
   Array1DReal ZCell;
   Array1DReal DvEdge;
   Array1DReal DcEdge;
   Array1DReal AngleEdge;
   Array1DReal AreaCell;
   Array2DReal EdgeSignOnCell;
   Array1DReal XVertex;
   Array1DReal YVertex;
   Array1DReal ZVertex;

   Array1DReal ThetaAbs;
   Array2DReal XPCell;
   Array2DReal YPCell;
   Array2DReal Angle2DCell;
   Array3DReal BCell;
   Array2DI4 CellListCell;

 protected:
   KOKKOS_INLINE_FUNCTION static void LeastSquaresFit(const Array1DReal XP,
                                                      const Array1DReal YP,
                                                      const int NEdges,
                                                      Array2DReal B) {
      constexpr int NA        = 6;
      Real P[MaxMaxEdges][NA] = {};

      // (polynomial_order == 2) is the only order supported
      // The first row is for the origin (0,0) since the data
      // is relative to the ICell centroid.
      P[0][0] = 1._Real;
      for (int I = 1; I <= NEdges; ++I) {
         P[I][0] = 1._Real;
         P[I][1] = XP[I - 1];
         P[I][2] = YP[I - 1];
         P[I][3] = XP[I - 1] * XP[I - 1];
         P[I][4] = XP[I - 1] * YP[I - 1];
         P[I][5] = YP[I - 1] * YP[I - 1];
      }
      poly_fit_2<NA>(P, B, 1 + NEdges);
   }

   KOKKOS_INLINE_FUNCTION static void DeterminePlanerPatchGeometry(
       const int ICell, const int NEdges, const Array1DI4 EdgesOnCell,
       const Array2DI4 CellsOnEdge, const Array1DReal AngleEdge,
       const Array1DReal DcEdge, Array1DReal XP, Array1DReal YP,
       Array1DReal Angle2D) {
      for (int I = 0; I < NEdges; ++I) {
         const auto IEdge = EdgesOnCell[I];
         Angle2D[I]       = AngleEdge(IEdge);
         if (ICell != CellsOnEdge(IEdge, 0))
            Angle2D[I] -= Pii;
         XP[I] = DcEdge(IEdge) * Kokkos::cos(Angle2D[I]);
         YP[I] = DcEdge(IEdge) * Kokkos::sin(Angle2D[I]);
      }
   }
   KOKKOS_INLINE_FUNCTION static void DetermineSphericalPatchGeometry(
       const int NEdges, const Array1DI4 EdgesOnCell,
       const Array2DI4 VerticesOnEdge, const Array1DReal XCell,
       const Array1DReal YCell, const Array1DReal ZCell,
       const Array1DReal XVertex, const Array1DReal YVertex,
       const Array1DReal ZVertex, const Array1DI4 CellList, Array1DReal XP,
       Array1DReal YP, Array1DReal Angle2D, Real &ThetaAbs) {
      const Real length_scale = 1._Real;
      const Real sphereRadius = distance(XCell(0), YCell(0), ZCell(0));
      Real XC[MaxMaxEdges]    = {};
      Real YC[MaxMaxEdges]    = {};
      Real ZC[MaxMaxEdges]    = {};

      Real Thetat_prev = 0;
      Real Thetav_prev = 0;

      for (int I = 0; I <= NEdges; ++I) {
         const auto J = CellList[I];
         XC[I]        = XCell(J) / sphereRadius;
         YC[I]        = YCell(J) / sphereRadius;
         ZC[I]        = ZCell(J) / sphereRadius;
      }
      if (ZC[0] == 1.0_Real)
         ThetaAbs = Pii / 2._Real;
      else if (1 - ZC[0] < 1.0e-6)
         ThetaAbs = Pii / 2._Real - (1 - ZC[0]) * std::atan2(YC[0], XC[0]);
      else
         ThetaAbs =
             Pii / 2._Real - sphere_angle(XC[0], YC[0], ZC[0], XC[1], YC[1],
                                          ZC[1], 0._Real, 0._Real, 1._Real);

      for (int I = 0; I < NEdges; ++I) {
         int Ip1 = I + 1, Ip2 = I + 2;
         if (NEdges < Ip2)
            Ip2 = 1;
         // angles from cell center to neighbor centers (thetav)
         const Real Thetav = sphere_angle(XC[0], YC[0], ZC[0], XC[Ip1], YC[Ip1],
                                          ZC[Ip1], XC[Ip2], YC[Ip2], ZC[Ip2]);
         Real Dl_sphere    = sphereRadius * arc_length(XC[0], YC[0], ZC[0],
                                                       XC[Ip1], YC[Ip1], ZC[Ip1]);

         Dl_sphere /= length_scale;
         // Thetat = 0.  this defines the x direction,
         // cell center 0 -> this defines the x direction, longitude line
         const Real Thetat = I ? Thetat_prev + Thetav_prev : ThetaAbs;
         Thetat_prev       = Thetat;
         Thetav_prev       = Thetav;

         XP[I] = Kokkos::cos(Thetat) * Dl_sphere;
         YP[I] = Kokkos::sin(Thetat) * Dl_sphere;

         const I4 IEdge = EdgesOnCell[I];
         Real XV[2] = {}, YV[2] = {}, ZV[2] = {}, EC[3] = {};
         for (int J = 0; J < 2; ++J) {
            const I4 vert = VerticesOnEdge(IEdge, J);
            XV[J]         = XVertex(vert) / sphereRadius;
            YV[J]         = YVertex(vert) / sphereRadius;
            ZV[J]         = ZVertex(vert) / sphereRadius;
         }
         arc_bisect(XV[0], YV[0], ZV[0], XV[1], YV[1], ZV[1], EC[0], EC[1],
                    EC[2]);
         Real ThTmp = sphere_angle(XC[0], YC[0], ZC[0], XC[Ip1], YC[Ip1],
                                   ZC[Ip1], EC[0], EC[1], EC[2]);
         ThTmp += Thetat;
         Angle2D[I] = ThTmp;
      }
   }
};

class MasksAndCoefficients {

   KOKKOS_INLINE_FUNCTION static void swap(Array2DI4 &vec, const int m,
                                           const int n) {
      for (int i : {0, 1}) {
         const I4 j = vec(i, m);
         vec(i, m)  = vec(i, n);
         vec(i, n)  = j;
      }
   }
   // Sort the second dimension (values) of vec based on the first (keys):
   // Array1DI4 keys   = Kokkos::subview(vec, 0, Kokkos::ALL);
   // Array1DI4 values = Kokkos::subview(vec, 1, Kokkos::ALL);
   KOKKOS_INLINE_FUNCTION static int partition(Array2DI4 &vec, const int low,
                                               const int high) {
      // Selecting last element as the pivot
      const I4 pivot = vec(0, high);
      int i          = (low - 1);
      for (int j = low; j < high; ++j) {
         if (vec(0, j) <= pivot) {
            i++;
            swap(vec, i, j);
         }
      }
      // Put pivot to its position
      swap(vec, i + 1, high);
      // Return the point of partition
      return (i + 1);
   }

   KOKKOS_INLINE_FUNCTION static void sort_by_key(Array2DI4 &vec, const I4 low,
                                                  const I4 high) {
      // Base case: This part will be executed till the starting
      // index low is lesser than the ending index high
      if (low < high) {
         // pi is Partitioning Index, vec[pi] is now at
         // right place
         const I4 pi = partition(vec, low, high);
         // Separately sort elements before and after the
         // Partition Index pi
         sort_by_key(vec, low, pi - 1);
         sort_by_key(vec, pi + 1, high);
      }
   }
   KOKKOS_INLINE_FUNCTION static bool is_sorted(Array2DI4 &vec) {
      bool sorted = true;
      const I4 N  = vec.extent(1) - 1;
      for (I4 I = 0; I < N && sorted; ++I)
         if (vec(0, I + 1) < vec(0, I))
            sorted = false;
      return sorted;
   }
   KOKKOS_INLINE_FUNCTION static I4 search(Array1DI4 &vec, const I4 x) {
      I4 low = 0, high = vec.extent(0) - 1;
      while (low <= high) {
         const I4 mid = low + (high - low) / 2;
         if (vec(mid) == x)
            return mid;
         else if (vec(mid) < x)
            low = mid + 1;
         else
            high = mid - 1;
      }
      // If we reach here, then element was not present
      return -1;
   }
   KOKKOS_INLINE_FUNCTION static bool found_in_list(const Array1DI4 &List,
                                                    const I4 N, const I4 X) {
      bool found = false;
      for (I4 I = 0; I < N && !found; ++I)
         if (X == List(I))
            found = true;
      return found;
   }

 public:
   MasksAndCoefficients(HorzMesh const *Mesh, const Array3DReal DerivTwo,
                        Array1DI4 NAdvCellsForEdge, Array2DI4 AdvCellsForEdge,
                        Array1DI4 AdvMaskHighOrder, Array2DReal AdvCoefs,
                        Array2DReal AdvCoefs3rd);

   KOKKOS_FUNCTION void operator()(const int IEdge) const {

      // Array1DI4 PatchCellList = Kokkos::subview(PatchCellLists, IEdge,
      // Kokkos::ALL); for (I4 I=0; I< PatchCellList.extent(0))
      //   PatchCellList[I] = -1;
      NAdvCellsForEdge(IEdge) = 0;
      Array1DI4 CellIndex     = Kokkos::subview(CellIndx, IEdge, Kokkos::ALL);
      Array2DI4 CellIndexSorted =
          Kokkos::subview(CellIndxSorted, IEdge, Kokkos::ALL, Kokkos::ALL);
      for (int I = 0; I < CellIndexSorted.extent(0); ++I)
         for (int K = 0; K < CellIndexSorted.extent(1); ++K)
            CellIndexSorted(I, K) = MaxI4;
      const int Cell1 = CellsOnEdge(IEdge, 0);
      const int Cell2 = CellsOnEdge(IEdge, 1);
      // at boundaries, must stay at low order
      AdvMaskHighOrder(IEdge) = 1;
      for (int K = 0; K < NEdgesOnCell(Cell1); ++K)
         if (CellsOnCell(Cell1, K) == NCellsGlobal)
            AdvMaskHighOrder(IEdge) = 0;

      for (int K = 0; K < NEdgesOnCell(Cell2); ++K)
         if (CellsOnCell(Cell2, K) == NCellsGlobal)
            AdvMaskHighOrder(IEdge) = 0;
      // do only if this edge flux is needed to update owned cells
      if (Cell1 < NCellsAll && Cell2 < NCellsAll) {
         // Insert cellsOnEdge to list of advection cells
         // insert_into_list(PatchCellList,Cell1);
         // insert_into_list(PatchCellList,Cell2);
         CellIndex(0)          = Cell1;
         CellIndex(1)          = Cell2;
         CellIndexSorted(0, 0) = CellID(Cell1);
         CellIndexSorted(1, 0) = Cell1;
         CellIndexSorted(0, 1) = CellID(Cell2);
         CellIndexSorted(1, 1) = Cell2;
         int N                 = 2;
         // Build unique list of cells used for advection on edge
         // by expanding to the extended neighbor cells
         for (int I = 0; I < NEdgesOnCell(Cell1); ++I) {
            const I4 CellOnCell = CellsOnCell(Cell1, I);
            if (!found_in_list(CellIndex, N, CellOnCell)) {
               CellIndex(N)          = CellOnCell;
               CellIndexSorted(0, N) = CellID(CellOnCell);
               CellIndexSorted(1, N) = CellOnCell;
               ++N;
            }
         }
         for (int I = 0; I < NEdgesOnCell(Cell2); ++I) {
            const I4 CellOnCell = CellsOnCell(Cell2, I);
            if (!found_in_list(CellIndex, N, CellOnCell)) {
               CellIndex(N)          = CellOnCell;
               CellIndexSorted(0, N) = CellID(CellOnCell);
               CellIndexSorted(1, N) = CellOnCell;
               ++N;
            }
         }
         // sort the cell indices by cellID
         sort_by_key(CellIndexSorted, 0, CellIndexSorted.extent(1) - 1);

         Array1DI4 keys   = Kokkos::subview(CellIndexSorted, 0, Kokkos::ALL);
         Array1DI4 values = Kokkos::subview(CellIndexSorted, 1, Kokkos::ALL);
         //  store local cell indices for high-order calculations
         NAdvCellsForEdge(IEdge) = N;
         for (int ICell = 0; ICell < N; ++ICell)
            AdvCellsForEdge(IEdge, ICell) = values(ICell);
         // equation 7 in Skamarock, W. C., & Gassmann, A. (2011):
         // F(u,psi)_{i+1/2} = u_{i+1/2} *
         //  [   1/2 (psi_{i+1} + psi_i)                       term 1
         //    - 1/12(dx^2psi_{i+1} + dx^2psi_i)               term 2
         //    + sign(u) beta/12 (dx^2psi_{i+1} - dx^2psi_i)]  term 3
         //                                         (note minus sign)
         //
         // advCoefs accounts for terms 1 and 2 in SG11 equation 7.
         // Term 1 is the 2nd-order flux-function term. advCoefs
         // accounts for this with the "+ 0.5" lines below. In the
         // advection routines that use these coefficients, the
         // 2nd-order flux loop is then skipped. Term 2 is the
         // 4th-order flux-function term. advCoefs accounts for
         // term 3, the beta term. beta > 0 corresponds to the
         // third-order flux function. The - sign in the derivTwo
         // accumulation is for the i+1 part of term 3, while
         // the + sign is for the i part.

         for (int I = 0; I < NAdvCellsMax; ++I) {
            AdvCoefs(I, IEdge)    = 0._Real;
            AdvCoefs3rd(I, IEdge) = 0._Real;
         }
         // pull together third and fourth order contributions to the flux first
         // from cell1
         Array1DI4 keys_edge = Kokkos::subview(
             keys, Kokkos::make_pair(0, NAdvCellsForEdge(IEdge)));
         if (const I4 I = search(keys_edge, CellID(Cell1)); - 1 != I) {
            AdvCoefs(I, IEdge) += DerivTwo(0, 0, IEdge);
            AdvCoefs3rd(I, IEdge) += DerivTwo(0, 0, IEdge);
         }
         for (int ICell = 0; ICell < NEdgesOnCell(Cell1); ++ICell) {
            if (const I4 I =
                    search(keys_edge, CellID(CellsOnCell(Cell1, ICell)));
                - 1 != I) {
               AdvCoefs(I, IEdge) += DerivTwo(ICell + 1, 0, IEdge);
               AdvCoefs3rd(I, IEdge) += DerivTwo(ICell + 1, 0, IEdge);
            }
         }
         // pull together third and fourth order contributions to the flux first
         // from cell2
         if (const I4 I = search(keys_edge, CellID(Cell2)); - 1 != I) {
            AdvCoefs(I, IEdge) += DerivTwo(0, 1, IEdge);
            AdvCoefs3rd(I, IEdge) -= DerivTwo(0, 1, IEdge);
         }
         for (int ICell = 0; ICell < NEdgesOnCell(Cell2); ++ICell) {
            if (const I4 I =
                    search(keys_edge, CellID(CellsOnCell(Cell2, ICell)));
                - 1 != I) {
               AdvCoefs(I, IEdge) += DerivTwo(ICell + 1, 1, IEdge);
               AdvCoefs3rd(I, IEdge) -= DerivTwo(ICell + 1, 1, IEdge);
            }
         }
         for (int ICell = 0; ICell < NAdvCellsForEdge(IEdge); ++ICell) {
            AdvCoefs(ICell, IEdge) *= -DcEdge(IEdge) * DcEdge(IEdge) / 12._Real;
            AdvCoefs3rd(ICell, IEdge) *=
                -DcEdge(IEdge) * DcEdge(IEdge) / 12._Real;
         }
         // 2nd order centered contribution place this in the main flux weights
         if (const I4 I = search(keys_edge, CellID(Cell1)); - 1 != I) {
            AdvCoefs(I, IEdge) += 0.5_Real;
         }
         if (const I4 I = search(keys_edge, CellID(Cell2)); - 1 != I) {
            AdvCoefs(I, IEdge) += 0.5_Real;
         }
         // multiply by edge length - thus the flux is just dt*ru times the
         // results of the vector-vector multiply
         for (int ICell = 0; ICell < NAdvCellsForEdge(IEdge); ++ICell) {
            AdvCoefs(ICell, IEdge) *= DvEdge(IEdge);
            AdvCoefs3rd(ICell, IEdge) *= DvEdge(IEdge);
         }
      }
   }

 private:
   const I4 MaxI4 = std::numeric_limits<I4>::max();
   const I4 NCellsGlobal;
   const I4 NCellsAll;
   const I4 NAdvCellsMax;
   Array2DI4 PatchCellLists;
   Array1DI4 NAdvCellsForEdge;
   Array2DI4 AdvCellsForEdge;
   Array1DI4 NEdgesOnEdge;
   Array1DI4 NEdgesOnCell;
   Array2DI4 CellIndx;
   Array3DI4 CellIndxSorted;
   Array1DI4 CellID;
   Array1DI4 AdvMaskHighOrder;
   Array2DI4 EdgesOnEdge;
   Array2DI4 CellsOnCell;
   Array2DI4 CellsOnEdge;
   Array1DReal DcEdge;
   Array1DReal DvEdge;
   Array2DReal AdvCoefs;
   Array2DReal AdvCoefs3rd;
   Array3DReal DerivTwo;
};
} // namespace OMEGA
#endif
