#ifndef OMEGA_OCEAN_TEST_COMMON_H
#define OMEGA_OCEAN_TEST_COMMON_H

#include "DataTypes.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

namespace OMEGA {

// check if two real numbers are equal with a given relative tolerance
KOKKOS_INLINE_FUNCTION
bool isApprox(Real X, Real Y, Real RTol, Real ATol = 0) {
   if (Kokkos::isnan(X) || Kokkos::isnan(Y) || Kokkos::isinf(X) ||
       Kokkos::isinf(Y)) {
      return false; // Treat NaN or Inf as failure
   }

   return Kokkos::abs(X - Y) <=
          Kokkos::max(ATol, RTol * Kokkos::max(Kokkos::abs(X), Kokkos::abs(Y)));
}

// convert spherical components of a vector to Cartesian
KOKKOS_INLINE_FUNCTION void sphereToCartVec(Real (&CartVec)[3],
                                            const Real (&SphereVec)[2],
                                            Real Lon, Real Lat) {
   using std::cos;
   using std::sin;
   CartVec[0] = -sin(Lon) * SphereVec[0] - sin(Lat) * cos(Lon) * SphereVec[1];
   CartVec[1] = cos(Lon) * SphereVec[0] - sin(Lat) * sin(Lon) * SphereVec[1];
   CartVec[2] = cos(Lat) * SphereVec[1];
}

// returns Cartesian components of unit vector tangent to the spherical arc
// between Cartesian points X1 and X2 parametrized with t
KOKKOS_INLINE_FUNCTION void tangentVector(Real (&TanVec)[3],
                                          const Real (&X1)[3],
                                          const Real (&X2)[3], Real t = 0) {
   const Real Radius =
       Kokkos::sqrt(X1[0] * X1[0] + X1[1] * X1[1] + X1[2] * X1[2]);
   Real XC[3];
   Real DX[3];
   for (int Dim = 0; Dim < 3; ++Dim) {
      XC[Dim] = (1 - t) * X1[Dim] + t * X2[Dim];
      DX[Dim] = X2[Dim] - X1[Dim];
   }
   const Real XCDotDX = XC[0] * DX[0] + XC[1] * DX[1] + XC[2] * DX[2];
   const Real NormXC =
       Kokkos::sqrt(XC[0] * XC[0] + XC[1] * XC[1] + XC[2] * XC[2]);

   for (int Dim = 0; Dim < 3; ++Dim) {
      const Real NormXC3 = NormXC * NormXC * NormXC;
      TanVec[Dim] =
          Radius / NormXC * DX[Dim] - (Radius * XCDotDX) / NormXC3 * XC[Dim];
   }

   const Real NormTanVec = Kokkos::sqrt(
       TanVec[0] * TanVec[0] + TanVec[1] * TanVec[1] + TanVec[2] * TanVec[2]);
   for (int Dim = 0; Dim < 3; ++Dim) {
      TanVec[Dim] /= NormTanVec;
   }
}

enum class EdgeComponent { Normal, Tangential };
enum class Geometry { Planar, Spherical };
enum class ExchangeHalos { Yes, No };

// helper to get vertical iteration bounds that can be provided
// either as an integer or an array
template <class T> int getVertBound(const T &VertBound, int I) {
   if constexpr (std::is_integral_v<T>) {
      return VertBound;
   } else {
      static_assert(Kokkos::is_view_v<T>);
      return VertBound(I);
   }
}

// set scalar field on chosen elements (cells/vertices/edges) based on
// analytical formula and optionally exchange halos
template <class Functor, class Array, class VertMin, class VertMax>
int setScalar(const Functor &Fun, const Array &ScalarElement, Geometry Geom,
              const HorzMesh *Mesh, MeshElement Element, const VertMin &VMin,
              const VertMax &VMax,
              ExchangeHalos ExchangeHalosOpt = ExchangeHalos::Yes) {

   int Err = 0;

   int NElementsOwned;
   Array1DReal XElement, YElement;
   Array1DReal LonElement, LatElement;

   switch (Element) {
   case OnCell:
      NElementsOwned = Mesh->NCellsOwned;
      XElement       = createDeviceMirrorCopy(Mesh->XCellH);
      YElement       = createDeviceMirrorCopy(Mesh->YCellH);
      LonElement     = createDeviceMirrorCopy(Mesh->LonCellH);
      LatElement     = createDeviceMirrorCopy(Mesh->LatCellH);
      break;
   case OnVertex:
      NElementsOwned = Mesh->NVerticesOwned;
      XElement       = createDeviceMirrorCopy(Mesh->XVertexH);
      YElement       = createDeviceMirrorCopy(Mesh->YVertexH);
      LonElement     = createDeviceMirrorCopy(Mesh->LonVertexH);
      LatElement     = createDeviceMirrorCopy(Mesh->LatVertexH);
      break;
   case OnEdge:
      NElementsOwned = Mesh->NEdgesOwned;
      XElement       = createDeviceMirrorCopy(Mesh->XEdgeH);
      YElement       = createDeviceMirrorCopy(Mesh->YEdgeH);
      LonElement     = createDeviceMirrorCopy(Mesh->LonEdgeH);
      LatElement     = createDeviceMirrorCopy(Mesh->LatEdgeH);
      break;
   default:
      LOG_ERROR("setScalar: element needs to be one of (OnCell, OnVertex, "
                "OnEdge)");
      return 1;
   }

   if constexpr (Array::rank == 1) {
      parallelFor(
          {NElementsOwned}, KOKKOS_LAMBDA(int IElement) {
             if (Geom == Geometry::Planar) {
                const Real X            = XElement(IElement);
                const Real Y            = YElement(IElement);
                ScalarElement(IElement) = Fun(X, Y);
             } else {
                const Real Lon          = LonElement(IElement);
                const Real Lat          = LatElement(IElement);
                ScalarElement(IElement) = Fun(Lon, Lat);
             }
          });
   }

   if constexpr (Array::rank == 2) {
      const int NVertLayers = ScalarElement.extent_int(1);

      parallelForOuter(
          {NElementsOwned},
          KOKKOS_LAMBDA(int IElement, const TeamMember &Team) {
             const int KMin   = getVertBound(VMin, IElement);
             const int KMax   = getVertBound(VMax, IElement);
             const int KRange = KMax - KMin + 1;
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KOff) {
                    const int K = KMin + KOff;
                    if (Geom == Geometry::Planar) {
                       const Real X               = XElement(IElement);
                       const Real Y               = YElement(IElement);
                       ScalarElement(IElement, K) = Fun(X, Y);
                    } else {
                       const Real Lon             = LonElement(IElement);
                       const Real Lat             = LatElement(IElement);
                       ScalarElement(IElement, K) = Fun(Lon, Lat);
                    }
                 });
          });
   }

   if constexpr (Array::rank == 3) {
      const int NTracers = ScalarElement.extent_int(0);
      parallelForOuter(
          {NTracers, NElementsOwned},
          KOKKOS_LAMBDA(int L, int IElement, const TeamMember &Team) {
             const int KMin   = getVertBound(VMin, IElement);
             const int KMax   = getVertBound(VMax, IElement);
             const int KRange = KMax - KMin + 1;
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KOff) {
                    const int K = KMin + KOff;
                    if (Geom == Geometry::Planar) {
                       const Real X                  = XElement(IElement);
                       const Real Y                  = YElement(IElement);
                       ScalarElement(L, IElement, K) = Fun(X, Y);
                    } else {
                       const Real Lon                = LonElement(IElement);
                       const Real Lat                = LatElement(IElement);
                       ScalarElement(L, IElement, K) = Fun(Lon, Lat);
                    }
                 });
          });
   }

   if (ExchangeHalosOpt == ExchangeHalos::Yes) {
      auto MyHalo = Halo::getDefault();
      Err         = MyHalo->exchangeFullArrayHalo(ScalarElement, Element);
      if (Err != 0)
         LOG_ERROR("setScalarElement: error in halo exchange");
   }
   return Err;
}

// This overload calls setScalar with vertical bounds based on the array size
template <class Functor, class Array>
int setScalar(const Functor &Fun, const Array &ScalarElement, Geometry Geom,
              const HorzMesh *Mesh, MeshElement Element,
              ExchangeHalos ExchangeHalosOpt = ExchangeHalos::Yes) {
   const int VMin = 0;
   const int VMax = ScalarElement.extent_int(Array::rank - 1) - 1;
   return setScalar(Fun, ScalarElement, Geom, Mesh, Element, VMin, VMax,
                    ExchangeHalosOpt);
}

enum class CartProjection { Yes, No };

// set vector field on edges based on analytical formula and optionally
// exchange halos
template <class Functor, class Array, class VertMin, class VertMax>
int setVectorEdge(const Functor &Fun, const Array &VectorFieldEdge,
                  EdgeComponent EdgeComp, Geometry Geom, const HorzMesh *Mesh,
                  const VertMin &VMin, const VertMax &VMax,
                  ExchangeHalos ExchangeHalosOpt   = ExchangeHalos::Yes,
                  CartProjection CartProjectionOpt = CartProjection::Yes) {

   int Err = 0;

   auto XEdge = createDeviceMirrorCopy(Mesh->XEdgeH);
   auto YEdge = createDeviceMirrorCopy(Mesh->YEdgeH);
   auto ZEdge = createDeviceMirrorCopy(Mesh->ZEdgeH);

   auto XCell = createDeviceMirrorCopy(Mesh->XCellH);
   auto YCell = createDeviceMirrorCopy(Mesh->YCellH);
   auto ZCell = createDeviceMirrorCopy(Mesh->ZCellH);

   auto XVertex = createDeviceMirrorCopy(Mesh->XVertexH);
   auto YVertex = createDeviceMirrorCopy(Mesh->YVertexH);
   auto ZVertex = createDeviceMirrorCopy(Mesh->ZVertexH);

   auto LonEdge = createDeviceMirrorCopy(Mesh->LonEdgeH);
   auto LatEdge = createDeviceMirrorCopy(Mesh->LatEdgeH);

   auto &AngleEdge      = Mesh->AngleEdge;
   auto &CellsOnEdge    = Mesh->CellsOnEdge;
   auto &VerticesOnEdge = Mesh->VerticesOnEdge;

   auto ProjectVector = KOKKOS_LAMBDA(int IEdge) {
      Real VecFieldEdge;
      if (Geom == Geometry::Planar) {
         const Real XE = XEdge(IEdge);
         const Real YE = YEdge(IEdge);

         Real VecField[2];
         Fun(VecField, XE, YE);

         if (EdgeComp == EdgeComponent::Normal) {
            const Real EdgeNormalX = std::cos(AngleEdge(IEdge));
            const Real EdgeNormalY = std::sin(AngleEdge(IEdge));
            VecFieldEdge =
                EdgeNormalX * VecField[0] + EdgeNormalY * VecField[1];
         }

         if (EdgeComp == EdgeComponent::Tangential) {
            const Real EdgeTangentX = -std::sin(AngleEdge(IEdge));
            const Real EdgeTangentY = std::cos(AngleEdge(IEdge));
            VecFieldEdge =
                EdgeTangentX * VecField[0] + EdgeTangentY * VecField[1];
         }
      } else {
         const Real LonE = LonEdge(IEdge);
         const Real LatE = LatEdge(IEdge);

         Real VecField[2];
         Fun(VecField, LonE, LatE);

         if (CartProjectionOpt == CartProjection::Yes) {
            Real VecFieldCart[3];
            sphereToCartVec(VecFieldCart, VecField, LonE, LatE);

            const Real EdgeCoords[3] = {XEdge[IEdge], YEdge[IEdge],
                                        ZEdge[IEdge]};

            if (EdgeComp == EdgeComponent::Normal) {
               const int JCell1         = CellsOnEdge(IEdge, 1);
               const Real CellCoords[3] = {XCell(JCell1), YCell(JCell1),
                                           ZCell(JCell1)};

               Real EdgeNormal[3];
               tangentVector(EdgeNormal, EdgeCoords, CellCoords);
               VecFieldEdge = EdgeNormal[0] * VecFieldCart[0] +
                              EdgeNormal[1] * VecFieldCart[1] +
                              EdgeNormal[2] * VecFieldCart[2];
            }

            if (EdgeComp == EdgeComponent::Tangential) {
               const int JVertex1         = VerticesOnEdge(IEdge, 1);
               const Real VertexCoords[3] = {
                   XVertex(JVertex1), YVertex(JVertex1), ZVertex(JVertex1)};

               Real EdgeTangent[3];
               tangentVector(EdgeTangent, EdgeCoords, VertexCoords);
               VecFieldEdge = EdgeTangent[0] * VecFieldCart[0] +
                              EdgeTangent[1] * VecFieldCart[1] +
                              EdgeTangent[2] * VecFieldCart[2];
            }
         } else {
            if (EdgeComp == EdgeComponent::Normal) {
               const Real EdgeNormalX = std::cos(AngleEdge(IEdge));
               const Real EdgeNormalY = std::sin(AngleEdge(IEdge));
               VecFieldEdge =
                   EdgeNormalX * VecField[0] + EdgeNormalY * VecField[1];
            }

            if (EdgeComp == EdgeComponent::Tangential) {
               const Real EdgeTangentX = -std::sin(AngleEdge(IEdge));
               const Real EdgeTangentY = std::cos(AngleEdge(IEdge));
               VecFieldEdge =
                   EdgeTangentX * VecField[0] + EdgeTangentY * VecField[1];
            }
         }
      }
      return VecFieldEdge;
   };

   if constexpr (Array::rank == 1) {
      parallelFor(
          {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
             VectorFieldEdge(IEdge) = ProjectVector(IEdge);
          });
   }

   if constexpr (Array::rank == 2) {
      parallelForOuter(
          {Mesh->NEdgesOwned},
          KOKKOS_LAMBDA(int IEdge, const TeamMember &Team) {
             const int KMin   = getVertBound(VMin, IEdge);
             const int KMax   = getVertBound(VMax, IEdge);
             const int KRange = KMax - KMin + 1;
             parallelForInner(
                 Team, KRange, INNER_LAMBDA(int KOff) {
                    const int K               = KMin + KOff;
                    VectorFieldEdge(IEdge, K) = ProjectVector(IEdge);
                 });
          });
   }

   if (ExchangeHalosOpt == ExchangeHalos::Yes) {
      auto MyHalo = Halo::getDefault();
      Err         = MyHalo->exchangeFullArrayHalo(VectorFieldEdge, OnEdge);
      if (Err != 0)
         LOG_ERROR("setVectorEdge: error in halo exchange");
   }
   return Err;
}

// This overload calls setVectorEdge with vertical bounds based on the array
// size
template <class Functor, class Array>
int setVectorEdge(const Functor &Fun, const Array &VectorFieldEdge,
                  EdgeComponent EdgeComp, Geometry Geom, const HorzMesh *Mesh,
                  ExchangeHalos ExchangeHalosOpt   = ExchangeHalos::Yes,
                  CartProjection CartProjectionOpt = CartProjection::Yes) {

   const int VMin = 0;
   const int VMax = VectorFieldEdge.extent_int(Array::rank - 1) - 1;
   return setVectorEdge(Fun, VectorFieldEdge, EdgeComp, Geom, Mesh, VMin, VMax,
                        ExchangeHalosOpt, CartProjectionOpt);
}

template <class Reducer> Real reduceArray(const Array1DReal &Arr, int Extent0) {
   Real Res;
   Reducer R(Res);
   parallelReduce(
       {Extent0}, KOKKOS_LAMBDA(int I, Real &Accum) { R.join(Accum, Arr(I)); },
       R);

   return Res;
}

template <class Reducer> Real reduceArray(const Array1DReal &Arr) {
   return reduceArray<Reducer>(Arr, Arr.extent_int(0));
}

template <class Reducer, class VertMin, class VertMax>
Real reduceArray(const Array2DReal &Arr, int Extent0, const VertMin &VMin,
                 const VertMax &VMax) {
   Real Res;
   Reducer ROuter(Res);

   parallelReduceOuter(
       {Extent0},
       KOKKOS_LAMBDA(int I, const TeamMember &Team, Real &AccumOuter) {
          Real ResInner;
          Reducer RInner(ResInner);

          const int KMin   = getVertBound(VMin, I);
          const int KMax   = getVertBound(VMax, I);
          const int KRange = KMax - KMin + 1;

          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, Real &AccumInner) {
                 const int K = KMin + KOff;
                 RInner.join(AccumInner, Arr(I, K));
              },
              RInner);

          Kokkos::single(PerTeam(Team),
                         [&]() { ROuter.join(AccumOuter, ResInner); });
       },
       ROuter);

   return Res;
}

template <class Reducer> Real reduceArray(const Array2DReal &Arr, int Extent0) {
   return reduceArray<Reducer>(Arr, Extent0, 0, Arr.extent_int(1) - 1);
}

template <class Reducer> Real reduceArray(const Array2DReal &Arr) {
   return reduceArray<Reducer>(Arr, Arr.extent_int(0), 0,
                               Arr.extent_int(1) - 1);
}

template <class Reducer, class VertMin, class VertMax>
Real reduceArray(const Array3DReal &Arr, int Extent0, int Extent1,
                 const VertMin &VMin, const VertMax &VMax) {
   Real Res;
   Reducer ROuter(Res);

   parallelReduceOuter(
       {Extent0, Extent1},
       KOKKOS_LAMBDA(int L, int I, const TeamMember &Team, Real &AccumOuter) {
          Real ResInner;
          Reducer RInner(ResInner);

          const int KMin   = getVertBound(VMin, I);
          const int KMax   = getVertBound(VMax, I);
          const int KRange = KMax - KMin + 1;

          parallelReduceInner(
              Team, KRange,
              INNER_LAMBDA(int KOff, Real &AccumInner) {
                 const int K = KMin + KOff;
                 RInner.join(AccumInner, Arr(L, I, K));
              },
              RInner);

          Kokkos::single(PerTeam(Team),
                         [&]() { ROuter.join(AccumOuter, ResInner); });
       },
       ROuter);

   return Res;
}

template <class Reducer>
Real reduceArray(const Array3DReal &Arr, int Extent0, int Extent1) {
   return reduceArray<Reducer>(Arr, Extent0, Extent1, 0, Arr.extent_int(2) - 1);
}

template <class Reducer> Real reduceArray(const Array3DReal &Arr, int Extent0) {
   return reduceArray<Reducer>(Arr, Extent0, Arr.extent_int(1), 0,
                               Arr.extent_int(2) - 1);
}

template <class Reducer> Real reduceArray(const Array3DReal &Arr) {
   return reduceArray<Reducer>(Arr, Arr.extent_int(0), Arr.extent_int(1), 0,
                               Arr.extent_int(2) - 1);
}

template <class... ArgTypes> Real maxVal(ArgTypes &&...Args) {
   return reduceArray<Kokkos::Max<Real>>(std::forward<ArgTypes>(Args)...);
}

template <class... ArgTypes> Real sum(ArgTypes &&...Args) {
   return reduceArray<Kokkos::Sum<Real>>(std::forward<ArgTypes>(Args)...);
}

struct ErrorMeasures {
   Real LInf;
   Real L2;
};

// compute global normalized error measures based on the difference
// between two fields
template <class Array>
int computeErrors(ErrorMeasures &ErrorMeasures, const Array &NumFieldElement,
                  const Array &ExactFieldElement, const HorzMesh *Mesh,
                  MeshElement Element) {

   int Err = 0;

   int NElementsOwned;
   Array1DReal AreaElement;

   switch (Element) {
   case OnCell:
      NElementsOwned = Mesh->NCellsOwned;
      AreaElement    = Mesh->AreaCell;
      break;
   case OnVertex:
      NElementsOwned = Mesh->NVerticesOwned;
      AreaElement    = Mesh->AreaTriangle;
      break;
   case OnEdge:
      NElementsOwned = Mesh->NEdgesOwned;
      // need to compute areas associated with edges since we don't store those
      // in the mesh class
      {
         auto &DcEdge = Mesh->DcEdge;
         auto &DvEdge = Mesh->DvEdge;
         AreaElement  = Array1DReal("AreaEdge", Mesh->NEdgesOwned);
         parallelFor(
             {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
                AreaElement(IEdge) = DcEdge(IEdge) * DvEdge(IEdge) / 2;
             });
      }
      break;
   default:
      LOG_ERROR("computeErrors: element needs to be one of (OnCell, OnVertex, "
                "OnEdge)");
      return 1;
   }

   Array LInfElement;
   Array L2Element;
   Array LInfScaleElement;
   Array L2ScaleElement;

   // Compute element-wise errors
   if constexpr (Array::rank == 1) {
      LInfElement = Array("LInfElement", NElementsOwned);
      L2Element   = Array("L2Element", NElementsOwned);

      LInfScaleElement = Array("LInfScaleElement", NElementsOwned);
      L2ScaleElement   = Array("L2ScaleElement", NElementsOwned);

      parallelFor(
          {NElementsOwned}, KOKKOS_LAMBDA(int IElement) {
             const Real NumValElement   = NumFieldElement(IElement);
             const Real ExactValElement = ExactFieldElement(IElement);

             // Errors
             LInfElement(IElement) = std::abs(NumValElement - ExactValElement);
             LInfScaleElement(IElement) = std::abs(ExactValElement);
             L2Element(IElement)        = AreaElement(IElement) *
                                   LInfElement(IElement) *
                                   LInfElement(IElement);
             L2ScaleElement(IElement) = AreaElement(IElement) *
                                        LInfScaleElement(IElement) *
                                        LInfScaleElement(IElement);
          });
   }
   if constexpr (Array::rank == 2) {
      const int NVertLayers = NumFieldElement.extent_int(1);

      LInfElement = Array("LInfElement", NElementsOwned, NVertLayers);
      L2Element   = Array("L2Element", NElementsOwned, NVertLayers);

      LInfScaleElement = Array("LInfScaleElement", NElementsOwned, NVertLayers);
      L2ScaleElement   = Array("L2ScaleElement", NElementsOwned, NVertLayers);

      parallelFor(
          {NElementsOwned, NVertLayers}, KOKKOS_LAMBDA(int IElement, int K) {
             const Real NumValElement   = NumFieldElement(IElement, K);
             const Real ExactValElement = ExactFieldElement(IElement, K);

             // Errors
             LInfElement(IElement, K) =
                 std::abs(NumValElement - ExactValElement);
             LInfScaleElement(IElement, K) = std::abs(ExactValElement);
             L2Element(IElement, K)        = AreaElement(IElement) *
                                      LInfElement(IElement, K) *
                                      LInfElement(IElement, K);
             L2ScaleElement(IElement, K) = AreaElement(IElement) *
                                           LInfScaleElement(IElement, K) *
                                           LInfScaleElement(IElement, K);
          });
   }

   if constexpr (Array::rank == 3) {
      const int NTracers    = NumFieldElement.extent_int(0);
      const int NVertLayers = NumFieldElement.extent_int(2);

      LInfElement = Array("LInfElement", NTracers, NElementsOwned, NVertLayers);
      L2Element   = Array("L2Element", NTracers, NElementsOwned, NVertLayers);

      LInfScaleElement =
          Array("LInfScaleElement", NTracers, NElementsOwned, NVertLayers);
      L2ScaleElement =
          Array("L2ScaleElement", NTracers, NElementsOwned, NVertLayers);

      parallelFor(
          {NTracers, NElementsOwned, NVertLayers},
          KOKKOS_LAMBDA(int L, int IElement, int K) {
             const Real NumValElement   = NumFieldElement(L, IElement, K);
             const Real ExactValElement = ExactFieldElement(L, IElement, K);

             // Errors
             LInfElement(L, IElement, K) =
                 std::abs(NumValElement - ExactValElement);
             LInfScaleElement(L, IElement, K) = std::abs(ExactValElement);
             L2Element(L, IElement, K)        = AreaElement(IElement) *
                                         LInfElement(L, IElement, K) *
                                         LInfElement(L, IElement, K);
             L2ScaleElement(L, IElement, K) = AreaElement(IElement) *
                                              LInfScaleElement(L, IElement, K) *
                                              LInfScaleElement(L, IElement, K);
          });
   }

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfElement);
   const Real L2ErrorLoc   = sum(L2Element);
   const Real LInfScaleLoc = maxVal(LInfScaleElement);
   const Real L2ScaleLoc   = sum(L2ScaleElement);

   MPI_Comm Comm = MachEnv::getDefault()->getComm();
   Real LInfError, LInfScale;
   Err +=
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err +=
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err += MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err += MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);

   if (Err != 0)
      LOG_ERROR("computeErrors: MPI Allreduce error");

   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   ErrorMeasures.L2   = L2Error;
   ErrorMeasures.LInf = LInfError;

   return Err;
}

inline int checkErrors(const std::string &TestSuite,
                       const std::string &Variable, const ErrorMeasures &Errors,
                       const ErrorMeasures &ExpectedErrors, Real RTol,
                       Real ATol = 0) {
   int Err = 0;
   if (!isApprox(Errors.LInf, ExpectedErrors.LInf, RTol, ATol)) {
      Err++;
      LOG_ERROR("{}: {} LInf FAIL, expected {}, got {}", TestSuite, Variable,
                ExpectedErrors.LInf, Errors.LInf);
   }
   if (!isApprox(Errors.L2, ExpectedErrors.L2, RTol, ATol)) {
      Err++;
      LOG_ERROR("{}: {} L2 FAIL, expected {}, got {}", TestSuite, Variable,
                ExpectedErrors.L2, Errors.L2);
   }
   return Err;
}

} // namespace OMEGA
#endif
