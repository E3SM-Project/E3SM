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
inline bool isApprox(Real X, Real Y, Real RTol) {
   return std::abs(X - Y) <= RTol * std::max(std::abs(X), std::abs(Y));
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

// set scalar field on chosen elements (cells/vertices/edges) based on
// analytical formula and optionally exchange halos
template <class Functor>
int setScalar(const Functor &Fun, const Array2DReal &ScalarElement,
              Geometry Geom, const HorzMesh *Mesh, MeshElement Element,
              int NVertLevels,
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

   parallelFor(
       {NElementsOwned, NVertLevels}, KOKKOS_LAMBDA(int IElement, int K) {
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

   if (ExchangeHalosOpt == ExchangeHalos::Yes) {
      auto MyHalo         = Halo::getDefault();
      auto ScalarElementH = createHostMirrorCopy(ScalarElement);
      Err = MyHalo->exchangeFullArrayHalo(ScalarElementH, Element);
      if (Err != 0)
         LOG_ERROR("setScalarElement: error in halo exchange");
      deepCopy(ScalarElement, ScalarElementH);
   }
   return Err;
}

enum class CartProjection { Yes, No };

// set vector field on edges based on analytical formula and optionally
// exchange halos
template <class Functor>
int setVectorEdge(const Functor &Fun, const Array2DReal &VectorFieldEdge,
                  EdgeComponent EdgeComp, Geometry Geom, const HorzMesh *Mesh,
                  int NVertLevels,
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

   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int K) {
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
          VectorFieldEdge(IEdge, K) = VecFieldEdge;
       });

   if (ExchangeHalosOpt == ExchangeHalos::Yes) {
      auto MyHalo           = Halo::getDefault();
      auto VectorFieldEdgeH = createHostMirrorCopy(VectorFieldEdge);
      Err = MyHalo->exchangeFullArrayHalo(VectorFieldEdgeH, OnEdge);
      if (Err != 0)
         LOG_ERROR("setVectorEdge: error in halo exchange");
      deepCopy(VectorFieldEdge, VectorFieldEdgeH);
   }
   return Err;
}

inline Real maxVal(const Array2DReal &Arr) {
   Real MaxVal;

   parallelReduce(
       {Arr.extent_int(0), Arr.extent_int(1)},
       KOKKOS_LAMBDA(int I, int J, Real &Accum) {
          Accum = Kokkos::max(Arr(I, J), Accum);
       },
       Kokkos::Max<Real>(MaxVal));

   return MaxVal;
}

inline Real sum(const Array2DReal &Arr, int Extent0, int Extent1) {
   Real Sum;

   parallelReduce(
       {Extent0, Extent1},
       KOKKOS_LAMBDA(int I, int J, Real &Accum) { Accum += Arr(I, J); }, Sum);

   return Sum;
}

inline Real sum(const Array2DReal &Arr) {
   return sum(Arr, Arr.extent_int(0), Arr.extent_int(1));
}

inline Real sum(const Array2DReal &Arr, int Extent0) {
   return sum(Arr, Extent0, Arr.extent_int(1));
}

struct ErrorMeasures {
   Real LInf;
   Real L2;
};

// compute global normalized error measures based on the difference
// between two fields
inline int computeErrors(ErrorMeasures &ErrorMeasures,
                         const Array2DReal &NumFieldElement,
                         const Array2DReal &ExactFieldElement,
                         const HorzMesh *Mesh, MeshElement Element,
                         int NVertLevels) {

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

   // Compute element-wise errors
   Array2DReal LInfElement("LInfElement", NElementsOwned, NVertLevels);
   Array2DReal L2Element("L2Element", NElementsOwned, NVertLevels);

   Array2DReal LInfScaleElement("LInfScaleElement", NElementsOwned,
                                NVertLevels);
   Array2DReal L2ScaleElement("L2ScaleElement", NElementsOwned, NVertLevels);

   parallelFor(
       {NElementsOwned, NVertLevels}, KOKKOS_LAMBDA(int IElement, int K) {
          const Real NumValElement   = NumFieldElement(IElement, K);
          const Real ExactValElement = ExactFieldElement(IElement, K);

          // Errors
          LInfElement(IElement, K) = std::abs(NumValElement - ExactValElement);
          LInfScaleElement(IElement, K) = std::abs(ExactValElement);
          L2Element(IElement, K)        = AreaElement(IElement) *
                                   LInfElement(IElement, K) *
                                   LInfElement(IElement, K);
          L2ScaleElement(IElement, K) = AreaElement(IElement) *
                                        LInfScaleElement(IElement, K) *
                                        LInfScaleElement(IElement, K);
       });

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
                       const ErrorMeasures &ExpectedErrors, Real RTol) {
   int Err = 0;
   if (!isApprox(Errors.LInf, ExpectedErrors.LInf, RTol)) {
      Err++;
      LOG_ERROR("{}: {} LInf FAIL, expected {}, got {}", TestSuite, Variable,
                ExpectedErrors.LInf, Errors.LInf);
   }
   if (!isApprox(Errors.L2, ExpectedErrors.L2, RTol)) {
      Err++;
      LOG_ERROR("{}: {} L2 FAIL, expected {}, got {}", TestSuite, Variable,
                ExpectedErrors.L2, Errors.L2);
   }
   return Err;
}

} // namespace OMEGA
#endif
