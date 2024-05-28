#ifndef OMEGA_OCEAN_TEST_COMMON_H
#define OMEGA_OCEAN_TEST_COMMON_H

#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
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

// set scalar field on cells based on analytical formula and optionally
// exchange halos
template <class Functor>
int setScalarCell(const Functor &Fun, const Array2DReal &ScalarCell,
                  Geometry Geom, const HorzMesh *Mesh, int NVertLevels,
                  bool ExchangeHalos = true) {
   int Err = 0;

   auto XCell = createDeviceMirrorCopy(Mesh->XCellH);
   auto YCell = createDeviceMirrorCopy(Mesh->YCellH);

   auto LonCell = createDeviceMirrorCopy(Mesh->LonCellH);
   auto LatCell = createDeviceMirrorCopy(Mesh->LatCellH);

   parallelFor(
       {Mesh->NCellsOwned, NVertLevels}, KOKKOS_LAMBDA(int ICell, int K) {
          if (Geom == Geometry::Planar) {
             const Real X         = XCell(ICell);
             const Real Y         = YCell(ICell);
             ScalarCell(ICell, K) = Fun(X, Y);
          } else {
             const Real Lon       = LonCell(ICell);
             const Real Lat       = LatCell(ICell);
             ScalarCell(ICell, K) = Fun(Lon, Lat);
          }
       });

   if (ExchangeHalos) {
      auto MyHalo      = Halo::getDefault();
      auto ScalarCellH = createHostMirrorCopy(ScalarCell);
      Err              = MyHalo->exchangeFullArrayHalo(ScalarCellH, OnCell);
      if (Err != 0)
         LOG_ERROR("setScalarCell: error in halo exchange");
      deepCopy(ScalarCell, ScalarCellH);
   }
   return Err;
}

// set scalar field on vertices based on analytical formula and optionally
// exchange halos
template <class Functor>
int setScalarVertex(const Functor &Fun, const Array2DReal &ScalarVertex,
                    Geometry Geom, const HorzMesh *Mesh, int NVertLevels,
                    bool ExchangeHalos = true) {

   int Err = 0;

   auto XVertex = createDeviceMirrorCopy(Mesh->XVertexH);
   auto YVertex = createDeviceMirrorCopy(Mesh->YVertexH);

   auto LonVertex = createDeviceMirrorCopy(Mesh->LonVertexH);
   auto LatVertex = createDeviceMirrorCopy(Mesh->LatVertexH);

   parallelFor(
       {Mesh->NVerticesOwned, NVertLevels}, KOKKOS_LAMBDA(int IVertex, int K) {
          if (Geom == Geometry::Planar) {
             const Real X             = XVertex(IVertex);
             const Real Y             = YVertex(IVertex);
             ScalarVertex(IVertex, K) = Fun(X, Y);
          } else {
             const Real Lon           = LonVertex(IVertex);
             const Real Lat           = LatVertex(IVertex);
             ScalarVertex(IVertex, K) = Fun(Lon, Lat);
          }
       });

   if (ExchangeHalos) {
      auto MyHalo        = Halo::getDefault();
      auto ScalarVertexH = createHostMirrorCopy(ScalarVertex);
      Err = MyHalo->exchangeFullArrayHalo(ScalarVertexH, OnVertex);
      if (Err != 0)
         LOG_ERROR("setScalarVertex: error in halo exchange");
      deepCopy(ScalarVertex, ScalarVertexH);
   }
   return Err;
}

// set scalar field on edges based on analytical formula and optionally
// exchange halos
template <class Functor>
int setScalarEdge(const Functor &Fun, const Array2DReal &ScalarFieldEdge,
                  Geometry Geom, const HorzMesh *Mesh, int NVertLevels,
                  bool ExchangeHalos = true) {

   int Err = 0;

   auto XEdge = createDeviceMirrorCopy(Mesh->XEdgeH);
   auto YEdge = createDeviceMirrorCopy(Mesh->YEdgeH);

   auto LonEdge = createDeviceMirrorCopy(Mesh->LonEdgeH);
   auto LatEdge = createDeviceMirrorCopy(Mesh->LatEdgeH);

   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int K) {
          Real VecFieldEdge;
          if (Geom == Geometry::Planar) {
             const Real XE = XEdge(IEdge);
             const Real YE = YEdge(IEdge);

             ScalarFieldEdge(IEdge, K) = Fun(XE, YE);
          } else {
             const Real LonE = LonEdge(IEdge);
             const Real LatE = LatEdge(IEdge);

             ScalarFieldEdge(IEdge, K) = Fun(LonE, LatE);
          }
       });

   if (ExchangeHalos) {
      auto MyHalo           = Halo::getDefault();
      auto ScalarFieldEdgeH = createHostMirrorCopy(ScalarFieldEdge);
      Err = MyHalo->exchangeFullArrayHalo(ScalarFieldEdgeH, OnEdge);
      if (Err != 0)
         LOG_ERROR("setScalarEdge: error in halo exchange");
      deepCopy(ScalarFieldEdge, ScalarFieldEdgeH);
   }

   return Err;
}

// set vector field on edges based on analytical formula and optionally
// exchange halos
template <class Functor>
int setVectorEdge(const Functor &Fun, const Array2DReal &VectorFieldEdge,
                  EdgeComponent EdgeComp, Geometry Geom, const HorzMesh *Mesh,
                  int NVertLevels, bool ExchangeHalos = true) {

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

             bool UseCartesianProjection = true;

             if (UseCartesianProjection) {
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

   if (ExchangeHalos) {
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

inline Real sum(const Array2DReal &Arr) {
   Real Sum;

   parallelReduce(
       {Arr.extent_int(0), Arr.extent_int(1)},
       KOKKOS_LAMBDA(int I, int J, Real &Accum) { Accum += Arr(I, J); }, Sum);

   return Sum;
}

struct ErrorMeasures {
   Real LInf;
   Real L2;
};

// compute global normalized error measures based on the difference
// between two cell fields
inline int computeErrorsCell(ErrorMeasures &ErrorMeasures,
                             const Array2DReal &NumFieldCell,
                             const Array2DReal &ExactFieldCell,
                             const HorzMesh *Mesh, int NVertLevels) {

   int Err = 0;

   auto &AreaCell = Mesh->AreaCell;

   // Compute element-wise errors
   Array2DReal LInfCell("LInfCell", Mesh->NCellsOwned, NVertLevels);
   Array2DReal L2Cell("L2Cell", Mesh->NCellsOwned, NVertLevels);

   Array2DReal LInfScaleCell("LInfScaleCell", Mesh->NCellsOwned, NVertLevels);
   Array2DReal L2ScaleCell("L2ScaleCell", Mesh->NCellsOwned, NVertLevels);

   parallelFor(
       {Mesh->NCellsOwned, NVertLevels}, KOKKOS_LAMBDA(int ICell, int K) {
          const Real NumValCell   = NumFieldCell(ICell, K);
          const Real ExactValCell = ExactFieldCell(ICell, K);

          // Errors
          LInfCell(ICell, K)      = std::abs(NumValCell - ExactValCell);
          LInfScaleCell(ICell, K) = std::abs(ExactValCell);
          L2Cell(ICell, K) =
              AreaCell(ICell) * LInfCell(ICell, K) * LInfCell(ICell, K);
          L2ScaleCell(ICell, K) = AreaCell(ICell) * LInfScaleCell(ICell, K) *
                                  LInfScaleCell(ICell, K);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfCell);
   const Real L2ErrorLoc   = sum(L2Cell);
   const Real LInfScaleLoc = maxVal(LInfScaleCell);
   const Real L2ScaleLoc   = sum(L2ScaleCell);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
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
      LOG_ERROR("computeErrorsCell: MPI Allreduce error");

   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   ErrorMeasures.L2   = L2Error;
   ErrorMeasures.LInf = LInfError;

   return Err;
}

// compute global normalized error measures based on the difference
// between two vertex fields
inline int computeErrorsVertex(ErrorMeasures &ErrorMeasures,
                               const Array2DReal &NumFieldVertex,
                               const Array2DReal &ExactFieldVertex,
                               const HorzMesh *Mesh, int NVertLevels) {

   int Err = 0;

   const auto &AreaTriangle = Mesh->AreaTriangle;

   // Compute element-wise errors
   Array2DReal LInfVertex("LInfVertex", Mesh->NVerticesOwned, NVertLevels);
   Array2DReal LInfScaleVertex("LInfScaleVertex", Mesh->NVerticesOwned,
                               NVertLevels);
   Array2DReal L2Vertex("L2Vertex", Mesh->NVerticesOwned, NVertLevels);
   Array2DReal L2ScaleVertex("L2ScaleVertex", Mesh->NVerticesOwned,
                             NVertLevels);
   parallelFor(
       {Mesh->NVerticesOwned, NVertLevels}, KOKKOS_LAMBDA(int IVertex, int K) {
          const Real NumValVertex   = NumFieldVertex(IVertex, K);
          const Real ExactValVertex = ExactFieldVertex(IVertex, K);

          // Errors
          LInfVertex(IVertex, K)      = std::abs(NumValVertex - ExactValVertex);
          LInfScaleVertex(IVertex, K) = std::abs(ExactValVertex);
          L2Vertex(IVertex, K)        = AreaTriangle(IVertex) *
                                 LInfVertex(IVertex, K) *
                                 LInfVertex(IVertex, K);
          L2ScaleVertex(IVertex, K) = AreaTriangle(IVertex) *
                                      LInfScaleVertex(IVertex, K) *
                                      LInfScaleVertex(IVertex, K);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfVertex);
   const Real LInfScaleLoc = maxVal(LInfScaleVertex);
   const Real L2ErrorLoc   = sum(L2Vertex);
   const Real L2ScaleLoc   = sum(L2ScaleVertex);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
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
      LOG_ERROR("computeErrorsCell: MPI Allreduce error");

   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   ErrorMeasures.LInf = LInfError;
   ErrorMeasures.L2   = L2Error;

   return Err;
}

// compute global normalized error measures based on the difference
// between two edge fields
inline int computeErrorsEdge(ErrorMeasures &ErrMeasures,
                             const Array2DReal &FieldEdge,
                             const Array2DReal &ExactFieldEdge,
                             const HorzMesh *Mesh, int NVertLevels) {

   int Err = 0;

   const auto &DcEdge = Mesh->DcEdge;
   const auto &DvEdge = Mesh->DvEdge;

   // Compute element-wise errors
   Array2DReal LInfEdge("LInfEdge", Mesh->NEdgesOwned, NVertLevels);
   Array2DReal L2Edge("L2Edge", Mesh->NEdgesOwned, NVertLevels);
   Array2DReal LInfScaleEdge("LInfScaleEdge", Mesh->NEdgesOwned, NVertLevels);
   Array2DReal L2ScaleEdge("L2ScaleEdge", Mesh->NEdgesOwned, NVertLevels);

   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int K) {
          const Real NumValEdge   = FieldEdge(IEdge, K);
          const Real ExactValEdge = ExactFieldEdge(IEdge, K);

          LInfEdge(IEdge, K)      = std::abs(NumValEdge - ExactValEdge);
          LInfScaleEdge(IEdge, K) = std::abs(ExactValEdge);
          const Real AreaEdge     = DcEdge(IEdge) * DvEdge(IEdge) / 2;
          L2Edge(IEdge, K) = AreaEdge * LInfEdge(IEdge, K) * LInfEdge(IEdge, K);
          L2ScaleEdge(IEdge, K) =
              AreaEdge * LInfScaleEdge(IEdge, K) * LInfScaleEdge(IEdge, K);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfEdge);
   const Real LInfScaleLoc = maxVal(LInfScaleEdge);
   const Real L2ErrorLoc   = sum(L2Edge);
   const Real L2ScaleLoc   = sum(L2ScaleEdge);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();

   Real LInfError, LInfScale;
   Err +=
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err +=
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);

   Real L2Error, L2Scale;
   Err += MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err += MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);

   if (Err != 0)
      LOG_ERROR("computeErrorsCell: MPI Allreduce error");

   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   ErrMeasures.L2   = L2Error;
   ErrMeasures.LInf = LInfError;
   return Err;
}

} // namespace OMEGA
#endif
