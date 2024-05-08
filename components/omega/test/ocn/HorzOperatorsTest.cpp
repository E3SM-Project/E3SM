#include "HorzOperators.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "mpi.h"

#include <cmath>

using namespace OMEGA;

// check if two real numbers are equal with a given relative tolerance
bool isApprox(Real X, Real Y, Real RTol) {
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

enum class EdgeOrientation { Normal, Tangential };

template <class Functor>
void computeVecFieldEdge(const Functor &Fun, const Array1DReal &VecFieldArr,
                         EdgeOrientation EdgeOrient, const HorzMesh *Mesh) {
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
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          Real VecFieldEdge;
#ifdef HORZOPERATORS_TEST_PLANE
          const Real XE = XEdge(IEdge);
          const Real YE = YEdge(IEdge);

          Real VecField[2];
          Fun(VecField, XE, YE);

          if (EdgeOrient == EdgeOrientation::Normal) {
             const Real EdgeNormalX = std::cos(AngleEdge(IEdge));
             const Real EdgeNormalY = std::sin(AngleEdge(IEdge));
             VecFieldEdge =
                 EdgeNormalX * VecField[0] + EdgeNormalY * VecField[1];
          }

          if (EdgeOrient == EdgeOrientation::Tangential) {
             const Real EdgeTangentX = -std::sin(AngleEdge(IEdge));
             const Real EdgeTangentY = std::cos(AngleEdge(IEdge));
             VecFieldEdge =
                 EdgeTangentX * VecField[0] + EdgeTangentY * VecField[1];
          }
#else
          const Real LonE                = LonEdge(IEdge);
          const Real LatE                = LatEdge(IEdge);

          Real VecField[2];
          Fun(VecField, LonE, LatE);

          bool UseCartesianProjection = true;

          if (UseCartesianProjection) {
            Real VecFieldCart[3];
            sphereToCartVec(VecFieldCart, VecField, LonE, LatE);

            const Real EdgeCoords[3] = {XEdge[IEdge], YEdge[IEdge], ZEdge[IEdge]};

            if (EdgeOrient == EdgeOrientation::Normal) {
              const int JCell1 = CellsOnEdge(IEdge, 1);
              const Real CellCoords[3] = {XCell(JCell1), YCell(JCell1), ZCell(JCell1)};

              Real EdgeNormal[3];
              tangentVector(EdgeNormal, EdgeCoords, CellCoords);
              VecFieldEdge = EdgeNormal[0] * VecFieldCart[0] +
                             EdgeNormal[1] * VecFieldCart[1] +
                             EdgeNormal[2] * VecFieldCart[2];
            }

            if (EdgeOrient == EdgeOrientation::Tangential) {
              const int JVertex1 = VerticesOnEdge(IEdge, 1);
              const Real VertexCoords[3] = {XVertex(JVertex1), YVertex(JVertex1), ZVertex(JVertex1)};

              Real EdgeTangent[3];
              tangentVector(EdgeTangent, EdgeCoords, VertexCoords);
              VecFieldEdge = EdgeTangent[0] * VecFieldCart[0] +
                             EdgeTangent[1] * VecFieldCart[1] +
                             EdgeTangent[2] * VecFieldCart[2];
            }
          } else {
            if (EdgeOrient == EdgeOrientation::Normal) {
              const Real EdgeNormalX      = std::cos(AngleEdge(IEdge));
              const Real EdgeNormalY      = std::sin(AngleEdge(IEdge));
              VecFieldEdge = EdgeNormalX * VecField[0] + EdgeNormalY * VecField[1];
            }

            if (EdgeOrient == EdgeOrientation::Tangential) {
              const Real EdgeTangentX  = -std::sin(AngleEdge(IEdge));
              const Real EdgeTangentY  = std::cos(AngleEdge(IEdge));
              VecFieldEdge = EdgeTangentX * VecField[0] + EdgeTangentY * VecField[1];
            }
          }
#endif
          VecFieldArr[IEdge] = VecFieldEdge;
       });
}

// temporary replacement for YAKL intrinsics
Real maxVal(const Array1DReal &Arr) {
   Real MaxVal;

   parallelReduce(
       {Arr.extent_int(0)},
       KOKKOS_LAMBDA(int I, Real &Accum) {
          Accum = Kokkos::max(Arr(I), Accum);
       },
       Kokkos::Max<Real>(MaxVal));

   return MaxVal;
}

Real sum(const Array1DReal &Arr) {
   Real Sum;

   parallelReduce(
       {Arr.extent_int(0)},
       KOKKOS_LAMBDA(int I, Real &Accum) { Accum += Arr(I); }, Sum);

   return Sum;
}

#ifdef HORZOPERATORS_TEST_PLANE
// analytical expressions for scalar and vector fields used as input for
// operator tests together with exact values of operators for computing errors
// and expected error values
struct TestSetup {
   Real Pi = M_PI;

   // lengths of periodic planar mesh
   // TODO: get this from the horizontal mesh once it supports periodic planar
   // meshes
   Real Lx = 1;
   Real Ly = std::sqrt(3) / 2;

   Real ExpectedDivErrorLInf   = 0.00124886886594427027;
   Real ExpectedDivErrorL2     = 0.00124886886590974385;
   Real ExpectedGradErrorLInf  = 0.00125026071878537952;
   Real ExpectedGradErrorL2    = 0.00134354611117262204;
   Real ExpectedCurlErrorLInf  = 0.161365663569699946;
   Real ExpectedCurlErrorL2    = 0.161348016897141039;
   Real ExpectedReconErrorLInf = 0.00450897496974901352;
   Real ExpectedReconErrorL2   = 0.00417367308684470691;

   KOKKOS_FUNCTION Real exactScalar(Real X, Real Y) const {
      return std::sin(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactGradScalarX(Real X, Real Y) const {
      return 2 * Pi / Lx * std::cos(2 * Pi * X / Lx) *
             std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactGradScalarY(Real X, Real Y) const {
      return 2 * Pi / Ly * std::sin(2 * Pi * X / Lx) *
             std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactVecX(Real X, Real Y) const {
      return std::sin(2 * Pi * X / Lx) * std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactVecY(Real X, Real Y) const {
      return std::cos(2 * Pi * X / Lx) * std::sin(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactDivVec(Real X, Real Y) const {
      return 2 * Pi * (1. / Lx + 1. / Ly) * std::cos(2 * Pi * X / Lx) *
             std::cos(2 * Pi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactCurlVec(Real X, Real Y) const {
      return 2 * Pi * (-1. / Lx + 1. / Ly) * std::sin(2 * Pi * X / Lx) *
             std::sin(2 * Pi * Y / Ly);
   }
};
#endif

#ifdef HORZOPERATORS_TEST_SPHERE_1
// analytical expressions for scalar and vector fields used as input for
// operator tests together with exact values of operators for computing errors
// and expected error values
struct TestSetup {
   // radius of spherical mesh
   // TODO: get this from the mesh
   Real Radius = 6371220;

   Real ExpectedDivErrorLInf   = 0.013659577398978353;
   Real ExpectedDivErrorL2     = 0.00367052484586382743;
   Real ExpectedGradErrorLInf  = 0.00187912292540628936;
   Real ExpectedGradErrorL2    = 0.00149841802817334306;
   Real ExpectedCurlErrorLInf  = 0.0271404735181308317;
   Real ExpectedCurlErrorL2    = 0.025202316610921989;
   Real ExpectedReconErrorLInf = 0.0206375134079833517;
   Real ExpectedReconErrorL2   = 0.00692590524910695858;

   KOKKOS_FUNCTION Real exactScalar(Real Lon, Real Lat) const {
      return Radius * std::cos(Lon) * std::pow(std::cos(Lat), 4);
   }

   KOKKOS_FUNCTION Real exactGradScalarX(Real Lon, Real Lat) const {
      return -std::sin(Lon) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_FUNCTION Real exactGradScalarY(Real Lon, Real Lat) const {
      return -4 * std::cos(Lon) * std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_FUNCTION Real exactVecX(Real Lon, Real Lat) const {
      return -Radius * std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_FUNCTION Real exactVecY(Real Lon, Real Lat) const {
      return -4 * Radius * std::sin(Lon) * std::cos(Lon) *
             std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_FUNCTION Real exactDivVec(Real Lon, Real Lat) const {
      return std::sin(Lon) * std::cos(Lon) * std::pow(std::cos(Lat), 2) *
             (20 * std::pow(std::sin(Lat), 2) - 6);
   }

   KOKKOS_FUNCTION Real exactCurlVec(Real Lon, Real Lat) const {
      return -4 * std::pow(std::cos(Lon), 2) * std::pow(std::cos(Lat), 2) *
             std::sin(Lat);
   }
};
#endif

#ifdef HORZOPERATORS_TEST_SPHERE_2
// analytical expressions for scalar and vector fields used as input for
// operator tests together with exact values of operators for computing errors
// and expected error values
struct TestSetup {
   // radius of spherical mesh
   // TODO: get this from the mesh
   Real Radius = 6371220;

   Real ExpectedDivErrorLInf   = 1.37734693033362766e-10;
   Real ExpectedDivErrorL2     = 0.000484370621558727582;
   Real ExpectedGradErrorLInf  = 0.000906351303388669991;
   Real ExpectedGradErrorL2    = 0.000949206041390823676;
   Real ExpectedCurlErrorLInf  = 0.00433205620592059647;
   Real ExpectedCurlErrorL2    = 0.00204725417666192042;
   Real ExpectedReconErrorLInf = 0.0254271921029878764;
   Real ExpectedReconErrorL2   = 0.00419630561428921064;

   KOKKOS_FUNCTION Real exactScalar(Real Lon, Real Lat) const {
      return -Radius * std::pow(std::sin(Lat), 2);
   }

   KOKKOS_FUNCTION Real exactGradScalarX(Real Lon, Real Lat) const { return 0; }

   KOKKOS_FUNCTION Real exactGradScalarY(Real Lon, Real Lat) const {
      return -2 * std::sin(Lat) * std::cos(Lat);
   }

   KOKKOS_FUNCTION Real exactVecX(Real Lon, Real Lat) const {
      return std::cos(Lat);
   }

   KOKKOS_FUNCTION Real exactVecY(Real Lon, Real Lat) const { return 0; }

   KOKKOS_FUNCTION Real exactDivVec(Real Lon, Real Lat) const { return 0; }

   KOKKOS_FUNCTION Real exactCurlVec(Real Lon, Real Lat) const {
      return 2 * std::sin(Lat) / Radius;
   }
};
#endif

int testDivergence(Real RTol) {
   int Err;
   TestSetup Setup;

   const auto &Mesh = HorzMesh::getDefault();

   // Prepare operator input
   Array1DReal VecEdge("VecEdge", Mesh->NEdgesSize);

   computeVecFieldEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeOrientation::Normal, Mesh);

   // Perform halo exchange
   auto MyHalo   = Halo::getDefault();
   auto VecEdgeH = createHostMirrorCopy(VecEdge);
   MyHalo->exchangeFullArrayHalo(VecEdgeH, OnEdge);
   deepCopy(VecEdge, VecEdgeH);

#ifdef HORZOPERATORS_TEST_PLANE
   auto XCell = createDeviceMirrorCopy(Mesh->XCellH);
   auto YCell = createDeviceMirrorCopy(Mesh->YCellH);
#else
   auto XCell = createDeviceMirrorCopy(Mesh->LonCellH);
   auto YCell = createDeviceMirrorCopy(Mesh->LatCellH);
#endif
   auto &AreaCell = Mesh->AreaCell;

   // Compute element-wise errors
   Array1DReal LInfCell("LInfCell", Mesh->NCellsOwned);
   Array1DReal L2Cell("L2Cell", Mesh->NCellsOwned);

   Array1DReal LInfScaleCell("LInfScaleCell", Mesh->NCellsOwned);
   Array1DReal L2ScaleCell("L2ScaleCell", Mesh->NCellsOwned);
   DivergenceOnCell DivergenceCell(Mesh);
   parallelFor(
       {Mesh->NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          // Numerical result
          const Real DivCellNum = DivergenceCell(ICell, VecEdge);

          // Exact result
          const Real X            = XCell(ICell);
          const Real Y            = YCell(ICell);
          const Real DivCellExact = Setup.exactDivVec(X, Y);

          // Errors
          LInfCell(ICell)      = std::abs(DivCellNum - DivCellExact);
          LInfScaleCell(ICell) = std::abs(DivCellExact);
          L2Cell(ICell) = AreaCell(ICell) * LInfCell(ICell) * LInfCell(ICell);
          L2ScaleCell(ICell) =
              AreaCell(ICell) * LInfScaleCell(ICell) * LInfScaleCell(ICell);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfCell);
   const Real L2ErrorLoc   = sum(L2Cell);
   const Real LInfScaleLoc = maxVal(LInfScaleCell);
   const Real L2ScaleLoc   = sum(L2ScaleCell);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
   Real LInfError, LInfScale;
   Err =
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err =
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err = MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err = MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);
   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   // Check error values
   if (Err == 0 && isApprox(LInfError, Setup.ExpectedDivErrorLInf, RTol) &&
       isApprox(L2Error, Setup.ExpectedDivErrorL2, RTol)) {
      return 0;
   } else {
      return 1;
   }
}

int testGradient(Real RTol) {
   int Err;
   TestSetup Setup;

   const auto &Mesh = HorzMesh::getDefault();
#ifdef HORZOPERATORS_TEST_PLANE
   const auto XCell = createDeviceMirrorCopy(Mesh->XCellH);
   const auto YCell = createDeviceMirrorCopy(Mesh->YCellH);
#else
   const auto XCell = createDeviceMirrorCopy(Mesh->LonCellH);
   const auto YCell = createDeviceMirrorCopy(Mesh->LatCellH);
#endif

   // Prepare operator input
   Array1DReal ScalarCell("ScalarCell", Mesh->NCellsSize);
   parallelFor(
       {Mesh->NCellsOwned}, KOKKOS_LAMBDA(int ICell) {
          const Real X      = XCell(ICell);
          const Real Y      = YCell(ICell);
          ScalarCell(ICell) = Setup.exactScalar(X, Y);
       });

   // Perform halo exchange
   auto MyHalo      = Halo::getDefault();
   auto ScalarCellH = createHostMirrorCopy(ScalarCell);
   MyHalo->exchangeFullArrayHalo(ScalarCellH, OnCell);
   deepCopy(ScalarCell, ScalarCellH);

   // Compute exact result
   Array1DReal ExactGradEdge("ExactGradEdge", Mesh->NEdgesOwned);
   computeVecFieldEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactGradScalarX(X, Y);
          VecField[1] = Setup.exactGradScalarY(X, Y);
       },
       ExactGradEdge, EdgeOrientation::Normal, Mesh);

   const auto &DcEdge = Mesh->DcEdge;
   const auto &DvEdge = Mesh->DvEdge;
   // Compute element-wise errors
   Array1DReal LInfEdge("LInfEdge", Mesh->NEdgesOwned);
   Array1DReal L2Edge("L2Edge", Mesh->NEdgesOwned);
   Array1DReal LInfScaleEdge("LInfScaleEdge", Mesh->NEdgesOwned);
   Array1DReal L2ScaleEdge("L2ScaleEdge", Mesh->NEdgesOwned);
   GradientOnEdge GradientEdge(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          // Numerical result
          const Real GradScalarNum = GradientEdge(IEdge, ScalarCell);

          // Exact result
          const Real GradScalarExact = ExactGradEdge(IEdge);

          LInfEdge(IEdge)      = std::abs(GradScalarNum - GradScalarExact);
          LInfScaleEdge(IEdge) = std::abs(GradScalarExact);
          const Real AreaEdge  = DcEdge(IEdge) * DvEdge(IEdge) / 2;
          L2Edge(IEdge)        = AreaEdge * LInfEdge(IEdge) * LInfEdge(IEdge);
          L2ScaleEdge(IEdge) =
              AreaEdge * LInfScaleEdge(IEdge) * LInfScaleEdge(IEdge);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfEdge);
   const Real LInfScaleLoc = maxVal(LInfScaleEdge);
   const Real L2ErrorLoc   = sum(L2Edge);
   const Real L2ScaleLoc   = sum(L2ScaleEdge);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
   Real LInfError, LInfScale;
   Err =
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err =
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err = MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err = MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);
   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   // Check error values
   if (Err == 0 && isApprox(LInfError, Setup.ExpectedGradErrorLInf, RTol) &&
       isApprox(L2Error, Setup.ExpectedGradErrorL2, RTol)) {
      return 0;
   } else {
      return 1;
   }
}

int testCurl(Real RTol) {
   int Err;
   TestSetup Setup;
   const auto &Mesh = HorzMesh::getDefault();

   // Prepare operator input
   Array1DReal VecEdge("VecEdge", Mesh->NEdgesSize);

   computeVecFieldEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeOrientation::Normal, Mesh);

   // Perform halo exchange
   auto MyHalo   = Halo::getDefault();
   auto VecEdgeH = createHostMirrorCopy(VecEdge);
   MyHalo->exchangeFullArrayHalo(VecEdgeH, OnEdge);
   deepCopy(VecEdge, VecEdgeH);

#ifdef HORZOPERATORS_TEST_PLANE
   const auto XVertex = createDeviceMirrorCopy(Mesh->XVertexH);
   const auto YVertex = createDeviceMirrorCopy(Mesh->YVertexH);
#else
   const auto XVertex = createDeviceMirrorCopy(Mesh->LonVertexH);
   const auto YVertex = createDeviceMirrorCopy(Mesh->LatVertexH);
#endif
   const auto &AreaTriangle = Mesh->AreaTriangle;

   // Compute element-wise errors
   Array1DReal LInfVertex("LInfVertex", Mesh->NVerticesOwned);
   Array1DReal LInfScaleVertex("LInfScaleVertex", Mesh->NVerticesOwned);
   Array1DReal L2Vertex("L2Vertex", Mesh->NVerticesOwned);
   Array1DReal L2ScaleVertex("L2ScaleVertex", Mesh->NVerticesOwned);
   CurlOnVertex CurlVertex(Mesh);
   parallelFor(
       {Mesh->NVerticesOwned}, KOKKOS_LAMBDA(int IVertex) {
          // Numerical result
          const Real CurlNum = CurlVertex(IVertex, VecEdge);

          // Exact result
          const Real X         = XVertex(IVertex);
          const Real Y         = YVertex(IVertex);
          const Real CurlExact = Setup.exactCurlVec(X, Y);

          // Errors
          LInfVertex(IVertex)      = std::abs(CurlNum - CurlExact);
          LInfScaleVertex(IVertex) = std::abs(CurlExact);
          L2Vertex(IVertex) =
              AreaTriangle(IVertex) * LInfVertex(IVertex) * LInfVertex(IVertex);
          L2ScaleVertex(IVertex) = AreaTriangle(IVertex) *
                                   LInfScaleVertex(IVertex) *
                                   LInfScaleVertex(IVertex);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfVertex);
   const Real LInfScaleLoc = maxVal(LInfScaleVertex);
   const Real L2ErrorLoc   = sum(L2Vertex);
   const Real L2ScaleLoc   = sum(L2ScaleVertex);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
   Real LInfError, LInfScale;
   Err =
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err =
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err = MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err = MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);
   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   // Check error values
   if (Err == 0 && isApprox(LInfError, Setup.ExpectedCurlErrorLInf, RTol) &&
       isApprox(L2Error, Setup.ExpectedCurlErrorL2, RTol)) {
      return 0;
   } else {
      return 1;
   }
}

int testRecon(Real RTol) {
   int Err;
   TestSetup Setup;

   const auto &Mesh = HorzMesh::getDefault();

   // Prepare operator input
   Array1DReal VecEdge("VecEdge", Mesh->NEdgesSize);

   computeVecFieldEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeOrientation::Normal, Mesh);

   // Perform halo exchange
   auto MyHalo   = Halo::getDefault();
   auto VecEdgeH = createHostMirrorCopy(VecEdge);
   MyHalo->exchangeFullArrayHalo(VecEdgeH, OnEdge);
   deepCopy(VecEdge, VecEdgeH);

   // Compute exact result
   Array1DReal ExactReconEdge("ExactReconEdge", Mesh->NEdgesOwned);

   computeVecFieldEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       ExactReconEdge, EdgeOrientation::Tangential, Mesh);

   const auto &DcEdge = Mesh->DcEdge;
   const auto &DvEdge = Mesh->DvEdge;

   // Compute element-wise errors
   Array1DReal LInfEdge("LInfEdge", Mesh->NEdgesOwned);
   Array1DReal LInfScaleEdge("LInfScaleEdge", Mesh->NEdgesOwned);
   Array1DReal L2Edge("L2Edge", Mesh->NEdgesOwned);
   Array1DReal L2ScaleEdge("L2ScaleEdge", Mesh->NEdgesOwned);
   TangentialReconOnEdge TanReconEdge(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          // Numerical result
          const Real VecReconNum = TanReconEdge(IEdge, VecEdge);

          // Exact result
          const Real VecReconExact = ExactReconEdge(IEdge);

          // Errors
          LInfEdge(IEdge)      = std::abs(VecReconNum - VecReconExact);
          LInfScaleEdge(IEdge) = std::abs(VecReconExact);
          const Real AreaEdge  = DcEdge(IEdge) * DvEdge(IEdge) / 2;
          L2Edge(IEdge)        = AreaEdge * LInfEdge(IEdge) * LInfEdge(IEdge);
          L2ScaleEdge(IEdge) =
              AreaEdge * LInfScaleEdge(IEdge) * LInfScaleEdge(IEdge);
       });

   // Compute global normalized error norms
   const Real LInfErrorLoc = maxVal(LInfEdge);
   const Real LInfScaleLoc = maxVal(LInfScaleEdge);
   const Real L2ErrorLoc   = sum(L2Edge);
   const Real L2ScaleLoc   = sum(L2ScaleEdge);

   MPI_Comm Comm = MachEnv::getDefaultEnv()->getComm();
   Real LInfError, LInfScale;
   Err =
       MPI_Allreduce(&LInfErrorLoc, &LInfError, 1, MPI_RealKind, MPI_MAX, Comm);
   Err =
       MPI_Allreduce(&LInfScaleLoc, &LInfScale, 1, MPI_RealKind, MPI_MAX, Comm);
   if (LInfScale > 0) {
      LInfError /= LInfScale;
   }

   Real L2Error, L2Scale;
   Err = MPI_Allreduce(&L2ErrorLoc, &L2Error, 1, MPI_RealKind, MPI_SUM, Comm);
   Err = MPI_Allreduce(&L2ScaleLoc, &L2Scale, 1, MPI_RealKind, MPI_SUM, Comm);
   if (L2Scale > 0) {
      L2Error = std::sqrt(L2Error / L2Scale);
   } else {
      L2Error = std::sqrt(L2Error);
   }

   // Check error values
   if (Err == 0 && isApprox(LInfError, Setup.ExpectedReconErrorLInf, RTol) &&
       isApprox(L2Error, Setup.ExpectedReconErrorL2, RTol)) {
      return 0;
   } else {
      return 1;
   }
}

//------------------------------------------------------------------------------
// The initialization routine for Operators testing
int initOperatorsTest(int argc, char *argv[]) {

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);

   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefaultEnv();
   MPI_Comm DefComm = DefEnv->getComm();

   Err = IO::init(DefComm);
   if (Err != 0) {
      LOG_ERROR("OperatorsTest: error initializing parallel IO");
   }

#ifdef HORZOPERATORS_TEST_PLANE
   Err = Decomp::init("OmegaPlanarMesh.nc");
#else
   Err = Decomp::init("OmegaSphereMesh.nc");
#endif
   if (Err != 0) {
      LOG_ERROR("OperatorsTest: error initializing default decomposition");
   }

   Err = Halo::init();
   if (Err != 0) {
      LOG_ERROR("OperatorsTest: error initializing default halo");
   }

   Err = HorzMesh::init();
   if (Err != 0) {
      LOG_ERROR("OperatorsTest: error initializing default mesh");
   }

   return Err;
}

void finalizeOperatorsTest() {
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
   Kokkos::finalize();
   MPI_Finalize();
}

int main(int argc, char *argv[]) {
   int Err = initOperatorsTest(argc, argv);
   if (Err != 0) {
      LOG_CRITICAL("OperatorsTest: Error initializing");
   }

   const Real RTol = sizeof(Real) == 4 ? 1e-2 : 1e-10;

   int DivErr = testDivergence(RTol);
   if (DivErr == 0) {
      LOG_INFO("OperatorsTest: Divergence PASS");
   } else {
      Err = DivErr;
      LOG_INFO("OperatorsTest: Divergence FAIL");
   }

   int GradErr = testGradient(RTol);
   if (GradErr == 0) {
      LOG_INFO("OperatorsTest: Gradient PASS");
   } else {
      Err = GradErr;
      LOG_INFO("OperatorsTest: Gradient FAIL");
   }

   int CurlErr = testCurl(RTol);
   if (CurlErr == 0) {
      LOG_INFO("OperatorsTest: Curl PASS");
   } else {
      Err = CurlErr;
      LOG_INFO("OperatorsTest: Curl FAIL");
   }

   int ReconErr = testRecon(RTol);
   if (Err == 0) {
      LOG_INFO("OperatorsTest: Recon PASS");
   } else {
      Err = ReconErr;
      LOG_INFO("OperatorsTest: Recon FAIL");
   }

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Successful completion");
   }

   finalizeOperatorsTest();
} // end of main
//===-----------------------------------------------------------------------===/
