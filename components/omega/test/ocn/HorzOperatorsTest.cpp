#include "HorzOperators.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "mpi.h"

#include <cmath>

using namespace OMEGA;

// analytical expressions for scalar and vector fields used as input for
// operator tests together with exact values of operators for computing errors
// and expected error values
struct TestSetupPlane {

   // lengths of periodic planar mesh
   // TODO: get this from the horizontal mesh once it supports periodic planar
   // meshes
   Real Lx = 1;
   Real Ly = SqrtThree / 2;

   ErrorMeasures ExpectedDivErrors         = {0.00124886886594427027,
                                              0.00124886886590974385};
   ErrorMeasures ExpectedGradErrors        = {0.00125026071878537952,
                                              0.00134354611117262204};
   ErrorMeasures ExpectedCurlErrors        = {0.161365663569699946,
                                              0.161348016897141039};
   ErrorMeasures ExpectedReconErrors       = {0.00450897496974901352,
                                              0.00417367308684470691};
   ErrorMeasures ExpectedAnisoInterpErrors = {0.0026762081503380526,
                                              0.003058198461518835};
   ErrorMeasures ExpectedIsoInterpErrors   = {0.004279097382993937,
                                              0.004200067675522098};

   KOKKOS_FUNCTION Real exactScalar(Real X, Real Y) const {
      return std::sin(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactGradScalarX(Real X, Real Y) const {
      return TwoPi / Lx * std::cos(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactGradScalarY(Real X, Real Y) const {
      return TwoPi / Ly * std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactVecX(Real X, Real Y) const {
      return std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactVecY(Real X, Real Y) const {
      return std::cos(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactDivVec(Real X, Real Y) const {
      return TwoPi * (1. / Lx + 1. / Ly) * std::cos(TwoPi * X / Lx) *
             std::cos(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real exactCurlVec(Real X, Real Y) const {
      return TwoPi * (-1. / Lx + 1. / Ly) * std::sin(TwoPi * X / Lx) *
             std::sin(TwoPi * Y / Ly);
   }
};

struct TestSetupSphere1 {
   // radius of spherical mesh
   // TODO: get this from the mesh
   Real Radius = REarth;

   ErrorMeasures ExpectedDivErrors         = {0.013652414501664885,
                                              0.0036904315983599676};
   ErrorMeasures ExpectedGradErrors        = {0.0019094381714837498,
                                              0.0015218320661105687};
   ErrorMeasures ExpectedCurlErrors        = {0.02713957370128636,
                                              0.025202095212756463};
   ErrorMeasures ExpectedReconErrors       = {0.0206375134079833517,
                                              0.00692590524910695858};
   ErrorMeasures ExpectedAnisoInterpErrors = {0.0024015775047603197,
                                              0.0018490649516209202};
   ErrorMeasures ExpectedIsoInterpErrors   = {0.007438367234983312,
                                              0.0029921955942401697};

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

struct TestSetupSphere2 {
   // radius of spherical mesh
   // TODO: get this from the mesh
   Real Radius = REarth;

   ErrorMeasures ExpectedDivErrors         = {1.37734693033362766e-10,
                                              0.000484370621558727582};
   ErrorMeasures ExpectedGradErrors        = {0.000906351303388669991,
                                              0.000949206041390823676};
   ErrorMeasures ExpectedCurlErrors        = {0.00433205620592059647,
                                              0.00204725417666192042};
   ErrorMeasures ExpectedReconErrors       = {0.0254271921029878764,
                                              0.00419630561428921064};
   ErrorMeasures ExpectedAnisoInterpErrors = {0.0014465229922953644,
                                              0.001643777653612931};
   ErrorMeasures ExpectedIsoInterpErrors   = {0.004755875091568591,
                                              0.0025556382734782538};

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

#ifdef HORZOPERATORS_TEST_PLANE
constexpr Geometry Geom          = Geometry::Planar;
constexpr char DefaultMeshFile[] = "OmegaPlanarMesh.nc";
#else
constexpr Geometry Geom          = Geometry::Spherical;
constexpr char DefaultMeshFile[] = "OmegaSphereMesh.nc";
#endif

#if defined HORZOPERATORS_TEST_PLANE
using TestSetup = TestSetupPlane;
#elif defined HORZOPERATORS_TEST_SPHERE_1
using TestSetup = TestSetupSphere1;
#elif defined HORZOPERATORS_TEST_SPHERE_2
using TestSetup = TestSetupSphere2;
#endif

int testDivergence(Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto &Mesh      = HorzMesh::getDefault();
   const int NVertLayers = 16;

   // Prepare operator input
   Array2DReal VecEdge("VecEdge", Mesh->NEdgesSize, NVertLayers);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeComponent::Normal, Geom, Mesh);

   // Compute exact result
   Array2DReal ExactDivCell("ExactDivCell", Mesh->NCellsOwned, NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.exactDivVec(X, Y); },
       ExactDivCell, Geom, Mesh, OnCell, ExchangeHalos::No);

   // Compute numerical result
   Array2DReal NumDivCell("NumDivCell", Mesh->NCellsOwned, NVertLayers);
   DivergenceOnCell DivergenceCell(Mesh);
   parallelFor(
       {Mesh->NCellsOwned, NVertLayers}, KOKKOS_LAMBDA(int ICell, int K) {
          DivergenceCell(NumDivCell, ICell, K, VecEdge);
       });

   // Compute error measures
   ErrorMeasures DivErrors;
   Err += computeErrors(DivErrors, NumDivCell, ExactDivCell, Mesh, OnCell);
   // Check error values
   Err += checkErrors("OperatorsTest", "Divergence", DivErrors,
                      Setup.ExpectedDivErrors, RTol);

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Divergence PASS");
   }

   return Err;
}

int testGradient(Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto &Mesh      = HorzMesh::getDefault();
   const int NVertLayers = 16;

   // Prepare operator input
   Array2DReal ScalarCell("ScalarCell", Mesh->NCellsSize, NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real Coord1, Real Coord2) {
          return Setup.exactScalar(Coord1, Coord2);
       },
       ScalarCell, Geom, Mesh, OnCell);

   // Compute exact result
   Array2DReal ExactGradEdge("ExactGradEdge", Mesh->NEdgesOwned, NVertLayers);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactGradScalarX(X, Y);
          VecField[1] = Setup.exactGradScalarY(X, Y);
       },
       ExactGradEdge, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::No);

   // Compute numerical result
   GradientOnEdge GradientEdge(Mesh);
   Array2DReal NumGradEdge("NumGradEdge", Mesh->NEdgesOwned, NVertLayers);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int K) {
          GradientEdge(NumGradEdge, IEdge, K, ScalarCell);
       });

   // Compute error measures
   ErrorMeasures GradErrors;
   Err += computeErrors(GradErrors, NumGradEdge, ExactGradEdge, Mesh, OnEdge);
   // Check error values
   Err += checkErrors("OperatorsTest", "Gradient", GradErrors,
                      Setup.ExpectedGradErrors, RTol);

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Gradient PASS");
   }

   return Err;
}

int testCurl(Real RTol) {
   int Err = 0;
   TestSetup Setup;
   const auto &Mesh      = HorzMesh::getDefault();
   const int NVertLayers = 16;

   // Prepare operator input
   Array2DReal VecEdge("VecEdge", Mesh->NEdgesSize, NVertLayers);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeComponent::Normal, Geom, Mesh);

   // Compute exact result
   Array2DReal ExactCurlVertex("ExactCurlVertex", Mesh->NVerticesOwned,
                               NVertLayers);
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.exactCurlVec(X, Y); },
       ExactCurlVertex, Geom, Mesh, OnVertex, ExchangeHalos::No);

   // Compute numerical result
   Array2DReal NumCurlVertex("NumCurlVertex", Mesh->NVerticesOwned,
                             NVertLayers);
   CurlOnVertex CurlVertex(Mesh);
   parallelFor(
       {Mesh->NVerticesOwned, NVertLayers}, KOKKOS_LAMBDA(int IVertex, int K) {
          CurlVertex(NumCurlVertex, IVertex, K, VecEdge);
       });

   // Compute error measures
   ErrorMeasures CurlErrors;
   Err += computeErrors(CurlErrors, NumCurlVertex, ExactCurlVertex, Mesh,
                        OnVertex);
   // Check error values
   Err += checkErrors("OperatorsTest", "Curl", CurlErrors,
                      Setup.ExpectedCurlErrors, RTol);

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Curl PASS");
   }

   return Err;
}

int testRecon(Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto &Mesh      = HorzMesh::getDefault();
   const int NVertLayers = 16;

   // Prepare operator input
   Array2DReal VecEdge("VecEdge", Mesh->NEdgesSize, NVertLayers);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeComponent::Normal, Geom, Mesh);

   // Compute exact result
   Array2DReal ExactReconEdge("ExactReconEdge", Mesh->NEdgesOwned, NVertLayers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       ExactReconEdge, EdgeComponent::Tangential, Geom, Mesh,
       ExchangeHalos::No);

   // Compute numerical result
   Array2DReal NumReconEdge("NumReconEdge", Mesh->NEdgesOwned, NVertLayers);
   TangentialReconOnEdge TanReconEdge(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLayers}, KOKKOS_LAMBDA(int IEdge, int K) {
          TanReconEdge(NumReconEdge, IEdge, K, VecEdge);
       });

   // Compute error measures
   ErrorMeasures ReconErrors;
   Err +=
       computeErrors(ReconErrors, NumReconEdge, ExactReconEdge, Mesh, OnEdge);
   // Check error values
   Err += checkErrors("OperatorsTest", "Recon", ReconErrors,
                      Setup.ExpectedReconErrors, RTol);

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Recon PASS");
   }

   return Err;
}

int testInterpCellToEdge(Real RTol) {
   int Err = 0;
   TestSetup Setup;
   const auto &Mesh = HorzMesh::getDefault();

   // Prepare operator input
   Array1DReal ScalarCell("ScalarCell", Mesh->NCellsSize);
   Err += setScalar(
       KOKKOS_LAMBDA(Real Coord1, Real Coord2) {
          return Setup.exactScalar(Coord1, Coord2);
       },
       ScalarCell, Geom, Mesh, OnCell);

   // Compute exact result
   Array1DReal ExactScalarEdge("ExactScalarEdge", Mesh->NEdgesOwned);
   Err += setScalar(
       KOKKOS_LAMBDA(Real Coord1, Real Coord2) {
          return Setup.exactScalar(Coord1, Coord2);
       },
       ExactScalarEdge, Geom, Mesh, OnEdge, ExchangeHalos::No);

   // Compute numerical result
   Array1DReal IsoNumScalarEdge("IsoNumScalarEdge", Mesh->NEdgesOwned);
   Array1DReal AnisoNumScalarEdge("AnisoNumScalarEdge", Mesh->NEdgesOwned);
   InterpCellToEdge Interp(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned}, KOKKOS_LAMBDA(int IEdge) {
          AnisoNumScalarEdge(IEdge) =
              Interp(IEdge, ScalarCell, InterpCellToEdgeOption::Anisotropic);
          IsoNumScalarEdge(IEdge) =
              Interp(IEdge, ScalarCell, InterpCellToEdgeOption::Isotropic);
       });

   // Compute error measures
   ErrorMeasures AnisoInterpErrors;
   Err += computeErrors(AnisoInterpErrors, AnisoNumScalarEdge, ExactScalarEdge,
                        Mesh, OnEdge);

   ErrorMeasures IsoInterpErrors;
   Err += computeErrors(IsoInterpErrors, IsoNumScalarEdge, ExactScalarEdge,
                        Mesh, OnEdge);

   // Check error values
   Err += checkErrors("OperatorsTest", "AnisoInterpCellToEdge",
                      AnisoInterpErrors, Setup.ExpectedAnisoInterpErrors, RTol);

   Err += checkErrors("OperatorsTest", "IsoInterpCellToEdge", IsoInterpErrors,
                      Setup.ExpectedIsoInterpErrors, RTol);

   if (Err == 0) {
      LOG_INFO("OperatorsTest: InterpCellToEdge PASS");
   }

   return Err;
}

int testsecondderivativeoncellDeterminePlanerPatchGeometry(Real RTol) {
   int Err            = 0;
   const int C[38][2] = {
       {2200000, 2000000}, {2800000, 2000000}, {1300000, 2519615},
       {1900000, 2519615}, {3100000, 2519615}, {3700000, 2519615},
       {1000000, 3039230}, {2200000, 3039230}, {2800000, 3039230},
       {4000000, 3039230}, {1300000, 3558846}, {1900000, 3558846},
       {3100000, 3558846}, {3700000, 3558846}, {1000000, 4078461},
       {2200000, 4078461}, {2800000, 4078461}, {4000000, 4078461},
       {1300000, 4598076}, {1900000, 4598076}, {3100000, 4598076},
       {3700000, 4598076}, {1000000, 5117691}, {2200000, 5117691},
       {2800000, 5117691}, {4000000, 5117691}, {1300000, 5637307},
       {1900000, 5637307}, {3100000, 5637307}, {3700000, 5637307},
       {1000000, 6156922}, {2200000, 6156922}, {2800000, 6156922},
       {4000000, 6156922}, {1300000, 6676537}, {1900000, 6676537},
       {3100000, 6676537}, {3700000, 6676537}};
   const int E[12][6] = {{0, 1, 4, 8, 7, 3},       {2, 3, 7, 11, 10, 6},
                         {4, 5, 9, 13, 12, 8},     {7, 8, 12, 16, 15, 11},
                         {10, 11, 15, 19, 18, 14}, {12, 13, 17, 21, 20, 16},
                         {15, 16, 20, 24, 23, 19}, {18, 19, 23, 27, 26, 22},
                         {20, 21, 25, 29, 28, 24}, {23, 24, 28, 32, 31, 27},
                         {26, 27, 31, 35, 34, 30}, {28, 29, 33, 37, 36, 32}};
   const int CtoC     = 1039230;
   const R8 Pii       = 3.141592653589793_Real;

   std::map<std::pair<int, int>, int> Edges;
   for (int c = 0; c < 12; ++c) {
      for (int i = 0; i < 6; ++i) {
         const int j = (i + 1) % 6;
         const std::pair<int, int> edge(std::min(E[c][i], E[c][j]),
                                        std::max(E[c][i], E[c][j]));
         if (!Edges.count(edge))
            Edges[edge] = Edges.size();
      }
   }
   const int NEdges = Edges.size();
   std::vector<std::array<int, 6>> edgesOnCell(12);
   std::vector<std::vector<int>> cellsOnEdge(NEdges);
   std::vector<double> angleEdge(NEdges);
   std::vector<double> dcEdge(NEdges, -1);
   for (int c = 0; c < 12; ++c) {
      for (int i = 0; i < 6; ++i) {
         const int j = (i + 1) % 6;
         const std::pair<int, int> edge(std::min(E[c][i], E[c][j]),
                                        std::max(E[c][i], E[c][j]));
         const int e       = Edges.at(edge);
         edgesOnCell[c][i] = e;
         if (cellsOnEdge[e].empty()) {
            const int v0[2]    = {C[E[c][i]][0], C[E[c][i]][1]};
            const int v1[2]    = {C[E[c][j]][0], C[E[c][j]][1]};
            const int v[2]     = {v1[0] - v0[0], v1[1] - v0[1]};
            const double theta = std::atan2(-v[0], v[1]);
            angleEdge[e]       = theta;
         }
         cellsOnEdge[e].push_back(c);
      }
   }
   for (int e = 0; e < cellsOnEdge.size(); ++e) {
      if (1 < cellsOnEdge[e].size())
         dcEdge[e] = CtoC;
   }
   Array2DI4 EdgesOnCell("EdgesOnCell", 12, 6);
   auto EdgesOnCellH = createHostMirrorCopy(EdgesOnCell);
   for (int c = 0; c < 12; ++c)
      for (int i = 0; i < 6; ++i)
         EdgesOnCellH(c, i) = edgesOnCell[c][i];
   OMEGA::deepCopy(EdgesOnCell, EdgesOnCellH);

   Array2DI4 CellsOnEdge("CellsOnEdge", cellsOnEdge.size(), 2);
   auto CellsOnEdgeH = createHostMirrorCopy(CellsOnEdge);
   for (int e = 0; e < cellsOnEdge.size(); ++e)
      for (int i = 0; i < cellsOnEdge[e].size(); ++i)
         CellsOnEdgeH(e, i) = cellsOnEdge[e][i];
   OMEGA::deepCopy(CellsOnEdge, CellsOnEdgeH);

   Array1DReal AngleEdge("AndleEdge", angleEdge.size());
   auto AngleEdgeH = createHostMirrorCopy(AngleEdge);
   for (int e = 0; e < angleEdge.size(); ++e)
      AngleEdgeH(e) = angleEdge[e];
   OMEGA::deepCopy(AngleEdge, AngleEdgeH);

   Array1DReal DcEdge("DcEdge", dcEdge.size());
   auto DcEdgeH = createHostMirrorCopy(DcEdge);
   for (int e = 0; e < dcEdge.size(); ++e)
      DcEdgeH(e) = dcEdge[e];
   OMEGA::deepCopy(DcEdge, DcEdgeH);

   class SecondDerivativeOnCellTest : public SecondDerivativeOnCell {
    public:
      virtual ~SecondDerivativeOnCellTest() {}
      KOKKOS_INLINE_FUNCTION static void
      Test(const int ICell, const int NEdges, const Array1DI4 EdgesOnCell,
           const Array2DI4 CellsOnEdge, const Array1DReal AngleEdge,
           const Array1DReal DcEdge, Array1DReal XP, Array1DReal YP,
           Array1DReal Angle2D) {
         SecondDerivativeOnCell::DeterminePlanerPatchGeometry(
             ICell, NEdges, EdgesOnCell, CellsOnEdge, AngleEdge, DcEdge, XP, YP,
             Angle2D);
      }
   };
   for (int ICell = 3; ICell < 7; ICell += 3) {
      const I4 NEdges   = 6;
      const I4 MaxEdges = 7;
      Array1DReal XP("XP", MaxEdges);
      Array1DReal YP("YP", MaxEdges);
      Array1DReal Angle2D("Array2D", MaxEdges);
      const Array1DI4 edgesOnCell =
          Kokkos::subview(EdgesOnCell, ICell, Kokkos::ALL);
      Kokkos::parallel_for(
          1, KOKKOS_LAMBDA(int /* dummy */) {
             SecondDerivativeOnCellTest::Test(ICell, NEdges, edgesOnCell,
                                              CellsOnEdge, AngleEdge, DcEdge,
                                              XP, YP, Angle2D);
          });
      auto XPH      = createHostMirrorCopy(XP);
      auto YPH      = createHostMirrorCopy(YP);
      auto Angle2DH = createHostMirrorCopy(Angle2D);
      OMEGA::deepCopy(XPH, XP);
      OMEGA::deepCopy(YPH, YP);
      OMEGA::deepCopy(Angle2DH, Angle2D);
      for (int i = 0; i < 6; ++i) {
         const Real RTol    = sizeof(Real) == 4 ? 1e-3 : 1e-5;
         const double theta = -Pii / 2 + (5 == i ? -Pii / 3 : i * Pii / 3);
         const double x     = CtoC * std::cos(theta);
         const double y     = CtoC * std::sin(theta);
         if (!isApprox(XPH[i], x, RTol)) {
            Err++;
            LOG_ERROR(
                "{}: FAIL, expected {}, got {}",
                "testsecondderivativeoncellDeterminePlanerPatchGeometry:x", x,
                XP[i]);
         }
         if (!isApprox(YPH[i], y, RTol)) {
            Err++;
            LOG_ERROR(
                "{}: FAIL, expected {}, got {}",
                "testsecondderivativeoncellDeterminePlanerPatchGeometry:y", y,
                YP[i]);
         }
         if (!isApprox(Angle2DH[i], theta, RTol)) {
            Err++;
            LOG_ERROR(
                "{}: FAIL, expected {}, got {}",
                "testsecondderivativeoncellDeterminePlanerPatchGeometry:theta",
                theta, Angle2D[i]);
         }
      }
   }
   if (Err == 0)
      LOG_INFO("testsecondderivativeoncellDeterminePlanerPatchGeometry: "
               "Successful completion");
   return Err;
}

int testsecondderivativeoncellLeastSquaresFit(Real RTol) {
   int Err = 0;

   const int CtoC       = 7;
   const R8 Pii         = 3.141592653589793_Real;
   const I4 MaxMaxEdges = 10;
   const int NEdges     = 6;
   class SecondDerivativeOnCellTest : public SecondDerivativeOnCell {
    public:
      virtual ~SecondDerivativeOnCellTest() {}
      KOKKOS_INLINE_FUNCTION static void Test(const Array1DReal XP,
                                              const Array1DReal YP,
                                              const int NEdges, Array2DReal B) {
         SecondDerivativeOnCell::LeastSquaresFit(XP, YP, NEdges, B);
      }
   };
   Array1DReal XP("XP", NEdges);
   Array1DReal YP("YP", NEdges);
   auto XPH = createHostMirrorCopy(XP);
   auto YPH = createHostMirrorCopy(YP);
   for (int i = 0; i < NEdges; ++i) {
      const double theta = -Pii / 2 + i * Pii / 3;
      XPH[i]             = CtoC * std::cos(theta);
      YPH[i]             = CtoC * std::sin(theta);
   }
   OMEGA::deepCopy(XP, XPH);
   OMEGA::deepCopy(YP, YPH);
   Array2DReal B("B", MaxMaxEdges, MaxMaxEdges);
   SecondDerivativeOnCellTest::Test(XP, YP, NEdges, B);

   // From a Mathematica script.
   const Real c0      = 0.0412393049421;
   const Real c1      = 0.0238095238095;
   const Real c2      = 0.00680272108844;
   const Real c3      = 0.0117826585549;
   const Real M[6][7] = {
       {1., 0, 0, 0, 0, 0, 0},
       {0, 0, c0, c0, 0, -c0, -c0},
       {0, -0.047619047619, -c1, c1, 0.047619047619, c1, -c1},
       {-0.0204081632653, -0.00340136054422, c2, c2, -0.00340136054422, c2, c2},
       {0, 0, -c3, c3, 0, -c3, c3},
       {-0.0204081632653, 0.0102040816327, 0, 0, 0.0102040816327, 0, 0}};

   auto BH = createHostMirrorCopy(B);
   OMEGA::deepCopy(BH, B);
   for (int i = 0; i < MaxMaxEdges; ++i) {
      for (int j = 0; j < MaxMaxEdges; ++j) {
         const Real RTol = sizeof(Real) == 4 ? 1e-3 : 1e-5;
         const Real m    = (i < 6 && j < 7) ? M[i][j] : 0;
         if (i < 6 && j < 7 && !m) {
            if (1e-12 < std::abs(BH(i, j))) {
               Err++;
               LOG_ERROR("{}: FAIL, expected {}, got {}",
                         "testsecondderivativeoncellLeastSquaresFit", m,
                         B(i, j));
            }
         } else {
            if (!isApprox(B(i, j), m, RTol)) {
               Err++;
               LOG_ERROR("{}: FAIL, expected {}, got {}",
                         "testsecondderivativeoncellLeastSquaresFit", m,
                         B(i, j));
            }
         }
      }
   }
   if (Err == 0)
      LOG_INFO(
          "testsecondderivativeoncellLeastSquaresFit: Successful completion");
   return Err;
}

int testsecondderivativeoncellDetermineSphericalPatchGeometry(Real RTol) {
   int Err            = 0;
   const int C[38][2] = {
       {2200000, 2000000}, {2800000, 2000000}, {1300000, 2519615},
       {1900000, 2519615}, {3100000, 2519615}, {3700000, 2519615},
       {1000000, 3039230}, {2200000, 3039230}, {2800000, 3039230},
       {4000000, 3039230}, {1300000, 3558846}, {1900000, 3558846},
       {3100000, 3558846}, {3700000, 3558846}, {1000000, 4078461},
       {2200000, 4078461}, {2800000, 4078461}, {4000000, 4078461},
       {1300000, 4598076}, {1900000, 4598076}, {3100000, 4598076},
       {3700000, 4598076}, {1000000, 5117691}, {2200000, 5117691},
       {2800000, 5117691}, {4000000, 5117691}, {1300000, 5637307},
       {1900000, 5637307}, {3100000, 5637307}, {3700000, 5637307},
       {1000000, 6156922}, {2200000, 6156922}, {2800000, 6156922},
       {4000000, 6156922}, {1300000, 6676537}, {1900000, 6676537},
       {3100000, 6676537}, {3700000, 6676537}};
   const int E[12][6] = {{0, 1, 4, 8, 7, 3},       {2, 3, 7, 11, 10, 6},
                         {4, 5, 9, 13, 12, 8},     {7, 8, 12, 16, 15, 11},
                         {10, 11, 15, 19, 18, 14}, {12, 13, 17, 21, 20, 16},
                         {15, 16, 20, 24, 23, 19}, {18, 19, 23, 27, 26, 22},
                         {20, 21, 25, 29, 28, 24}, {23, 24, 28, 32, 31, 27},
                         {26, 27, 31, 35, 34, 30}, {28, 29, 33, 37, 36, 32}};
   const int CtoC     = 1039230;
   const int R        = 100 * CtoC;
   const R8 Pii       = 3.141592653589793_Real;

   // Project coordinates to sphere.

   std::map<std::pair<int, int>, int> Edges;
   for (int c = 0; c < 12; ++c) {
      for (int i = 0; i < 6; ++i) {
         const int j = (i + 1) % 6;
         const std::pair<int, int> edge(std::min(E[c][i], E[c][j]),
                                        std::max(E[c][i], E[c][j]));
         if (!Edges.count(edge))
            Edges[edge] = Edges.size();
      }
   }
   const int NEdges = Edges.size();
   std::vector<std::array<int, 6>> edgesOnCell(12);
   std::vector<std::vector<int>> cellsOnEdge(NEdges);
   std::vector<double> angleEdge(NEdges);
   std::vector<std::array<int, 2>> verticesOnEdge(NEdges);
   std::vector<double> dcEdge(NEdges, -1);
   for (int c = 0; c < 12; ++c) {
      for (int i = 0; i < 6; ++i) {
         const int j = (i + 1) % 6;
         const std::pair<int, int> edge(std::min(E[c][i], E[c][j]),
                                        std::max(E[c][i], E[c][j]));
         const int e       = Edges.at(edge);
         edgesOnCell[c][i] = e;
         if (cellsOnEdge[e].empty()) {
            const int v0[2]      = {C[E[c][i]][0], C[E[c][i]][1]};
            const int v1[2]      = {C[E[c][j]][0], C[E[c][j]][1]};
            const int v[2]       = {v1[0] - v0[0], v1[1] - v0[1]};
            const double theta   = std::atan2(-v[0], v[1]);
            angleEdge[e]         = theta;
            verticesOnEdge[e][0] = E[c][i];
            verticesOnEdge[e][1] = E[c][j];
         }
         cellsOnEdge[e].push_back(c);
      }
   }
   for (int e = 0; e < cellsOnEdge.size(); ++e) {
      if (1 < cellsOnEdge[e].size())
         dcEdge[e] = CtoC;
   }

   double T[2] = {};
   for (int I = 0; I < 6; ++I)
      for (int J = 0; J < 2; ++J)
         T[J] += C[E[3][I]][J];

   T[0] /= 6;
   T[1] /= 6;
   double X[38][2] = {};
   for (int I = 0; I < 38; ++I)
      for (int J = 0; J < 2; ++J)
         X[I][J] = C[I][J] - T[J];

   Array1DReal XCell("XCell", 12);
   Array1DReal YCell("YCell", 12);
   Array1DReal ZCell("ZCell", 12);
   auto XCellH = createHostMirrorCopy(XCell);
   auto YCellH = createHostMirrorCopy(YCell);
   auto ZCellH = createHostMirrorCopy(ZCell);
   std::vector<double> xcellh(12, 0);
   std::vector<double> ycellh(12, 0);
   std::vector<double> zcellh(12, 0);

   for (int I = 0; I < 12; ++I) {
      T[0] = T[1] = 0;
      for (int J = 0; J < 6; ++J) {
         for (int K = 0; K < 2; ++K)
            T[K] += X[E[I][J]][K];
      }
      T[0] /= 6;
      T[1] /= 6;
      const double d = std::sqrt(T[0] * T[0] + T[1] * T[1]);
      xcellh[I]      = R * std::sin(d / R) * T[0] / d;
      ycellh[I]      = R * std::sin(d / R) * T[1] / d;
      zcellh[I]      = R * std::cos(d / R);
   }
   for (int I = 0; I < 12; ++I) {
      XCellH(I) = xcellh[I];
      YCellH(I) = ycellh[I];
      ZCellH(I) = zcellh[I];
   }
   OMEGA::deepCopy(XCell, XCellH);
   OMEGA::deepCopy(YCell, YCellH);
   OMEGA::deepCopy(ZCell, ZCellH);

   Array1DReal XVertex("XVertex", 38);
   Array1DReal YVertex("YVertex", 38);
   Array1DReal ZVertex("ZVertex", 38);
   auto XVertexH = createHostMirrorCopy(XVertex);
   auto YVertexH = createHostMirrorCopy(YVertex);
   auto ZVertexH = createHostMirrorCopy(ZVertex);

   for (int I = 0; I < 38; ++I) {
      const double x = X[I][0], y = X[I][1];
      const double d = std::sqrt(x * x + y * y);
      XVertexH(I)    = R * std::sin(d / R) * X[I][0] / d;
      YVertexH(I)    = R * std::sin(d / R) * X[I][1] / d;
      ZVertexH(I)    = R * std::cos(d / R);
   }
   OMEGA::deepCopy(XVertex, XVertexH);
   OMEGA::deepCopy(YVertex, YVertexH);
   OMEGA::deepCopy(ZVertex, ZVertexH);

   Array2DI4 EdgesOnCell("EdgesOnCell", 12, 6);
   auto EdgesOnCellH = createHostMirrorCopy(EdgesOnCell);
   for (int c = 0; c < 12; ++c)
      for (int i = 0; i < 6; ++i)
         EdgesOnCellH(c, i) = edgesOnCell[c][i];
   OMEGA::deepCopy(EdgesOnCell, EdgesOnCellH);

   Array2DI4 VerticesOnEdge("VerticesOnEdge", NEdges, 2);
   auto VerticesOnEdgeH = createHostMirrorCopy(VerticesOnEdge);
   for (int c = 0; c < NEdges; ++c)
      for (int i = 0; i < 2; ++i)
         VerticesOnEdgeH(c, i) = verticesOnEdge[c][i];
   OMEGA::deepCopy(VerticesOnEdge, VerticesOnEdgeH);

   Array2DI4 CellsOnEdge("CellsOnEdge", cellsOnEdge.size(), 2);
   auto CellsOnEdgeH = createHostMirrorCopy(CellsOnEdge);
   for (int e = 0; e < cellsOnEdge.size(); ++e)
      for (int i = 0; i < cellsOnEdge[e].size(); ++i)
         CellsOnEdgeH(e, i) = cellsOnEdge[e][i];
   OMEGA::deepCopy(CellsOnEdge, CellsOnEdgeH);

   Array1DReal AngleEdge("AndleEdge", angleEdge.size());
   auto AngleEdgeH = createHostMirrorCopy(AngleEdge);
   for (int e = 0; e < angleEdge.size(); ++e)
      AngleEdgeH(e) = angleEdge[e];
   OMEGA::deepCopy(AngleEdge, AngleEdgeH);

   Array1DReal DcEdge("DcEdge", dcEdge.size());
   auto DcEdgeH = createHostMirrorCopy(DcEdge);
   for (int e = 0; e < dcEdge.size(); ++e)
      DcEdgeH(e) = dcEdge[e];
   OMEGA::deepCopy(DcEdge, DcEdgeH);

   class SecondDerivativeOnCellTest : public SecondDerivativeOnCell {
    public:
      virtual ~SecondDerivativeOnCellTest() {}
      KOKKOS_INLINE_FUNCTION static void
      Test(const int NEdges, const Array1DI4 EdgesOnCell,
           const Array2DI4 VerticiesOnEdge, const Array1DReal XCell,
           const Array1DReal YCell, const Array1DReal ZCell,
           const Array1DReal XVertex, const Array1DReal YVertex,
           const Array1DReal ZVertex, const Array1DI4 CellList, Array1DReal XP,
           Array1DReal YP, Array1DReal Angle2D, Real &ThetaAbs) {
         SecondDerivativeOnCell::DetermineSphericalPatchGeometry(
             NEdges, EdgesOnCell, VerticiesOnEdge, XCell, YCell, ZCell, XVertex,
             YVertex, ZVertex, CellList, XP, YP, Angle2D, ThetaAbs);
      }
   };
   {
      const I4 NEdges = 6;
      Array1DReal XP("XP", NEdges);
      Array1DReal YP("YP", NEdges);
      Array1DReal Angle2D("Angle2D", NEdges);
      ;
      Array1DI4 CellList("CellList", NEdges);
      {
         const I4 List[NEdges + 1] = {3, 0, 2, 5, 6, 4, 1};
         auto CellListH            = createHostMirrorCopy(CellList);
         for (int I = 0; I < NEdges + 1; ++I)
            CellListH[I] = List[I];
         OMEGA::deepCopy(CellList, CellListH);
      }
      Array1DI4 edgesOnCell = Kokkos::subview(EdgesOnCell, 3, Kokkos::ALL);
      Kokkos::parallel_for(
          1, KOKKOS_LAMBDA(int /* dummy */) {
             Real ThetaAbs = {};
             SecondDerivativeOnCellTest::Test(
                 NEdges, edgesOnCell, VerticesOnEdge, XCell, YCell, ZCell,
                 XVertex, YVertex, ZVertex, CellList, XP, YP, Angle2D,
                 ThetaAbs);
          });
      auto XPH      = createHostMirrorCopy(XP);
      auto YPH      = createHostMirrorCopy(YP);
      auto Angle2DH = createHostMirrorCopy(Angle2D);
      for (int i = 0; i < 6; ++i) {
         const Real RTol = sizeof(Real) == 4 ? 1e-3 : 1e-5;
         T[0] = T[1] = 0;
         for (int J = 0; J < 6; ++J)
            for (int K = 0; K < 2; ++K)
               T[K] += X[E[CellList[i + 1]][J]][K];
         T[0] /= 6;
         T[1] /= 6;

         const double phi = Pii / 2 + i * 2 * Pii / 6;
         // Rotate the mesh by Pi.  This is what is done in
         // DetermineSphericalPatchGeometry
         const double x = -T[0]; // same as CtoC*std::cos(phi);
         const double y = -T[1]; // same as CtoC*std::sin(phi);

         if (!isApprox(1 + XPH[i], 1 + x, RTol)) {
            Err++;
            LOG_ERROR(
                "{}: FAIL, expected {}, got {}",
                "testsecondderivativeoncellDetermineSphericalPatchGeometry:x",
                x, XPH[i]);
         }
         if (!isApprox(1 + YPH[i], 1 + y, RTol)) {
            Err++;
            LOG_ERROR(
                "{}: FAIL, expected {}, got {}",
                "testsecondderivativeoncellDetermineSphericalPatchGeometry:y",
                y, YPH[i]);
         }
         if (!isApprox(Angle2DH[i], phi, RTol)) {
            Err++;
            LOG_ERROR(
                "{}: FAIL, expected {}, got {}",
                "testsecondderivativeoncellDetermineSphericalPatchGeometry:phi",
                phi, Angle2DH[i]);
         }
      }
   }
   if (Err == 0)
      LOG_INFO("testsecondderivativeoncellDetermineSphericalPatchGeometry: "
               "Successful completion");
   return Err;
}

int testsecondderivativeoncellconstructor(Real RTol) {
   int Err          = 0;
   const auto &Mesh = HorzMesh::getDefault();
   // Compute numerical result
   const auto MaxEdges2   = Mesh->MaxEdges2;
   const auto NEdgesAll   = Mesh->NEdgesAll;
   const auto NCellsOwned = Mesh->NCellsOwned;

   Array3DReal DerivTwo("DerivTwo", MaxEdges2, 2, NEdgesAll);
   SecondDerivativeOnCell DerivativeOnCell(Mesh);
   parallelFor(
       {NCellsOwned},
       KOKKOS_LAMBDA(int ICell) { DerivativeOnCell(DerivTwo, ICell); });
   if (Err == 0)
      LOG_INFO("OperatorsTest: TestSecondDerivativeOnCell PASS");
   return Err;
}
//------------------------------------------------------------------------------
// The initialization routine for Operators testing
int initOperatorsTest(const std::string &MeshFile) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   IO::init(DefComm);
   Decomp::init(MeshFile);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("OperatorsTest: error initializing default halo");
   }

   HorzMesh::init();

   return Err;
}

void finalizeOperatorsTest() {
   HorzMesh::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int operatorsTest(const std::string &MeshFile = DefaultMeshFile) {
   int Err = initOperatorsTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("OperatorsTest: Error initializing");
   }

   const Real RTol = sizeof(Real) == 4 ? 1e-2 : 1e-10;

   Err += testDivergence(RTol);
   Err += testGradient(RTol);
   Err += testCurl(RTol);
   Err += testRecon(RTol);
   Err += testInterpCellToEdge(RTol);
   Err += testsecondderivativeoncellconstructor(RTol);
   Err += testsecondderivativeoncellDeterminePlanerPatchGeometry(RTol);
   Err += testsecondderivativeoncellLeastSquaresFit(RTol);
   Err += testsecondderivativeoncellDetermineSphericalPatchGeometry(RTol);

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Successful completion");
   }
   finalizeOperatorsTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += operatorsTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
