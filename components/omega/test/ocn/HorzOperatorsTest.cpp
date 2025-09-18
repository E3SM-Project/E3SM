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
#include "VertCoord.h"
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
   Real Ly = std::sqrt(3) / 2;

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

   ErrorMeasures ExpectedDivErrors         = {0.013659577398978353,
                                              0.00367052484586382743};
   ErrorMeasures ExpectedGradErrors        = {0.00187912292540628936,
                                              0.00149841802817334306};
   ErrorMeasures ExpectedCurlErrors        = {0.0271404735181308317,
                                              0.025202316610921989};
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

   VertCoord::init1();

   HorzMesh::init();

   return Err;
}

void finalizeOperatorsTest() {
   HorzMesh::clear();
   VertCoord::clear();
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
