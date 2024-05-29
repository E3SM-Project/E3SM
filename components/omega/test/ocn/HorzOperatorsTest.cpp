#include "HorzOperators.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "mpi.h"

#include <cmath>

using namespace OMEGA;

static_assert(VecLength == 1, "HorzOperatorsTest needs vector length to be 1");

// analytical expressions for scalar and vector fields used as input for
// operator tests together with exact values of operators for computing errors
// and expected error values
struct TestSetupPlane {
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

struct TestSetupSphere1 {
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

struct TestSetupSphere2 {
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
   const int NVertLevels = 1;

   // Prepare operator input
   Array2DReal VecEdge("VecEdge", Mesh->NEdgesSize, NVertLevels);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeComponent::Normal, Geom, Mesh, NVertLevels);

   // Compute exact result
   Array2DReal ExactDivCell("ExactDivCell", Mesh->NCellsOwned, NVertLevels);
   Err += setScalarCell(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.exactDivVec(X, Y); },
       ExactDivCell, Geom, Mesh, NVertLevels, false);

   // Compute numerical result
   Array2DReal NumDivCell("NumDivCell", Mesh->NCellsOwned, NVertLevels);
   DivergenceOnCell DivergenceCell(Mesh);
   parallelFor(
       {Mesh->NCellsOwned, NVertLevels}, KOKKOS_LAMBDA(int ICell, int K) {
          DivergenceCell(NumDivCell, ICell, K, VecEdge);
       });

   // Compute error measures
   ErrorMeasures DivErrors;
   Err += computeErrorsCell(DivErrors, NumDivCell, ExactDivCell, Mesh,
                            NVertLevels);

   // Check error values
   if (!isApprox(DivErrors.LInf, Setup.ExpectedDivErrorLInf, RTol)) {
      Err++;
      LOG_ERROR("OperatorsTest: Divergence LInf FAIL");
   }

   if (!isApprox(DivErrors.L2, Setup.ExpectedDivErrorL2, RTol)) {
      Err++;
      LOG_ERROR("OperatorsTest: Divergence L2 FAIL");
   }

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Divergence PASS");
   }

   return Err;
}

int testGradient(Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto &Mesh      = HorzMesh::getDefault();
   const int NVertLevels = 1;

   // Prepare operator input
   Array2DReal ScalarCell("ScalarCell", Mesh->NCellsSize, NVertLevels);
   Err += setScalarCell(
       KOKKOS_LAMBDA(Real Coord1, Real Coord2) {
          return Setup.exactScalar(Coord1, Coord2);
       },
       ScalarCell, Geom, Mesh, NVertLevels);

   // Compute exact result
   Array2DReal ExactGradEdge("ExactGradEdge", Mesh->NEdgesOwned, NVertLevels);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactGradScalarX(X, Y);
          VecField[1] = Setup.exactGradScalarY(X, Y);
       },
       ExactGradEdge, EdgeComponent::Normal, Geom, Mesh, NVertLevels, false);

   // Compute numerical result
   GradientOnEdge GradientEdge(Mesh);
   Array2DReal NumGradEdge("NumGradEdge", Mesh->NEdgesOwned, NVertLevels);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int K) {
          GradientEdge(NumGradEdge, IEdge, K, ScalarCell);
       });

   // Compute error measures
   ErrorMeasures GradErrors;
   Err += computeErrorsEdge(GradErrors, NumGradEdge, ExactGradEdge, Mesh,
                            NVertLevels);

   // Check error values
   if (!isApprox(GradErrors.LInf, Setup.ExpectedGradErrorLInf, RTol)) {
      Err++;
      LOG_ERROR("OperatorsTest: Gradient LInf FAIL");
   }

   if (!isApprox(GradErrors.L2, Setup.ExpectedGradErrorL2, RTol)) {
      Err++;
      LOG_ERROR("OperatorsTest: Gradient L2 FAIL");
   }

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Gradient PASS");
   }

   return Err;
}

int testCurl(Real RTol) {
   int Err = 0;
   TestSetup Setup;
   const auto &Mesh      = HorzMesh::getDefault();
   const int NVertLevels = 1;

   // Prepare operator input
   Array2DReal VecEdge("VecEdge", Mesh->NEdgesSize, NVertLevels);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeComponent::Normal, Geom, Mesh, NVertLevels);

   // Compute exact result
   Array2DReal ExactCurlVertex("ExactCurlVertex", Mesh->NVerticesOwned,
                               NVertLevels);
   Err += setScalarVertex(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.exactCurlVec(X, Y); },
       ExactCurlVertex, Geom, Mesh, NVertLevels, false);

   // Compute numerical result
   Array2DReal NumCurlVertex("NumCurlVertex", Mesh->NVerticesOwned,
                             NVertLevels);
   CurlOnVertex CurlVertex(Mesh);
   parallelFor(
       {Mesh->NVerticesOwned, NVertLevels}, KOKKOS_LAMBDA(int IVertex, int K) {
          CurlVertex(NumCurlVertex, IVertex, K, VecEdge);
       });

   // Compute error measures
   ErrorMeasures CurlErrors;
   Err += computeErrorsVertex(CurlErrors, NumCurlVertex, ExactCurlVertex, Mesh,
                              NVertLevels);

   // Check error values
   if (!isApprox(CurlErrors.LInf, Setup.ExpectedCurlErrorLInf, RTol)) {
      Err++;
      LOG_ERROR("OperatorsTest: Curl LInf FAIL");
   }

   if (!isApprox(CurlErrors.L2, Setup.ExpectedCurlErrorL2, RTol)) {
      Err++;
      LOG_ERROR("OperatorsTest: Curl L2 FAIL");
   }

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Curl PASS");
   }

   return Err;
}

int testRecon(Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto &Mesh      = HorzMesh::getDefault();
   const int NVertLevels = 1;

   // Prepare operator input
   Array2DReal VecEdge("VecEdge", Mesh->NEdgesSize, NVertLevels);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       VecEdge, EdgeComponent::Normal, Geom, Mesh, NVertLevels);

   // Compute exact result
   Array2DReal ExactReconEdge("ExactReconEdge", Mesh->NEdgesOwned, NVertLevels);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.exactVecX(X, Y);
          VecField[1] = Setup.exactVecY(X, Y);
       },
       ExactReconEdge, EdgeComponent::Tangential, Geom, Mesh, NVertLevels,
       false);

   // Compute numerical result
   Array2DReal NumReconEdge("NumReconEdge", Mesh->NEdgesOwned, NVertLevels);
   TangentialReconOnEdge TanReconEdge(Mesh);
   parallelFor(
       {Mesh->NEdgesOwned, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int K) {
          TanReconEdge(NumReconEdge, IEdge, K, VecEdge);
       });

   // Compute error measures
   ErrorMeasures ReconErrors;
   Err += computeErrorsEdge(ReconErrors, NumReconEdge, ExactReconEdge, Mesh,
                            NVertLevels);

   if (!isApprox(ReconErrors.LInf, Setup.ExpectedReconErrorLInf, RTol)) {
      Err++;
      LOG_ERROR("OperatorsTest: Recon LInf FAIL");
   }

   if (!isApprox(ReconErrors.L2, Setup.ExpectedReconErrorL2, RTol)) {
      Err++;
      LOG_ERROR("OperatorsTest: Recon L2 FAIL");
   }

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Recon PASS");
   }

   return Err;
}

//------------------------------------------------------------------------------
// The initialization routine for Operators testing
int initOperatorsTest(const std::string &MeshFile) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefaultEnv();
   MPI_Comm DefComm = DefEnv->getComm();

   int IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("OperatorsTest: error initializing parallel IO");
   }

   int DecompErr = Decomp::init(MeshFile);
   if (DecompErr != 0) {
      Err++;
      LOG_ERROR("OperatorsTest: error initializing default decomposition");
   }

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("OperatorsTest: error initializing default halo");
   }

   int MeshErr = HorzMesh::init();
   if (MeshErr != 0) {
      Err++;
      LOG_ERROR("OperatorsTest: error initializing default mesh");
   }

   return Err;
}

void finalizeOperatorsTest() {
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

void operatorsTest(const std::string &MeshFile = DefaultMeshFile) {
   int Err = initOperatorsTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("OperatorsTest: Error initializing");
   }

   const Real RTol = sizeof(Real) == 4 ? 1e-2 : 1e-10;

   Err += testDivergence(RTol);
   Err += testGradient(RTol);
   Err += testCurl(RTol);
   Err += testRecon(RTol);

   if (Err == 0) {
      LOG_INFO("OperatorsTest: Successful completion");
   }
   finalizeOperatorsTest();
}

int main(int argc, char *argv[]) {
   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);

   operatorsTest();

   Kokkos::finalize();
   MPI_Finalize();
} // end of main
//===-----------------------------------------------------------------------===/
