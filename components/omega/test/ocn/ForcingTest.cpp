//===-- ocn/ForcingTest.cpp - Forcing Unit Test ---------------*- C++ -*-===//
//
/// \file
/// \brief Scoped tests for the Forcing class (SfcStress + SfcStressForcing)
//
//===----------------------------------------------------------------------===//

#include "Forcing.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "forcingVars/SfcStressForcingVars.h"
#include "mpi.h"

#include <limits>
#include <string>

using namespace OMEGA;

namespace {

struct TestSetupPlane {

   Real Lx = 1;
   Real Ly = SqrtThree / 2;

   ErrorMeasures ExpectedNormalStressErrors = {0.0033910709836867704,
                                               0.0039954090464502795};
   KOKKOS_FUNCTION Real sfcStressX(Real X, Real Y) const {
      return std::cos(TwoPi * X / Lx) * std::sin(TwoPi * Y / Ly);
   }

   KOKKOS_FUNCTION Real sfcStressY(Real X, Real Y) const {
      return std::sin(TwoPi * X / Lx) * std::cos(TwoPi * Y / Ly);
   }
};

struct TestSetupSphere {
   ErrorMeasures ExpectedNormalStressErrors = {0.0038588958862868362,
                                               0.003813760171030077};
   KOKKOS_FUNCTION Real sfcStressX(Real Lon, Real Lat) const {
      return -4 * std::sin(Lon) * std::cos(Lon) * std::pow(std::cos(Lat), 3) *
             std::sin(Lat);
   }

   KOKKOS_FUNCTION Real sfcStressY(Real Lon, Real Lat) const {
      return -std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 3);
   }
};

#ifdef FORCING_TEST_PLANE
constexpr Geometry Geom          = Geometry::Planar;
constexpr char DefaultMeshFile[] = "OmegaPlanarMesh.nc";
using TestSetup                  = TestSetupPlane;
#else
constexpr Geometry Geom          = Geometry::Spherical;
constexpr char DefaultMeshFile[] = "OmegaSphereMesh.nc";
using TestSetup                  = TestSetupSphere;
#endif

Forcing *TestForcing = nullptr;

int testSfcStressForcingVars(Real RTol) {
   int Err = 0;
   TestSetup Setup;

   const auto Mesh = HorzMesh::getDefault();

   // Compute exact result

   Array1DReal ExactNormalStressEdge("ExactNormalStressEdge",
                                     Mesh->NEdgesOwned);
   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real X, Real Y) {
          VecField[0] = Setup.sfcStressX(X, Y);
          VecField[1] = Setup.sfcStressY(X, Y);
       },
       ExactNormalStressEdge, EdgeComponent::Normal, Geom, Mesh,
       ExchangeHalos::No);

   SfcStressForcingVars SfcStressForcing("", Mesh);
   SfcStressForcing.InterpChoice = InterpCellToEdgeOption::Anisotropic;

   // Set inputs
   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.sfcStressX(X, Y); },
       SfcStressForcing.ZonalStressCell, Geom, Mesh, OnCell);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.sfcStressY(X, Y); },
       SfcStressForcing.MeridStressCell, Geom, Mesh, OnCell);

   // Compute numerical result
   parallelFor(
       {Mesh->NEdgesOwned},
       KOKKOS_LAMBDA(int IEdge) { SfcStressForcing.computeVarsOnEdge(IEdge); });
   const auto &NumNormalStressEdge = SfcStressForcing.NormalStressEdge;

   // Compute error measures and check error values
   ErrorMeasures NormalStressErrors;
   Err += computeErrors(NormalStressErrors, NumNormalStressEdge,
                        ExactNormalStressEdge, Mesh, OnEdge);

   Err += checkErrors("ForcingVarsTest", "NormalStress", NormalStressErrors,
                      Setup.ExpectedNormalStressErrors, RTol);

   if (Err == 0) {
      LOG_INFO("ForcingVarsTest: testSfcStressForcingVars PASS");
   }

   return Err;
}

int initForcingTest(const std::string &MeshFile) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);

   Config("Omega");
   Config::readAll("omega.yml");

   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);
   Decomp::init(MeshFile);

   Field::init(ModelClock);
   IOStream::init(ModelClock);

   Err = Halo::init();
   if (Err != 0) {
      ABORT_ERROR("ForcingTest: error initializing default halo");
   }

   HorzMesh::init(ModelClock);

   const HorzMesh *DefMesh = HorzMesh::getDefault();
   Halo *DefHalo           = Halo::getDefault();
   TestForcing             = Forcing::create("Default", DefMesh, DefHalo);
   if (TestForcing == nullptr) {
      ABORT_ERROR("ForcingTest: failed creating default forcing instance");
   }

   Config *OmegaConfig = Config::getOmegaConfig();
   TestForcing->readConfigOptions(OmegaConfig);

   return 0;
}

void finalizeForcingTest() {
   Forcing::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   IOStream::finalize();
   TimeStepper::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int testForcingInitAndConfig() {
   // Verify forcing setup consumes Omega.SfcStress config and maps InterpType.
   int Err = 0;

   Forcing *DefForcing = TestForcing;
   if (DefForcing == nullptr) {
      LOG_ERROR("ForcingTest: test forcing instance is null");
      return 1;
   }

   Config *OmegaConfig = Config::getOmegaConfig();
   if (OmegaConfig == nullptr) {
      LOG_ERROR("ForcingTest: Omega config unavailable");
      return 1;
   }

   Error CfgErr;
   Config SfcStressConfig("SfcStress");
   CfgErr = OmegaConfig->get(SfcStressConfig);
   CHECK_ERROR_ABORT(CfgErr, "ForcingTest: missing Omega.SfcStress config");

   std::string InterpType;
   CfgErr = SfcStressConfig.get("InterpType", InterpType);
   CHECK_ERROR_ABORT(CfgErr, "ForcingTest: missing SfcStress.InterpType");

   InterpCellToEdgeOption ExpectedChoice;
   if (InterpType == "Isotropic") {
      ExpectedChoice = InterpCellToEdgeOption::Isotropic;
   } else if (InterpType == "Anisotropic") {
      ExpectedChoice = InterpCellToEdgeOption::Anisotropic;
   } else {
      LOG_ERROR("ForcingTest: unknown InterpType in config: {}", InterpType);
      return 1;
   }

   if (DefForcing->SfcStressForcing.InterpChoice != ExpectedChoice) {
      LOG_ERROR("ForcingTest: InterpChoice mismatch after forcing setup");
      Err++;
   }

   if (Err == 0) {
      LOG_INFO("ForcingTest: Init/config PASS");
   }

   return Err;
}

int testForcingComputeAll() {
   // Verify edge-normal stress is the expected projection of zonal/meridional
   // stress components on edge orientation.
   int Err = 0;

   const HorzMesh *Mesh = HorzMesh::getDefault();
   Forcing *DefForcing  = TestForcing;
   if (Mesh == nullptr || DefForcing == nullptr) {
      LOG_ERROR("ForcingTest: missing mesh or forcing for compute test");
      return 1;
   }

   auto &ZonalStressCell = DefForcing->SfcStressForcing.ZonalStressCell;
   auto &MeridStressCell = DefForcing->SfcStressForcing.MeridStressCell;
   auto &NormalStress    = DefForcing->SfcStressForcing.NormalStressEdge;
   const auto AngleEdge  = Mesh->AngleEdge;

   // First pass: pure zonal stress should project to cos(AngleEdge).
   deepCopy(ZonalStressCell, 1._Real);
   deepCopy(MeridStressCell, 0._Real);
   deepCopy(NormalStress, 0._Real);

   DefForcing->computeAll();

   Real MaxErrCos = 0._Real;
   parallelReduce(
       {Mesh->NEdgesOwned},
       KOKKOS_LAMBDA(int IEdge, Real &LocalMax) {
          const Real Expected = Kokkos::cos(AngleEdge(IEdge));
          const Real AbsErr   = Kokkos::abs(NormalStress(IEdge) - Expected);
          if (AbsErr > LocalMax) {
             LocalMax = AbsErr;
          }
       },
       Kokkos::Max<Real>(MaxErrCos));

   Real GlobalMaxErrCos = 0._Real;
   MPI_Datatype MpiReal =
       sizeof(Real) == sizeof(double) ? MPI_DOUBLE : MPI_FLOAT;
   MPI_Allreduce(&MaxErrCos, &GlobalMaxErrCos, 1, MpiReal, MPI_MAX,
                 MachEnv::getDefault()->getComm());

   // Second pass: pure meridional stress should project to sin(AngleEdge).
   deepCopy(ZonalStressCell, 0._Real);
   deepCopy(MeridStressCell, 1._Real);
   deepCopy(NormalStress, 0._Real);

   DefForcing->computeAll();

   Real MaxErrSin = 0._Real;
   parallelReduce(
       {Mesh->NEdgesOwned},
       KOKKOS_LAMBDA(int IEdge, Real &LocalMax) {
          const Real Expected = Kokkos::sin(AngleEdge(IEdge));
          const Real AbsErr   = Kokkos::abs(NormalStress(IEdge) - Expected);
          if (AbsErr > LocalMax) {
             LocalMax = AbsErr;
          }
       },
       Kokkos::Max<Real>(MaxErrSin));

   Real GlobalMaxErrSin = 0._Real;
   MPI_Allreduce(&MaxErrSin, &GlobalMaxErrSin, 1, MpiReal, MPI_MAX,
                 MachEnv::getDefault()->getComm());

   const Real Tol       = 1e-11;
   const Real GlobalErr = Kokkos::max(GlobalMaxErrCos, GlobalMaxErrSin);
   // Expected outcome: both projection passes remain below strict tolerance.
   if (GlobalErr > Tol) {
      LOG_ERROR(
          "ForcingTest: normal stress mismatch in cos/sin checks, max error "
          "{} > tol {}",
          GlobalErr, Tol);
      Err++;
   }

   if (Err == 0) {
      LOG_INFO("ForcingTest: computeAll PASS");
   }

   return Err;
}

int testForcingAPISmoke() {
   // Verify named Forcing lifecycle semantics: create/get/exists/erase.
   int Err = 0;

   const std::string Name = "UnitTestForcing";

   if (Forcing::exists(Name)) {
      Forcing::erase(Name);
   }

   Forcing *Named =
       Forcing::create(Name, HorzMesh::getDefault(), Halo::getDefault());
   if (Named == nullptr) {
      LOG_ERROR("ForcingTest: failed creating named forcing instance");
      Err++;
   }

   if (!Forcing::exists(Name)) {
      LOG_ERROR("ForcingTest: created forcing instance not found in map");
      Err++;
   }

   if (Forcing::get(Name) == nullptr) {
      LOG_ERROR("ForcingTest: get() failed for named forcing instance");
      Err++;
   }

   Forcing::erase(Name);
   if (Forcing::exists(Name)) {
      LOG_ERROR("ForcingTest: erase() failed for named forcing instance");
      Err++;
   }

   if (Err == 0) {
      LOG_INFO("ForcingTest: API smoke PASS");
   }

   return Err;
}

int forcingTest() {
   // Execute Forcing checks in order: init/config, compute projection,
   // direct vars analytic check, then API lifecycle smoke test.
   int Err         = 0;
   const Real RTol = sizeof(Real) == 4 ? 1e-2 : 2e-4;

   Err += initForcingTest(DefaultMeshFile);
   Err += testForcingInitAndConfig();
   Err += testForcingComputeAll();
   Err += testSfcStressForcingVars(RTol);
   Err += testForcingAPISmoke();

   if (Err == 0) {
      LOG_INFO("ForcingTest: Successful completion");
   }

   finalizeForcingTest();
   return Err;
}

} // namespace

int main(int argc, char *argv[]) {
   int RetErr = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetErr = forcingTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return RetErr;
}
