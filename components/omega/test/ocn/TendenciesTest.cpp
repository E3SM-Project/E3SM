#include "Tendencies.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "TimeStepper.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>

using namespace OMEGA;

struct TestSetup {
   Real Radius = 6371220;

   KOKKOS_FUNCTION Real layerThickness(Real Lon, Real Lat) const {
      return (2 + std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }

   KOKKOS_FUNCTION Real velocityX(Real Lon, Real Lat) const {
      return -Radius * std::pow(std::sin(Lon), 2) * std::pow(std::cos(Lat), 3);
   }

   KOKKOS_FUNCTION Real velocityY(Real Lon, Real Lat) const {
      return -4 * Radius * std::sin(Lon) * std::cos(Lon) *
             std::pow(std::cos(Lat), 3) * std::sin(Lat);
   }

   KOKKOS_FUNCTION Real tracer(Real Lon, Real Lat) const {
      return (2 - std::cos(Lon) * std::pow(std::cos(Lat), 4));
   }
};

constexpr Geometry Geom   = Geometry::Spherical;
constexpr int NVertLevels = 60;

int initState() {
   int Err = 0;

   TestSetup Setup;
   auto *Mesh  = HorzMesh::getDefault();
   auto *State = OceanState::getDefault();

   Array2DReal LayerThickCell;
   Array2DReal NormalVelEdge;
   Err += State->getLayerThickness(LayerThickCell, 0);
   Err += State->getNormalVelocity(NormalVelEdge, 0);

   Array3DReal TracersArray;
   Err += Tracers::getAll(TracersArray, 0);
   const auto &TracersCell = TracersArray;

   int NTracers = Tracers::getNumTracers();

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.layerThickness(X, Y); },
       LayerThickCell, Geom, Mesh, OnCell, NVertLevels);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracer(X, Y); },
       TracersCell, Geom, Mesh, OnCell, NVertLevels, NTracers);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real Lon, Real Lat) {
          VecField[0] = Setup.velocityX(Lon, Lat);
          VecField[1] = Setup.velocityY(Lon, Lat);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh, NVertLevels,
       ExchangeHalos::Yes, CartProjection::No);

   return Err;
}

//------------------------------------------------------------------------------
// The initialization routine for tendencies testing
int initTendenciesTest(const std::string &mesh) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("TendenciesTest: Error reading config file");
      return Err;
   }

   int TimeStepperErr = TimeStepper::init1();
   if (TimeStepperErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default time stepper");
   }

   int IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing parallel IO");
   }

   int DecompErr = Decomp::init(mesh);
   if (DecompErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default decomposition");
   }

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default halo");
   }

   int MeshErr = HorzMesh::init();
   if (MeshErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default mesh");
   }

   int TracerErr = Tracers::init();
   if (TracerErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing tracer infrastructure");
   }

   const auto &Mesh = HorzMesh::getDefault();
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLevels", NVertLevels);

   int StateErr = OceanState::init();
   if (StateErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default state");
   }

   int AuxStateErr = AuxiliaryState::init();
   if (AuxStateErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default aux state");
   }

   return Err;
}

int testTendencies() {
   int Err = 0;

   // test initialization
   int TendenciesErr = Tendencies::init();
   if (TendenciesErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default tendencies");
   }

   // test retrievel of default
   Tendencies *DefTendencies = Tendencies::getDefault();

   if (DefTendencies) {
      LOG_INFO("TendenciesTest: Default tendencies retrieval PASS");
   } else {
      Err++;
      LOG_INFO("TendenciesTest: Default tendencies retrieval FAIL");
      return -1;
   }

   const auto *Mesh = HorzMesh::getDefault();
   // test creation of another tendencies
   Config *Options = Config::getOmegaConfig();
   Tendencies::create("TestTendencies", Mesh, 12, 3, Options);

   // test retrievel of another tendencies
   if (Tendencies::get("TestTendencies")) {
      LOG_INFO("TendenciesTest: Non-default tendencies retrieval PASS");
   } else {
      Err++;
      LOG_INFO("TendenciesTest: Non-default tendencies retrieval FAIL");
   }

   // test erase
   Tendencies::erase("TestTendencies");

   if (Tendencies::get("TestTendencies")) {
      Err++;
      LOG_INFO("TendenciesTest: Non-default tendencies erase FAIL");
   } else {
      LOG_INFO("TendenciesTest: Non-default tendencies erase PASS");
   }

   // put NANs in every tendency variables
   deepCopy(DefTendencies->LayerThicknessTend, NAN);
   deepCopy(DefTendencies->NormalVelocityTend, NAN);
   deepCopy(DefTendencies->TracerTend, NAN);

   // compute tendencies
   const auto *State    = OceanState::getDefault();
   const auto *AuxState = AuxiliaryState::getDefault();
   Array3DReal TracerArray;
   Err += Tracers::getAll(TracerArray, 0);
   int ThickTimeLevel = 0;
   int VelTimeLevel   = 0;
   TimeInstant Time;
   DefTendencies->computeAllTendencies(State, AuxState, TracerArray,
                                       ThickTimeLevel, VelTimeLevel, Time);

   // check that everything got computed correctly
   int NCellsOwned    = Mesh->NCellsOwned;
   int NEdgesOwned    = Mesh->NEdgesOwned;
   int NVerticesOwned = Mesh->NVerticesOwned;
   int NTracers       = Tracers::getNumTracers();

   const Real LayerThickTendSum =
       sum(DefTendencies->LayerThicknessTend, NCellsOwned);
   if (!Kokkos::isfinite(LayerThickTendSum) || LayerThickTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: LayerThickTend FAIL");
   }

   const Real NormVelTendSum =
       sum(DefTendencies->NormalVelocityTend, NEdgesOwned);
   if (!Kokkos::isfinite(NormVelTendSum) || NormVelTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: NormVelTendSum FAIL");
   }

   const Real TraceTendSum =
       sum(DefTendencies->TracerTend, NTracers, NCellsOwned);
   if (!Kokkos::isfinite(TraceTendSum) || TraceTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: TraceTendSum FAIL");
   }

   Tendencies::clear();

   return Err;
}

void finalizeTendenciesTest() {
   Tracers::clear();
   AuxiliaryState::clear();
   OceanState::clear();
   Field::clear();
   Dimension::clear();
   TimeStepper::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int tendenciesTest(const std::string &MeshFile = "OmegaMesh.nc") {
   int Err = initTendenciesTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("TendenciesTest: Error initializing");
   }

   const auto &Mesh = HorzMesh::getDefault();

   Err += initState();

   Err += testTendencies();

   if (Err == 0) {
      LOG_INFO("TendenciesTest: Successful completion");
   }
   finalizeTendenciesTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);

   RetVal += tendenciesTest();

   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
