#include "Tendencies.h"
#include "AuxiliaryState.h"
#include "Config.h"
#include "CustomTendencyTerms.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
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
#include "OmegaKokkos.h"
#include "PGrad.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "VertCoord.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>

using namespace OMEGA;

struct TestSetup {
   Real Radius = REarth;

   KOKKOS_FUNCTION Real pseudoThickness(Real Lon, Real Lat) const {
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
constexpr int NVertLayers = 60;

int initState() {
   int Err = 0;

   TestSetup Setup;
   auto *Mesh   = HorzMesh::getDefault();
   auto *VCoord = VertCoord::getDefault();
   auto *State  = OceanState::getDefault();

   // Define tendency fields
   int NDims = 2;
   std::vector<std::string> DimNamesThickness(NDims);
   DimNamesThickness[0] = "NCells";

   Array2DReal PseudoThickCell = State->getPseudoThickness(0);
   Array2DReal NormalVelEdge   = State->getNormalVelocity(0);

   Array3DReal TracersArray = Tracers::getAll(0);
   const auto &TracersCell  = TracersArray;

   int NTracers = Tracers::getNumTracers();

   deepCopy(PseudoThickCell, NAN);
   deepCopy(NormalVelEdge, NAN);
   deepCopy(TracersCell, NAN);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.pseudoThickness(X, Y); },
       PseudoThickCell, Geom, Mesh, OnCell, VCoord->MinLayerCell,
       VCoord->MaxLayerCell, ExchangeHalos::Yes, SetBoundary::Yes);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracer(X, Y); },
       TracersCell, Geom, Mesh, OnCell, VCoord->MinLayerCell,
       VCoord->MaxLayerCell, ExchangeHalos::Yes, SetBoundary::Yes);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real Lon, Real Lat) {
          VecField[0] = Setup.velocityX(Lon, Lat);
          VecField[1] = Setup.velocityY(Lon, Lat);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh,
       VCoord->MinLayerEdgeTop, VCoord->MaxLayerEdgeBot, ExchangeHalos::Yes,
       CartProjection::No, SetBoundary::Yes);

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
   Config::readAll("omega.yml");

   TimeStepper::init1();

   IO::init(DefComm);
   Decomp::init(mesh);

   // Initialize streams
   IOStream::init();

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default halo");
   }

   HorzMesh::init();
   VertCoord::init();
   Tracers::init();
   VertAdv::init();
   PressureGrad::init();
   Eos::init();

   int StateErr = OceanState::init();
   if (StateErr != 0) {
      Err++;
      LOG_ERROR("TendenciesTest: error initializing default state");
   }

   AuxiliaryState::init();

   return Err;
}

int testTendencies() {
   int Err = 0;
   Error Err1;

   // test initialization
   Tendencies::init();

   // test retrievel of default
   Tendencies *DefTendencies = Tendencies::getDefault();

   if (DefTendencies) {
      LOG_INFO("TendenciesTest: Default tendencies retrieval PASS");
   } else {
      Err++;
      LOG_INFO("TendenciesTest: Default tendencies retrieval FAIL");
      return -1;
   }

   const auto Mesh     = HorzMesh::getDefault();
   const auto VCoord   = VertCoord::getDefault();
   const auto VAdv     = VertAdv::getDefault();
   const auto PGrad    = PressureGrad::getDefault();
   const auto EqState  = Eos::getInstance();
   VCoord->NVertLayers = 12;

   // test creation of another tendencies

   TimeInterval ZeroTimeStep; // Zero-length time step placeholder
   Config *Options = Config::getOmegaConfig();
   Config TendConfig("Tendencies");
   Err1             = Options->get(TendConfig);
   int NTracersTest = 3;

   Tendencies::create("TestTendencies", Mesh, VCoord, VAdv, PGrad, EqState,
                      NTracersTest, ZeroTimeStep, &TendConfig);

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

   VCoord->NVertLayers = NVertLayers;

   // put NANs in every tendency variables
   deepCopy(DefTendencies->PseudoThicknessTend, NAN);
   deepCopy(DefTendencies->NormalVelocityTend, NAN);
   deepCopy(DefTendencies->TracerTend, NAN);

   // compute tendencies
   const auto *State       = OceanState::getDefault();
   const auto *AuxState    = AuxiliaryState::getDefault();
   Array3DReal TracerArray = Tracers::getAll(0);
   int ThickTimeLevel      = 0;
   int VelTimeLevel        = 0;
   int TracerTimeLevel     = 0;
   TimeInstant Time;
   TimeInterval Interval(1., TimeUnits::Seconds);
   DefTendencies->computeAllTendencies(State, AuxState, TracerArray,
                                       ThickTimeLevel, VelTimeLevel,
                                       TracerTimeLevel, Time, Interval);

   // check that everything got computed correctly
   int NCellsOwned = Mesh->NCellsOwned;
   int NEdgesOwned = Mesh->NEdgesOwned;
   int NTracers    = Tracers::getNumTracers();

   const Real PseudoThickTendSum =
       sum(DefTendencies->PseudoThicknessTend, NCellsOwned,
           VCoord->MinLayerCell, VCoord->MaxLayerCell);
   if (!Kokkos::isfinite(PseudoThickTendSum) || PseudoThickTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: PseudoThickTend FAIL");
   }

   const Real NormVelTendSum =
       sum(DefTendencies->NormalVelocityTend, NEdgesOwned,
           VCoord->MinLayerEdgeBot, VCoord->MaxLayerEdgeTop);
   if (!Kokkos::isfinite(NormVelTendSum) || NormVelTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: NormVelTendSum FAIL");
   }

   const Real TraceTendSum =
       sum(DefTendencies->TracerTend, NTracers, NCellsOwned,
           VCoord->MinLayerCell, VCoord->MaxLayerCell);
   if (!Kokkos::isfinite(TraceTendSum) || TraceTendSum == 0) {
      Err++;
      LOG_ERROR("TendenciesTest: TraceTendSum FAIL");
   }

   Tendencies::clear();

   return Err;
}

void finalizeTendenciesTest() {
   Tracers::clear();
   PressureGrad::clear();
   Eos::destroyInstance();
   AuxiliaryState::clear();
   OceanState::clear();
   VertAdv::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   TimeStepper::clear();
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
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += tendenciesTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
