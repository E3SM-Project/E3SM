#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>

using namespace OMEGA;

struct TestSetup {
   Real Radius = REarth;

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
constexpr int NVertLayers = 60;

int initState() {
   int Err = 0;

   TestSetup Setup;
   auto *Mesh  = HorzMesh::getDefault();
   auto *State = OceanState::getDefault();

   Array2DReal LayerThickCell;
   Array2DReal NormalVelEdge;
   Err = State->getLayerThickness(LayerThickCell, 0);
   Err += State->getNormalVelocity(NormalVelEdge, 0);
   Array3DReal TracerArray;
   Err += Tracers::getAll(TracerArray, 0);

   int NTracers = Tracers::getNumTracers();

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.layerThickness(X, Y); },
       LayerThickCell, Geom, Mesh, OnCell);

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.tracer(X, Y); },
       TracerArray, Geom, Mesh, OnCell);

   Err += setVectorEdge(
       KOKKOS_LAMBDA(Real(&VecField)[2], Real Lon, Real Lat) {
          VecField[0] = Setup.velocityX(Lon, Lat);
          VecField[1] = Setup.velocityY(Lon, Lat);
       },
       NormalVelEdge, EdgeComponent::Normal, Geom, Mesh, ExchangeHalos::Yes,
       CartProjection::No);

   return Err;
}

//------------------------------------------------------------------------------
// The initialization routine for aux vars testing
int initAuxStateTest(const std::string &mesh) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("------ Auxiliary State unit tests ------");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   TimeStepper::init1();

   IO::init(DefComm);
   Decomp::init(mesh);

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("AuxStateTest: error initializing default halo");
   }

   VertCoord::init1();
   HorzMesh::init();

   Tracers::init();
   int StateErr = OceanState::init();
   if (StateErr != 0) {
      Err++;
      LOG_ERROR("AuxStateTest: error initializing default state");
   }

   return Err;
}

int testAuxState() {
   int Err = 0;

   // test initialization
   AuxiliaryState::init();

   // test retrievel of default
   AuxiliaryState *DefAuxState = AuxiliaryState::getDefault();

   if (DefAuxState) {
      LOG_INFO("AuxStateTest: Default auxstate retrieval PASS");
   } else {
      Err++;
      LOG_INFO("AuxStateTest: Default auxstate retrieval FAIL");
      return -1;
   }

   const auto *Mesh = HorzMesh::getDefault();
   auto *MeshHalo   = Halo::getDefault();
   // test creation of another auxiliary state
   AuxiliaryState::create("AnotherAuxState", Mesh, MeshHalo, 12, 3);

   // test retrievel of another
   if (AuxiliaryState::get("AnotherAuxState")) {
      LOG_INFO("AuxStateTest: Non-default auxstate retrieval PASS");
   } else {
      Err++;
      LOG_INFO("AuxStateTest: Non-default auxstate retrieval FAIL");
   }

   // test erase
   AuxiliaryState::erase("AnotherAuxState");

   if (AuxiliaryState::get("AnotherAuxState")) {
      Err++;
      LOG_INFO("AuxStateTest: Non-default auxstate retrieval FAIL");
   } else {
      LOG_INFO("AuxStateTest: Non-default auxstate retrieval PASS");
   }

   // put NANs in every auxiliary variable
   deepCopy(DefAuxState->KineticAux.KineticEnergyCell, NAN);
   deepCopy(DefAuxState->KineticAux.VelocityDivCell, NAN);

   deepCopy(DefAuxState->LayerThicknessAux.FluxLayerThickEdge, NAN);
   deepCopy(DefAuxState->LayerThicknessAux.MeanLayerThickEdge, NAN);

   deepCopy(DefAuxState->VorticityAux.RelVortVertex, NAN);
   deepCopy(DefAuxState->VorticityAux.NormRelVortVertex, NAN);
   deepCopy(DefAuxState->VorticityAux.NormPlanetVortVertex, NAN);
   deepCopy(DefAuxState->VorticityAux.NormRelVortEdge, NAN);
   deepCopy(DefAuxState->VorticityAux.NormPlanetVortEdge, NAN);

   deepCopy(DefAuxState->VelocityDel2Aux.Del2Edge, NAN);
   deepCopy(DefAuxState->VelocityDel2Aux.Del2DivCell, NAN);
   deepCopy(DefAuxState->VelocityDel2Aux.Del2RelVortVertex, NAN);

   deepCopy(DefAuxState->TracerAux.HTracersEdge, NAN);
   deepCopy(DefAuxState->TracerAux.Del2TracersCell, NAN);

   // compute auxiliary variables
   const auto *State = OceanState::getDefault();
   Array3DReal TracerArray;
   Err += Tracers::getAll(TracerArray, 0);
   DefAuxState->computeAll(State, TracerArray, 0);

   // check that everything got computed correctly
   int NCellsOwned    = Mesh->NCellsOwned;
   int NEdgesOwned    = Mesh->NEdgesOwned;
   int NVerticesOwned = Mesh->NVerticesOwned;
   int NTracers       = Tracers::getNumTracers();

   const Real KineticEnergySum =
       sum(DefAuxState->KineticAux.KineticEnergyCell, NCellsOwned);
   if (!Kokkos::isfinite(KineticEnergySum)) {
      Err++;
      LOG_ERROR("AuxStateTest: KineticEnergy FAIL");
   }

   const Real VelocityDivSum =
       sum(DefAuxState->KineticAux.VelocityDivCell, NCellsOwned);
   if (!Kokkos::isfinite(VelocityDivSum)) {
      Err++;
      LOG_ERROR("AuxStateTest: VelocityDivCell FAIL");
   }

   const Real FluxLayerThickSum =
       sum(DefAuxState->LayerThicknessAux.FluxLayerThickEdge, NEdgesOwned);
   if (!Kokkos::isfinite(FluxLayerThickSum)) {
      Err++;
      LOG_ERROR("AuxStateTest: FluxLayerThickEdge FAIL");
   }

   const Real MeanLayerThickSum =
       sum(DefAuxState->LayerThicknessAux.MeanLayerThickEdge, NEdgesOwned);
   if (!Kokkos::isfinite(MeanLayerThickSum)) {
      Err++;
      LOG_ERROR("AuxStateTest: MeanLayerThickEdge FAIL");
   }

   const Real RelVortVSum =
       sum(DefAuxState->VorticityAux.RelVortVertex, NVerticesOwned);
   if (!Kokkos::isfinite(RelVortVSum)) {
      Err++;
      LOG_ERROR("AuxStateTest: RelVortVertex FAIL");
   }

   const Real NormRelVortVSum =
       sum(DefAuxState->VorticityAux.NormRelVortVertex, NVerticesOwned);
   if (!Kokkos::isfinite(NormRelVortVSum)) {
      Err++;
      LOG_ERROR("AuxStateTest: NormRelVortVertex FAIL");
   }

   const Real NormPlanetVortVSum =
       sum(DefAuxState->VorticityAux.NormPlanetVortVertex, NVerticesOwned);
   if (!Kokkos::isfinite(NormPlanetVortVSum)) {
      Err++;
      LOG_ERROR("AuxStateTest: NormPlanetVortVertex FAIL");
   }

   const Real NormRelVortESum =
       sum(DefAuxState->VorticityAux.NormRelVortEdge, NEdgesOwned);
   if (!Kokkos::isfinite(NormRelVortESum)) {
      Err++;
      LOG_ERROR("AuxStateTest: NormRelVortEdge FAIL");
   }

   const Real NormPlanetVortESum =
       sum(DefAuxState->VorticityAux.NormPlanetVortEdge, NEdgesOwned);
   if (!Kokkos::isfinite(NormPlanetVortESum)) {
      Err++;
      LOG_ERROR("AuxStateTest: NormPlanetVortEdge FAIL");
   }

   const Real HTracersESum =
       sum(DefAuxState->TracerAux.HTracersEdge, NTracers, NEdgesOwned);
   if (!Kokkos::isfinite(HTracersESum)) {
      Err++;
      LOG_ERROR("AuxStateTest: HTracersOnEdge FAIL");
   }

   const Real Del2TracersCSum =
       sum(DefAuxState->TracerAux.Del2TracersCell, NTracers, NCellsOwned);
   if (!Kokkos::isfinite(Del2TracersCSum)) {
      Err++;
      LOG_ERROR("AuxStateTest: Del2TracersOnCell FAIL");
   }

   AuxiliaryState::clear();

   return Err;
}

void finalizeAuxStateTest() {
   Tracers::clear();
   OceanState::clear();
   VertCoord::clear();
   Field::clear();
   Dimension::clear();
   TimeStepper::clear();
   HorzMesh::clear();
   VertCoord::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int auxStateTest(const std::string &MeshFile = "OmegaMesh.nc") {
   int Err = initAuxStateTest(MeshFile);
   if (Err != 0) {
      LOG_CRITICAL("AuxStateTest: Error initializing");
   }

   const auto &Mesh = HorzMesh::getDefault();

   Err += initState();

   Err += testAuxState();

   if (Err == 0) {
      LOG_INFO("AuxStateTest: Successful completion");
   }
   finalizeAuxStateTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += auxStateTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;
   if (RetVal == 0)
      LOG_INFO("------ Auxiliary State unit tests successful ------");

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
