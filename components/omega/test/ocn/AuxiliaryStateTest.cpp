#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOField.h"
#include "Logging.h"
#include "MachEnv.h"
#include "MetaData.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
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
};

constexpr Geometry Geom   = Geometry::Spherical;
constexpr int NVertLevels = 60;

int initState() {
   int Err = 0;

   TestSetup Setup;
   auto *Mesh  = HorzMesh::getDefault();
   auto *State = OceanState::getDefault();

   const auto &LayerThickCell = State->LayerThickness[0];
   const auto &NormalVelEdge  = State->NormalVelocity[0];

   Err += setScalar(
       KOKKOS_LAMBDA(Real X, Real Y) { return Setup.layerThickness(X, Y); },
       LayerThickCell, Geom, Mesh, OnCell, NVertLevels);

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
// The initialization routine for aux vars testing
int initAuxStateTest(const std::string &mesh) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   int IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("AuxStateTest: error initializing parallel IO");
   }

   int DecompErr = Decomp::init(mesh);
   if (DecompErr != 0) {
      Err++;
      LOG_ERROR("AuxStateTest: error initializing default decomposition");
   }

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("AuxStateTest: error initializing default halo");
   }

   int MeshErr = HorzMesh::init();
   if (MeshErr != 0) {
      Err++;
      LOG_ERROR("AuxStateTest: error initializing default mesh");
   }

   const auto &Mesh = HorzMesh::getDefault();
   MetaDim::create("NCells", Mesh->NCellsSize);
   MetaDim::create("NVertices", Mesh->NVerticesSize);
   MetaDim::create("NEdges", Mesh->NEdgesSize);
   MetaDim::create("NVertLevels", NVertLevels);

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
   int AuxStateErr = AuxiliaryState::init();
   if (AuxStateErr != 0) {
      Err++;
      LOG_ERROR("AuxStateTest: error initializing default aux state");
   }

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
   // test creation of another auxiliary state
   AuxiliaryState::create("AnotherAuxState", Mesh, 12);

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

   // compute auxiliary variables
   const auto *State = OceanState::getDefault();
   DefAuxState->computeAll(State, 0);

   // check that everything got computed correctly
   int NCellsOwned    = Mesh->NCellsOwned;
   int NEdgesOwned    = Mesh->NEdgesOwned;
   int NVerticesOwned = Mesh->NVerticesOwned;

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

   AuxiliaryState::clear();

   return Err;
}

void finalizeAuxStateTest() {
   OceanState::clear();
   IOField::clear();
   HorzMesh::clear();
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

   RetVal += auxStateTest();

   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
