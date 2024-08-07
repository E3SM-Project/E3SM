#include "TimeStepper.h"
#include "../ocn/OceanTestCommon.h"
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
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "TendencyTerms.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace OMEGA;

constexpr Geometry Geom   = Geometry::Planar;
constexpr int NVertLevels = 1;

struct DecayThicknessTendency {
   void operator()(Array2DReal ThicknessTend, OceanState *State,
                   AuxiliaryState *AuxState, int TimeLevel, Real Time) const {}
};

struct DecayVelocityTendency {
   Real Coeff = 0.5;

   Real exactSolution(Real Time) { return std::exp(-Coeff * Time); }

   void operator()(Array2DReal NormalVelTend, OceanState *State,
                   AuxiliaryState *AuxState, int TimeLevel, Real Time) const {

      auto *Mesh                = HorzMesh::getDefault();
      auto NVertLevels          = NormalVelTend.extent_int(1);
      const auto &NormalVelEdge = State->NormalVelocity[TimeLevel];

      OMEGA_SCOPE(LocCoeff, Coeff);

      parallelFor(
          {Mesh->NEdgesAll, NVertLevels}, KOKKOS_LAMBDA(int IEdge, int K) {
             NormalVelTend(IEdge, K) -= LocCoeff * NormalVelEdge(IEdge, K);
          });
   }
};

int initState() {
   int Err = 0;

   auto *Mesh  = HorzMesh::getDefault();
   auto *State = OceanState::get("TestState");

   const auto &LayerThickCell = State->LayerThickness[0];
   const auto &NormalVelEdge  = State->NormalVelocity[0];

   deepCopy(LayerThickCell, 1);
   deepCopy(NormalVelEdge, 1);

   return Err;
}

int createExactSolution(Real TimeEnd) {
   int Err = 0;

   auto *DefHalo = Halo::getDefault();
   auto *DefMesh = HorzMesh::getDefault();

   auto *ExactState =
       OceanState::create("Exact", DefMesh, DefHalo, NVertLevels, 1);

   const auto &LayerThickCell = ExactState->LayerThickness[0];
   const auto &NormalVelEdge  = ExactState->NormalVelocity[0];

   deepCopy(LayerThickCell, 1);
   deepCopy(NormalVelEdge, DecayVelocityTendency{}.exactSolution(TimeEnd));

   return Err;
}

ErrorMeasures computeErrors() {
   const auto *DefMesh = HorzMesh::getDefault();

   const auto *State      = OceanState::get("TestState");
   const auto *ExactState = OceanState::get("Exact");

   const auto &NormalVelEdge      = State->NormalVelocity[0];
   const auto &ExactNormalVelEdge = ExactState->NormalVelocity[0];

   ErrorMeasures VelErrors;
   computeErrors(VelErrors, NormalVelEdge, ExactNormalVelEdge, DefMesh, OnEdge,
                 NVertLevels);

   return VelErrors;
}

//------------------------------------------------------------------------------
// The initialization routine for time stepper testing
int initTimeStepperTest(const std::string &mesh) {
   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // default init

   int IOErr = IO::init(DefComm);
   if (IOErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing parallel IO");
   }

   int DecompErr = Decomp::init(mesh);
   if (DecompErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing default decomposition");
   }

   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing default halo");
   }

   int MeshErr = HorzMesh::init();
   if (MeshErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing default mesh");
   }

   // non-default init

   auto *DefMesh = HorzMesh::getDefault();
   auto *DefHalo = Halo::getDefault();

   MetaDim::create("NCells", DefMesh->NCellsSize);
   MetaDim::create("NVertices", DefMesh->NVerticesSize);
   MetaDim::create("NEdges", DefMesh->NEdgesSize);
   MetaDim::create("NVertLevels", NVertLevels);

   const int NTimeLevels = 2;
   auto *TestOceanState  = OceanState::create("TestState", DefMesh, DefHalo,
                                              NVertLevels, NTimeLevels);
   if (!TestOceanState) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test state");
   }

   auto *TestAuxState =
       AuxiliaryState::create("TestAuxState", DefMesh, NVertLevels);
   if (!TestAuxState) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test auxiliary state");
   }

   Config Options;
   auto *TestTendencies =
       Tendencies::create("TestTendencies", DefMesh, NVertLevels, &Options,
                          DecayThicknessTendency{}, DecayVelocityTendency{});
   if (!TestTendencies) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test tendencies");
   }

   return Err;
}

void timeLoop(Real TimeEnd, Real TimeStep) {
   const auto *Stepper = TimeStepper::get("TestTimeStepper");
   auto *State         = OceanState::get("TestState");

   const int NSteps = std::ceil(TimeEnd / TimeStep);
   TimeStep         = TimeEnd / NSteps;

   for (int Step = 0; Step < NSteps; ++Step) {
      const Real Time = Step * TimeStep;
      Stepper->doStep(State, Time, TimeStep);
   }
}

void finalizeTimeStepperTest() {

   MetaDim::destroy("NCells");
   MetaDim::destroy("NVertices");
   MetaDim::destroy("NEdges");
   MetaDim::destroy("NVertLevels");

   TimeStepper::clear();
   Tendencies::clear();
   AuxiliaryState::clear();
   OceanState::clear();
   IOField::clear();
   HorzMesh::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int testTimeStepper(const std::string &Name, TimeStepperType Type,
                    Real ExpectedOrder, Real ATol) {
   int Err = 0;

   auto *DefMesh        = HorzMesh::getDefault();
   auto *DefHalo        = Halo::getDefault();
   auto *TestAuxState   = AuxiliaryState::get("TestAuxState");
   auto *TestTendencies = Tendencies::get("TestTendencies");

   auto *TestTimeStepper = TimeStepper::create(
       "TestTimeStepper", Type, TestTendencies, TestAuxState, DefMesh, DefHalo);

   if (!TestTimeStepper) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test time stepper {}", Name);
   }

   int NRefinements = 2;
   std::vector<ErrorMeasures> Errors(NRefinements);

   const Real TimeEnd      = 2;
   const Real BaseTimeStep = 0.2;

   const static bool CallOnlyOnce = [=]() {
      createExactSolution(TimeEnd);
      return true;
   }();

   Real TimeStep = BaseTimeStep;
   for (int RefLevel = 0; RefLevel < NRefinements; ++RefLevel) {
      Err += initState();

      timeLoop(TimeEnd, TimeStep);

      Errors[RefLevel] = computeErrors();

      TimeStep /= 2;
   }

   std::vector<Real> ConvRates(NRefinements - 1);
   for (int RefLevel = 0; RefLevel < NRefinements - 1; ++RefLevel) {
      ConvRates[RefLevel] =
          std::log2(Errors[RefLevel].LInf / Errors[RefLevel + 1].LInf);
   }

   if (std::abs(ConvRates.back() - ExpectedOrder) > ATol) {
      Err++;
      LOG_ERROR(
          "Wrong convergence rate for time stepper {}, got {:.3f}, expected {}",
          Name, ConvRates.back(), ExpectedOrder);
   }

   TimeStepper::erase("TestTimeStepper");

   return Err;
}

int timeStepperTest(const std::string &MeshFile = "OmegaMesh.nc") {

   int Err = initTimeStepperTest(MeshFile);

   if (Err != 0) {
      LOG_CRITICAL("TimeStepperTest: Error initializing");
   }

   Real ExpectedOrder = 4;
   Real ATol          = 0.1;
   Err += testTimeStepper("RungeKutta4", TimeStepperType::RungeKutta4,
                          ExpectedOrder, ATol);

   ExpectedOrder = 1;
   ATol          = 0.1;
   Err += testTimeStepper("ForwardBackward", TimeStepperType::ForwardBackward,
                          ExpectedOrder, ATol);

   if (Err == 0) {
      LOG_INFO("TimeStepperTest: Successful completion");
   }

   finalizeTimeStepperTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);

   RetVal += timeStepperTest();

   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
