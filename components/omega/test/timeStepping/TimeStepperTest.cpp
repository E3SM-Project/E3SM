//===-- Test driver for OMEGA time steppers -----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA time steppers
///
/// This driver tests the OMEGA time stepping module. It checks that every
/// implemented time stepping scheme converges to an exact solution at
/// the expected theoretical rate. To only measure time convergence,
/// the initial conditions and tendency terms are set up such that for
/// every mesh element an independent ODE is solved.
///
//
//===-----------------------------------------------------------------------===/

#include "TimeStepper.h"
#include "../ocn/OceanTestCommon.h"
#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "TendencyTerms.h"
#include "TimeMgr.h"
#include "mpi.h"

#include <cmath>
#include <iomanip>
#include <iostream>

using namespace OMEGA;

// Geometry doesn't matter for this test
constexpr Geometry Geom = Geometry::Planar;
// Only one vertical level is needed
constexpr int NVertLevels = 1;

// Custom tendency for normal velocity
// du/dt = -coeff * u
struct DecayVelocityTendency {
   Real Coeff = 0.5;

   // exact solution assumes that this is the only tendency active
   // the solution is exponential decay
   Real exactSolution(Real Time) { return std::exp(-Coeff * Time); }

   void operator()(Array2DReal NormalVelTend, const OceanState *State,
                   const AuxiliaryState *AuxState, int ThickTimeLevel,
                   int VelTimeLevel, TimeInstant Time) const {

      auto *Mesh                = HorzMesh::getDefault();
      auto NVertLevels          = NormalVelTend.extent_int(1);
      const auto &NormalVelEdge = State->NormalVelocity[VelTimeLevel];

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

   // Initially set thickness and velocity to 1
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

   // There are no thickness tendencies in this test, so exact thickness ==
   // initial thickness
   deepCopy(LayerThickCell, 1);
   // Normal velocity decays exponentially
   deepCopy(NormalVelEdge, DecayVelocityTendency{}.exactSolution(TimeEnd));

   return Err;
}

ErrorMeasures computeErrors() {
   const auto *DefMesh = HorzMesh::getDefault();

   const auto *State      = OceanState::get("TestState");
   const auto *ExactState = OceanState::get("Exact");

   const auto &NormalVelEdge      = State->NormalVelocity[0];
   const auto &ExactNormalVelEdge = ExactState->NormalVelocity[0];

   // Only velocity errors matters, because thickness remains constant
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

   // Default init

   initLogging(DefEnv);

   // Open config file
   OMEGA::Config("Omega");
   Err = OMEGA::Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("TimeStepperTest: Error reading config file");
      return Err;
   }

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

   // Non-default init
   // Creating non-default state and auxiliary state to use only one vertical
   // level

   auto *DefMesh = HorzMesh::getDefault();
   auto *DefHalo = Halo::getDefault();

   // Horz dimensions created in HorzMesh
   auto VertDim = Dimension::create("NVertLevels", NVertLevels);

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

   // Creating non-default tendencies with custom velocity tendencies
   auto *TestTendencies = Tendencies::create(
       "TestTendencies", DefMesh, NVertLevels, &Options,
       Tendencies::CustomTendencyType{}, DecayVelocityTendency{});
   if (!TestTendencies) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test tendencies");
   }

   // Disable all other tendencies
   TestTendencies->ThicknessFluxDiv.Enabled   = false;
   TestTendencies->PotientialVortHAdv.Enabled = false;
   TestTendencies->KEGrad.Enabled             = false;
   TestTendencies->SSHGrad.Enabled            = false;
   TestTendencies->VelocityDiffusion.Enabled  = false;
   TestTendencies->VelocityHyperDiff.Enabled  = false;

   return Err;
}

// Slightly adjust time step so that it evenly divides TimeEnd and return number
// of steps
int adjustTimeStep(TimeStepper *Stepper, Real TimeEnd) {
   TimeInterval TimeStep = Stepper->getTimeStep();
   Real TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);

   const int NSteps = std::ceil(TimeEnd / TimeStepSeconds);

   TimeStepSeconds = TimeEnd / NSteps;
   TimeStep.set(TimeStepSeconds, TimeUnits::Seconds);
   Stepper->setTimeStep(TimeInterval(TimeStepSeconds, TimeUnits::Seconds));

   return NSteps;
}

void timeLoop(TimeInstant TimeStart, Real TimeEnd) {
   auto *Stepper = TimeStepper::get("TestTimeStepper");
   auto *State   = OceanState::get("TestState");

   const int NSteps            = adjustTimeStep(Stepper, TimeEnd);
   const TimeInterval TimeStep = Stepper->getTimeStep();

   // Time loop
   for (int Step = 0; Step < NSteps; ++Step) {
      const TimeInstant Time = TimeStart + Step * TimeStep;
      Stepper->doStep(State, Time);
   }
}

void finalizeTimeStepperTest() {

   TimeStepper::clear();
   Tendencies::clear();
   AuxiliaryState::clear();
   OceanState::clear();
   Dimension::clear();
   Field::clear();
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

   Calendar TestCalendar("TestCalendar", CalendarNoCalendar);

   auto *TestTimeStepper = TimeStepper::create(
       "TestTimeStepper", Type, TestTendencies, TestAuxState, DefMesh, DefHalo);

   if (!TestTimeStepper) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test time stepper {}", Name);
   }

   int NRefinements = 2;
   std::vector<ErrorMeasures> Errors(NRefinements);

   // Start time = 0
   const TimeInstant TimeStart(&TestCalendar, 0, 0, 0, 0, 0, 0);

   const Real TimeEnd = 1;

   const Real BaseTimeStepSeconds = 0.2;

   // This creates global exact solution and needs to be done only once
   const static bool CallOnlyOnce = [=]() {
      createExactSolution(TimeEnd);
      return true;
   }();

   Real TimeStepSeconds = BaseTimeStepSeconds;

   // Convergence loop
   for (int RefLevel = 0; RefLevel < NRefinements; ++RefLevel) {
      TestTimeStepper->setTimeStep(
          TimeInterval(TimeStepSeconds, TimeUnits::Seconds));

      Err += initState();

      timeLoop(TimeStart, TimeEnd);

      Errors[RefLevel] = computeErrors();

      TimeStepSeconds /= 2;
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

   // Test convergence rate of different time steppers

   Real ExpectedOrder = 4;
   Real ATol          = 0.1;
   Err += testTimeStepper("RungeKutta4", TimeStepperType::RungeKutta4,
                          ExpectedOrder, ATol);

   ExpectedOrder = 1;
   ATol          = 0.1;
   Err += testTimeStepper("ForwardBackward", TimeStepperType::ForwardBackward,
                          ExpectedOrder, ATol);

   ExpectedOrder = 2;
   ATol          = 0.1;
   Err += testTimeStepper("RungeKutta2", TimeStepperType::RungeKutta2,
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
