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
#include "Tracers.h"
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

      auto *Mesh       = HorzMesh::getDefault();
      auto NVertLevels = NormalVelTend.extent_int(1);
      Array2DReal NormalVelEdge;
      State->getNormalVelocity(NormalVelEdge, VelTimeLevel);

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
   Array3DReal TracerArray;
   Err = Tracers::getAll(TracerArray, 0);

   Array2DReal LayerThickCell;
   Array2DReal NormalVelEdge;
   State->getLayerThickness(LayerThickCell, 0);
   State->getNormalVelocity(NormalVelEdge, 0);

   // Initially set thickness and velocity and tracers to 1
   deepCopy(LayerThickCell, 1);
   deepCopy(NormalVelEdge, 1);
   deepCopy(TracerArray, 1);

   return Err;
}

int createExactSolution(Real TimeEnd) {
   int Err = 0;

   auto *DefHalo = Halo::getDefault();
   auto *DefMesh = HorzMesh::getDefault();
   Array3DReal TracerArray;
   Err = Tracers::getAll(TracerArray, 0);

   auto *ExactState =
       OceanState::create("Exact", DefMesh, DefHalo, NVertLevels, 1);

   Array2DReal LayerThickCell;
   Array2DReal NormalVelEdge;
   ExactState->getLayerThickness(LayerThickCell, 0);
   ExactState->getNormalVelocity(NormalVelEdge, 0);

   // There are no thickness tendencies in this test, so exact thickness ==
   // initial thickness
   deepCopy(LayerThickCell, 1);
   // Normal velocity decays exponentially
   deepCopy(NormalVelEdge, DecayVelocityTendency{}.exactSolution(TimeEnd));
   // No tracer tendenciesk, final tracers == initial tracers
   deepCopy(TracerArray, 1);

   return Err;
}

ErrorMeasures computeErrors() {
   const auto *DefMesh = HorzMesh::getDefault();

   const auto *State      = OceanState::get("TestState");
   const auto *ExactState = OceanState::get("Exact");

   Array2DReal NormalVelEdge;
   Array2DReal ExactNormalVelEdge;
   State->getNormalVelocity(NormalVelEdge, 0);
   ExactState->getNormalVelocity(ExactNormalVelEdge, 0);

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
   Config("Omega");
   Err = Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("TimeStepperTest: Error reading config file");
      return Err;
   }

   // Reset NVertLevels to 1 regardless of config value
   Config *OmegaConfig = Config::getOmegaConfig();
   Config DimConfig("Dimension");
   Err = OmegaConfig->get(DimConfig);
   if (Err != 0) {
      LOG_CRITICAL("TimeStepperTest: Dimension group not found in Config");
      return Err;
   }
   Err = DimConfig.set("NVertLevels", NVertLevels);
   if (Err != 0) {
      LOG_CRITICAL("TimeStepperTest: Unable to reset NVertLevels in Config");
      return Err;
   }

   // Horz dimensions will be created in HorzMesh
   auto VertDim = Dimension::create("NVertLevels", NVertLevels);

   // Note that the default time stepper is not used in subsequent tests
   // but is initialized here because the number of time levels is needed
   // to initialize the Tracers. If a later timestepper test uses more time
   // levels than the default, this unit test may fail.
   int TSErr = TimeStepper::init1();
   if (TSErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing default time stepper");
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

   int TracerErr = Tracers::init();
   if (TracerErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing tracers infrastructure");
   }

   int AuxStateErr = AuxiliaryState::init();
   if (AuxStateErr != 0) {
      Err++;
      LOG_ERROR("TimeStepperTest: error initializing default aux state");
   }

   Err = Tendencies::init();
   if (Err != 0) {
      LOG_CRITICAL("Error initializing default tendencies");
      return Err;
   }

   // finish initializing default time stepper
   TSErr = TimeStepper::init2();
   if (TSErr != 0) {
      Err++;
      LOG_ERROR("error initializing default time stepper");
   }

   // Default time stepper never used and time stepper tests require
   // the no calendar option, so we reset calendar here
   Calendar::reset();
   Calendar::init("No Calendar");

   // Non-default init
   // Creating non-default state and auxiliary state to use only one vertical
   // level

   auto *DefMesh = HorzMesh::getDefault();
   auto *DefHalo = Halo::getDefault();

   int NTracers          = Tracers::getNumTracers();
   const int NTimeLevels = 2;
   auto *TestOceanState  = OceanState::create("TestState", DefMesh, DefHalo,
                                              NVertLevels, NTimeLevels);
   if (!TestOceanState) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test state");
   }

   auto *TestAuxState =
       AuxiliaryState::create("TestAuxState", DefMesh, NVertLevels, NTracers);
   if (!TestAuxState) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test auxiliary state");
   }

   Config Options;

   // Creating non-default tendencies with custom velocity tendencies
   auto *TestTendencies = Tendencies::create(
       "TestTendencies", DefMesh, NVertLevels, NTracers, &Options,
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
   TestTendencies->TracerHorzAdv.Enabled      = false;
   TestTendencies->TracerDiffusion.Enabled    = false;
   TestTendencies->TracerHyperDiff.Enabled    = false;

   return Err;
}

// Slightly adjust time step so that it evenly divides TimeEnd and return number
// of steps
int adjustTimeStep(TimeStepper *Stepper, Real TimeEnd) {
   TimeInterval TimeStep = Stepper->getTimeStep();
   R8 TimeStepSeconds;
   TimeStep.get(TimeStepSeconds, TimeUnits::Seconds);

   const int NSteps = std::ceil(TimeEnd / TimeStepSeconds);

   TimeStepSeconds = TimeEnd / NSteps;
   TimeStep.set(TimeStepSeconds, TimeUnits::Seconds);
   Stepper->changeTimeStep(TimeInterval(TimeStepSeconds, TimeUnits::Seconds));

   return NSteps;
}

void timeLoop(TimeInstant TimeStart, Real TimeEnd) {
   auto *Stepper = TimeStepper::get("TestTimeStepper");
   auto *State   = OceanState::get("TestState");

   const int NSteps            = adjustTimeStep(Stepper, TimeEnd);
   const TimeInterval TimeStep = Stepper->getTimeStep();
   TimeInstant CurTime         = Stepper->getStartTime();

   // Time loop
   for (int Step = 0; Step < NSteps; ++Step) {
      Stepper->doStep(State, CurTime);
   }
}

void finalizeTimeStepperTest() {

   Tracers::clear();
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

   // Set pointers to data
   auto *DefMesh        = HorzMesh::getDefault();
   auto *DefHalo        = Halo::getDefault();
   auto *TestAuxState   = AuxiliaryState::get("TestAuxState");
   auto *TestTendencies = Tendencies::get("TestTendencies");

   // Set time information
   const TimeInstant TimeStart(0, 0, 0, 0, 0, 0);

   const Real TimeEnd = 1;
   TimeInstant TimeEndTI(0, 0, 0, 0, 0, 1);

   const Real BaseTimeStepSeconds = 0.2;
   TimeInterval TimeStepTI(BaseTimeStepSeconds, TimeUnits::Seconds);

   auto *TestTimeStepper = TimeStepper::create(
       "TestTimeStepper", Type, TimeStart, TimeEndTI, TimeStepTI,
       TestTendencies, TestAuxState, DefMesh, DefHalo);

   if (!TestTimeStepper) {
      Err++;
      LOG_ERROR("TimeStepperTest: error creating test time stepper {}", Name);
   }

   int NRefinements = 2;
   std::vector<ErrorMeasures> Errors(NRefinements);

   // This creates global exact solution and needs to be done only once
   const static bool CallOnlyOnce = [=]() {
      createExactSolution(TimeEnd);
      return true;
   }();

   R8 TimeStepSeconds = BaseTimeStepSeconds;

   // Convergence loop
   for (int RefLevel = 0; RefLevel < NRefinements; ++RefLevel) {
      TestTimeStepper->changeTimeStep(
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

   LOG_INFO("----- Time Stepper Unit Test -----");

   RetVal += timeStepperTest();

   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
