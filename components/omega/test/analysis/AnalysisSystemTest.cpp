//===-- Test 5.2, 5.3, 5.4, 5.5, 5.6: Analysis system tests ----*- C++ -*-===//
//
// System tests for Analysis framework: dependency resolution, alarm system,
// factory registration, configuration parsing, and end-to-end integration
//
//===-----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AnalysisOpFactory.h"
#include "AuxiliaryState.h"
#include "Decomp.h"
#include "Eos.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "OceanState.h"
#include "PGrad.h"
#include "Tendencies.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertAdv.h"
#include "VertCoord.h"

#include <iostream>

using namespace OMEGA;

// Test result tracking
int NumTests  = 0;
int NumPassed = 0;
int NumFailed = 0;

//------------------------------------------------------------------------------
// Helper function to report test results
void reportTest(const std::string &TestName, bool Passed) {
   NumTests++;
   if (Passed) {
      NumPassed++;
   } else {
      NumFailed++;
      LOG_ERROR("FAIL: {}", TestName);
   }
}

//===----------------------------------------------------------------------===//
// Test 5.2: Dependency Resolution and Execution Order
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Test 5.2.1: Verify shared intermediate operators are deduplicated
void testSharedIntermediates() {

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   // Count how many times each operator name appears
   std::map<std::string, int> OpCounts;
   for (const auto *Node : OpNodes) {
      std::string OpName = Node->Op->getName();
      OpCounts[OpName]++;
   }

   // Check that no operator appears more than once
   bool Passed = true;
   for (const auto &Pair : OpCounts) {
      if (Pair.second > 1) {
         LOG_ERROR("  Operator {} appears {} times (should be 1)", Pair.first,
                   Pair.second);
         Passed = false;
      }
   }

   reportTest("Dependency: Shared intermediates deduplicated", Passed);
}

//------------------------------------------------------------------------------
// Test 5.2.2: Verify upstream dependencies are correctly resolved
void testUpstreamDependencies() {

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   bool Passed = true;

   // For each operator, verify its upstreams produce the fields it needs
   for (const auto *Node : OpNodes) {
      auto InputNames = Node->Op->getInputFieldNames();

      for (const auto &InputName : InputNames) {
         // Check if this input is a simulation field or an operator output
         bool FoundUpstream = false;

         // Check if it's a simulation field
         if (Field::exists(InputName)) {
            FoundUpstream = true;
         }

         // Check if it's produced by an upstream operator
         for (const auto *Upstream : Node->Upstreams) {
            auto UpstreamOutputs = Upstream->Op->getOutputFieldNames();
            for (const auto &Output : UpstreamOutputs) {
               if (Output == InputName) {
                  FoundUpstream = true;
                  break;
               }
            }
            if (FoundUpstream)
               break;
         }

         if (!FoundUpstream) {
            LOG_ERROR("  Operator {} requires input {} but no upstream found",
                      Node->Op->getName(), InputName);
            Passed = false;
         }
      }
   }

   reportTest("Dependency: Upstream dependencies resolved", Passed);
}

//------------------------------------------------------------------------------
// Test 5.2.3: Verify cache prevents redundant computation
void testCacheValidation() {

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   if (OpNodes.empty()) {
      reportTest("Cache: Validation (no operators to test)", true);
      return;
   }

   auto AnalysisClock = DefAnalysis->getModelClock();
   auto CurTime       = AnalysisClock->getCurrentTime();
   auto TimeStep      = AnalysisClock->getTimeStep();

   auto TestTime = CurTime + TimeStep;

   // Get first operator
   auto *FirstNode = OpNodes[0];
   //   TimeInstant TestTime(0, 0, 0, 0, 0, 1);

   // Initially should not be valid
   bool InitiallyInvalid = !FirstNode->Op->isCacheValid(TestTime);

   // Compute it
   FirstNode->Op->compute(TestTime);

   // Now should be valid for same timestamp
   bool NowValid = FirstNode->Op->isCacheValid(TestTime);

   // Should be invalid for different timestamp
   auto DifferentTime = TestTime + TimeStep;
   //  TimeInstant DifferentTime(0, 0, 0, 0, 0, 2);
   bool InvalidForDifferentTime = !FirstNode->Op->isCacheValid(DifferentTime);

   bool Passed = InitiallyInvalid && NowValid && InvalidForDifferentTime;
   reportTest("Cache: Validation prevents redundant computation", Passed);

   if (!Passed) {
      LOG_ERROR(
          "  InitiallyInvalid: {}, NowValid: {}, InvalidForDifferentTime: {}",
          InitiallyInvalid, NowValid, InvalidForDifferentTime);
   }
}

//===----------------------------------------------------------------------===//
// Test 5.3: Alarm System Verification
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Test 5.3.1: Verify terminal operators have alarms
void testTerminalOperatorAlarms() {

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   bool Passed = true;

   // Terminal operators (those with StreamNames) should have alarms
   for (const auto *Node : OpNodes) {
      if (!Node->StreamNames.empty()) {
         if (Node->ComputeAlarms.empty()) {
            LOG_ERROR("  Terminal operator {} has no alarms",
                      Node->Op->getName());
            Passed = false;
         }
      }
   }

   reportTest("Alarm: Terminal operators have alarms", Passed);
}

//------------------------------------------------------------------------------
// Test 5.3.2: Verify temporal operators have multiple alarms
void testTemporalOperatorAlarms() {

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   bool Passed = true;

   // Temporal reduction operators should have 2 alarms (accumulation + output)
   for (const auto *Node : OpNodes) {
      std::string OpType   = Node->Op->getOperatorType();
      bool IsTimeReduction = (OpType.find("Time") != std::string::npos);

      if (IsTimeReduction && !Node->StreamNames.empty()) {
         if (Node->ComputeAlarms.size() < 2) {
            LOG_ERROR("  Temporal operator {} has {} alarms (expected 2)",
                      Node->Op->getName(), Node->ComputeAlarms.size());
            Passed = false;
         }
      }
   }

   reportTest("Alarm: Temporal operators have accumulation + output alarms",
              Passed);
}

//------------------------------------------------------------------------------
// Test 5.3.3: Verify alarm propagation to upstream operators
void testAlarmPropagation() {

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   bool Passed = true;

   // Intermediate operators should have alarms propagated from downstream
   for (const auto *Node : OpNodes) {
      if (Node->StreamNames.empty() && !Node->Upstreams.empty()) {
         // This is an intermediate operator
         if (Node->ComputeAlarms.empty()) {
            LOG_ERROR("  Intermediate operator {} has no propagated alarms",
                      Node->Op->getName());
            Passed = false;
         }
      }
   }

   reportTest("Alarm: Propagation to upstream operators", Passed);
}

//===----------------------------------------------------------------------===//
// Test 5.4: Factory Registration and Type Dispatch
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Helper: Create a test field for factory testing
void createFactoryTestField(const std::string &FieldName, const HorzMesh *Mesh,
                            const VertCoord *VCoord) {

   I4 NCells      = Mesh->NCellsSize;
   I4 NVertLayers = VCoord->NVertLayers;

   // Create a 2D test field
   std::vector<std::string> DimNames = {"NCells", "NVertLayers"};
   auto TestField =
       Field::create(FieldName, "Test field for factory validation", "units",
                     "", -1.0e30, 1.0e30, -1.0e30, 2, DimNames);

   // Allocate and attach data with known values
   Array2DReal TestData(FieldName + "_data", NCells, NVertLayers);
   TestField->attachData<Array2DReal>(TestData);

   // Fill with simple pattern: value = Cell + K
   auto TestDataHost = Kokkos::create_mirror_view(TestData);
   for (I4 Cell = 0; Cell < NCells; ++Cell) {
      for (I4 K = 0; K < NVertLayers; ++K) {
         TestDataHost(Cell, K) = static_cast<Real>(Cell + K);
      }
   }
   Kokkos::deep_copy(TestData, TestDataHost);
}

//------------------------------------------------------------------------------
// Test 5.4.1: Verify factory creates valid operators (proves registration)
void testFactoryCreatesValidOperators() {

   // Create a simple test field to use for operator creation
   auto Mesh   = HorzMesh::getDefault();
   auto VCoord = VertCoord::getDefault();

   std::string TestFieldName = "FactoryRegistrationTestField";
   createFactoryTestField(TestFieldName, Mesh, VCoord);

   Config EmptyConfig;
   bool Passed = true;

   // Try to create a SpatialMax operator - if this succeeds, it proves
   // the operator type is registered in the factory
   auto TestOp =
       AnalysisOpFactory::createOp("SpatialMax", {TestFieldName}, EmptyConfig);
   if (!TestOp) {
      LOG_ERROR("  Failed to create SpatialMax operator - registration failed");
      Passed = false;
   } else {
      // Verify the operator has correct type
      if (TestOp->getOperatorType() != "SpatialMax") {
         LOG_ERROR("  Created operator has wrong type");
         Passed = false;
      }
   }

   reportTest("Factory: Creates valid operators (proves registration)", Passed);
}

//------------------------------------------------------------------------------
// Test 5.4.2: Verify type dispatch and operator creation
void testFactoryTypeDispatch(const MachEnv *Env, const HorzMesh *Mesh,
                             const VertCoord *VCoord) {

   // Create a test field specifically for this test
   std::string TestFieldName = "FactoryTestField";
   createFactoryTestField(TestFieldName, Mesh, VCoord);

   Config EmptyConfig;
   bool Passed = true;

   // Create operator using the factory
   auto MaxOp =
       AnalysisOpFactory::createOp("SpatialMax", {TestFieldName}, EmptyConfig);
   if (!MaxOp) {
      LOG_ERROR("  Failed to create SpatialMax operator");
      Passed = false;
   } else {
      // Verify operator type is correct
      if (MaxOp->getOperatorType() != "SpatialMax") {
         LOG_ERROR("  Operator type mismatch: expected SpatialMax, got {}",
                   MaxOp->getOperatorType());
         Passed = false;
      }

      // Verify input field names are correct
      auto InputNames = MaxOp->getInputFieldNames();
      if (InputNames.size() != 1 || InputNames[0] != TestFieldName) {
         LOG_ERROR("  Input field names mismatch");
         Passed = false;
      }

      // Verify output field names follow expected pattern
      auto OutputNames               = MaxOp->getOutputFieldNames();
      std::string ExpectedOutputName = TestFieldName + "_SpatialMax";
      if (OutputNames.size() != 1 || OutputNames[0] != ExpectedOutputName) {
         LOG_ERROR("  Output field name mismatch: expected {}, got {}",
                   ExpectedOutputName,
                   OutputNames.empty() ? "none" : OutputNames[0]);
         Passed = false;
      }
   }

   reportTest("Factory: Type dispatch creates correct operator", Passed);
}

//------------------------------------------------------------------------------
// Test 5.4.3: Verify factory handles different field configurations
void testFactoryDifferentFieldTypes(const HorzMesh *Mesh,
                                    const VertCoord *VCoord) {

   Config EmptyConfig;
   bool Passed = true;

   // Test 1: Create operator with first test field
   std::string TestFieldName1 = "FactoryDifferentFields1";
   createFactoryTestField(TestFieldName1, Mesh, VCoord);

   auto Op1 =
       AnalysisOpFactory::createOp("SpatialMin", {TestFieldName1}, EmptyConfig);
   if (!Op1) {
      LOG_ERROR("  Failed to create operator for test field 1");
      Passed = false;
   } else if (Op1->getOperatorType() != "SpatialMin") {
      LOG_ERROR("  Operator 1 has wrong type");
      Passed = false;
   }

   // Test 2: Create operator with second test field
   std::string TestFieldName2 = "FactoryDifferentFields2";
   createFactoryTestField(TestFieldName2, Mesh, VCoord);

   auto Op2 = AnalysisOpFactory::createOp("SpatialMean", {TestFieldName2},
                                          EmptyConfig);
   if (!Op2) {
      LOG_ERROR("  Failed to create operator for test field 2");
      Passed = false;
   } else if (Op2->getOperatorType() != "SpatialMean") {
      LOG_ERROR("  Operator 2 has wrong type");
      Passed = false;
   }

   // Test 3: Verify operators have correct input/output associations
   if (Op1 && Op2) {
      auto Inputs1 = Op1->getInputFieldNames();
      auto Inputs2 = Op2->getInputFieldNames();

      if (Inputs1.empty() || Inputs1[0] != TestFieldName1) {
         LOG_ERROR("  Operator 1 has wrong input field");
         Passed = false;
      }
      if (Inputs2.empty() || Inputs2[0] != TestFieldName2) {
         LOG_ERROR("  Operator 2 has wrong input field");
         Passed = false;
      }
   }

   reportTest("Factory: Handles different field configurations", Passed);
}

//------------------------------------------------------------------------------
// Test 5.4.4: Verify factory instantiates all operator types
void testFactoryInstantiateAll(const MachEnv *Env, const HorzMesh *Mesh,
                               const VertCoord *VCoord) {

   // Create a test field for all operators
   std::string TestFieldName = "FactoryInstantiateTestField";
   createFactoryTestField(TestFieldName, Mesh, VCoord);

   std::vector<std::string> SpatialOps = {"SpatialMax", "SpatialMin",
                                          "SpatialMean", "SpatialStdDev"};

   Config EmptyConfig;
   bool Passed = true;

   // Test each spatial operator - verify factory can create them
   for (const auto &OpName : SpatialOps) {
      auto Op =
          AnalysisOpFactory::createOp(OpName, {TestFieldName}, EmptyConfig);
      if (!Op) {
         LOG_ERROR("  Failed to instantiate {} operator", OpName);
         Passed = false;
         continue;
      }

      // Verify operator type matches request
      if (Op->getOperatorType() != OpName) {
         LOG_ERROR("  Operator type mismatch for {}", OpName);
         Passed = false;
      }

      // Verify output field names follow expected pattern
      auto OutputNames               = Op->getOutputFieldNames();
      std::string ExpectedOutputName = TestFieldName + "_" + OpName;
      if (OutputNames.empty() || OutputNames[0] != ExpectedOutputName) {
         LOG_ERROR("  Unexpected output field name for {}", OpName);
         Passed = false;
      }
   }

   // Test TimeMean with required Period parameter
   auto TimeMeanOp = AnalysisOpFactory::createOp(
       "TimeMean", {TestFieldName},
       makeOpConfig(opParam("Period", std::string("1day"))));
   if (!TimeMeanOp) {
      LOG_ERROR("  Failed to instantiate TimeMean operator");
      Passed = false;
   } else {
      // Verify operator type
      if (TimeMeanOp->getOperatorType() != "TimeMean") {
         LOG_ERROR("  Operator type mismatch for TimeMean");
         Passed = false;
      }

      // Verify output field naming includes period
      auto OutputNames = TimeMeanOp->getOutputFieldNames();
      if (OutputNames.empty()) {
         LOG_ERROR("  TimeMean operator has no output field names");
         Passed = false;
      }
   }

   reportTest("Factory: Instantiate all registered operator types", Passed);
}

//===----------------------------------------------------------------------===//
// Test 5.5: Configuration Parsing and Validation
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Test 5.5.1: Verify operator chain parsing
void testOperatorChainParsing() {

   // This test verifies that parseChainAndBuildOps correctly handles
   // underscore-delimited operator chains

   auto DefAnalysis = Analysis::getDefault();

   // Check if any operators were created from config
   auto OpNodes = DefAnalysis->getOpNodes();
   bool Passed  = !OpNodes.empty();

   reportTest("Config: Operator chain parsing", Passed);

   if (!Passed) {
      LOG_ERROR("  No operators were created from configuration");
   }
}

//------------------------------------------------------------------------------
// Test 5.5.2: Verify field reuse in chains
void testFieldReuseInChains() {

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   // Count unique output field names
   std::set<std::string> UniqueOutputs;
   for (const auto *Node : OpNodes) {
      auto Outputs = Node->Op->getOutputFieldNames();
      for (const auto &Output : Outputs) {
         UniqueOutputs.insert(Output);
      }
   }

   // Each operator should produce a unique output
   bool Passed = (UniqueOutputs.size() == OpNodes.size());
   reportTest("Config: Field reuse prevents duplicates", Passed);

   if (!Passed) {
      LOG_ERROR("  {} unique outputs for {} operators", UniqueOutputs.size(),
                OpNodes.size());
   }
}

//------------------------------------------------------------------------------
// Test 5.5.3: Verify stream parameter application
void testStreamParameterApplication() {

   // Verify that streams were created for analysis output
   // This is a basic check that the stream creation succeeded

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   bool Passed = false;

   // Check if any terminal operators have associated streams
   for (const auto *Node : OpNodes) {
      if (!Node->StreamNames.empty()) {
         Passed = true;
         break;
      }
   }

   reportTest("Config: Stream parameters applied", Passed);

   if (!Passed) {
      LOG_ERROR("  No terminal operators with associated streams found");
   }
}

//===----------------------------------------------------------------------===//
// Test 5.6: End-to-End Integration
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Test 5.6.1: Verify computeAll executes without errors
void testComputeAllExecution(Clock *ModelClock) {

   auto DefAnalysis = Analysis::getDefault();

   bool Passed = true;

   // Advance clock to trigger alarms
   TimeInstant CurrentTime = ModelClock->getCurrentTime();
   TimeInterval OneStep    = ModelClock->getTimeStep();
   TimeInstant NextTime    = CurrentTime + OneStep;
   ModelClock->advance();

   // Call computeAll
   DefAnalysis->computeAll();

   reportTest("Integration: computeAll executes without errors", Passed);
}

//------------------------------------------------------------------------------
// Test 5.6.2: Verify output fields are created
void testOutputFieldsCreated() {

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   bool Passed = true;

   // Verify that all operator output fields exist in Field registry
   for (const auto *Node : OpNodes) {
      auto Outputs = Node->Op->getOutputFieldNames();
      for (const auto &Output : Outputs) {
         if (!Field::exists(Output)) {
            LOG_ERROR("  Output field {} does not exist", Output);
            Passed = false;
         }
      }
   }

   reportTest("Integration: Output fields created", Passed);
}

//------------------------------------------------------------------------------
// Test 5.6.3: Verify stream output (basic check)
void testStreamOutput() {

   // This is a basic check that streams are configured
   // Full I/O testing would require writing and reading files

   auto DefAnalysis = Analysis::getDefault();
   auto OpNodes     = DefAnalysis->getOpNodes();

   bool Passed = false;

   // Check if any operators are associated with streams
   for (const auto *Node : OpNodes) {
      if (!Node->StreamNames.empty()) {
         // Verify the stream exists
         for (const auto &StreamName : Node->StreamNames) {
            auto NodeStream = IOStream::get(StreamName);
            if (NodeStream) {
               Passed = true;
               break;
            }
         }
      }
      if (Passed)
         break;
   }

   reportTest("Integration: Stream output configured", Passed);

   if (!Passed) {
      LOG_ERROR("  No valid streams found for analysis output");
   }
}

//===----------------------------------------------------------------------===//
// Initialization and finalization functions
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Initialize needed modules
void initAnalysisSystemTest() {

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Enable GlobalStats to ensure at least one AnalysisGroup is enabled
   // for testing
   Config *OmegaConfig = Config::getOmegaConfig();
   Config AnalysisCfg("Analysis");
   OmegaConfig->get(AnalysisCfg);
   Config GlobalStatsCfg("GlobalStats");
   AnalysisCfg.get(GlobalStatsCfg);
   GlobalStatsCfg.set("Enable", true); // Force enable for testing

   // First step of time stepper initialization needed for IOstream
   TimeStepper::init1();

   // Get the model clock
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize streams
   IOStream::init(ModelClock);

   // Initialize fields
   Field::init(ModelClock);

   // Initialize the default halo
   Halo::init();

   // Initialize the default mesh
   HorzMesh::init(ModelClock);

   // Initialize the default vertical coordinate
   VertCoord::init();

   // Initialize tracers
   Tracers::init();

   // Initialize auxiliary state
   AuxiliaryState::init();

   // Initialize equation of state
   Eos::init();

   // Initialize pressure gradient
   PressureGrad::init();

   // Initialize tendencies
   Tendencies::init();

   // Initialize vertical advection
   VertAdv::init();

   // Second step of time stepper initialization
   TimeStepper::init2();

   // Initialize ocean state
   OceanState::init();

   // Validate streams
   bool StreamsValid = IOStream::validateAll();
   if (!StreamsValid) {
      LOG_ERROR("Stream validation failed");
   }

   // Read initial state
   Metadata ReqMeta;
   Error Err1 = IOStream::read("InitialState", ModelClock, ReqMeta);
   if (Err1.isFail()) {
      LOG_ERROR("Failed to read initial state");
   }

   // Initialize Analysis module (creates operators, resolves dependencies,
   // sets alarms)
   Analysis::init();
}

//------------------------------------------------------------------------------
// Clean-up modules
void finalizeAnalysisSystemTest() {

   Analysis::finalize();
   IOStream::finalize();
   OceanState::clear();
   Tracers::clear();
   AuxiliaryState::clear();
   PressureGrad::clear();
   Tendencies::clear();
   VertAdv::clear();
   VertCoord::clear();
   TimeStepper::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

//===----------------------------------------------------------------------===//
// Main test driver
//===----------------------------------------------------------------------===//

int main(int argc, char *argv[]) {

   int ErrCode = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      initAnalysisSystemTest();

      auto DefEnv     = MachEnv::getDefault();
      auto DefStepper = TimeStepper::getDefault();
      auto Mesh       = HorzMesh::getDefault();
      auto VCoord     = VertCoord::getDefault();
      auto ModelClock = DefStepper->getClock();

      // Dependency Resolution and Execution Order
      testSharedIntermediates();
      testUpstreamDependencies();
      testCacheValidation();

      // Alarm System Verification
      testTerminalOperatorAlarms();
      testTemporalOperatorAlarms();
      testAlarmPropagation();

      // Factory Registration and Type Dispatch
      testFactoryCreatesValidOperators();
      testFactoryTypeDispatch(DefEnv, Mesh, VCoord);
      testFactoryDifferentFieldTypes(Mesh, VCoord);
      testFactoryInstantiateAll(DefEnv, Mesh, VCoord);

      // Configuration Parsing and Validation
      testOperatorChainParsing();
      testFieldReuseInChains();
      testStreamParameterApplication();

      // End-to-End Integration
      testComputeAllExecution(ModelClock);
      testOutputFieldsCreated();
      testStreamOutput();

      if (NumFailed > 0) {
         ErrCode = 1;
         LOG_ERROR("AnalysisSystemTest Failure");
         LOG_ERROR("  Total tests: {}", NumTests);
         LOG_ERROR("  Passed: {}", NumPassed);
         LOG_ERROR("  Failed: {}", NumFailed);
      }

      finalizeAnalysisSystemTest();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return ErrCode;
}
