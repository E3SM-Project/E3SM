#include "SfcCoupling.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Eos.h"
#include "Field.h"
#include "Forcing.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "VertCoord.h"
#include "mpi.h"

using namespace OMEGA;

struct TestSetup {

   std::map<std::string, int> ImportIdx = {{"Foxx_taux", 3}, {"Foxx_tauy", 8}};
   std::map<std::string, int> ExportIdx = {
       {"So_t", 2}, {"So_s", 4}, {"So_u", 6}, {"So_v", 9}, {"So_ssh", 1}};
};

CouplingInitParams mockCouplingInitParams(
    CouplingLayout Layout                        = CouplingLayout::MCT,
    std::optional<TimeInterval> CouplingTimeStep = std::nullopt) {

   TestSetup Setup;

   auto *DefTimeStepper = TimeStepper::getDefault();
   TimeInterval CouplingTimeStep_ =
       CouplingTimeStep.value_or(DefTimeStepper->getTimeStep());

   CouplingInitParams CouplingParams{.NImportFields    = 10,
                                     .NExportFields    = 10,
                                     .ImportIdx        = Setup.ImportIdx,
                                     .ExportIdx        = Setup.ExportIdx,
                                     .CouplingTimeStep = CouplingTimeStep_,
                                     .Layout           = Layout};

   return CouplingParams;
}

std::string toString(const CouplingLayout &Layout) {
   switch (Layout) {
   case CouplingLayout::MCT:
      return "MCT";
   case CouplingLayout::MOAB:
      return "MOAB";
   default:
      return "Unknown";
   }
}

HostArray1DReal makeCellVarryingArray(const std::string &Name, const int NCells,
                                      const Real Base) {
   HostArray1DReal Array(Name, NCells);
   for (int Cell = 0; Cell < NCells; Cell++) {
      Array(Cell) = Base + Cell;
   }
   return Array;
}

// Shaed index formula for packing/unpacking from rae coupler buffer
int flatIdx(const CouplingLayout Layout, const int Cell, const int Field,
            const int NCells, const int NFields) {
   if (Layout == CouplingLayout::MCT)
      return Cell * NFields + Field;
   else // MOAB
      return Field * NCells + Cell;
}

int initSfcCouplingTest(const std::string &MeshFile) {

   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("------ Surface Coupling unit tests ------");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Initialize time stepper and get model clock
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);
   Decomp::init(MeshFile);

   // Initialize streams
   Field::init(ModelClock);
   IOStream::init(ModelClock);

   // TODO: Need to initialize halo? SfcCoupling has no lateral exchanges
   int HaloErr = Halo::init();
   if (HaloErr != 0) {
      Err++;
      LOG_ERROR("SfcCouplingTest: Error initializing defualt halo");
   }

   HorzMesh::init(ModelClock);

   VertCoord::init();

   int StateErr = OceanState::init();
   if (StateErr != 0) {
      Err++;
      LOG_ERROR("SfcCouplingTest: Error initializing defualt ocean state");
   }

   Eos::init();
   Forcing::init();
   Tracers::init();

   return Err;
}

int testImportFromCoupler(const CouplingLayout Layout) {

   int Err = 0;

   auto CouplingParams = mockCouplingInitParams(Layout);
   Err += SfcCoupling::init(CouplingParams);

   SfcCoupling *DefCoupling = SfcCoupling::getDefault();

   int NCells   = DefCoupling->NCellsOwned;
   int NImports = DefCoupling->NImportFields;
   int NExports = DefCoupling->NExportFields;

   int TauxIdx = CouplingParams.ImportIdx.at("Foxx_taux");
   int TauyIdx = CouplingParams.ImportIdx.at("Foxx_tauy");

   std::vector<Real> CplToOcnData(NCells * NImports, 0.0);
   std::vector<Real> OcnToCplData(NCells * NExports, 0.0);

   HostArray1DReal ExpectedSfcStressZonal =
       makeCellVarryingArray("ExpectedSfcStressZonal", NCells, Real(TauxIdx));
   HostArray1DReal ExpectedSfcStressMerid =
       makeCellVarryingArray("ExpectedSfcStressMerid", NCells, Real(TauyIdx));

   for (int Cell = 0; Cell < NCells; Cell++) {
      CplToOcnData[flatIdx(Layout, Cell, TauxIdx, NCells, NImports)] =
          ExpectedSfcStressZonal(Cell);
      CplToOcnData[flatIdx(Layout, Cell, TauyIdx, NCells, NImports)] =
          ExpectedSfcStressMerid(Cell);
   }

   DefCoupling->attachData(CplToOcnData.data(), OcnToCplData.data());
   DefCoupling->importFromCoupler();

   auto ImportPass = arraysEqual(DefCoupling->CplToOcn.SfcStressZonal,
                                 ExpectedSfcStressZonal) &&
                     arraysEqual(DefCoupling->CplToOcn.SfcStressMerid,
                                 ExpectedSfcStressMerid);

   if (ImportPass) {
      LOG_INFO("SfcCouplingTest: importFromCoupler with {} layout PASS",
               toString(Layout));
   } else {
      Err++;
      LOG_ERROR("SfcCouplingTest: importFromCoupler with {} layout FAIL",
                toString(Layout));
   }

   SfcCoupling::clear();

   return Err;
}

int testApplyImportFields() {

   int Err = 0;

   auto CouplingParams = mockCouplingInitParams();
   Err += SfcCoupling::init(CouplingParams);

   Forcing *DefForcing      = Forcing::getDefault();
   SfcCoupling *DefCoupling = SfcCoupling::getDefault();

   int NCells   = DefCoupling->NCellsOwned;
   int NImports = DefCoupling->NImportFields;
   int NExports = DefCoupling->NExportFields;

   int TauxIdx = CouplingParams.ImportIdx.at("Foxx_taux");
   int TauyIdx = CouplingParams.ImportIdx.at("Foxx_tauy");

   std::vector<Real> CplToOcnData(NCells * NImports, 0.0);
   std::vector<Real> OcnToCplData(NCells * NExports, 0.0);

   int Offset                             = 27;
   HostArray1DReal ExpectedSfcStressZonal = makeCellVarryingArray(
       "ExpectedSfcStressZonal", NCells, Real(TauxIdx + Offset));
   HostArray1DReal ExpectedSfcStressMerid = makeCellVarryingArray(
       "ExpectedSfcStressMerid", NCells, Real(TauyIdx + Offset));

   // Copy the expected values into the CplToOcn fields directly
   deepCopy(DefCoupling->CplToOcn.SfcStressZonal, ExpectedSfcStressZonal);
   deepCopy(DefCoupling->CplToOcn.SfcStressMerid, ExpectedSfcStressMerid);

   DefCoupling->applyImportFields(DefForcing);

   auto SfcStressZonalOwned = Kokkos::subview(
       DefForcing->SfcStressForcing.ZonalStressCell, std::pair(0, NCells));
   auto SfcStressMeridOwned = Kokkos::subview(
       DefForcing->SfcStressForcing.MeridStressCell, std::pair(0, NCells));

   auto ApplyPass = arraysEqual(SfcStressZonalOwned, ExpectedSfcStressZonal) &&
                    arraysEqual(SfcStressMeridOwned, ExpectedSfcStressMerid);

   if (ApplyPass) {
      LOG_INFO("SfcCouplingTest: applyImportFields PASS");
   } else {
      Err++;
      LOG_ERROR("SfcCouplingTest: applyImportFields FAIL");
   }

   SfcCoupling::clear();

   return Err;
}

int testUpdateExportFields(const I4 NSteps) {

   int Err = 0;

   auto *DefStepper  = TimeStepper::getDefault();
   Clock *ModelClock = DefStepper->getClock();

   // Reset the shared clock
   ModelClock->setCurrentTime(DefStepper->getStartTime());

   // Coupling interval spans NSteps ocean timesteps
   auto CouplingParams = mockCouplingInitParams(
       CouplingLayout::MCT, DefStepper->getTimeStep() * NSteps);

   Err += SfcCoupling::init(CouplingParams);

   SfcCoupling *DefCoupling = SfcCoupling::getDefault();
   OceanState *DefState     = OceanState::getDefault();
   VertCoord *DefVertCoord  = VertCoord::getDefault();

   int NCells = DefCoupling->NCellsOwned;

   Real TempBase  = static_cast<Real>(CouplingParams.ExportIdx.at("So_t"));
   Real SalinBase = static_cast<Real>(CouplingParams.ExportIdx.at("So_s"));

   I4 TempIdx, SalinIdx;
   Tracers::getIndex(TempIdx, "Temperature");
   Tracers::getIndex(SalinIdx, "Salinity");

   while (!DefCoupling->CouplingAlarm.isRinging()) {
      Real CurrStep = static_cast<Real>(DefCoupling->getNAccumSteps());

      HostArray2DReal TempH  = Tracers::getHostByIndex(0, TempIdx);
      HostArray2DReal SalinH = Tracers::getHostByIndex(0, SalinIdx);

      for (int Cell = 0; Cell < NCells; Cell++) {
         int KSfc = DefVertCoord->MinLayerCell(Cell);

         TempH(Cell, KSfc)  = TempBase + Cell + CurrStep;
         SalinH(Cell, KSfc) = SalinBase + Cell + CurrStep;
      }
      Tracers::copyToDevice(0);

      DefCoupling->updateExportFields(DefState, Tracers::getAll(0));

      ModelClock->advance();
   }

   // Sanity check: alarm should ring after NSteps
   if (DefCoupling->getNAccumSteps() != NSteps) {
      Err++;
      LOG_ERROR("SfcCouplingTest: updateExportFields FAIL - "
                "NAccumSteps = {}, expected {}",
                DefCoupling->getNAccumSteps(), NSteps);
   }

   auto AvgTempH =
       createHostMirrorCopy(DefCoupling->OcnToCpl.AvgSfcTemperature);
   auto AvgSalinH = createHostMirrorCopy(DefCoupling->OcnToCpl.AvgSfcSalinity);

   Real Tol        = 1e-10;
   Real StepOffset = static_cast<Real>(NSteps - 1) / 2.0;
   HostArray1DReal ExpectedTemp =
       makeCellVarryingArray("ExpectedTemp", NCells, TempBase + StepOffset);
   HostArray1DReal ExpectedSalin =
       makeCellVarryingArray("ExpectedSalin", NCells, SalinBase + StepOffset);

   I4 TempErr  = 0;
   I4 SalinErr = 0;

   // Will this fail for single precision, given the tolerance?
   for (int Cell = 0; Cell < NCells; Cell++) {
      if (std::abs(AvgTempH(Cell) - ExpectedTemp(Cell)) > Tol) {
         TempErr++;
      }

      if (std::abs(AvgSalinH(Cell) - ExpectedSalin(Cell)) > Tol) {
         SalinErr++;
      }
   }

   if (TempErr == 0) {
      LOG_INFO("SfcCouplingTest: updateExportFields PASS - "
               "AvgSfcTemperature within tolerance of {}",
               Tol);
   } else {
      Err += TempErr;
      LOG_ERROR("SfcCouplingTest: updateExportFields FAIL - "
                "AvgSfcTemperature outside tolerance of {} for {} cells",
                Tol, TempErr);
   }

   if (SalinErr == 0) {
      LOG_INFO("SfcCouplingTest: updateExportFields PASS - "
               "AvgSfcSalinity within tolerance of {}",
               Tol);
   } else {
      Err += SalinErr;
      LOG_ERROR("SfcCouplingTest: updateExportFields FAIL - "
                "AvgSfcSalinity outside tolerance of {} for {} cells",
                Tol, SalinErr);
   }

   // reset model clock to the start time for any subsequent tests
   ModelClock->setCurrentTime(DefStepper->getStartTime());

   SfcCoupling::clear();
   return Err;
}

int testExportToCoupler(const CouplingLayout Layout) {

   int Err = 0;

   auto CouplingParams = mockCouplingInitParams(Layout);
   Err += SfcCoupling::init(CouplingParams);

   SfcCoupling *DefCoupling = SfcCoupling::getDefault();
   VertCoord *DefVertCoord  = VertCoord::getDefault();

   int NCells   = DefCoupling->NCellsOwned;
   int NImports = DefCoupling->NImportFields;
   int NExports = DefCoupling->NExportFields;

   std::vector<Real> CplToOcnData(NCells * NImports, 0.0);
   std::vector<Real> OcnToCplData(NCells * NExports, 0.0);

   int TempIdx  = CouplingParams.ExportIdx.at("So_t");
   int SalinIdx = CouplingParams.ExportIdx.at("So_s");
   int VelUIdx  = CouplingParams.ExportIdx.at("So_u");
   int VelVIdx  = CouplingParams.ExportIdx.at("So_v");
   int SshIdx   = CouplingParams.ExportIdx.at("So_ssh");

   DefCoupling->attachData(CplToOcnData.data(), OcnToCplData.data());

   HostArray1DReal ExpectedTemp =
       makeCellVarryingArray("ExpectedTemp", NCells, Real(TempIdx));
   HostArray1DReal ExpectedSalin =
       makeCellVarryingArray("ExpectedSalin", NCells, Real(SalinIdx));
   HostArray1DReal ExpectedVelU =
       makeCellVarryingArray("ExpectedVelU", NCells, Real(VelUIdx));
   HostArray1DReal ExpectedVelV =
       makeCellVarryingArray("ExpectedVelV", NCells, Real(VelVIdx));
   HostArray1DReal ExpectedSsh =
       makeCellVarryingArray("ExpectedSsh", NCells, Real(SshIdx));

   deepCopy(DefCoupling->OcnToCpl.AvgSfcTemperature, ExpectedTemp);
   deepCopy(DefCoupling->OcnToCpl.AvgSfcSalinity, ExpectedSalin);
   deepCopy(DefCoupling->OcnToCpl.AvgSfcVelocityZonal, ExpectedVelU);
   deepCopy(DefCoupling->OcnToCpl.AvgSfcVelocityMerid, ExpectedVelV);

   auto SshCellOwned =
       Kokkos::subview(DefVertCoord->SshCell, std::pair(0, NCells));
   deepCopy(SshCellOwned, ExpectedSsh);

   DefCoupling->exportToCoupler();

   // Check 1: exportToCoupler properly packs into OcnToCplView
   int PackErr = 0;
   for (int Cell = 0; Cell < NCells; Cell++) {
      if (OcnToCplData[flatIdx(Layout, Cell, TempIdx, NCells, NExports)] !=
          ExpectedTemp(Cell)) {
         PackErr++;
      }
      if (OcnToCplData[flatIdx(Layout, Cell, SalinIdx, NCells, NExports)] !=
          ExpectedSalin(Cell)) {
         PackErr++;
      }
      if (OcnToCplData[flatIdx(Layout, Cell, VelUIdx, NCells, NExports)] !=
          ExpectedVelU(Cell)) {
         PackErr++;
      }
      if (OcnToCplData[flatIdx(Layout, Cell, VelVIdx, NCells, NExports)] !=
          ExpectedVelV(Cell)) {
         PackErr++;
      }
      if (OcnToCplData[flatIdx(Layout, Cell, SshIdx, NCells, NExports)] !=
          ExpectedSsh(Cell)) {
         PackErr++;
      }
   }

   if (PackErr == 0) {
      LOG_INFO("SfcCouplingTest: exportToCoupler with {} layout PASS",
               toString(Layout));
   } else {
      Err += PackErr;
      LOG_ERROR("SfcCouplingTest: exportToCoupler with {} layout FAIL - "
                "{} packing errors",
                toString(Layout), PackErr);
   }

   // Check 2: resetFields() zeroed the running-average accumulators
   HostArray1DReal ZeroArray("ZeroArray", NCells);
   deepCopy(ZeroArray, 0.0);

   // Refresh OcnToCpl's own host mirrors post-reset, rather than creating
   // separate test-only mirrors
   DefCoupling->OcnToCpl.copyToHost();

   int ResetErr = 0;
   if (!arraysEqual(DefCoupling->OcnToCpl.AvgSfcTemperatureH, ZeroArray))
      ResetErr++;
   if (!arraysEqual(DefCoupling->OcnToCpl.AvgSfcSalinityH, ZeroArray))
      ResetErr++;
   if (!arraysEqual(DefCoupling->OcnToCpl.AvgSfcVelocityZonalH, ZeroArray))
      ResetErr++;
   if (!arraysEqual(DefCoupling->OcnToCpl.AvgSfcVelocityMeridH, ZeroArray))
      ResetErr++;

   if (ResetErr == 0) {
      LOG_INFO("SfcCouplingTest: exportToCoupler resetFields PASS");
   } else {
      Err += ResetErr;
      LOG_ERROR("SfcCouplingTest: exportToCoupler resetFields FAIL");
   }

   // Check 3: NAccumSteps was reset to 0
   if (DefCoupling->getNAccumSteps() != 0) {
      Err++;
      LOG_ERROR("SfcCouplingTest: exportToCoupler reset counter FAIL");
   } else {
      LOG_INFO("SfcCouplingTest: exportToCoupler reset counter PASS");
   }

   SfcCoupling::clear();

   return Err;
}
void finalizeSfcCouplingTest() {

   Tracers::clear();
   Forcing::clear();
   OceanState::clear();
   VertCoord::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   TimeStepper::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int sfcCouplingTest(const std::string &MeshFile = "OmegaMesh.nc") {

   int Err = initSfcCouplingTest(MeshFile);

   if (Err != 0) {
      LOG_CRITICAL("SfcCouplingTest: Error initializing");
   }

   Err += testImportFromCoupler(CouplingLayout::MCT);
   Err += testImportFromCoupler(CouplingLayout::MOAB);

   Err += testApplyImportFields();

   Err += testUpdateExportFields(1);
   Err += testUpdateExportFields(5);

   Err += testExportToCoupler(CouplingLayout::MCT);
   Err += testExportToCoupler(CouplingLayout::MOAB);

   if (Err == 0) {
      LOG_INFO("SfcCouplingTest: Successful completion");
   }

   finalizeSfcCouplingTest();

   return Err;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   RetVal += sfcCouplingTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;
   if (RetVal == 0)
      LOG_INFO("------ SufaceCoupling unit tests successful ------");

   return RetVal;
} // end of main
