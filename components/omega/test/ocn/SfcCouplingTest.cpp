#include "SfcCoupling.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Eos.h"
#include "Field.h"
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

   std::map<std::string, int> ImportIdx = {{"Foxx_taux", 0}, {"Foxx_tauy", 1}};
   std::map<std::string, int> ExportIdx = {
       {"So_t", 0}, {"So_u", 1}, {"So_v", 2}};
};

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

int initRawData() {
   int Err = 0;

   return Err;
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

   return Err;
}

int testSfcCoupling(const CouplingLayout Layout,
                    const TimeInterval CouplingTimeStep) {

   int Err = 0;

   TestSetup Setup;

   CouplingInitParams CouplingParams{.ImportIdx        = Setup.ImportIdx,
                                     .ExportIdx        = Setup.ExportIdx,
                                     .CouplingTimeStep = CouplingTimeStep,
                                     .Layout           = Layout};

   Err += SfcCoupling::init(CouplingParams);

   // test retrival of default
   SfcCoupling *DefCoupling = SfcCoupling::getDefault();

   if (DefCoupling) {
      LOG_INFO("SfcCouplingTest: Default retrival PASS");
   } else {
      Err++;
      LOG_ERROR("SfcCouplingTest: Default retrival FAIL");
      return -1;
   }

   int NCells   = DefCoupling->NCellsOwned;
   int NImports = DefCoupling->NImportFields;
   int NExports = DefCoupling->NExportFields;

   std::vector<Real> CplToOcnData(NCells * NImports, 0.0);
   std::vector<Real> OcnToCplData(NCells * NExports, 0.0);

   // Index formula depend on the Layout
   auto flatIdx = [&](int Cell, int Field) -> int {
      if (Layout == CouplingLayout::MCT)
         return Cell * NImports + Field;
      else // MOAB
         return Field * NCells + Cell;
   };

   for (int Cell = 0; Cell < NCells; Cell++)
      for (int Field = 0; Field < NImports; Field++)
         CplToOcnData[flatIdx(Cell, Field)] = static_cast<Real>(Field);

   DefCoupling->attachData(CplToOcnData.data(), OcnToCplData.data());
   DefCoupling->importFromCoupler();

   bool ImportPass = true;
   for (int Cell = 0; Cell < NCells; Cell++) {
      if (DefCoupling->CplToOcn.SurfaceStressZonal(Cell) != Real(0))
         ImportPass = false;
      if (DefCoupling->CplToOcn.SurfaceStressMeridional(Cell) != Real(1))
         ImportPass = false;
   }

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

void finalizeSfcCouplingTest() {

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

int surfaceCouplingTest(const std::string &MeshFile = "OmegaMesh.nc") {

   int Err = initSfcCouplingTest(MeshFile);

   if (Err != 0) {
      LOG_CRITICAL("SfcCouplingTest: Error initializing");
   }

   auto *DefTimeStepper     = TimeStepper::getDefault();
   TimeInterval OcnTimeStep = DefTimeStepper->getTimeStep();

   Err += testSfcCoupling(CouplingLayout::MCT, OcnTimeStep);
   Err += testSfcCoupling(CouplingLayout::MOAB, OcnTimeStep);

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

   RetVal += surfaceCouplingTest();

   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;
   if (RetVal == 0)
      LOG_INFO("------ SufaceCoupling unit tests successful ------");

   return RetVal;
} // end of main
