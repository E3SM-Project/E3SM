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

   std::map<std::string, int> ImportIdx = {
       {"Foxx_taux", 0},  {"Foxx_tauy", 1},   {"Foxx_swnet", 2},
       {"Foxx_lwnet", 3}, {"Foxx_lat", 4},    {"Foxx_sen", 5},
       {"Foxx_lwup", 6},  {"Faxa_lwdn", 7},   {"Fioi_melth", 8},
       {"Fioi_bergh", 9}, {"Faxa_snow", 10},  {"Faxa_rain", 11},
       {"Foxx_evap", 12}, {"Fioi_meltw", 13}, {"Fioi_bergw", 14},
       {"Fioi_salt", 15}, {"Foxx_rofl", 16},  {"Foxx_rofi", 17},
       {"Si_ifrac", 18},  {"Si_bpress", 19},  {"Sa_pslv", 20}};

   std::map<std::string, int> ExportIdx = {
       {"So_t", 0},        {"So_s", 1},        {"So_u", 2},    {"So_v", 3},
       {"So_ssh", 4},      {"So_dhdx", 5},     {"So_dhdy", 6}, {"Fioo_q", 7},
       {"Fioo_frazil", 8}, {"Faoo_h2otemp", 9}};
};

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

int testSfcCoupling() {

   int Err = 0;

   TestSetup Setup;

   CouplingInitParams CouplingParams;
   CouplingParams.ImportIdx        = Setup.ImportIdx;
   CouplingParams.ExportIdx        = Setup.ExportIdx;
   CouplingParams.CouplingTimeStep = TimeInterval(3600.0, TimeUnits::Seconds);
   CouplingParams.Layout           = CouplingLayout::MCT;

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

   Err += testSfcCoupling();

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
