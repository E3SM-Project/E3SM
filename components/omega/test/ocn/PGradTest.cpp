//===-- Test driver for OMEGA Pressure Gradient (PGrad) --------------*- C++-*-===/
//
/// \file
/// \brief Test driver for PressureGrad module
//
//===----------------------------------------------------------------------===/

#include "PGrad.h"

#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Eos.h"
#include "Field.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "PGrad.h"
#include "TimeStepper.h"
#include "VertCoord.h"
#include "mpi.h"

using namespace OMEGA;

void initPGradTest() {

   Error Err;
   int Err1;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Read default config if present
   Config("Omega");
   Config::readAll("omega.yml");

   // First step of time stepper initialization needed for IOstream
   TimeStepper::init1();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize the default halo
   Err1 = Halo::init();
   if (Err1 != 0) {
      LOG_ERROR("PGrad: error initializing default halo");
      Err += Error(ErrorCode::Fail,
                   "PGrad: error initializing default halo");
   }

   // Begin initialization of the default vertical coordinate
   VertCoord::init1();

   // Initialize the default mesh
   HorzMesh::init();

   // Complete initialization of the default vertical coordinate
   VertCoord::init2();

   // Initialize the equation of state
   Eos::init();

   // Initialize ocean state
   OceanState::init();


   CHECK_ERROR_ABORT(Err, "PGrad: error during initialization");
}

int main(int argc, char *argv[]) {
   int RetVal = 0;
   int Err;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      initPGradTest();

      // Initialize PressureGrad manager
      PressureGrad::init();

      MachEnv *DefEnv = MachEnv::getDefault();
      HorzMesh *DefMesh = HorzMesh::getDefault();
      VertCoord *DefVCoord = VertCoord::getDefault();
      OceanState *DefState = OceanState::getDefault();
      Eos *DefEos = Eos::getInstance();
      Config *Options = Config::getOmegaConfig();

      I4 NEdgesOwned = DefMesh->NEdgesOwned;
      I4 NVertLayers = VertCoord::getDefault()->NVertLayers;
      I4 NChunks = NVertLayers / VecLength;

      // create arrays: Tend on edges, Pressure/Geopotential/SpecVol on cells
      Array2DReal Tend("Tend", DefMesh->NEdgesAll, NVertLayers);
      Array2DReal Pressure("Pressure", DefMesh->NCellsAll, NVertLayers);
      Array2DReal Geopotential("Geopotential", DefMesh->NCellsAll, NVertLayers);
      Array2DReal SpecVol("SpecVol", DefMesh->NCellsAll, NVertLayers);

      // initialize arrays to zero
      parallelFor({DefMesh->NEdgesAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
         Tend(i,k) = 0.0_Real;
      });
      parallelFor({DefMesh->NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
         Pressure(i,k) = 0.0_Real;
         Geopotential(i,k) = 0.0_Real;
         SpecVol(i,k) = 1.0_Real; // unit specific volume
      });

      // set a simple pressure difference between two neighboring cells
      if (DefMesh->NCellsAll >= 2) {
         Pressure(0,0) = 1.0_Real;
         Pressure(1,0) = 2.0_Real;
      }

      // call computePressureGrad
      PressureGrad *Pg = PressureGrad::getDefault();
      if (!Pg) {
         LOG_INFO("PGrad: default instance not present, creating via create");
         Pg = PressureGrad::create("default", DefMesh, DefVCoord, Options);
      }

      int TimeLevel = 0;
      if (Pg) {
         Pg->computePressureGrad(Tend, DefState, DefVCoord, DefEos, TimeLevel);

         // simple check: ensure Tend contains some finite values (not all zero)
         I4 CountNonZero = 0;
         parallelReduce("countNonZero", {NEdgesOwned, NVertLayers},
                        KOKKOS_LAMBDA(int e, int k, I4 &acc) {
                           if (Tend(e,k) != 0.0_Real) acc++;
                        }, CountNonZero);

         if (CountNonZero > 0) {
            LOG_INFO("PGrad: computePressureGrad produced non-zero tendencies PASS");
         } else {
            RetVal += 1;
            LOG_ERROR("PGrad: computePressureGrad produced all-zero tendencies FAIL");
         }
      } else {
         RetVal += 1;
         LOG_ERROR("PGrad: failed to obtain PressureGrad instance FAIL");
      }

      // cleanup
      PressureGrad::clear();
      VertCoord::clear();
      OceanState::clear();
      HorzMesh::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
      FieldGroup::clear();
      Field::clear();
      Dimension::clear();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256) RetVal = 255;
   return RetVal;
}
