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
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "PGrad.h"
#include "Tracers.h"
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

   // Initialize streams
   IOStream::init();

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

   // Initialize tracers
   Tracers::init();

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
      //VertCoord *DefVCoord = VertCoord::getDefault();
      VertCoord *VCoord = VertCoord::getDefault();
      OceanState *DefState = OceanState::getDefault();
      Decomp *DefDecomp = Decomp::getDefault();
      Eos *DefEos = Eos::getInstance();
      Config *Options = Config::getOmegaConfig();


      // create arrays: Tend on edges, Pressure/Geopotential/SpecVol on cells
      Array2DReal Tend("Tend", DefMesh->NEdgesSize, VCoord->NVertLayers);
      Array2DReal SpecVolOld("SpecVolOld", DefMesh->NCellsSize, VCoord->NVertLayers);
      Array2DReal PressureMidOld("PressureMidOld", DefMesh->NCellsSize, VCoord->NVertLayers);
      Array1DReal SurfacePressure("SurfacePressure", DefMesh->NCellsSize);

      I4 NVertLayers = 60;
      Real dC = 30000.0_Real;
      for (int refinement = 0; refinement < 4; ++refinement) {

         LOG_INFO("PGradTest: Starting refinement level {}", refinement);
         VCoord->NVertLayers = NVertLayers;
         VCoord->NVertLayersP1 = NVertLayers + 1;

         auto &MinLayerCell = VCoord->MinLayerCell;
         auto &MaxLayerCell = VCoord->MaxLayerCell;
         parallelFor({DefMesh->NCellsAll}, KOKKOS_LAMBDA(int i) {
            MinLayerCell(i) = 0;
            MaxLayerCell(i) = NVertLayers - 1;
         });

         auto &CellsOnEdge = DefMesh->CellsOnEdge;
         auto &DcEdge = DefMesh->DcEdge;
         auto &EdgeMask = DefMesh->EdgeMask;
         parallelFor({DefMesh->NEdgesAll}, KOKKOS_LAMBDA(int e) {
            CellsOnEdge(e, 0)= 0;
            CellsOnEdge(e, 1)= 1;
            DcEdge(e) = dC;
            //for (int k = 0; k < NVertLayers; ++k) {
            //   EdgeMask(e, k) = 1.0_Real;
            //}
         });     

         // Fetch reference desnity from Config
         Real Density0;
         Error ErrorCode;
         Config TendConfig("Tendencies");
         ErrorCode.reset();
         ErrorCode += Options->get(TendConfig);
         CHECK_ERROR_ABORT(ErrorCode, "VertCoord: Tendencies group not found in Config");
         ErrorCode += TendConfig.get("Density0", Density0);
         CHECK_ERROR_ABORT(ErrorCode, "VertCoord: Density0 not found in TendConfig");
      
         I4 TimeLevel = 0;

         // get state and tracer arrays
         Array2DReal LayerThick;
         DefState->getLayerThickness(LayerThick, TimeLevel);
         Array2DReal Temp;
         Array2DReal Salinity;
         Err = Tracers::getByName(Temp, TimeLevel, "Temperature");
         Err = Tracers::getByName(Salinity, TimeLevel, "Salinity");



         // set Z interface and mid-point locations
         Real ZBottom = -1000.0_Real;
         Real dZ = 2.0_Real * (-ZBottom / NVertLayers);
         auto &BottomDepth = VCoord->BottomDepth;
         auto &ZInterface = VCoord->ZInterface;
         auto &ZMid = VCoord->ZMid;
         Real tilt_factor = 0.495_Real;
         //Real tilt_factor = 0.45_Real;
         parallelFor({DefMesh->NCellsAll}, KOKKOS_LAMBDA(int i) {
            ZInterface(i, NVertLayers) = ZBottom;
            SurfacePressure(i) = 0.0_Real;
            BottomDepth(i) = 0.0_Real;
            for (int k = NVertLayers - 1; k >= 0; --k) {
               Real x = (k + i) % 2;
               Real dz = ( 2.0_Real * tilt_factor - 1.0_Real ) * x * dZ + (1.0_Real - tilt_factor) * dZ; // staggered layer thickness
               ZInterface(i, k) = ZInterface(i, k + 1) + dz;
               LayerThick(i, k) = ZInterface(i, k) - ZInterface(i, k + 1);
               ZMid(i, k) = 0.5_Real * (ZInterface(i, k) + ZInterface(i, k + 1));
               BottomDepth(i) += dz;
            }
         });

         LOG_INFO("NVertLayers = {}", NVertLayers);
         LOG_INFO("dC = {}", dC);
         //for (int i = 0; i < 2; ++i) {
         //   for (int k = 0; k <= NVertLayers; ++k) {
         //      LOG_INFO("ZInterface({}, {}) = {}", i, k, ZInterface(i, k));
         //   }
         //}
         //LOG_INFO("NVertLayers = {}", NVertLayers);
         for (int i = 0; i < 2; ++i) {
            //for (int k = 0; k <= NVertLayers; ++k) {
            for (int k = 0; k < 2; ++k) {
               LOG_INFO("LayerThick({}, {}) = {}", i, k, LayerThick(i, k));
            }
         }
         //LOG_INFO("NVertLayers = {}", NVertLayers);
         //for (int i = 0; i < 2; ++i) {
         //   for (int k = 0; k <= NVertLayers; ++k) {
         //      LOG_INFO("ZMid({}, {}) = {}", i, k, ZMid(i, k));
         //   }
         //}

         // set simple temperature and salinity profiles
         auto &SpecVol = DefEos->SpecVol;
         parallelFor({DefMesh->NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
            Real T0 = 30.0;
            Real TB = 5.0;
            Real S0 = 30.0;
            Real SB = 40.0;


            Real phi0  = (ZMid(i, k) - ZBottom) / (-ZBottom);
            Real phiB  = 1.0_Real - phi0;

            Temp(i, k) =  T0 * phi0 + TB * phiB;
            Salinity(i, k) = S0 * phi0 + SB * phiB;
            SpecVol(i, k) = 1.0_Real / Density0;
            SpecVolOld(i, k) = SpecVol(i, k);
         });

         //for (int i = 0; i < 2; ++i) {
         //   for (int k = 0; k < NVertLayers; ++k) {
         //      LOG_INFO("Temp({}, {}) = {}", i, k, Temp(i, k));
         //   }
         //}
         //for (int i = 0; i < 2; ++i) {
         //   for (int k = 0; k < NVertLayers; ++k) {
         //      LOG_INFO("Salinity({}, {}) = {}", i, k, Salinity(i, k));
         //   }
         //}

         // iterate to converge SpecVol
         auto &PressureMid = VCoord->PressureMid;
         VCoord->computePressure(LayerThick, SurfacePressure);
         deepCopy(PressureMidOld, PressureMid);
         //for (int i = 0; i < 2; ++i) {
         //   for (int k = 0; k < NVertLayers; ++k) {
         //      LOG_INFO("PressureMid({}, {}) = {}, PressureMidOld({}, {}) = {}", i, k, PressureMid(i, k), i, k, PressureMidOld(i, k));
         //   }
         //}
         for (int iteration = 0; iteration < 15; ++iteration) {

            // compute specific volume from EOS
            VCoord->computePressure(LayerThick, SurfacePressure);
            DefEos->computeSpecVol(Temp, Salinity, PressureMid);

            // compute psuedo thickness from specific volume
            parallelFor({DefMesh->NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
               LayerThick(i, k) = 1.0_Real / (SpecVol(i, k) * Density0) * (ZInterface(i, k) - ZInterface(i, k + 1));
            });

            // compute difference from previous iteration
            Real max_value = 0.0_Real;
            parallelReduce({DefMesh->NCellsAll, NVertLayers},
                           KOKKOS_LAMBDA(int i, int k, Real &max) {
                              Real diff = Kokkos::abs(SpecVol(i, k) - SpecVolOld(i, k));
                              if (diff > max) max = diff;
                           }, Kokkos::Max<Real>(max_value)  );

            // check convergence               
            if (max_value < 1e-12_Real) {
               LOG_INFO("converged: max diff = {}", max_value);
               break;
            } else {
               //LOG_INFO("max diff = {}", max_value);
               parallelFor({DefMesh->NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
                  SpecVolOld(i, k) = SpecVol(i, k);
               });
            }

         }
         //for (int i = 0; i < 2; ++i) {
         //   for (int k = 0; k < NVertLayers; ++k) {
         //      LOG_INFO("SpecVol({}, {}) = {}", i, k, SpecVol(i, k));
         //   }
         //}

         // compute geopotential
         //VCoord->computeZHeight(LayerThick, SpecVol);
         Array1DReal SelfAttractionLoading("SelfAttractionLoading", DefMesh->NCellsAll);
         Array1DReal TidalPotential("TidalPotential", DefMesh->NCellsAll);
         deepCopy(TidalPotential, 0.0_Real);
         deepCopy(SelfAttractionLoading, 0.0_Real);
         VCoord->computeGeopotential(TidalPotential, SelfAttractionLoading);

         // create PressureGrad instance 
         PressureGrad *DefPGrad = PressureGrad::getDefault();
         if (!DefPGrad) {
            LOG_INFO("PGrad: default instance not created by init");
         }
         //auto PGradTest = PressureGrad::create("Test", DefMesh, VCoord, Options);

         //for (int i = 0; i < 2; ++i) {
         //  for (int k = 0; k < NVertLayers; ++k) {
         //     Tend(i, k) = 0.0_Real;
         //     LOG_INFO("Tend({}, {}) = {}", i, k, Tend(i, k));
         //  }
         //} 

         // compute pressure gradient
         //PGradTest->computePressureGrad(Tend, DefState, VCoord, DefEos, TimeLevel);
         deepCopy(Tend, 0.0_Real);
         DefPGrad->computePressureGrad(Tend, DefState, VCoord, DefEos, TimeLevel);
         //for (int i = 0; i < 2; ++i) {
         //   for (int k = 0; k < NVertLayers; ++k) {
         //      LOG_INFO("Tend({}, {}) = {}", i, k, Tend(i, k));
         //   }
         //}
         Real max_value = 0.0_Real;
         parallelReduce({DefMesh->NEdgesAll, NVertLayers - 1},
                        KOKKOS_LAMBDA(int i, int k, Real &max) {
                           Real val = Kokkos::abs(Tend(i + 1, k));
                           if (val > max) max = val;
                        }, Kokkos::Max<Real>(max_value)  );
         Real sum_value = 0.0_Real;
         parallelReduce({DefMesh->NEdgesAll, NVertLayers - 1},
                        KOKKOS_LAMBDA(int i, int k, Real &lsum) {
                           lsum += Tend(i + 1, k) * Tend(i + 1, k);
                        }, Kokkos::Sum<Real>(sum_value)  );
         LOG_INFO("refinement level {}: max |Tend| = {}, average Tend = {}", refinement, max_value, std::sqrt(sum_value) / (DefMesh->NEdgesAll * (NVertLayers - 1)));

         dC = dC * 2.0_Real;
         NVertLayers = NVertLayers / 2;

      } // refinement loop

      // cleanup
      PressureGrad::clear();
      IOStream::finalize();
      TimeStepper::clear();
      Tracers::clear();
      Eos::destroyInstance();
      OceanState::clear();
      VertCoord::clear();
      Field::clear();
      Dimension::clear();
      HorzMesh::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256) RetVal = 255;
   return RetVal;
}
