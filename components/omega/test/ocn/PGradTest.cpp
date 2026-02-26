//===-- Test driver for OMEGA Pressure Gradient (PGrad) --------------*-
// C++-*-===/
//
/// \file
/// \brief Test driver for PressureGrad module
//
//===----------------------------------------------------------------------===/

#include "PGrad.h"

#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Eos.h"
#include "Error.h"
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "PGrad.h"
#include "Pacer.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "mpi.h"

using namespace OMEGA;

void initPGradTest() {

   Error Err;
   int Err1;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
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
      Err += Error(ErrorCode::Fail, "PGrad: error initializing default halo");
   }

   // Initialize the default mesh
   HorzMesh::init();

   // Initialize the default vertical coordinate
   VertCoord::init();

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
      // Initialize default PressureGrad
      PressureGrad::init();

      MachEnv *DefEnv      = MachEnv::getDefault();
      HorzMesh *DefMesh    = HorzMesh::getDefault();
      VertCoord *VCoord    = VertCoord::getDefault();
      OceanState *DefState = OceanState::getDefault();
      Decomp *DefDecomp    = Decomp::getDefault();
      Eos *DefEos          = Eos::getInstance();
      Config *Options      = Config::getOmegaConfig();

      // create arrays: Tend on edges, Pressure/Geopotential/SpecVol on cells
      Array2DReal Tend("Tend", DefMesh->NEdgesSize, VCoord->NVertLayers);
      Array2DReal SpecVolOld("SpecVolOld", DefMesh->NCellsSize,
                             VCoord->NVertLayers);
      Array2DReal PressureMidOld("PressureMidOld", DefMesh->NCellsSize,
                                 VCoord->NVertLayers);
      Array1DReal SurfacePressure("SurfacePressure", DefMesh->NCellsSize);

      I4 NEdgesAll = DefMesh->NEdgesAll;
      I4 NCellsAll = DefMesh->NCellsAll;

      I4 NVertLayers  = 60;
      Real dC         = 30000.0_Real;
      I4 NRefinements = 4;
      HostArray1DReal RMSE("RMSE", NRefinements);
      for (int refinement = 0; refinement < NRefinements; ++refinement) {

         LOG_INFO("PGradTest: Starting refinement level {}", refinement);
         VCoord->NVertLayers   = NVertLayers;
         VCoord->NVertLayersP1 = NVertLayers + 1;

         auto &MinLayerCell = VCoord->MinLayerCell;
         auto &MaxLayerCell = VCoord->MaxLayerCell;
         parallelFor(
             {NCellsAll}, KOKKOS_LAMBDA(int i) {
                MinLayerCell(i) = 0;
                MaxLayerCell(i) = NVertLayers - 1;
             });

         auto &MinLayerEdgeBot = VCoord->MinLayerEdgeBot;
         auto &MaxLayerEdgeTop = VCoord->MaxLayerEdgeTop;
         parallelFor(
             {NEdgesAll}, KOKKOS_LAMBDA(int i) {
                MinLayerEdgeBot(i) = 0;
                MaxLayerEdgeTop(i) = NVertLayers - 1;
                // MaxLayerEdgeTop(i) = 4;
             });

         auto &CellsOnEdge = DefMesh->CellsOnEdge;
         auto &DcEdge      = DefMesh->DcEdge;
         auto &EdgeMask    = VCoord->EdgeMask;
         parallelFor(
             {NEdgesAll}, KOKKOS_LAMBDA(int i) {
                CellsOnEdge(i, 0) = 0;
                CellsOnEdge(i, 1) = 1;
                DcEdge(i)         = dC;
             });

         // Fetch reference desnity from Config
         Real Density0;
         Density0 = RhoSw;

         I4 TimeLevel = 0;

         // get state and tracer arrays
         Array2DReal LayerThick = DefState->getLayerThickness(TimeLevel);
         Array2DReal Temp       = Tracers::getByName(TimeLevel, "Temperature");
         Array2DReal Salinity   = Tracers::getByName(TimeLevel, "Salinity");

         // set Z interface and mid-point locations
         Real ZBottom      = -1000.0_Real;
         Real dZ           = 2.0_Real * (-ZBottom / NVertLayers);
         auto &BottomDepth = VCoord->BottomDepth;
         auto &ZInterface  = VCoord->ZInterface;
         auto &ZMid        = VCoord->ZMid;
         Real tilt_factor  = 0.495_Real;
         // Real tilt_factor = 0.45_Real;
         parallelFor(
             {NCellsAll}, KOKKOS_LAMBDA(int i) {
                ZInterface(i, NVertLayers) = ZBottom;
                SurfacePressure(i)         = 0.0_Real;
                BottomDepth(i)             = 0.0_Real;
                for (int k = NVertLayers - 1; k >= 0; --k) {
                   Real x  = (k + i) % 2;
                   Real dz = (2.0_Real * tilt_factor - 1.0_Real) * x * dZ +
                             (1.0_Real - tilt_factor) *
                                 dZ; // staggered layer thickness
                   ZInterface(i, k) = ZInterface(i, k + 1) + dz;
                   LayerThick(i, k) = ZInterface(i, k) - ZInterface(i, k + 1);
                   ZMid(i, k) =
                       0.5_Real * (ZInterface(i, k) + ZInterface(i, k + 1));
                   BottomDepth(i) += dz;
                }
             });

         LOG_INFO("NVertLayers = {}", NVertLayers);
         LOG_INFO("dC = {}", dC);
         DefState->copyToHost(0);
         HostArray2DReal LayerThickH = DefState->getLayerThicknessH(TimeLevel);
         for (int i = 0; i < 2; ++i) {
            for (int k = 0; k < 2; ++k) {
               LOG_INFO("LayerThick({}, {}) = {}", i, k, LayerThickH(i, k));
            }
         }

         // set simple temperature and salinity profiles
         auto &SpecVol = DefEos->SpecVol;
         parallelFor(
             {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
                Real T0 = 30.0;
                Real TB = 5.0;
                Real S0 = 30.0;
                Real SB = 40.0;

                Real phi0 = (ZMid(i, k) - ZBottom) / (-ZBottom);
                Real phiB = 1.0_Real - phi0;

                Temp(i, k)       = T0 * phi0 + TB * phiB;
                Salinity(i, k)   = S0 * phi0 + SB * phiB;
                SpecVol(i, k)    = 1.0_Real / Density0;
                SpecVolOld(i, k) = SpecVol(i, k);
             });

         // Iterate to converge LayerThick, SpecVol, PressureMid
         auto &PressureMid = VCoord->PressureMid;
         VCoord->computePressure(LayerThick, SurfacePressure);
         deepCopy(PressureMidOld, PressureMid);
         for (int iteration = 0; iteration < 15; ++iteration) {

            // compute specific volume from EOS
            VCoord->computePressure(LayerThick, SurfacePressure);
            DefEos->computeSpecVol(Temp, Salinity, PressureMid);

            // compute psuedo thickness from specific volume
            parallelFor(
                {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
                   LayerThick(i, k) = 1.0_Real / (SpecVol(i, k) * Density0) *
                                      (ZInterface(i, k) - ZInterface(i, k + 1));
                });

            // compute difference from previous iteration
            Real max_value = 0.0_Real;
            parallelReduce(
                {NCellsAll, NVertLayers},
                KOKKOS_LAMBDA(int i, int k, Real &max) {
                   Real diff = Kokkos::abs(SpecVol(i, k) - SpecVolOld(i, k));
                   if (diff > max)
                      max = diff;
                },
                Kokkos::Max<Real>(max_value));

            // check convergence
            if (max_value < 1e-12_Real) {
               LOG_INFO("converged: max diff = {}", max_value);
               break;
            } else {
               parallelFor(
                   {NCellsAll, NVertLayers}, KOKKOS_LAMBDA(int i, int k) {
                      SpecVolOld(i, k) = SpecVol(i, k);
                   });
            }
         }

         // compute pressure once more with converged LayerThick
         VCoord->computePressure(LayerThick, SurfacePressure);

         // compute z levels
         VCoord->computeZHeight(LayerThick, SpecVol);

         // get PressureGrad instance
         PressureGrad *DefPGrad = PressureGrad::getDefault();
         if (!DefPGrad) {
            LOG_INFO("PGrad: default instance not created by init");
         }

         // compute pressure gradient
         deepCopy(Tend, 0.0_Real);
         DefPGrad->computePressureGrad(Tend, DefState, VCoord, DefEos,
                                       TimeLevel);

         // compute errors
         Real max_value = 0.0_Real;
         parallelReduce(
             {NEdgesAll, NVertLayers - 2},
             KOKKOS_LAMBDA(int i, int k, Real &max) {
                Real val = Kokkos::abs(Tend(i, k + 1));
                if (val > max)
                   max = val;
             },
             Kokkos::Max<Real>(max_value));
         Real sum_value = 0.0_Real;
         parallelReduce(
             {NEdgesAll, NVertLayers - 2},
             KOKKOS_LAMBDA(int i, int k, Real &lsum) {
                lsum += Tend(i, k + 1) * Tend(i, k + 1);
             },
             Kokkos::Sum<Real>(sum_value));
         Real rmse = std::sqrt(sum_value / (NEdgesAll * (NVertLayers - 2)));
         RMSE(refinement) = rmse;

         LOG_INFO("refinement level {}: max |Tend| = {}, average Tend = {}",
                  refinement, max_value, rmse);

         // coarsen for next iteration
         dC          = dC * 2.0_Real;
         NVertLayers = NVertLayers / 2;

      } // refinement loop

      // Test for second order convergence
      if (RMSE(NRefinements) < RMSE(0) / pow(4.0_Real, NRefinements - 1)) {
         RetVal = 0;
      } else {
         RetVal = 1;
      }

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

   if (RetVal >= 256)
      RetVal = 255;
   return RetVal;
}
