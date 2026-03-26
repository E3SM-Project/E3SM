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

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");

   {
      initPGradTest();
      // Initialize default PressureGrad
      PressureGrad::init();

      HorzMesh *DefMesh    = HorzMesh::getDefault();
      VertCoord *VCoord    = VertCoord::getDefault();
      OceanState *DefState = OceanState::getDefault();
      Eos *DefEos          = Eos::getInstance();

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
      Real DC         = 30000.0_Real;
      I4 NRefinements = 4;
      HostArray1DReal Rmse("Rmse", NRefinements);
      for (int Refinement = 0; Refinement < NRefinements; ++Refinement) {

         LOG_INFO("PGradTest: Starting refinement level {}", Refinement);
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
             });

         auto &CellsOnEdge = DefMesh->CellsOnEdge;
         auto &DcEdge      = DefMesh->DcEdge;
         parallelFor(
             {NEdgesAll}, KOKKOS_LAMBDA(int i) {
                CellsOnEdge(i, 0) = 0;
                CellsOnEdge(i, 1) = 1;
                DcEdge(i)         = DC;
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
         Real DZ           = 2.0_Real * (-ZBottom / NVertLayers);
         auto &BottomDepth = VCoord->BottomDepth;
         auto &ZInterface  = VCoord->ZInterface;
         auto &ZMid        = VCoord->ZMid;
         Real TiltFactor   = 0.495_Real;
         parallelFor(
             {NCellsAll}, KOKKOS_LAMBDA(int i) {
                ZInterface(i, NVertLayers) = ZBottom;
                SurfacePressure(i)         = 0.0_Real;
                BottomDepth(i)             = 0.0_Real;
                for (int k = NVertLayers - 1; k >= 0; --k) {
                   Real X  = (k + i) % 2;
                   Real Dz = (2.0_Real * TiltFactor - 1.0_Real) * X * DZ +
                             (1.0_Real - TiltFactor) *
                                 DZ; // staggered layer thickness
                   ZInterface(i, k) = ZInterface(i, k + 1) + Dz;
                   LayerThick(i, k) = ZInterface(i, k) - ZInterface(i, k + 1);
                   ZMid(i, k) =
                       0.5_Real * (ZInterface(i, k) + ZInterface(i, k + 1));
                   BottomDepth(i) += Dz;
                }
             });

         LOG_INFO("NVertLayers = {}", NVertLayers);
         LOG_INFO("dC = {}", DC);
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

                Real Phi0 = (ZMid(i, k) - ZBottom) / (-ZBottom);
                Real PhiB = 1.0_Real - Phi0;

                Temp(i, k)       = T0 * Phi0 + TB * PhiB;
                Salinity(i, k)   = S0 * Phi0 + SB * PhiB;
                SpecVol(i, k)    = 1.0_Real / Density0;
                SpecVolOld(i, k) = SpecVol(i, k);
             });

         // Iterate to converge LayerThick, SpecVol, PressureMid
         auto &PressureMid = VCoord->PressureMid;
         VCoord->computePressure(LayerThick, SurfacePressure);
         deepCopy(PressureMidOld, PressureMid);
         for (int Iteration = 0; Iteration < 15; ++Iteration) {

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
            Real MaxValue = 0.0_Real;
            parallelReduce(
                {NCellsAll, NVertLayers},
                KOKKOS_LAMBDA(int i, int k, Real &max) {
                   Real Diff = Kokkos::abs(SpecVol(i, k) - SpecVolOld(i, k));
                   if (Diff > max)
                      max = Diff;
                },
                Kokkos::Max<Real>(MaxValue));

            // check convergence
            if (MaxValue < 1e-12_Real) {
               LOG_INFO("converged: max diff = {}", MaxValue);
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
         
         const auto &PressureInterface = VCoord->PressureInterface;
         DefPGrad->computePressureGrad(Tend, PressureMid, PressureInterface,
                                       SpecVol, ZInterface, 
                                       LayerThick);

         // compute errors
         Real MaxValue = 0.0_Real;
         parallelReduce(
             {NEdgesAll, NVertLayers - 2},
             KOKKOS_LAMBDA(int i, int k, Real &max) {
                Real Val = Kokkos::abs(Tend(i, k + 1));
                if (Val > max)
                   max = Val;
             },
             Kokkos::Max<Real>(MaxValue));
         Real SumValue = 0.0_Real;
         parallelReduce(
             {NEdgesAll, NVertLayers - 2},
             KOKKOS_LAMBDA(int i, int k, Real &LSum) {
                LSum += Tend(i, k + 1) * Tend(i, k + 1);
             },
             Kokkos::Sum<Real>(SumValue));
         Real RmseVal = std::sqrt(SumValue / (NEdgesAll * (NVertLayers - 2)));
         Rmse(Refinement) = RmseVal;

         LOG_INFO("refinement level {}: max |Tend| = {}, average Tend = {}",
                  Refinement, MaxValue, RmseVal);

         // coarsen for next iteration
         DC          = DC * 2.0_Real;
         NVertLayers = NVertLayers / 2;

      } // refinement loop

      // Test for second order convergence
      // resolution (dC) increases in refimenent loop
      if (Rmse(0) < Rmse(NRefinements - 1) / pow(4.0_Real, NRefinements - 1)) {
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
