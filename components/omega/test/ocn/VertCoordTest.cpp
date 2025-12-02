//===-- Test driver for OMEGA Vertical Coordinate ----------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA VertCoord class
///
///
//
//===-----------------------------------------------------------------------===/

#include "VertCoord.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "mpi.h"

#include <algorithm>

using namespace OMEGA;

int initVertCoordTest() {

   int Err = 0;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Open config file
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
   Err = Halo::init();
   if (Err != 0)
      LOG_ERROR("VertCoordTest: error initializing default halo");

   // Begin initialization of the default vertical coordinate
   VertCoord::init1();

   // Initialize the default mesh
   HorzMesh::init();

   // Complete initialization of the default vertical coordinate
   VertCoord::init2();

   return Err;
} // end initVertCoordTest

//------------------------------------------------------------------------------
// The test driver for VertCoord test
//
int main(int argc, char *argv[]) {

   // Initialize error code
   Error ErrAll;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      int Err = initVertCoordTest();
      if (Err != 0)
         ABORT_ERROR("VertCoordTest: Error initializing");

      auto *DefVertCoord = VertCoord::getDefault();
      auto *DefMesh      = HorzMesh::getDefault();

      I4 NCellsSize   = DefMesh->NCellsSize;
      I4 NCellsAll    = DefMesh->NCellsAll;
      I4 NEdgesAll    = DefMesh->NEdgesAll;
      I4 NVerticesAll = DefMesh->NVerticesAll;
      I4 VertexDegree = DefMesh->VertexDegree;
      I4 NVertLayers  = DefVertCoord->NVertLayers;

      // Tests for computePressure

      Array2DReal LayerThickness("LayerThickness", NCellsSize, NVertLayers);
      Array1DReal SurfacePressure("SurfacePressure", NCellsSize);

      /// Initialize layer thickness and surface pressure so that resulting
      /// interface pressure is the number of layers above plus one
      Real Rho0 = DefVertCoord->Rho0;
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             SurfacePressure(ICell) = 1.0_Real;
             for (int K = 0; K < NVertLayers; K++) {
                LayerThickness(ICell, K) = 1.0_Real / (Gravity * Rho0);
             }
          });
      Kokkos::fence();

      /// Call function and get host copies of outputs
      DefVertCoord->computePressure(LayerThickness, SurfacePressure);
      auto PressInterfH = createHostMirrorCopy(DefVertCoord->PressureInterface);
      auto PressMidH    = createHostMirrorCopy(DefVertCoord->PressureMid);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            // Interface pressure at layer K should be K+1
            Real Expected = K + 1;
            Real Diff     = std::abs(PressInterfH(ICell, K) - Expected);
            if (Diff > 1e-10) {
               Err += 1;
            }
            // Mid pressures at layer K should be K+1.5
            Expected = K + 1.5;
            Diff     = std::abs(PressMidH(ICell, K) - Expected);
            if (Diff > 1e-10) {
               Err += Err + 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO(
             "VertCoordTest: computePressure with uniform LayerThickness PASS");
      } else {
         ErrAll += Error(
             ErrorCode::Fail,
             "VertCoordTest: computePressure with uniform LayerThickness FAIL");
      }

      /// Initialize layer thickness and surface pressure so that the resulting
      /// interface pressure is (K+1)*K/2 + the cell number
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             SurfacePressure(ICell) = 1.0_Real * ICell;
             for (int K = 0; K < NVertLayers; K++) {
                LayerThickness(ICell, K) = (K + 1.0_Real) / (Gravity * Rho0);
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computePressure(LayerThickness, SurfacePressure);
      auto PressInterfH2 =
          createHostMirrorCopy(DefVertCoord->PressureInterface);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// Interface pressure should be (K+1)*K/2 + the cell number
            Real Expected = ((K + 1.0_Real) * K) / 2.0_Real + ICell;
            Real Diff     = std::abs(PressInterfH2(ICell, K) - Expected);
            if (Diff > 1e-10) {
               Err += 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computePressure with non-uniform "
                  "LayerThickness PASS");
      } else {
         ErrAll += Error(
             ErrorCode::Fail,
             "VertCoordTest: computePressure with non-uniform LayerThickness "
             "FAIL");
      }

      // Tests for computeZHeight

      Array2DReal SpecVol("SpecVol", NCellsSize, NVertLayers);
      Array1DReal BottomDepth("BottomDepth", NCellsSize);
      Array1DReal MaxLyrCellReal("MaxLyrCellReal", NCellsSize);
      deepCopy(MaxLyrCellReal, DefVertCoord->MaxLayerCell);

      auto &BotDepth = DefVertCoord->BottomDepth;

      /// Initialize bottom depth, layer thickness and specific volume so that
      /// the resulting interface z value is the negative layer number
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             BotDepth(ICell) = MaxLyrCellReal(ICell) + 1.0_Real;
             for (int K = 0; K < NVertLayers; K++) {
                LayerThickness(ICell, K) = (ICell + 1.0_Real) / Rho0;
                SpecVol(ICell, K)        = 1.0_Real / (ICell + 1.0_Real);
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computeZHeight(LayerThickness, SpecVol);
      auto ZInterfH = createHostMirrorCopy(DefVertCoord->ZInterface);
      auto ZMidH    = createHostMirrorCopy(DefVertCoord->ZMid);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// Z value at interface K should be -K
            Real Expected = -K;
            Real Diff     = std::abs(ZInterfH(ICell, K) - Expected);
            if (Diff > 1e-10) {
               Err += 1;
            }
            /// Z value at mid point of layer K should be -(K + .5)
            Expected = -K - 0.5;
            Diff     = std::abs(ZMidH(ICell, K) - Expected);
            if (Diff > 1e-10) {
               Err += 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO(
             "VertCoordTest: computeZHeight with uniform LayerThickness PASS");
      } else {
         ErrAll += Error(
             ErrorCode::Fail,
             "VertCoordTest: computeZHeight with uniform LayerThickness FAIL");
      }

      /// Initialize bottom depth, layer thickness and specific volume so that
      /// the resulting interface z value is -(K+1)*K/2
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             BotDepth(ICell) = (MaxLyrCellReal(ICell) + 2) *
                               (MaxLyrCellReal(ICell) + 1) / 2.0_Real;
             for (int K = 0; K < NVertLayers; K++) {
                LayerThickness(ICell, K) = (K + 1) / Rho0;
                SpecVol(ICell, K)        = 1.0_Real;
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computeZHeight(LayerThickness, SpecVol);
      auto ZInterfH2 = createHostMirrorCopy(DefVertCoord->ZInterface);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// Z value at interface should be -(K+1)*K/2
            Real Expected = -((K + 1.0_Real) * K) / 2.0_Real;
            Real Diff     = std::abs(ZInterfH2(ICell, K) - Expected);
            if (Diff > 1e-10) {
               Err += 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computeZHeight with non-uniform "
                  "LayerThickness PASS");
      } else {
         ErrAll += Error(
             ErrorCode::Fail,
             "VertCoordTest: computeZHeight with non-uniform LayerThickness "
             "FAIL");
      }

      // Tests for computeGeopotential
      Array1DReal TidalPotential("TidalPotential", NCellsSize);
      Array1DReal SelfAttractionLoading("SelfAttractionLoading", NCellsSize);

      auto &ZMid = DefVertCoord->ZMid;

      /// Initialize z mid, tidal potential and SAL so that the resulting
      /// geopotential is the cell number + layer number
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             TidalPotential(ICell)        = ICell;
             SelfAttractionLoading(ICell) = -ICell;
             for (int K = 0; K < NVertLayers; K++) {
                ZMid(ICell, K) = (ICell + K) / Gravity;
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computeGeopotential(TidalPotential, SelfAttractionLoading);
      auto GeopotentialMidH =
          createHostMirrorCopy(DefVertCoord->GeopotentialMid);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// Geopotential should be cell number + layer number
            Real Expected = ICell + K;
            Real Diff     = std::abs(GeopotentialMidH(ICell, K) - Expected);
            if (Diff > 1e-10) {
               Err += 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computeGeopotential PASS");
      } else {
         ErrAll +=
             Error(ErrorCode::Fail, "VertCoordTest: computeGeopotential FAIL");
      }

      // Tests for computeTargetThickness
      auto &RefLayerThick = DefVertCoord->RefLayerThickness;

      /// Initialize surface pressure, vertical coord weights, ref layer
      /// thickness, and layer thickness so that the resulting target thickness
      /// is 2 (perturbation is evenly distributed amoung layers)
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             SurfacePressure(ICell) = 0.0;
             for (int K = 0; K < NVertLayers; K++) {
                RefLayerThick(ICell, K)  = 1.0;
                LayerThickness(ICell, K) = 2.0;
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computePressure(LayerThickness, SurfacePressure);
      DefVertCoord->computeTargetThickness();
      auto LayerThicknessTargetH =
          createHostMirrorCopy(DefVertCoord->LayerThicknessTarget);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// target thickness should be 2
            Real Expected = 2.0;
            Real Diff = std::abs(LayerThicknessTargetH(ICell, K) - Expected);
            if (Diff > 1e-10) {
               Err += 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computeTargetThickness with uniform "
                  "distribution PASS");
      } else {
         ErrAll += Error(
             ErrorCode::Fail,
             "VertCoordTest: computeTargetThickness with uniform distribution "
             "FAIL");
      }

      /// Intialize surface pressure, vertical coord weights, ref layer
      /// thickness, and layer thickness so that the resulting target thickness
      /// is the max number of layers + 2 in the top layer and 1 elsewhere
      /// (perturbation is distributed to top layer only)
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             SurfacePressure(ICell) = 0.0;
             for (int K = 0; K < NVertLayers; K++) {
                RefLayerThick(ICell, K)  = 1.0;
                LayerThickness(ICell, K) = 2.0;
             }
          });

      auto &MovementWgts = DefVertCoord->VertCoordMovementWeights;
      parallelFor(
          {NVertLayers}, KOKKOS_LAMBDA(int K) { MovementWgts(K) = 0.0; });
      parallelFor({1}, KOKKOS_LAMBDA(const int &) { MovementWgts(0) = 1.0; });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computePressure(LayerThickness, SurfacePressure);
      DefVertCoord->computeTargetThickness();
      auto LayerThicknessTargetH2 =
          createHostMirrorCopy(DefVertCoord->LayerThicknessTarget);
      Err = 0;

      /// Check results
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            Real Expected;
            if (K == 0) {
               /// target thickness is number of layers + 2 in top layer
               Expected = DefVertCoord->MaxLayerCellH(ICell) + 2;
            } else {
               /// target thickness is 1 in all other layer
               Expected = 1.0;
            }
            Real Diff = std::abs(LayerThicknessTargetH2(ICell, K) - Expected);
            if (Diff > 1e-10) {
               LOG_INFO("LayerThicknessTargetH({},{}) = {}, {}", ICell, K,
                        LayerThicknessTargetH2(ICell, K), Expected);
               Err += 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computeTargetThickness with top only "
                  "distribution PASS");
      } else {
         ErrAll += Error(
             ErrorCode::Fail,
             "VertCoordTest: computeTargetThickness with top only distribution "
             "FAIL");
      }

      // Tests for minMaxLayerEdge

      /// Initialize min/max number of cell layers such that
      /// the cellsOnEdge information can be used to determine min/max layer
      /// edge
      const auto &LocMinLayerCell = DefVertCoord->MinLayerCell;
      const auto &LocMaxLayerCell = DefVertCoord->MaxLayerCell;
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             LocMinLayerCell(ICell) = -2 * ICell;
             LocMaxLayerCell(ICell) = 2 * ICell;
          });
      Kokkos::fence();

      /// Call function, outputs are member variables of class
      DefVertCoord->minMaxLayerEdge();

      /// Check results
      Err = 0;
      for (int IEdge = 0; IEdge < NEdgesAll; IEdge++) {
         I4 Expected;

         /// Skip edges on boundary
         if ((DefMesh->CellsOnEdgeH(IEdge, 1) == NCellsAll) ||
             (DefMesh->CellsOnEdgeH(IEdge, 0) == NCellsAll)) {
            continue;
         }

         /// MinLayerEdgeTop is the min of the min cell values on edge
         Expected  = std::min(-2 * DefMesh->CellsOnEdgeH(IEdge, 0),
                              -2 * DefMesh->CellsOnEdgeH(IEdge, 1));
         Real Diff = std::abs(DefVertCoord->MinLayerEdgeTopH(IEdge) - Expected);
         if (Diff > 1e-10) {
            Err += 1;
         }
         /// MinLayerEdgeBot is the max of the min cell values on edge
         Expected = std::max(-2 * DefMesh->CellsOnEdgeH(IEdge, 0),
                             -2 * DefMesh->CellsOnEdgeH(IEdge, 1));
         Diff     = std::abs(DefVertCoord->MinLayerEdgeBotH(IEdge) - Expected);
         if (Diff > 1e-10) {
            Err += 1;
         }
         /// MaxLayerEdgeTop is the min of the max cell values on edge
         Expected = std::min(2 * DefMesh->CellsOnEdgeH(IEdge, 0),
                             2 * DefMesh->CellsOnEdgeH(IEdge, 1));
         Diff     = std::abs(DefVertCoord->MaxLayerEdgeTopH(IEdge) - Expected);
         if (Diff > 1e-10) {
            Err += 1;
         }
         /// MaxLayerEdgeBot is the max of the max cell values on edge
         Expected = std::max(2 * DefMesh->CellsOnEdgeH(IEdge, 0),
                             2 * DefMesh->CellsOnEdgeH(IEdge, 1));
         Diff     = std::abs(DefVertCoord->MaxLayerEdgeBotH(IEdge) - Expected);
         if (Diff > 1e-10) {
            Err += 1;
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: minMaxLayerEdge PASS");
      } else {
         ErrAll +=
             Error(ErrorCode::Fail, "VertCoordTest: minMaxLayerEdge FAIL");
      }

      // Tests for minMaxLayerVertex

      /// Use MinLayerCell, MaxLayerCell values initialized in previous test
      /// CellsOnVertex information can be used to determine min/max layer
      /// vertex

      /// Call function, outputs are member variables of class
      DefVertCoord->minMaxLayerVertex();

      /// Check results
      Err = 0;
      for (int IVertex = 0; IVertex < NVerticesAll; IVertex++) {

         /// Skip vertices on boundary
         I4 Boundary = 0;
         for (int I = 0; I < VertexDegree; I++) {
            if (DefMesh->CellsOnVertexH(IVertex, I) == NCellsAll) {
               Boundary += 1;
            }
         }
         if (Boundary > 0) {
            continue;
         }

         /// MinLayerVertexTop is the min of the min cell values on vertex
         I4 Expected = 1e7;
         for (int I = 0; I < VertexDegree; I++) {
            Expected =
                std::min(Expected, -2 * DefMesh->CellsOnVertexH(IVertex, I));
         }
         Real Diff =
             std::abs(DefVertCoord->MinLayerVertexTopH(IVertex) - Expected);
         if (Diff > 1e-10) {
            Err += 1;
         }

         /// MinLayerVertexBot is the max of the min cell values on vertex
         Expected = -1e7;
         for (int I = 0; I < VertexDegree; I++) {
            Expected =
                std::max(Expected, -2 * DefMesh->CellsOnVertexH(IVertex, I));
         }
         Diff = std::abs(DefVertCoord->MinLayerVertexBotH(IVertex) - Expected);
         if (Diff > 1e-10) {
            Err += 1;
         }

         /// MaxLayerVertexTop is the min of the max cell values on vertex
         Expected = 1e7;
         for (int I = 0; I < VertexDegree; I++) {
            Expected =
                std::min(Expected, 2 * DefMesh->CellsOnVertexH(IVertex, I));
         }
         Diff = std::abs(DefVertCoord->MaxLayerVertexTopH(IVertex) - Expected);
         if (Diff > 1e-10) {
            Err += 1;
         }

         /// MaxLayerVertexBot is the max of the max cell values on vertex
         Expected = -1e7;
         for (int I = 0; I < VertexDegree; I++) {
            Expected =
                std::max(Expected, 2 * DefMesh->CellsOnVertexH(IVertex, I));
         }
         Diff = std::abs(DefVertCoord->MaxLayerVertexBotH(IVertex) - Expected);
         if (Diff > 1e-10) {
            Err += 1;
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: minMaxLayerVertex PASS");
      } else {
         ErrAll +=
             Error(ErrorCode::Fail, "VertCoordTest: minMaxLayerVertex FAIL");
      }

      // Finalize Omega objects
      IOStream::finalize();
      TimeStepper::clear();
      HorzMesh::clear();
      VertCoord::clear();
      Dimension::clear();
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   CHECK_ERROR_ABORT(ErrAll, "VertCoord unit tests FAIL");

   return 0;

} // end of main
//===-----------------------------------------------------------------------===/
