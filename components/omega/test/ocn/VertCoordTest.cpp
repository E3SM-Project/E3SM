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
#include "Field.h"
#include "GlobalConstants.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanTestCommon.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "mpi.h"

#include <algorithm>
#include <random>

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

   // First step of time stepper initialization needed for IOstream and clock
   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize streams
   Field::init(ModelClock);
   IOStream::init(ModelClock);

   // Initialize the default halo
   Err = Halo::init();
   if (Err != 0)
      LOG_ERROR("VertCoordTest: error initializing default halo");

   // Initialize the default mesh
   HorzMesh::init(ModelClock);

   // Initialize the default vertical coordinate
   VertCoord::init();

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
      auto *DefHalo      = Halo::getDefault();
      auto *DefDecomp    = Decomp::getDefault();

      I4 NCellsSize     = DefMesh->NCellsSize;
      I4 NCellsOwned    = DefMesh->NCellsOwned;
      I4 NCellsAll      = DefMesh->NCellsAll;
      I4 NEdgesOwned    = DefMesh->NEdgesOwned;
      I4 NEdgesAll      = DefMesh->NEdgesAll;
      I4 NVerticesOwned = DefMesh->NVerticesOwned;
      I4 NVerticesAll   = DefMesh->NVerticesAll;
      I4 VertexDegree   = DefMesh->VertexDegree;
      I4 NVertLayers    = DefVertCoord->NVertLayers;

      // Tolerances for the comparisons below. In double precision these
      // reproduce the fixed 1e-10 absolute tolerance this test has always
      // used. Single precision cannot hold that: the expected values reach a
      // few thousand once cell and layer indices accumulate down a column, and
      // several of them (the sea surface height, the target thicknesses) are
      // arrived at by cancelling contributions from every layer, so the error
      // is set by the magnitude of the whole column rather than by the size of
      // the result. Hence a relative tolerance with an absolute floor. Both
      // stay far below the smallest separation any of these checks has to
      // resolve, which is the 0.5 between a layer midpoint and its interface.
      const Real RTol = sizeof(Real) == 4 ? 1.0e-4_Real : 0.0_Real;
      const Real ATol = sizeof(Real) == 4 ? 1.0e-2_Real : 1.0e-10_Real;

      // Rest bottom depth successful read
      R8 MaxBathy = -1e10;
      R8 MinBathy = 1e10;

      for (int ICell = 0; ICell < NCellsOwned; ++ICell) {
         if (DefVertCoord->BottomGeomDepthH(ICell) < MinBathy) {
            MinBathy = DefVertCoord->BottomGeomDepthH(ICell);
         }
         if (DefVertCoord->BottomGeomDepthH(ICell) > MaxBathy) {
            MaxBathy = DefVertCoord->BottomGeomDepthH(ICell);
         }
      }

      if ((MinBathy > 0) && (MaxBathy < 11000.)) {
         LOG_INFO("VertCoordTest: Bathy min/max test PASS");
      } else {
         ErrAll +=
             Error(ErrorCode::Fail, "VertCoordTest: Bathy min/max test FAIL");
      }

      // SurfacePressure is owned by VertCoord and its host mirror is sized on
      // cells. Its value is read from the InitialState stream (exercised in
      // StateTest); here we just verify the host mirror allocation.
      if (DefVertCoord->SurfacePressureH.extent_int(0) == NCellsSize) {
         LOG_INFO("VertCoordTest: SurfacePressureH size PASS");
      } else {
         ErrAll += Error(ErrorCode::Fail,
                         "VertCoordTest: SurfacePressureH size FAIL");
      }

      // Tests for computePressure

      Array2DReal PseudoThickness("PseudoThickness", NCellsSize, NVertLayers);
      Array1DReal SurfacePressure("SurfacePressure", NCellsSize);

      /// Initialize pseudo-thickness and surface pressure so that resulting
      /// interface pressure is the number of layers above plus one
      Real Rho0 = RhoSw;
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             SurfacePressure(ICell) = 1.0_Real;
             for (int K = 0; K < NVertLayers; K++) {
                PseudoThickness(ICell, K) = 1.0_Real / (Gravity * Rho0);
             }
          });
      Kokkos::fence();

      /// Call function and get host copies of outputs
      DefVertCoord->computePressure(PseudoThickness, SurfacePressure);
      auto PressInterfH = createHostMirrorCopy(DefVertCoord->PressureInterface);
      auto PressMidH    = createHostMirrorCopy(DefVertCoord->PressureMid);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            // Interface pressure at layer K should be K+1
            Real Expected = K + 1;
            if (!isApprox(PressInterfH(ICell, K), Expected, RTol, ATol)) {
               Err += 1;
            }
            // Mid pressures at layer K should be K+1.5
            Expected = K + 1.5;
            if (!isApprox(PressMidH(ICell, K), Expected, RTol, ATol)) {
               Err += 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computePressure with uniform PseudoThickness "
                  "PASS");
      } else {
         ErrAll += Error(ErrorCode::Fail, "VertCoordTest: computePressure with "
                                          "uniform PseudoThickness FAIL");
      }

      /// Initialize pseudo-thickness and surface pressure so that the resulting
      /// interface pressure is (K+1)*K/2 + the cell number
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             SurfacePressure(ICell) = 1.0_Real * ICell;
             for (int K = 0; K < NVertLayers; K++) {
                PseudoThickness(ICell, K) = (K + 1.0_Real) / (Gravity * Rho0);
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computePressure(PseudoThickness, SurfacePressure);
      auto PressInterfH2 =
          createHostMirrorCopy(DefVertCoord->PressureInterface);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// Interface pressure should be (K+1)*K/2 + the cell number
            Real Expected = ((K + 1.0_Real) * K) / 2.0_Real + ICell;
            if (!isApprox(PressInterfH2(ICell, K), Expected, RTol, ATol)) {
               Err += 1;
            }
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computePressure with non-uniform "
                  "PseudoThickness PASS");
      } else {
         ErrAll += Error(
             ErrorCode::Fail,
             "VertCoordTest: computePressure with non-uniform PseudoThickness "
             "FAIL");
      }

      // Tests for computeGeomZHeight

      Array2DReal SpecVol("SpecVol", NCellsSize, NVertLayers);
      Array1DReal BottomGeomDepth("BottomGeomDepth", NCellsSize);
      Array1DReal MaxLyrCellReal("MaxLyrCellReal", NCellsSize);
      deepCopy(MaxLyrCellReal, DefVertCoord->MaxLayerCell);

      auto &BotDepth = DefVertCoord->BottomGeomDepth;

      /// Initialize bottom depth, pseudo-thickness and specific volume so that
      /// the resulting interface z value is the negative layer number
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             BotDepth(ICell) = MaxLyrCellReal(ICell) + 1.0_Real;
             for (int K = 0; K < NVertLayers; K++) {
                PseudoThickness(ICell, K) = (ICell + 1.0_Real) / Rho0;
                SpecVol(ICell, K)         = 1.0_Real / (ICell + 1.0_Real);
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computeGeomZHeight(PseudoThickness, SpecVol);
      auto GeomZInterfH = createHostMirrorCopy(DefVertCoord->GeomZInterface);
      auto GeomZMidH    = createHostMirrorCopy(DefVertCoord->GeomZMid);
      auto SshCellH     = createHostMirrorCopy(DefVertCoord->SshCell);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// Z value at interface K should be -K
            Real Expected = -K;
            if (!isApprox(GeomZInterfH(ICell, K), Expected, RTol, ATol)) {
               Err += 1;
            }
            /// Z value at mid point of layer K should be -(K + .5)
            Expected = -K - 0.5;
            if (!isApprox(GeomZMidH(ICell, K), Expected, RTol, ATol)) {
               Err += 1;
            }
         }
         /// SshCell should equal ZInterface at the top of the active column
         Real Expected = -DefVertCoord->MinLayerCellH(ICell);
         if (!isApprox(SshCellH(ICell), Expected, RTol, ATol)) {
            Err += 1;
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computeGeomZHeight with uniform "
                  "PseudoThickness PASS");
      } else {
         ErrAll += Error(ErrorCode::Fail, "VertCoordTest: computeGeomZHeight "
                                          "with uniform PseudoThickness FAIL");
      }

      /// Initialize bottom depth, pseudo-thickness and specific volume so that
      /// the resulting interface z value is -(K+1)*K/2
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             BotDepth(ICell) = (MaxLyrCellReal(ICell) + 2) *
                               (MaxLyrCellReal(ICell) + 1) / 2.0_Real;
             for (int K = 0; K < NVertLayers; K++) {
                PseudoThickness(ICell, K) = (K + 1) / Rho0;
                SpecVol(ICell, K)         = 1.0_Real;
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computeGeomZHeight(PseudoThickness, SpecVol);
      auto ZInterfH2 = createHostMirrorCopy(DefVertCoord->GeomZInterface);
      auto SshCellH2 = createHostMirrorCopy(DefVertCoord->SshCell);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// Z value at interface should be -(K+1)*K/2
            Real Expected = -((K + 1.0_Real) * K) / 2.0_Real;
            if (!isApprox(ZInterfH2(ICell, K), Expected, RTol, ATol)) {
               Err += 1;
            }
         }
         /// SshCell should equal ZInterface at the top of the active column,
         /// which is -((MinLayer+1)*MinLayer)/2
         I4 MinK       = DefVertCoord->MinLayerCellH(ICell);
         Real Expected = -((MinK + 1.0_Real) * MinK) / 2.0_Real;
         if (!isApprox(SshCellH2(ICell), Expected, RTol, ATol)) {
            Err += 1;
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: computeGeomZHeight with non-uniform "
                  "PseudoThickness PASS");
      } else {
         ErrAll += Error(ErrorCode::Fail, "VertCoordTest: computeGeomZHeight "
                                          "with non-uniform PseudoThickness "
                                          "FAIL");
      }

      // Tests for computeGeopotential
      Array1DReal TidalPotential("TidalPotential", NCellsSize);
      Array1DReal SelfAttractionLoading("SelfAttractionLoading", NCellsSize);

      auto &GeomZMid = DefVertCoord->GeomZMid;

      /// Initialize z mid, tidal potential and SAL so that the resulting
      /// geopotential is the cell number + layer number
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             TidalPotential(ICell)        = ICell;
             SelfAttractionLoading(ICell) = -ICell;
             for (int K = 0; K < NVertLayers; K++) {
                GeomZMid(ICell, K) = (ICell + K) / Gravity;
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
            if (!isApprox(GeopotentialMidH(ICell, K), Expected, RTol, ATol)) {
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
      auto &RefPseudoThick = DefVertCoord->RefPseudoThickness;

      /// Initialize surface pressure, vertical coord weights, ref layer
      /// thickness, and pseudo-thickness so that the resulting target thickness
      /// is 2 (perturbation is evenly distributed amoung layers)
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             SurfacePressure(ICell) = 0.0;
             for (int K = 0; K < NVertLayers; K++) {
                RefPseudoThick(ICell, K)  = 1.0;
                PseudoThickness(ICell, K) = 2.0;
             }
          });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computePressure(PseudoThickness, SurfacePressure);
      DefVertCoord->computeTargetThickness();
      auto PseudoThicknessTargetH =
          createHostMirrorCopy(DefVertCoord->PseudoThicknessTarget);

      /// Check results
      Err = 0;
      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         for (int K = DefVertCoord->MinLayerCellH(ICell);
              K < DefVertCoord->MaxLayerCellH(ICell) + 1; K++) {
            /// target thickness should be 2
            Real Expected = 2.0;
            if (!isApprox(PseudoThicknessTargetH(ICell, K), Expected, RTol,
                          ATol)) {
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
      /// thickness, and pseudo-thickness so that the resulting target thickness
      /// is the max number of layers + 2 in the top layer and 1 elsewhere
      /// (perturbation is distributed to top layer only)
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             SurfacePressure(ICell) = 0.0;
             for (int K = 0; K < NVertLayers; K++) {
                RefPseudoThick(ICell, K)  = 1.0;
                PseudoThickness(ICell, K) = 2.0;
             }
          });

      auto &MovementWgts = DefVertCoord->VertCoordMovementWeights;
      parallelFor(
          {NVertLayers}, KOKKOS_LAMBDA(int K) { MovementWgts(K) = 0.0; });
      parallelFor({1}, KOKKOS_LAMBDA(const int &) { MovementWgts(0) = 1.0; });
      Kokkos::fence();

      /// Call functions and get host copy of output
      DefVertCoord->computePressure(PseudoThickness, SurfacePressure);
      DefVertCoord->computeTargetThickness();
      auto PseudoThicknessTargetH2 =
          createHostMirrorCopy(DefVertCoord->PseudoThicknessTarget);
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
            if (!isApprox(PseudoThicknessTargetH2(ICell, K), Expected, RTol,
                          ATol)) {
               LOG_INFO("PseudoThicknessTargetH({},{}) = {}, {}", ICell, K,
                        PseudoThicknessTargetH2(ICell, K), Expected);
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
      const auto &LocCellID       = DefDecomp->CellID;
      parallelFor(
          {NCellsAll}, KOKKOS_LAMBDA(int ICell) {
             LocMinLayerCell(ICell) = -2 * LocCellID(ICell);
             LocMaxLayerCell(ICell) = 2 * LocCellID(ICell);
          });
      Kokkos::fence();

      /// Call function, outputs are member variables of class
      DefVertCoord->minMaxLayerEdge(DefHalo);

      /// Check results
      Err = 0;
      for (int IEdge = 0; IEdge < NEdgesAll; IEdge++) {
         I4 Expected;

         /// Skip edges on boundary
         if ((DefMesh->CellsOnEdgeH(IEdge, 1) == NCellsAll) ||
             (DefMesh->CellsOnEdgeH(IEdge, 0) == NCellsAll)) {
            continue;
         }

         I4 CellID1 = DefDecomp->CellIDH(DefMesh->CellsOnEdgeH(IEdge, 0));
         I4 CellID2 = DefDecomp->CellIDH(DefMesh->CellsOnEdgeH(IEdge, 1));
         /// MinLayerEdgeTop is the min of the min cell values on edge
         Expected = std::min(-2 * CellID1, -2 * CellID2);
         if (!isApprox(DefVertCoord->MinLayerEdgeTopH(IEdge), Expected, RTol,
                       ATol)) {
            Err += 1;
         }
         /// MinLayerEdgeBot is the max of the min cell values on edge
         Expected = std::max(-2 * CellID1, -2 * CellID2);
         if (!isApprox(DefVertCoord->MinLayerEdgeBotH(IEdge), Expected, RTol,
                       ATol)) {
            Err += 1;
         }
         /// MaxLayerEdgeTop is the min of the max cell values on edge
         Expected = std::min(2 * CellID1, 2 * CellID2);
         if (!isApprox(DefVertCoord->MaxLayerEdgeTopH(IEdge), Expected, RTol,
                       ATol)) {
            Err += 1;
         }
         /// MaxLayerEdgeBot is the max of the max cell values on edge
         Expected = std::max(2 * CellID1, 2 * CellID2);
         if (!isApprox(DefVertCoord->MaxLayerEdgeBotH(IEdge), Expected, RTol,
                       ATol)) {
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
      DefVertCoord->minMaxLayerVertex(DefHalo);

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

         std::vector<I4> CellIDs = {
             DefDecomp->CellIDH(DefMesh->CellsOnVertexH(IVertex, 0)),
             DefDecomp->CellIDH(DefMesh->CellsOnVertexH(IVertex, 1)),
             DefDecomp->CellIDH(DefMesh->CellsOnVertexH(IVertex, 2))};
         /// MinLayerVertexTop is the min of the min cell values on vertex
         I4 Expected = 1e7;
         for (int I = 0; I < VertexDegree; I++) {
            Expected = std::min(Expected, -2 * CellIDs[I]);
         }
         if (!isApprox(DefVertCoord->MinLayerVertexTopH(IVertex), Expected,
                       RTol, ATol)) {
            Err += 1;
         }

         /// MinLayerVertexBot is the max of the min cell values on vertex
         Expected = -1e7;
         for (int I = 0; I < VertexDegree; I++) {
            Expected = std::max(Expected, -2 * CellIDs[I]);
         }
         if (!isApprox(DefVertCoord->MinLayerVertexBotH(IVertex), Expected,
                       RTol, ATol)) {
            Err += 1;
         }

         /// MaxLayerVertexTop is the min of the max cell values on vertex
         Expected = 1e7;
         for (int I = 0; I < VertexDegree; I++) {
            Expected = std::min(Expected, 2 * CellIDs[I]);
         }
         if (!isApprox(DefVertCoord->MaxLayerVertexTopH(IVertex), Expected,
                       RTol, ATol)) {
            Err += 1;
         }

         /// MaxLayerVertexBot is the max of the max cell values on vertex
         Expected = -1e7;
         for (int I = 0; I < VertexDegree; I++) {
            Expected = std::max(Expected, 2 * CellIDs[I]);
         }
         if (!isApprox(DefVertCoord->MaxLayerVertexBotH(IVertex), Expected,
                       RTol, ATol)) {
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

      // Tests for Masks

      /// Setup random values for MinLayerCell and MaxLayerCell near the top
      /// and bottom of each column
      I4 Seed = 12345;
      std::mt19937 Generator(Seed);

      I4 MinMin = 0;
      I4 MaxMin = 10;
      I4 MinMax = NVertLayers - 11;
      I4 MaxMax = NVertLayers - 1;

      std::uniform_int_distribution<int> MinDist(MinMin, MaxMin);
      std::uniform_int_distribution<int> MaxDist(MinMax, MaxMax);

      for (int ICell = 0; ICell < NCellsAll; ICell++) {
         DefVertCoord->MinLayerCellH(ICell) = MinDist(Generator);
         DefVertCoord->MaxLayerCellH(ICell) = MaxDist(Generator);
      }

      deepCopy(DefVertCoord->MinLayerCell, DefVertCoord->MinLayerCellH);
      deepCopy(DefVertCoord->MaxLayerCell, DefVertCoord->MaxLayerCellH);
      Kokkos::fence();

      /// Reset min/max layer arrays and set masks
      DefVertCoord->minMaxLayerEdge(DefHalo);
      DefVertCoord->minMaxLayerVertex(DefHalo);
      DefVertCoord->setMasks();

      Err = 0;
      for (int ICell = 0; ICell < NCellsOwned; ++ICell) {
         /// The sum of the masks in a column is the number of active
         /// layers in that column
         Real Expected = DefVertCoord->MaxLayerCellH(ICell) -
                         DefVertCoord->MinLayerCellH(ICell) + 1._Real;
         Real Sum = 0.;
         for (int K = 0; K < NVertLayers; ++K) {
            Sum += DefVertCoord->CellMaskH(ICell, K);
         }

         if (!isApprox(Sum, Expected, RTol, ATol)) {
            Err += 1;
         }
      }

      for (int IEdge = 0; IEdge < NEdgesOwned; ++IEdge) {
         Real Expected;

         /// The sum of the masks along an edge is the number of layers with
         /// active cells on both sides.
         if ((DefMesh->CellsOnEdgeH(IEdge, 1) == NCellsAll) ||
             (DefMesh->CellsOnEdgeH(IEdge, 0) == NCellsAll)) {
            Expected = 0.;
         } else {
            Expected = DefVertCoord->MaxLayerEdgeTopH(IEdge) -
                       DefVertCoord->MinLayerEdgeBotH(IEdge) + 1._Real;
         }
         Real Sum = 0.;
         for (int K = 0; K < NVertLayers; ++K) {
            Sum += DefVertCoord->EdgeMaskH(IEdge, K);
         }
         if (!isApprox(Sum, Expected, RTol, ATol)) {
            Err += 1;
         }
      }

      for (int IVertex = 0; IVertex < NVerticesOwned; ++IVertex) {
         Real Expected = DefVertCoord->MaxLayerVertexBotH(IVertex) -
                         DefVertCoord->MinLayerVertexTopH(IVertex) + 1._Real;

         Real Sum = 0.;
         for (int K = 0; K < NVertLayers; ++K) {
            Sum += DefVertCoord->VertexMaskH(IVertex, K);
         }
         if (!isApprox(Sum, Expected, RTol, ATol)) {
            Err += 1;
         }
      }

      /// Determine test pass/fail
      if (Err == 0) {
         LOG_INFO("VertCoordTest: setMasks PASS");
      } else {
         ErrAll += Error(ErrorCode::Fail, "VertCoordTest: setMasks FAIL");
      }

      // Finalize Omega objects
      IOStream::finalize();
      TimeStepper::clear();
      HorzMesh::clear();
      VertCoord::clear();
      Field::clear();
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
