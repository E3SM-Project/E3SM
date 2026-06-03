//===-- Test driver for Omega fill value standardization --------*- C++ -*-===//
//
/// \file
/// \brief Test driver for Omega fill value constants and auto-fill behavior
///
/// Tests:
///   1. FillValueI4/I8/R4/R8 constants match the NetCDF-C NC_FILL_* values
///   2. Field::attachData() auto-fills the array with the declared fill value
///   3. Inactive layers (k >= MaxLayerCell) contain FillValueReal after
///      VertCoord initialization and compute
///   4. NormalVelocity has the correct 3-zone pattern after IC read and
///   masking:
///      fully-inactive layers = FillValueReal, boundary layers = 0, active
///      layers
///      != FillValueReal (quiescent IC gives 0, which is valid)
//
//===----------------------------------------------------------------------===//

#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "Error.h"
#include "Field.h"
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
#include "Tendencies.h"
#include "TimeMgr.h"
#include "TimeStepper.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "mpi.h"

#include <netcdf.h>

using namespace OMEGA;

//------------------------------------------------------------------------------
void logNormalVelocityFillDiagnostics(const char *Stage) {

   auto *DefEnv       = MachEnv::getDefault();
   auto *DefMesh      = HorzMesh::getDefault();
   auto *DefVertCoord = VertCoord::getDefault();
   auto *DefState     = OceanState::getDefault();

   I4 NEdgesOwned = DefMesh->NEdgesOwned;
   I4 NVertLayers = DefVertCoord->NVertLayers;

   auto NVH = createHostMirrorCopy(DefState->NormalVelocity[0]);
   const auto &MinLayerEdgeTopH = DefVertCoord->MinLayerEdgeTopH;
   const auto &MaxLayerEdgeBotH = DefVertCoord->MaxLayerEdgeBotH;
   const auto &MinLayerEdgeBotH = DefVertCoord->MinLayerEdgeBotH;
   const auto &MaxLayerEdgeTopH = DefVertCoord->MaxLayerEdgeTopH;

   int ErrActiveLocal = 0;
   int NLogged        = 0;
   for (int IEdge = 0; IEdge < NEdgesOwned; ++IEdge) {
      int KTop = MinLayerEdgeTopH(IEdge);
      int KBot = std::max(0, (int)MaxLayerEdgeBotH(IEdge));
      int KMin = MinLayerEdgeBotH(IEdge);
      int KMax = MaxLayerEdgeTopH(IEdge);

      for (int K = 0; K < NVertLayers; ++K) {
         Real Val    = NVH(IEdge, K);
         bool Active = (KMax >= 0 && K >= KMin && K <= KMax);

         if (Active && Val == FillValueReal) {
            ++ErrActiveLocal;
            if (NLogged < 5) {
               LOG_INFO("FillValueTest: {} active fill sample: IEdge={} K={} "
                        "KTop={} KBot={} KMin={} KMax={} Val={}",
                        Stage, IEdge, K, KTop, KBot, KMin, KMax, Val);
               ++NLogged;
            }
         }
      }
   }

   int ErrActiveGlobal = 0;
   MPI_Allreduce(&ErrActiveLocal, &ErrActiveGlobal, 1, MPI_INT, MPI_SUM,
                 DefEnv->getComm());
   LOG_INFO("FillValueTest: {} active NormalVelocity fill count: local={} "
            "global={}",
            Stage, ErrActiveLocal, ErrActiveGlobal);
}

//------------------------------------------------------------------------------
// Full initialization needed for OceanState (matches StateTest.cpp pattern)
void initFillValueTest() {

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);

   Config("Omega");
   Config::readAll("omega.yml");

   TimeStepper::init1();
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   IO::init(DefComm);

   Field::init(ModelClock);
   IOStream::init(ModelClock);

   Decomp::init();

   int Err = Halo::init();
   if (Err != 0)
      ABORT_ERROR("FillValueTest: error initializing default halo");

   HorzMesh::init(ModelClock);
   VertCoord::init();
   Tracers::init();
   AuxiliaryState::init();
   PressureGrad::init();
   Tendencies::init();
   TimeStepper::init2();

   Err = OceanState::init();
   if (Err != 0)
      ABORT_ERROR("FillValueTest: error initializing OceanState");

   bool StreamsValid = IOStream::validateAll();
   if (!StreamsValid)
      ABORT_ERROR("FillValueTest: error validating IO streams");
}

//------------------------------------------------------------------------------
int main(int argc, char *argv[]) {

   Error ErrAll;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      initFillValueTest();

      // Read initial state so NormalVelocity has real values (quiescent IC:
      // NV=0 for all active layers). Mirrors the pattern in ocnInit /
      // StateTest.
      {
         Clock *ModelClock = TimeStepper::getDefault()->getClock();
         Metadata ReqMeta;
         Error Err1 = IOStream::read("InitialState", ModelClock, ReqMeta);
         CHECK_ERROR_ABORT(Err1, "FillValueTest: error reading initial state");
         logNormalVelocityFillDiagnostics("after InitialState read");

         OceanState::getDefault()->exchangeHalo(0);
         logNormalVelocityFillDiagnostics("after NormalVelocity halo exchange");

         // Apply 3-zone mask after IC read: boundary layers -> 0,
         // inactive layers -> FillValueReal, active layers unchanged.
         VertCoord::getDefault()->applyEdgeLayerMask(
             OceanState::getDefault()->NormalVelocity[0],
             HorzMesh::getDefault()->NEdgesAll);
         logNormalVelocityFillDiagnostics("after edge layer mask");

         OceanState::getDefault()->copyToHost(0);
      }

      auto *DefMesh      = HorzMesh::getDefault();
      auto *DefVertCoord = VertCoord::getDefault();
      auto *DefState     = OceanState::getDefault();

      // ------------------------------------------------------------------
      // Test 1: FillValue constants match NetCDF-C NC_FILL_* values
      // ------------------------------------------------------------------
      {
         int Err = 0;
         if (FillValueI4 != static_cast<I4>(NC_FILL_INT))
            ++Err;
         if (FillValueI8 != static_cast<I8>(NC_FILL_INT64))
            ++Err;
         if (static_cast<double>(FillValueR4) !=
             static_cast<double>(NC_FILL_FLOAT))
            ++Err;
         if (FillValueR8 != static_cast<R8>(NC_FILL_DOUBLE))
            ++Err;

         if (Err == 0) {
            LOG_INFO("FillValueTest: fill constant values PASS");
         } else {
            ErrAll += Error(ErrorCode::Fail,
                            "FillValueTest: fill constant values FAIL");
         }
      }

      // ------------------------------------------------------------------
      // Test 2: attachData() auto-fills the array with the declared fill value
      // ------------------------------------------------------------------
      {
         // Create a small non-distributed test dimension and field
         I4 NTest = 10;
         Dimension::create("FVTestDim", NTest);
         auto TestField =
             Field::create("FVTestField", "fill value test", "1", "", -1.0e40,
                           1.0e40, 1, {"FVTestDim"}, false);

         // Allocate a host array set to 0.0 (not the fill value)
         HostArray1DR8 TestArr("FVTestArr", NTest);
         Kokkos::deep_copy(TestArr, 0.0);

         // attachData() should automatically fill with FillValueR8
         TestField->attachData<HostArray1DR8>(TestArr);

         int Err = 0;
         for (int I = 0; I < NTest; ++I) {
            if (TestArr(I) != FillValueR8)
               ++Err;
         }

         if (Err == 0) {
            LOG_INFO("FillValueTest: attachData auto-fill PASS");
         } else {
            ErrAll += Error(ErrorCode::Fail,
                            "FillValueTest: attachData auto-fill FAIL");
         }
      }

      // ------------------------------------------------------------------
      // Test 3: Inactive layers (k >= MaxLayerCell) contain FillValueReal
      //         in a VertCoord field (GeomZMid) after initialization
      // ------------------------------------------------------------------
      {
         I4 NCellsOwned = DefMesh->NCellsOwned;
         I4 NVertLayers = DefVertCoord->NVertLayers;

         auto GeomZMidH = createHostMirrorCopy(DefVertCoord->GeomZMid);

         int Err = 0;
         for (int ICell = 0; ICell < NCellsOwned; ++ICell) {
            I4 KMax = DefVertCoord->MaxLayerCellH(ICell);
            for (int K = KMax; K < NVertLayers; ++K) {
               if (GeomZMidH(ICell, K) != FillValueReal)
                  ++Err;
            }
         }

         if (Err == 0) {
            LOG_INFO("FillValueTest: inactive layers contain fill value PASS");
         } else {
            ErrAll +=
                Error(ErrorCode::Fail,
                      "FillValueTest: inactive layers contain fill value FAIL");
         }
      }

      // ------------------------------------------------------------------
      // Test 4: NormalVelocity has the correct 3-zone pattern after IC read
      //         and applyEdgeLayerMask():
      //   - Fully inactive layers (outside [MinLayerEdgeTop, MaxLayerEdgeBot]):
      //     FillValueReal
      //   - Boundary layers ([MinLayerEdgeTop, MaxLayerEdgeBot] but outside
      //     [MinLayerEdgeBot, MaxLayerEdgeTop]): 0
      //   - Active layers ([MinLayerEdgeBot, MaxLayerEdgeTop]):
      //     != FillValueReal (quiescent IC gives 0, which is valid)
      // ------------------------------------------------------------------
      {
         I4 NEdgesOwned = DefMesh->NEdgesOwned;
         I4 NVertLayers = DefVertCoord->NVertLayers;

         auto NVH = createHostMirrorCopy(DefState->NormalVelocity[0]);
         const auto &MinLayerEdgeTopH = DefVertCoord->MinLayerEdgeTopH;
         const auto &MaxLayerEdgeBotH = DefVertCoord->MaxLayerEdgeBotH;
         const auto &MinLayerEdgeBotH = DefVertCoord->MinLayerEdgeBotH;
         const auto &MaxLayerEdgeTopH = DefVertCoord->MaxLayerEdgeTopH;

         int ErrInactive = 0; // non-fill in fully-inactive zone (bad)
         int ErrBoundary = 0; // non-zero in boundary zone (bad)
         int ErrActive   = 0; // fill value in active zone (bad)

         for (int IEdge = 0; IEdge < NEdgesOwned; ++IEdge) {
            int KTop = MinLayerEdgeTopH(IEdge);
            int KBot = std::max(0, (int)MaxLayerEdgeBotH(IEdge));
            int KMin = MinLayerEdgeBotH(IEdge);
            int KMax = MaxLayerEdgeTopH(IEdge); // -1 for land-adjacent edges

            for (int K = 0; K < NVertLayers; ++K) {
               Real Val      = NVH(IEdge, K);
               bool valid    = (K >= KTop && K <= KBot);
               bool active   = (KMax >= 0 && K >= KMin && K <= KMax);
               bool boundary = valid && !active;

               if (!valid && Val != FillValueReal)
                  ++ErrInactive;
               if (boundary && Val != 0._Real)
                  ++ErrBoundary;
               if (active && Val == FillValueReal)
                  ++ErrActive;
            }
         }

         if (ErrInactive == 0 && ErrBoundary == 0 && ErrActive == 0) {
            LOG_INFO(
                "FillValueTest: NormalVelocity fill/zero/real pattern PASS");
         } else {
            ErrAll += Error(
                ErrorCode::Fail,
                "FillValueTest: NormalVelocity fill/zero/real pattern FAIL"
                " (ErrInactive={}, ErrBoundary={}, ErrActive={})",
                ErrInactive, ErrBoundary, ErrActive);
         }
      }

      // ------------------------------------------------------------------
      // Finalize
      // ------------------------------------------------------------------
      OceanState::clear();
      Tendencies::clear();
      PressureGrad::clear();
      AuxiliaryState::clear();
      Tracers::clear();
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
   int RetVal = ErrAll.isFail() ? 1 : 0;
   MPI_Finalize();

   return RetVal;
}
//===----------------------------------------------------------------------===//
