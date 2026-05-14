//===-- Test driver for OMEGA Vertical Advection  ----------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA VertAdv class
///
///
//
//===-----------------------------------------------------------------------===/

#include "VertAdv.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "GlobalConstants.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "TimeStepper.h"
#include "VertCoord.h"

#include "Pacer.h"
#include <cmath>
#include <mpi.h>
#include <vector>

using namespace OMEGA;

constexpr Real Alpha  = 0.197_Real; // random phase angle
constexpr Real Length = 1024._Real;

// Test function used for tendTest
KOKKOS_FUNCTION Real testFunc(const Real X) {
   Real Arg = Pi * (X + Alpha) / Length;
   return std::sin(Arg) * std::sin(Arg) * std::cos(Arg) * std::cos(Arg);
}

// Function used for analytical derivative in tendTest
KOKKOS_FUNCTION Real testDeriv(const Real X) {
   Real C = Pi / Length;
   return (C / 2._Real) * std::sin(4._Real * C * (X + Alpha));
}

// For a given resolution, use computeStdVAdvTend to compute tracer
// tendencies, then calculate L2 error
Real tendTest(const int NLayers, VertAdv *VAdv) {

   // Setup test for this resolution
   Real Delta = Length / NLayers;

   VAdv->NVertLayers = NLayers;
   VAdv->TotalVerticalVelocity =
       Array2DReal("TotalVerticaVelocity", 1, NLayers + 1);
   VAdv->VertFlux = Array3DReal("VertFlux", 1, 1, NLayers + 1);
   Array1DReal XLayer("XLayer", NLayers);
   OMEGA_SCOPE(LocTotVertVel, VAdv->TotalVerticalVelocity);
   Array2DReal LayerThick("LayerThickness", 1, NLayers);
   Array3DReal TracerArray("TracerArray", 1, 1, NLayers);
   Array3DReal Tend("TracerArray", 1, 1, NLayers);
   TimeInterval ZeroTimeStep;

   // Set uniform velocity throughout domain
   deepCopy(VAdv->TotalVerticalVelocity, 1._Real);
   // Set velocities at top and bottom interfaces to 0.
   parallelFor(
       {1}, KOKKOS_LAMBDA(const int &) {
          LocTotVertVel(0, 0)       = 0._Real;
          LocTotVertVel(0, NLayers) = 0._Real;
       });
   // Set X values, LayerThickness, and Tracer values throughout domain
   parallelFor(
       {NLayers}, KOKKOS_LAMBDA(int K) {
          XLayer(K)            = Delta * (K + 0.5_Real);
          LayerThick(0, K)     = Delta;
          TracerArray(0, 0, K) = testFunc(XLayer(K));
       });

   // Compute tracer tendencies
   VAdv->computeTracerVAdvTend(Tend, TracerArray, LayerThick, ZeroTimeStep);

   HostArray3DReal TendH   = createHostMirrorCopy(Tend);
   HostArray1DReal XLayerH = createHostMirrorCopy(XLayer);
   Real SumSq              = 0._Real;
   int Count               = 0;

   // Compute errors at interior of range, ignoring layers near top and bottom
   // where lower order is used
   for (int K = 2; K < NLayers - 2; ++K) {
      Real Expected = testDeriv(XLayerH(K));
      Real Diff     = TendH(0, 0, K) / Delta - Expected;
      SumSq += Diff * Diff;
      Count++;
   }
   Real L2 = std::sqrt(SumSq / Count);

   return L2;
}

// Initialize needed modules
void initVertAdvTest() {

   I4 Err;

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
      ABORT_ERROR("VertAdvTest: error initializing default halo");

   // Initialize the default mesh
   HorzMesh::init();

   // Initialize the default vertical coordinate
   VertCoord::init(false);

   // Initialize tracers
   Tracers::init();

   // Initialize the default vertical advection
   VertAdv::init();
}

// Clean-up modules
void finalizeVertAdvTest() {

   IOStream::finalize();
   Tracers::clear();
   VertAdv::clear();
   VertCoord::clear();
   TimeStepper::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

int main(int argc, char *argv[]) {

   // Initialize error code
   Error ErrAll;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      initVertAdvTest();

      auto DefMesh    = HorzMesh::getDefault();
      auto DefDecomp  = Decomp::getDefault();
      auto DefVertAdv = VertAdv::getDefault();
      auto DefVCoord  = VertCoord::getDefault();

      // Only need to test within one column or along one edge and with one
      // tracer
      DefVertAdv->NCellsOwned = 1;
      DefVertAdv->NCellsAll   = 1;
      DefVertAdv->NEdgesOwned = 1;
      DefVertAdv->NEdgesAll   = 1;
      DefVertAdv->NTracers    = 1;

      I4 NVertLayers   = DefVertAdv->NVertLayers;
      I4 NVertLayersP1 = DefVertAdv->NVertLayersP1;

      Array2DReal FluxLayerThickEdge("FluxLayerThickEdge",
                                     DefDecomp->NEdgesSize, NVertLayers);
      Array2DReal NormalVelEdge("NormalVelEdge", DefDecomp->NEdgesSize,
                                NVertLayers);
      Array2DReal LayerThickCell("LayerThickCell", DefDecomp->NCellsSize,
                                 NVertLayers);

      const I4 ICell0        = 0;
      const I4 NEdgesOnCell0 = DefDecomp->NEdgesOnCellH(ICell0);

      // Test for computeVerticalVelocity
      I4 Err = 0;

      OMEGA_SCOPE(LocEOnC, DefDecomp->EdgesOnCell);
      OMEGA_SCOPE(LocESOnC, DefMesh->EdgeSignOnCell);
      OMEGA_SCOPE(LocTargetThick, DefVCoord->LayerThicknessTarget);

      // Set uniform values for layer thickness and velocity on edges of column
      parallelFor(
          {NVertLayers}, KOKKOS_LAMBDA(int K) {
             for (int J = 0; J < NEdgesOnCell0; ++J) {
                const I4 JEdge               = LocEOnC(ICell0, J);
                FluxLayerThickEdge(JEdge, K) = 1._Real;
                NormalVelEdge(JEdge, K)      = 1._Real * LocESOnC(ICell0, J);
                LayerThickCell(ICell0, K)    = 1._Real;
                LocTargetThick(ICell0, K)    = 1._Real;
             }
          });

      Real ProjDt = 1._Real;
      // Compute vertical velocity
      DefVertAdv->computeVerticalVelocity(NormalVelEdge, FluxLayerThickEdge,
                                          LayerThickCell, ProjDt);

      HostArray2DReal VertVelH =
          createHostMirrorCopy(DefVertAdv->VerticalVelocity);

      Real Tol = 1e-10;

      // Divergence is equal in each layer. Expected vertical velocity through
      // each interface is sum over negative divergence starting from bottom
      // layer
      Real Perim0 = 0.;
      for (int J = 0; J < NEdgesOnCell0; ++J) {
         Perim0 += DefMesh->DvEdgeH(DefMesh->EdgesOnCellH(ICell0, J));
      }
      Real InvAreaCell0 = 1._Real / DefMesh->AreaCellH(ICell0);

      for (int K = 1; K < NVertLayersP1; ++K) {
         Real Expected = (NVertLayers - K) * Perim0 * InvAreaCell0;
         Real Diff     = std::abs(VertVelH(ICell0, K) - Expected);
         if (Diff > Tol) {
            ++Err;
         }
      }

      if (Err != 0) {
         ErrAll += Error(ErrorCode::Fail,
                         "VertAdvTest: computeVerticalVelocity FAIL");
      }

      // Test for computeThicknessVAdvTend
      Err = 0;

      Array2DReal TendCell2D("TendCell2D", DefMesh->NCellsSize, NVertLayers);
      // Compute thickness tendencies for each layer using velocities from
      // previous tests
      DefVertAdv->computeThicknessVAdvTend(TendCell2D);
      HostArray2DReal TendCell2DH = createHostMirrorCopy(TendCell2D);

      // Expected thickness tendency for each layer is difference between
      // velocity at bottom and top interface, i.e. the divergence
      for (int K = 1; K < NVertLayers; ++K) {
         Real Expected = -Perim0 * InvAreaCell0;
         Real Diff     = std::abs(TendCell2DH(ICell0, K) - Expected);
         if (Diff > Tol) {
            ++Err;
         }
      }

      if (Err != 0) {
         ErrAll += Error(ErrorCode::Fail,
                         "VertAdvTest: computeThicknessVAdvTend FAIL");
      }

      // Test for computeVelocityVAdvTend
      Err = 0;

      // Set uniform vertical velocity in each column along the edge and
      // uniform layer thickness at the edge. Horizontal velocity is set
      // to a periodic function
      Array2DReal TendEdge2D("TendEdge2D", DefMesh->NEdgesSize, NVertLayers);
      const I4 IEdge0 = 0;
      OMEGA_SCOPE(LocCOnE, DefDecomp->CellsOnEdge);
      OMEGA_SCOPE(LocVertVel, DefVertAdv->TotalVerticalVelocity);
      parallelFor(
          {NVertLayers}, KOKKOS_LAMBDA(int K) {
             NormalVelEdge(IEdge0, K)      = sin(Pi * K / NVertLayers);
             FluxLayerThickEdge(IEdge0, K) = 1._Real;
             const I4 Cell1                = LocCOnE(IEdge0, 0);
             const I4 Cell2                = LocCOnE(IEdge0, 1);
             LocVertVel(Cell1, K)          = 1._Real;
             LocVertVel(Cell2, K)          = 1._Real;
          });

      // Compute velocity tendencies for each layer
      DefVertAdv->computeVelocityVAdvTend(TendEdge2D, NormalVelEdge,
                                          FluxLayerThickEdge);

      Tol = 1._Real / (NVertLayers * NVertLayers);

      HostArray2DReal TendEdge2DH = createHostMirrorCopy(TendEdge2D);

      // Expected velocity tendency is derivative of initial distribution
      for (int K = 0; K < NVertLayers; ++K) {
         Real Expected = (Pi / NVertLayers) * cos(Pi * K / NVertLayers);
         if (K == 0 or K == NVertLayers - 1)
            Expected /= 2._Real;
         Real Diff = std::abs(TendEdge2DH(IEdge0, K) - Expected);
         if (Diff > Tol) {
            ++Err;
         }
      }

      if (Err != 0) {
         ErrAll += Error(ErrorCode::Fail,
                         "VertAdvTest: computeVelocityVAdvTend FAIL");
      }

      // Tests for computeTracerVAdvTend

      DefVertAdv->VertAdvChoice = VertAdvOption::Standard;

      // Setup for convergence tests with computeStdVAdvTend
      constexpr I4 NLayersArray[]          = {128, 256, 512, 1024};
      constexpr VertFluxOption AllOrders[] = {VertFluxOption::Second,
                                              VertFluxOption::Third,
                                              VertFluxOption::Fourth};
      static const std::map<VertFluxOption, I4> OrderOfAcc = {
          {VertFluxOption::Second, 2},
          {VertFluxOption::Third, 3},
          {VertFluxOption::Fourth, 4}};
      static const std::map<VertFluxOption, std::string> OrderStr = {
          {VertFluxOption::Second, "2nd"},
          {VertFluxOption::Third, "3rd"},
          {VertFluxOption::Fourth, "4th"}};

      // Set Coef3rdOrder to 1 so VertFluxOption::Third is calculated with
      // purely 3rd-order and not a blend of 3rd- and 4th-order fluxes
      DefVertAdv->Coef3rdOrder = 1.;

      Tol = 0.1_Real;
      OMEGA_SCOPE(LocMaxLyrCell, DefVCoord->MaxLayerCell);

      // Compute L2 error for computeStdVAdvTend with each VertFluxOption over a
      // set of resolutions. Compare successive errors for each Order to confirm
      // expected order of accuracy.
      for (VertFluxOption Order : AllOrders) {
         Err                        = 0;
         DefVertAdv->VertFluxChoice = Order;
         std::vector<Real> L2Errors;
         for (I4 NLayers : NLayersArray) {
            parallelFor(
                {1},
                KOKKOS_LAMBDA(const int &) { LocMaxLyrCell(0) = NLayers - 1; });
            Real L2Err = tendTest(NLayers, DefVertAdv);
            L2Errors.push_back(L2Err);
         }
         Real ExpectedRat = std::pow(2._Real, OrderOfAcc.at(Order));
         for (I4 I = 0; I < L2Errors.size() - 1; ++I) {
            Real Rat    = L2Errors[I] / L2Errors[I + 1];
            Real RelErr = std::abs(Rat - ExpectedRat) / ExpectedRat;
            if (RelErr > Tol) {
               ++Err;
            }
         }
         if (Err != 0) {
            ErrAll +=
                Error(ErrorCode::Fail,
                      "VertAdvTest: computeStdVAdvTend with {} order FAIL",
                      OrderStr.at(Order));
         }
      }

      // Monotonicity check for computeFCTVAdvTend

      NVertLayers = 256;
      DefVertAdv->LowOrderVertFlux =
          Array3DReal("LowOrderVertFlux", 1, 1, NVertLayers + 1);
      DefVertAdv->VertAdvChoice  = VertAdvOption::FCT;
      DefVertAdv->VertFluxChoice = VertFluxOption::Third;
      DefVertAdv->Coef3rdOrder   = 0.25;
      std::vector<Real> TestVel  = {1._Real, -1._Real};
      std::vector<I4> StartIdx   = {NVertLayers / 2, NVertLayers / 4};
      parallelFor(
          {1},
          KOKKOS_LAMBDA(const int &) { LocMaxLyrCell(0) = NVertLayers - 1; });

      // Test a top-hat function for a uniform velocity in both directions
      for (int IMono = 0; IMono < 2; ++IMono) {
         Err = 0;
         Array3DReal TracerArray("TracerArray", 1, 1, NVertLayers);
         Array2DReal LayerThick("LayerThick", 1, NVertLayers);
         Array3DReal TendCell3D("TendCell3D", 1, 1, NVertLayers);
         deepCopy(LayerThick, 1._Real);
         DefVertAdv->NVertLayers = NVertLayers;
         LocVertVel = Array2DReal("TotalVerticaVelocity", 1, NVertLayers + 1);
         deepCopy(LocVertVel, TestVel[IMono]);
         parallelFor(
             {1}, KOKKOS_LAMBDA(const int &) {
                LocVertVel(0, 0)           = 0._Real;
                LocVertVel(0, NVertLayers) = 1._Real;
             });
         deepCopy(TracerArray, 0._Real);
         const I4 KStart = StartIdx[IMono];
         const I4 KRange = NVertLayers / 4;
         parallelFor(
             {KRange}, KOKKOS_LAMBDA(const int K) {
                TracerArray(0, 0, KStart + K) = 1._Real;
             });
         Real Dt = 0.5_Real;
         TimeInterval TimeStep(Dt, TimeUnits::Seconds);

         const I4 NTSteps = 100;

         // Integrate 100 time steps and confirm monotonicity
         for (int N = 0; N < NTSteps; ++N) {
            deepCopy(TendCell3D, 0._Real);
            DefVertAdv->computeTracerVAdvTend(TendCell3D, TracerArray,
                                              LayerThick, TimeStep);
            parallelFor(
                {NVertLayers}, KOKKOS_LAMBDA(const int K) {
                   TracerArray(0, 0, K) += Dt * TendCell3D(0, 0, K);
                });
         }
         HostArray3DReal TracerH = createHostMirrorCopy(TracerArray);
         HostArray3DReal TendH   = createHostMirrorCopy(TendCell3D);
         for (int K = 0; K < NVertLayers; ++K) {
            if (TracerH(0, 0, K) < 0._Real or TracerH(0, 0, K) > 1._Real) {
               ++Err;
            }
         }
         if (Err != 0) {
            ErrAll +=
                Error(ErrorCode::Fail,
                      "VertAdvTest: computeFCTVAdvTend with velocity = {} "
                      "monotonicity FAIL",
                      TestVel[IMono]);
         }
      }

      finalizeVertAdvTest();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   CHECK_ERROR_ABORT(ErrAll, "VertAdv unit tests FAIL");

   return 0;
}
