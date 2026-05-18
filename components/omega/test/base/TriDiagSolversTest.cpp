#include "TriDiagSolvers.h"
#include "../ocn/OceanTestCommon.h"

using namespace OMEGA;

void initTridiagonalTest() {

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
}

int testCorrectness(int NBatch, int NRow) {

   Array2DReal DL("DL", NBatch, NRow);
   Array2DReal D("D", NBatch, NRow);
   Array2DReal DU("DU", NBatch, NRow);
   Array2DReal X("X", NBatch, NRow);

   // Create prescribed tridiagonal matrix A = (DL, D, DU) and prescribed vector
   // X
   parallelFor(
       {NBatch, NRow}, KOKKOS_LAMBDA(int I, int K) {
          X(I, K) = 0.1 * K + 0.3 * (I % 12);

          DL(I, K) = K == 0 ? 0 : 1 + 0.1 * (I % 3) + 0.2 * (K % 7);
          D(I, K)  = 4 + 0.2 * (I % 11) + 0.1 * (K % 5);
          DU(I, K) = K == NRow - 1 ? 0 : 1 + 0.05 * (I % 5) - 0.1 * (K % 11);
       });

   // compute A X
   Array2DReal AX("AX", NBatch, NRow);
   parallelFor(
       {NBatch, NRow}, KOKKOS_LAMBDA(int I, int K) {
          AX(I, K) = D(I, K) * X(I, K);
          if (K > 0) {
             AX(I, K) += DL(I, K) * X(I, K - 1);
          }
          if (K < NRow - 1) {
             AX(I, K) += DU(I, K) * X(I, K + 1);
          }
       });

   // Solve the system A Y = AX
   // The solution Y is stored in AX
   TriDiagSolver::solve(DL, D, DU, AX);

   // Check that Y is very close to the initial vector X
   Real Error;
   parallelReduce(
       {NBatch, NRow},
       KOKKOS_LAMBDA(int I, int K, Real &Accum) {
          Accum = Kokkos::max(Accum, Kokkos::abs(X(I, K) - AX(I, K)));
       },
       Kokkos::Max<Real>(Error));

   if (Error > 1e-12) {
      LOG_ERROR(
          "TridiagonalSolver: Correctness for (NBatch, NRow) = ({}, {}) FAIL",
          NBatch, NRow);
      return 1;
   }
   return 0;
}

int testDiffusionCorrectness(int NBatch, int NRow) {

   Array2DReal G("G", NBatch, NRow);
   Array2DReal H("H", NBatch, NRow);
   Array2DReal X("X", NBatch, NRow);

   // Create prescribed tridiagonal matrix for a diffusion problem
   // using the parametrization A = (G, H) and prescribed vector X
   parallelFor(
       {NBatch, NRow}, KOKKOS_LAMBDA(int I, int K) {
          X(I, K) = 0.1 * K + 0.3 * (I % 12);
          G(I, K) = K == NRow - 1 ? 0 : 1 + 0.1 * (I % 3) + 0.2 * (K % 7);
          H(I, K) = 4 + 0.2 * (I % 11) + 0.1 * (K % 5);
       });

   // compute A X
   Array2DReal AX("AX", NBatch, NRow);
   parallelFor(
       {NBatch, NRow}, KOKKOS_LAMBDA(int I, int K) {
          const Real DL = K == 0 ? 0 : -G(I, K - 1);
          const Real DU = -G(I, K);
          const Real D  = H(I, K) - DL - DU;

          AX(I, K) = D * X(I, K);
          if (K > 0) {
             AX(I, K) += DL * X(I, K - 1);
          }
          if (K < NRow - 1) {
             AX(I, K) += DU * X(I, K + 1);
          }
       });

   // Solve the system A Y = AX
   // The solution Y is stored in AX
   TriDiagDiffSolver::solve(G, H, AX);

   // Check that Y is very close to the initial vector X
   Real Error;
   parallelReduce(
       {NBatch, NRow},
       KOKKOS_LAMBDA(int I, int K, Real &Accum) {
          Accum = Kokkos::max(Accum, Kokkos::abs(X(I, K) - AX(I, K)));
       },
       Kokkos::Max<Real>(Error));

   if (Error > 1e-12) {
      LOG_ERROR("TridiagonalSolver: Diffusion correctness for (NBatch, NRow) = "
                "({}, {}) FAIL",
                NBatch, NRow);
      return 1;
   }
   return 0;
}

KOKKOS_FUNCTION Real manufacturedSolution(Real X, Real T) {
   return Kokkos::cos(X) * Kokkos::sin(T);
}

KOKKOS_FUNCTION Real manufacturedDiffusivity(Real X, Real T) {
   return 2 + Kokkos::sin(X);
}

KOKKOS_FUNCTION Real manufacturedForcing(Real X, Real T) {
   using Kokkos::cos;
   using Kokkos::sin;
   return (2 * sin(T) * sin(X) + 2 * sin(T) + cos(T)) * cos(X);
}

Real runDiffManufactured(int NCells) {
   // Setup a 1D manufactured solution diffusion problem using NCells

   const int NVertices = NCells + 1;

   const Real TimeEnd  = 1;
   const Real TimeStep = 0.001 / (NCells / 100);

   const int NSteps = std::ceil(TimeEnd / TimeStep);

   Array1DReal XVertex("XVertex", NVertices);
   Array1DReal Diffusivity("Diffusivity", NVertices);
   parallelFor(
       {NVertices}, KOKKOS_LAMBDA(int IVertex) {
          const Real XVertexUni = IVertex * (1._Real / NCells);
          // A simple transformation to make grid spacing non-uniform
          XVertex(IVertex)     = Kokkos::tanh(5 * XVertexUni);
          Diffusivity(IVertex) = manufacturedDiffusivity(XVertex(IVertex), 0);
       });

   Array1DReal XCell("XCell", NCells);
   Array1DReal PseudoThick("PseudoThick", NCells);
   Array1DReal U("U", NCells);

   // Create initial condition
   parallelFor(
       {NCells}, KOKKOS_LAMBDA(int ICell) {
          XCell(ICell)       = (XVertex(ICell + 1) + XVertex(ICell)) / 2;
          PseudoThick(ICell) = XVertex(ICell + 1) - XVertex(ICell);
          U(ICell)           = manufacturedSolution(XCell(ICell), 0);
       });

   auto LConfig = TriDiagDiffSolver::makeLaunchConfig(1, NCells);

   // Integrate in time with backward Euler
   for (int Step = 0; Step < NSteps; ++Step) {
      const Real Time     = Step * TimeStep;
      const Real TimeNext = (Step + 1) * TimeStep;

      parallelForOuter(
          LConfig, KOKKOS_LAMBDA(int, const TeamMember &Team) {
             TriDiagDiffScratch Scratch(Team, NCells);

             // Setup the system to be solved
             parallelForInner(Team, NCells, [=](int ICell) {
                for (int IVec = 0; IVec < VecLength; ++IVec) {

                   // Forcing term from the manufactured solution
                   const Real F = manufacturedForcing(XCell(ICell), TimeNext);

                   Scratch.H(ICell, IVec) = PseudoThick(ICell);

                   if (ICell == NCells - 1) {
                      // Boundary condition
                      const Real XBnd = XVertex(ICell + 1);
                      const Real BoundaryCoeff =
                          -(2 + Kokkos::sin(XBnd)) * Kokkos::tan(XBnd);
                      Scratch.H(ICell, IVec) -= TimeStep * BoundaryCoeff;
                      Scratch.G(ICell, IVec) = 0;
                   } else {
                      const Real AvgPseudoThick =
                          (PseudoThick(ICell + 1) + PseudoThick(ICell)) / 2;
                      Scratch.G(ICell, IVec) =
                          Diffusivity(ICell + 1) * TimeStep / AvgPseudoThick;
                   }
                   // RHS
                   Scratch.X(ICell, IVec) =
                       PseudoThick(ICell) * (U(ICell) + TimeStep * F);
                }
             });

             // Solve the system
             teamBarrier(Team);
             TriDiagDiffSolver::solve(Team, Scratch);
             teamBarrier(Team);

             // Store the solution
             parallelForInner(Team, NCells, [=](int ICell) {
                U(ICell) = Scratch.X(ICell, 0);
             });
          });
   }

   // Compute L2 error
   Real L2Error;
   parallelReduce(
       {NCells},
       KOKKOS_LAMBDA(int ICell, Real &Accum) {
          const Real UExact = manufacturedSolution(XCell(ICell), TimeEnd);
          const Real DU     = U(ICell) - UExact;
          Accum += PseudoThick(ICell) * DU * DU;
       },
       L2Error);

   L2Error = Kokkos::sqrt(L2Error);

   return L2Error;
}

int testDiffusionManufactured() {
   int Err = 0;

   // Compute L2 error with 100 cells
   int NCells          = 100;
   const Real L2Err100 = runDiffManufactured(NCells);

   // Compute L2 error with 200 cells
   NCells *= 2;
   const Real L2Err200 = runDiffManufactured(NCells);

   const Real L2Rate = std::log2(L2Err100 / L2Err200);

   // Check convergence rate
   if (std::abs(L2Rate - 2) > 0.1) {
      Err += 1;
      LOG_ERROR("TridiagonalSolver: Wrong conv rate for manufactured solution, "
                "rate = {}",
                L2Rate);
   }

   // Check error magnitude
   const Real MaxL2Err = 2e-5;
   if (L2Err200 > MaxL2Err) {
      Err += 1;
      LOG_ERROR("TridiagonalSolver: Too large error for manufactured solution, "
                "error = {}, max error = {}",
                L2Err200, MaxL2Err);
   }

   return Err;
}

Real runDiffusionStability(bool UseGeneralSolver, Real DiffValue) {
   // Setup a 1D diffusion problem with discontinuous diffusivity

   const int NCells    = 100;
   const int NVertices = NCells + 1;

   const Real TimeEnd  = 100;
   const Real TimeStep = 1;

   const int NSteps = std::ceil(TimeEnd / TimeStep);

   // Problem domain is [0, 1]
   const Real DX = 1.0 / NCells;

   // Create discontinuous diffusivity with DiffValue in [0.3, 0.7] else 0
   Array1DReal XVertex("XVertex", NVertices);
   Array1DReal Diffusivity("Diffusivity", NVertices);
   parallelFor(
       {NVertices}, KOKKOS_LAMBDA(int IVertex) {
          XVertex(IVertex) = IVertex * DX;
          Diffusivity(IVertex) =
              Kokkos::abs(XVertex(IVertex) - 0.5) < 0.2 ? DiffValue : 0;
       });

   Array1DReal XCell("XCell", NCells);
   Array1DReal PseudoThick("PseudoThick", NCells);
   Array1DReal U("U", NCells);

   // Create initial condition
   parallelFor(
       {NCells}, KOKKOS_LAMBDA(int ICell) {
          XCell(ICell)       = ICell * DX + DX / 2;
          PseudoThick(ICell) = DX;
          const Real Tmp     = XCell(ICell) - 0.5_Real;
          U(ICell)           = Kokkos::exp(-Tmp * Tmp);
       });

   // Compute initial condition norm
   Real NormInit;
   parallelReduce(
       {NCells},
       KOKKOS_LAMBDA(int ICell, Real &Accum) {
          Accum += PseudoThick(ICell) * U(ICell) * U(ICell);
       },
       NormInit);
   NormInit = std::sqrt(NormInit);

   // Time integration using backward Euler
   for (int Step = 0; Step < NSteps; ++Step) {

      if (UseGeneralSolver) {
         auto LConfig = TriDiagSolver::makeLaunchConfig(1, NCells);

         parallelForOuter(
             LConfig, KOKKOS_LAMBDA(int, const TeamMember &Team) {
                TriDiagScratch Scratch(Team, NCells);

                // Setup the system to be solved in the form expected by the
                // general tridiagonal solver
                parallelForInner(Team, NCells, [=](int ICell) {
                   for (int IVec = 0; IVec < VecLength; ++IVec) {

                      if (ICell < NCells - 1) {
                         const Real AvgPseudoThick =
                             (PseudoThick(ICell + 1) + PseudoThick(ICell)) / 2;
                         Scratch.DU(ICell, IVec) = -Diffusivity(ICell + 1) *
                                                   TimeStep / AvgPseudoThick;
                      } else {
                         Scratch.DU(ICell, IVec) = 0;
                      }

                      if (ICell > 0) {
                         const Real AvgPseudoThick =
                             (PseudoThick(ICell) + PseudoThick(ICell - 1)) / 2;
                         Scratch.DL(ICell, IVec) =
                             -Diffusivity(ICell) * TimeStep / AvgPseudoThick;
                      } else {
                         Scratch.DL(ICell, IVec) = 0;
                      }

                      Scratch.D(ICell, IVec) = PseudoThick(ICell) -
                                               Scratch.DU(ICell, IVec) -
                                               Scratch.DL(ICell, IVec);

                      Scratch.X(ICell, IVec) = PseudoThick(ICell) * U(ICell);
                   }
                });

                // Solve the system
                teamBarrier(Team);
                TriDiagSolver::solve(Team, Scratch);
                teamBarrier(Team);

                // Save the solution
                parallelForInner(Team, NCells, [=](int ICell) {
                   U(ICell) = Scratch.X(ICell, 0);
                });
             });
      } else {
         auto LConfig = TriDiagDiffSolver::makeLaunchConfig(1, NCells);

         parallelForOuter(
             LConfig, KOKKOS_LAMBDA(int, const TeamMember &Team) {
                TriDiagDiffScratch Scratch(Team, NCells);

                // Setup the system to be solved in the form expected by the
                // specialized diffusion tridiagonal solver
                parallelForInner(Team, NCells, [=](int ICell) {
                   for (int IVec = 0; IVec < VecLength; ++IVec) {

                      Scratch.H(ICell, IVec) = PseudoThick(ICell);

                      if (ICell < NCells - 1) {
                         const Real AvgPseudoThick =
                             (PseudoThick(ICell + 1) + PseudoThick(ICell)) / 2;
                         Scratch.G(ICell, IVec) =
                             Diffusivity(ICell + 1) * TimeStep / AvgPseudoThick;
                      } else {
                         Scratch.G(ICell, IVec) = 0;
                      }

                      Scratch.X(ICell, IVec) = PseudoThick(ICell) * U(ICell);
                   }
                });

                // Solve the system
                teamBarrier(Team);
                TriDiagDiffSolver::solve(Team, Scratch);
                teamBarrier(Team);

                // Store the solution
                parallelForInner(Team, NCells, [=](int ICell) {
                   U(ICell) = Scratch.X(ICell, 0);
                });
             });
      }
   }

   // Compute the solution norm
   Real Norm;
   parallelReduce(
       {NCells},
       KOKKOS_LAMBDA(int ICell, Real &Accum) {
          Accum += PseudoThick(ICell) * U(ICell) * U(ICell);
       },
       Norm);
   Norm = std::sqrt(Norm);

   // Return normalized change in the norm
   return (Norm - NormInit) / NormInit;
}

int testDiffusionStability() {
   int Err = 0;

   bool UseGeneralSolver;

   // First check small change in diffusivity
   // Both general and specialized solver should work
   // and result in similar change in norms
   const Real SmallDiffValue = 1e2;

   UseGeneralSolver = true;
   const Real NormGeneralSmallDiff =
       runDiffusionStability(UseGeneralSolver, SmallDiffValue);

   UseGeneralSolver = false;
   const Real NormCustomSmallDiff =
       runDiffusionStability(UseGeneralSolver, SmallDiffValue);

   if (!isApprox(NormGeneralSmallDiff, NormCustomSmallDiff, 1e-3)) {
      Err += 1;
      LOG_ERROR("TridiagonalSolver: Different norms");
   }

   // Now check abrupt change in diffusivity
   // General solver is expected to fail (produce NaNs)
   // Specialized solver is expected to work
   const Real LargeDiffValue = 1e14;

   UseGeneralSolver = true;
   const Real NormGeneralLargeDiff =
       runDiffusionStability(UseGeneralSolver, LargeDiffValue);

   if (!std::isnan(NormGeneralLargeDiff)) {
      Err += 1;
      LOG_ERROR("TridiagonalSolver: Expected general solver to fail");
   }

   UseGeneralSolver = false;
   const Real NormCustomLargeDiff =
       runDiffusionStability(UseGeneralSolver, LargeDiffValue);

   if (std::isnan(NormCustomLargeDiff)) {
      Err += 1;
      LOG_ERROR("TridiagonalSolver: Expected custom solver to pass");
   }

   return Err;
}

int tridiagonalTest() {
   int Err = 0;

   initTridiagonalTest();

   LOG_INFO("----- Tridiagonal Unit Test -----");

   for (int NBatch : {1, 2, 4, 5, 11, 33, 102}) {
      for (int NRow : {3, 4, 5, 6, 11, 17, 64, 100}) {
         Err += testCorrectness(NBatch, NRow);
         Err += testDiffusionCorrectness(NBatch, NRow);
      }
   }

   Err += testDiffusionManufactured();

   Err += testDiffusionStability();

   if (Err == 0) {
      LOG_INFO("TridiagonalTest: Successful completion");
   }

   return Err;
}

void testPerformance(int NBatch, int NRow, int NRep) {

   Array2DReal DL("DL", NBatch, NRow);
   Array2DReal D("D", NBatch, NRow);
   Array2DReal DU("DU", NBatch, NRow);
   Array2DReal X("X", NBatch, NRow);

   parallelFor(
       {NBatch, NRow}, KOKKOS_LAMBDA(int I, int K) {
          X(I, K) = 0.1 * K + 0.3 * (I % 12);

          DL(I, K) = K == 0 ? 0 : 1 + 0.1 * (I % 3) + 0.2 * (K % 7);
          D(I, K)  = 4 + 0.2 * (I % 11) + 0.1 * (K % 5);
          DU(I, K) = K == NRow - 1 ? 0 : 1 + 0.05 * (I % 5) - 0.1 * (K % 11);
       });

   // Warmup
   TriDiagSolver::solve(DL, D, DU, X);

   Kokkos::fence();
   Kokkos::Timer Timer;
   for (int Rep = 0; Rep < NRep; ++Rep) {
      TriDiagSolver::solve(DL, D, DU, X);
   }
   Kokkos::fence();
   double TimeSeconds = Timer.seconds();

   double NSolvesPerSec = NRep * NBatch / TimeSeconds;
   double NBytes        = 5 * NRep * NBatch * NRow * sizeof(Real);
   double Bandwidth     = NBytes / 1e9 / TimeSeconds;

   printf("%12d %12d %18.2f %18e\n", NRow, NBatch, Bandwidth, NSolvesPerSec);
}

void testDiffusionPerformance(int NBatch, int NRow, int NRep) {

   Array2DReal G("G", NBatch, NRow);
   Array2DReal H("H", NBatch, NRow);
   Array2DReal X("X", NBatch, NRow);

   parallelFor(
       {NBatch, NRow}, KOKKOS_LAMBDA(int I, int K) {
          X(I, K) = 0.1 * K + 0.3 * (I % 12);
          G(I, K) = K == NRow - 1 ? 0 : 1 + 0.1 * (I % 3) + 0.2 * (K % 7);
          H(I, K) = 4 + 0.2 * (I % 11) + 0.1 * (K % 5);
       });

   // Warmup
   TriDiagDiffSolver::solve(G, H, X);

   Kokkos::fence();
   Kokkos::Timer Timer;
   for (int Rep = 0; Rep < NRep; ++Rep) {
      TriDiagDiffSolver::solve(G, H, X);
   }
   Kokkos::fence();
   double TimeSeconds = Timer.seconds();

   double NSolvesPerSec = NRep * NBatch / TimeSeconds;
   double NBytes        = 4 * NRep * NBatch * NRow * sizeof(Real);
   double Bandwidth     = NBytes / 1e9 / TimeSeconds;

   printf("%12d %12d %18.2f %18e\n", NRow, NBatch, Bandwidth, NSolvesPerSec);
}

int tridiagonalPerfTest() {
   initTridiagonalTest();

   const int NRows[] = {64};
#ifdef OMEGA_TARGET_DEVICE
   const int NBatches[] = {500, 1000, 5000, 10000};
   const int NWork      = 120 * 500 * 64;
#else
   const int NBatches[] = {100, 200, 400};
   const int NWork      = 40 * 100 * 64;
#endif

   printf("General solver performance\n");
   printf("%12s %12s %18s %18s\n", "NRow", "NBatch", "Bandwidth [GB/s]",
          "NSolvesPerSec");
   for (int NRow : NRows) {
      for (int NBatch : NBatches) {
         const int NRep = NWork / (NBatch * NRow);
         testPerformance(NBatch, NRow, NRep);
      }
   }

   printf("Diffusion solver performance\n");
   printf("%12s %12s %18s %18s\n", "NRow", "NBatch", "Bandwidth [GB/s]",
          "NSolvesPerSec");
   for (int NRow : NRows) {
      for (int NBatch : NBatches) {
         const int NRep = NWork / (NBatch * NRow);
         testDiffusionPerformance(NBatch, NRow, NRep);
      }
   }

   return 0;
}

int main(int argc, char *argv[]) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize(argc, argv);

   if (argc > 1 && std::string(argv[1]) == std::string("--perf")) {
      RetVal += tridiagonalPerfTest();
   } else {
      RetVal += tridiagonalTest();
   }

   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
