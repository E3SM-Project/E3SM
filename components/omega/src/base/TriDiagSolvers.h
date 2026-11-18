#ifndef OMEGA_TRIDIAGONAL_H
#define OMEGA_TRIDIAGONAL_H
//===-- base/TriDiagSolvers.h - tridiagonal solvers --------*- C++ -*-===//
//
/// \file
/// \brief Contains solvers for tridiagonal systems of equations
///
/// Batched Thomas and Parallel Cyclic Reduction (PCR) algorithms are
/// implemented for use on CPU and GPU platforms, respectively. Specialized
/// variants of both algorithms for diffusion-type systems, with better
/// stability, are also available.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

namespace OMEGA {

// Forward declaration of different solver types
struct ThomasSolver;
struct PCRSolver;
struct ThomasDiffusionSolver;
struct PCRDiffusionSolver;

#ifdef OMEGA_TARGET_DEVICE
// On GPUs we use the Parallel Cyclic Reduction algorithm
using TriDiagSolver     = PCRSolver;
using TriDiagDiffSolver = PCRDiffusionSolver;
#else
// On CPUs we use the Thomas algorithm
using TriDiagSolver     = ThomasSolver;
using TriDiagDiffSolver = ThomasDiffusionSolver;
#endif

// Type of real array of size (NRow, VecLength) in the scratch memory space
using TriDiagScratchArray =
    Kokkos::View<Real *[VecLength], MemLayout, ScratchMemSpace,
                 Kokkos::MemoryUnmanaged>;

// Scratch data for general tridiagonal solver
struct TriDiagScratch {
   TriDiagScratchArray DL; // lower diagonal
   TriDiagScratchArray D;  // main diagonal
   TriDiagScratchArray DU; // upper diagonal
   TriDiagScratchArray X; // rhs on input, contains solution after calling solve

   // Constructor takes team member and system size
   KOKKOS_FUNCTION TriDiagScratch(const TeamMember &Team, int NRow)
       : DL(teamScratch(Team), NRow), D(teamScratch(Team), NRow),
         DU(teamScratch(Team), NRow), X(teamScratch(Team), NRow) {}
};

// Thomas algorithm solver for general tridiagonal systems
struct ThomasSolver {

   // Create a Kokkos team policy for solving NBatch systems of size NRow
   // and set scratch size
   static LaunchConfig<1> makeLaunchConfig(int NBatch, int NRow) {
      const int NTeams   = (NBatch + VecLength - 1) / VecLength;
      const int NScratch = 4 * NRow * VecLength;
      return LaunchConfig({NTeams}, 1, TeamScratch<Real>(NScratch));
   }

   // Solve the system defined in the scratch data argument `Scratch`
   // This a team-level function that needs to be called inside a
   // parallel loop using TeamPolicy, hence it has a team member argument
   static void KOKKOS_FUNCTION solve(const TeamMember &Team,
                                     const TriDiagScratch &Scratch) {
      const int NRow = Scratch.X.extent_int(0);

      for (int K = 1; K < NRow; ++K) {
         for (int IVec = 0; IVec < VecLength; ++IVec) {
            const Real W = Scratch.DL(K, IVec) / Scratch.D(K - 1, IVec);
            Scratch.D(K, IVec) -= W * Scratch.DU(K - 1, IVec);
            Scratch.X(K, IVec) -= W * Scratch.X(K - 1, IVec);
         }
      }

      for (int IVec = 0; IVec < VecLength; ++IVec) {
         Scratch.X(NRow - 1, IVec) /= Scratch.D(NRow - 1, IVec);
      }

      for (int K = NRow - 2; K >= 0; --K) {
         for (int IVec = 0; IVec < VecLength; ++IVec) {
            Scratch.X(K, IVec) =
                (Scratch.X(K, IVec) -
                 Scratch.DU(K, IVec) * Scratch.X(K + 1, IVec)) /
                Scratch.D(K, IVec);
         }
      }
   }

   // Solve the system stored in the arrays DL, D, DU, and X
   // This function should be called outside of a parallel region
   static void solve(const Array2DReal &DL, const Array2DReal &D,
                     const Array2DReal &DU, const Array2DReal &X) {

      const int NBatch = X.extent_int(0);
      const int NRow   = X.extent_int(1);

      auto LConfig = makeLaunchConfig(NBatch, NRow);

      parallelForOuter(
          LConfig, KOKKOS_LAMBDA(const int IChunk, const TeamMember &Team) {
             const int IStart = IChunk * VecLength;

             TriDiagScratch Scratch(Team, NRow);

             for (int K = 0; K < NRow; ++K) {
                for (int IVec = 0; IVec < VecLength; ++IVec) {
                   const int I = IStart + IVec;
                   if (I < NBatch) {
                      Scratch.DL(K, IVec) = DL(I, K);
                      Scratch.D(K, IVec)  = D(I, K);
                      Scratch.DU(K, IVec) = DU(I, K);
                      Scratch.X(K, IVec)  = X(I, K);
                   }
                }
             }

             solve(Team, Scratch);

             for (int IVec = 0; IVec < VecLength; ++IVec) {
                for (int K = 0; K < NRow; ++K) {
                   const int I = IStart + IVec;
                   if (I < NBatch) {
                      X(I, K) = Scratch.X(K, IVec);
                   }
                }
             }
          });
   }
};

// Parallel Cyclic Reduction solver for general tridiagonal systems
struct PCRSolver {

   // Create a Kokkos team policy for solving NBatch systems of size NRow
   // and set scratch size
   static LaunchConfig<1> makeLaunchConfig(int NBatch, int NRow) {
      const int NScratch = 4 * NRow * VecLength;
      return LaunchConfig({NBatch}, NRow, TeamScratch<Real>(NScratch));
   }

   // Solve the system defined in the scratch data argument `Scratch`
   // This a team-level function that needs to be called inside a
   // parallel loop using TeamPolicy, hence it has a team member argument
   static void KOKKOS_FUNCTION solve(const TeamMember &Team,
                                     const TriDiagScratch &Scratch) {
      const int NRow = Scratch.X.extent_int(0);

      // Row index = Thread index
      const int K = Team.team_rank();

      // Number of reduction levels
      const int NLevels = Kokkos::ceil(Kokkos::log2(NRow));

      // Perform NLevels of parallel cyclic reduction
      for (int Lev = 1; Lev < NLevels; ++Lev) {

         const int HalfStride = 1 << (Lev - 1);

         int Kmh = K - HalfStride;
         Kmh     = Kmh < 0 ? 0 : Kmh;
         int Kph = K + HalfStride;
         Kph     = Kph >= NRow ? NRow - 1 : Kph;

         const Real alpha = -Scratch.DL(K, 0) / Scratch.D(Kmh, 0);
         const Real gamma = -Scratch.DU(K, 0) / Scratch.D(Kph, 0);

         // Compute new system coefficients
         const Real NewD = Scratch.D(K, 0) + alpha * Scratch.DU(Kmh, 0) +
                           gamma * Scratch.DL(Kph, 0);
         const Real NewX = Scratch.X(K, 0) + alpha * Scratch.X(Kmh, 0) +
                           gamma * Scratch.X(Kph, 0);
         const Real NewDL = alpha * Scratch.DL(Kmh, 0);
         const Real NewDU = gamma * Scratch.DU(Kph, 0);

         teamBarrier(Team);

         // Store new system coefficients
         Scratch.D(K, 0)  = NewD;
         Scratch.X(K, 0)  = NewX;
         Scratch.DL(K, 0) = NewDL;
         Scratch.DU(K, 0) = NewDU;

         teamBarrier(Team);
      }

      const int Stride = 1 << (NLevels - 1);

      // Need to solve 2x2 system
      if (K + Stride < NRow || K - Stride >= 0) {
         // The result is two values so only half of threads do work
         if (K < NRow / 2) {
            const Real Det = Scratch.D(K, 0) * Scratch.D(K + Stride, 0) -
                             Scratch.DL(K + Stride, 0) * Scratch.DU(K, 0);
            const Real Xk   = Scratch.X(K, 0);
            const Real Xkps = Scratch.X(K + Stride, 0);
            Scratch.X(K, 0) =
                (Scratch.D(K + Stride, 0) * Xk - Scratch.DU(K, 0) * Xkps) / Det;
            Scratch.X(K + Stride, 0) =
                (-Scratch.DL(K + Stride, 0) * Xk + Scratch.D(K, 0) * Xkps) /
                Det;
         }
      } else { // 1x1 system -> direct solution
         Scratch.X(K, 0) /= Scratch.D(K, 0);
      }
   }

   // Solve the system stored in the arrays DL, D, DU, and X
   // This function should be called outside of a parallel region
   static void solve(const Array2DReal &DL, const Array2DReal &D,
                     const Array2DReal &DU, const Array2DReal &X) {

      const int NBatch = X.extent_int(0);
      const int NRow   = X.extent_int(1);
      auto LConfig     = makeLaunchConfig(NBatch, NRow);

      parallelForOuter(
          LConfig, KOKKOS_LAMBDA(int I, const TeamMember &Team) {
             const int K = Team.team_rank();

             TriDiagScratch Scratch(Team, NRow);

             Scratch.DL(K, 0) = DL(I, K);
             Scratch.D(K, 0)  = D(I, K);
             Scratch.DU(K, 0) = DU(I, K);
             Scratch.X(K, 0)  = X(I, K);

             teamBarrier(Team);

             solve(Team, Scratch);

             teamBarrier(Team);

             X(I, K) = Scratch.X(K, 0);
          });
   }
};

// Scratch data for diffuson-type tridiagonal solver
// System is asumed to be of the form:
// -G(i) * x_{i-1} + (H(i) + G(i) + G(i+1)) * x_i - G(i+1) * x_{i+1} = y_i
struct TriDiagDiffScratch {
   TriDiagScratchArray G; // G above
   TriDiagScratchArray H; // H above
   TriDiagScratchArray X; // rhs on input, contains solution after calling solve
   TriDiagScratchArray Alpha; // internal workspace

   KOKKOS_FUNCTION TriDiagDiffScratch(const TeamMember &Team, int NRow)
       : G(teamScratch(Team), NRow), H(teamScratch(Team), NRow),
         X(teamScratch(Team), NRow), Alpha(teamScratch(Team), NRow) {}
};

// Thomas algorithm solver for diffusion-type tridiagonal systems
struct ThomasDiffusionSolver {

   // Create a Kokkos team policy for solving NBatch systems of size NRow
   // and set scratch size
   static LaunchConfig<1> makeLaunchConfig(int NBatch, int NRow) {
      const int NTeams   = (NBatch + VecLength - 1) / VecLength;
      const int NScratch = 4 * NRow * VecLength;
      return LaunchConfig({NTeams}, 1, TeamScratch<Real>(NScratch));
   }

   // Solve the system defined in the scratch data argument `Scratch`
   // This a team-level function that needs to be called inside a
   // parallel loop using TeamPolicy, hence it has a team member argument
   static void KOKKOS_FUNCTION solve(const TeamMember &Team,
                                     const TriDiagDiffScratch &Scratch) {
      const int NRow = Scratch.X.extent_int(0);

      for (int IVec = 0; IVec < VecLength; ++IVec) {
         Scratch.Alpha(0, IVec) = 0;
      }

      for (int K = 1; K < NRow; ++K) {
         for (int IVec = 0; IVec < VecLength; ++IVec) {
            Scratch.Alpha(K, IVec) =
                Scratch.G(K - 1, IVec) *
                (Scratch.H(K - 1, IVec) + Scratch.Alpha(K - 1, IVec)) /
                (Scratch.H(K - 1, IVec) + Scratch.Alpha(K - 1, IVec) +
                 Scratch.G(K - 1, IVec));
         }
      }

      for (int IVec = 0; IVec < VecLength; ++IVec) {
         Scratch.H(0, IVec) += Scratch.G(0, IVec);
      }

      for (int K = 1; K < NRow; ++K) {
         for (int IVec = 0; IVec < VecLength; ++IVec) {
            const Real AddH = Scratch.Alpha(K, IVec) + Scratch.G(K, IVec);

            Scratch.H(K, IVec) += AddH;
            Scratch.X(K, IVec) += Scratch.G(K - 1, IVec) /
                                  Scratch.H(K - 1, IVec) *
                                  Scratch.X(K - 1, IVec);
         }
      }

      for (int IVec = 0; IVec < VecLength; ++IVec) {
         Scratch.X(NRow - 1, IVec) /= Scratch.H(NRow - 1, IVec);
      }

      for (int K = NRow - 2; K >= 0; --K) {
         for (int IVec = 0; IVec < VecLength; ++IVec) {
            Scratch.X(K, IVec) = (Scratch.X(K, IVec) +
                                  Scratch.G(K, IVec) * Scratch.X(K + 1, IVec)) /
                                 Scratch.H(K, IVec);
         }
      }
   }

   // Solve the system stored in the arrays G, H, and X
   // This function should be called outside of a parallel region
   static void solve(const Array2DReal &G, const Array2DReal &H,
                     const Array2DReal &X) {
      const int NBatch = X.extent_int(0);
      const int NRow   = X.extent_int(1);

      auto LConfig = makeLaunchConfig(NBatch, NRow);

      parallelForOuter(
          LConfig, KOKKOS_LAMBDA(int IChunk, const TeamMember &Team) {
             const int IStart = IChunk * VecLength;

             TriDiagDiffScratch Scratch(Team, NRow);

             for (int K = 0; K < NRow; ++K) {
                for (int IVec = 0; IVec < VecLength; ++IVec) {
                   const int I = IStart + IVec;
                   if (I < NBatch) {
                      Scratch.G(K, IVec) = G(I, K);
                      Scratch.H(K, IVec) = H(I, K);
                      Scratch.X(K, IVec) = X(I, K);
                   }
                }
             }

             solve(Team, Scratch);

             for (int IVec = 0; IVec < VecLength; ++IVec) {
                for (int K = 0; K < NRow; ++K) {
                   const int I = IStart + IVec;
                   if (I < NBatch) {
                      X(I, K) = Scratch.X(K, IVec);
                   }
                }
             }
          });
   }
};

// Parallel Cyclic Reduction solver for diffusion-type tridiagonal systems
struct PCRDiffusionSolver {

   // Create a Kokkos team policy for solving NBatch systems of size NRow
   // and set scratch size
   static LaunchConfig<1> makeLaunchConfig(int NBatch, int NRow) {
      const int NScratch = 4 * NRow * VecLength;
      return LaunchConfig({NBatch}, NRow, TeamScratch<Real>(NScratch));
   }

   // Solve the system defined in the scratch data argument `Scratch`
   // This a team-level function that needs to be called inside a
   // parallel loop using TeamPolicy, hence it has a team member argument
   static void KOKKOS_FUNCTION solve(const TeamMember &Team,
                                     const TriDiagDiffScratch &Scratch) {
      const int NRow = Scratch.X.extent_int(0);

      // Row index = Thread index
      const int K = Team.team_rank();

      // Number of reduction levels
      const int NLevels = Kokkos::ceil(Kokkos::log2(NRow));

      // Perform NLevels of parallel cyclic reduction
      for (int Lev = 1; Lev < NLevels; ++Lev) {

         const int Stride     = 1 << Lev;
         const int HalfStride = 1 << (Lev - 1);

         int Kmh         = K - HalfStride;
         const Real Gkmh = Kmh < 0 ? 0 : Scratch.G(Kmh, 0);
         Kmh             = Kmh < 0 ? 0 : Kmh;

         const int Kms   = K - Stride;
         const Real Gkms = Kms < 0 ? 0 : Scratch.G(Kms, 0);

         int Kph = K + HalfStride;
         Kph     = Kph >= NRow ? NRow - 1 : Kph;

         const Real Alpha = Gkmh / (Scratch.H(Kmh, 0) + Gkms + Gkmh);
         const Real Beta =
             Scratch.G(K, 0) /
             (Scratch.H(Kph, 0) + Scratch.G(K, 0) + Scratch.G(Kph, 0));

         // Compute new system coefficients
         const Real NewG = Scratch.G(Kph, 0) * Beta;
         const Real NewX = Scratch.X(K, 0) + Alpha * Scratch.X(Kmh, 0) +
                           Beta * Scratch.X(Kph, 0);
         const Real NewH = Scratch.H(K, 0) + Alpha * Scratch.H(Kmh, 0) +
                           Beta * Scratch.H(Kph, 0);

         teamBarrier(Team);

         // Store new system coefficients
         Scratch.H(K, 0) = NewH;
         Scratch.G(K, 0) = NewG;
         Scratch.X(K, 0) = NewX;

         teamBarrier(Team);
      }

      const int Stride = 1 << (NLevels - 1);

      // Need to solve 2x2 system
      if (K + Stride < NRow || K - Stride >= 0) {
         // The result is two values so only half of threads do work
         if (K < NRow / 2) {
            const int Kms   = K - Stride;
            const Real Gkms = Kms < 0 ? 0 : Scratch.G(Kms, 0);

            const Real Dk   = Scratch.H(K, 0) + Gkms + Scratch.G(K, 0);
            const Real Dkps = Scratch.H(K + Stride, 0) + Scratch.G(K, 0) +
                              Scratch.G(K + Stride, 0);
            const Real DUk   = -Scratch.G(K, 0);
            const Real DLkps = -Scratch.G(K, 0);

            const Real Det = Dk * Dkps - DLkps * DUk;

            const Real Xk   = Scratch.X(K, 0);
            const Real Xkps = Scratch.X(K + Stride, 0);

            Scratch.X(K, 0)          = (Dkps * Xk - DUk * Xkps) / Det;
            Scratch.X(K + Stride, 0) = (-DLkps * Xk + Dk * Xkps) / Det;
         }

      } else { // 1x1 system -> direct solution
         int Kms         = K - Stride;
         const Real Gkms = Kms < 0 ? 0 : Scratch.G(Kms, 0);
         Scratch.X(K, 0) /= (Scratch.H(K, 0) + Gkms + Scratch.G(K, 0));
      }
   }

   // Solve the system stored in the arrays G, H, and X
   // This function should be called outside of a parallel region
   static void solve(const Array2DReal &G, const Array2DReal &H,
                     const Array2DReal &X) {
      const int NBatch = X.extent_int(0);
      const int NRow   = X.extent_int(1);

      auto LConfig = makeLaunchConfig(NBatch, NRow);
      parallelForOuter(
          LConfig, KOKKOS_LAMBDA(int I, const TeamMember &Team) {
             const int K = Team.team_rank();

             TriDiagDiffScratch Scratch(Team, NRow);

             Scratch.G(K, 0) = G(I, K);
             Scratch.H(K, 0) = H(I, K);
             Scratch.X(K, 0) = X(I, K);

             teamBarrier(Team);

             solve(Team, Scratch);

             teamBarrier(Team);

             X(I, K) = Scratch.X(K, 0);
          });
   }
};

} // namespace OMEGA
#endif
