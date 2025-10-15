//===-- Test driver for OMEGA Kokkos -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Kokkos
///
/// This driver tests the Omega wrappers of Kokkos capabilities.
///
//
//===-----------------------------------------------------------------------===/

#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"

#include "mpi.h"
#include <limits>
#include <sstream>

using namespace OMEGA;

template <class Array> bool hostArraysEqual(const Array &A, const Array &B) {
   bool Equal = true;
   for (size_t I = 0; I < A.span(); I++) {
      if (A.data()[I] != B.data()[I]) {
         Equal = false;
         break;
      }
   }
   return Equal;
}

static KOKKOS_FUNCTION int f2(int J1, int J2, int N1, int N2) {
   return -(N1 * N2) / 4 + J1 + N1 * J2;
}

static KOKKOS_FUNCTION int f3(int J1, int J2, int J3, int N1, int N2, int N3) {
   return -(N1 * N2 * N3) / 4 + J1 + N1 * (J2 + N2 * J3);
}

// Helpers to add test arguments to error messages
static std::string errorMsg(const std::string &Msg, int N1) {
   std::stringstream Stream;
   Stream << Msg << " N1 = " << N1;
   return Stream.str();
}

static std::string errorMsg(const std::string &Msg, int N1, int N2) {
   std::stringstream Stream;
   Stream << Msg << " N1 = " << N1 << " N2 = " << N2;
   return Stream.str();
}

static std::string errorMsg(const std::string &Msg, int N1, int N2, int N3) {
   std::stringstream Stream;
   Stream << Msg << " N1 = " << N1 << " N2 = " << N2 << " N3 = " << N3;
   return Stream.str();
}

Error testHiparFor1DFor1D(int N1) {
   Error Err;

   const int N2 = N1;

   HostArray2DI4 RefAH("RefA2H", N1, N2);
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < J1; ++J2) {
         RefAH(J1, J2) = f2(J1, J2, N1, N2);
      }
   }

   Array2DI4 A("A2", N1, N2);
   parallelForOuter(
       {N1}, KOKKOS_LAMBDA(int J1, const TeamMember &Team) {
          parallelForInner(
              Team, J1,
              INNER_LAMBDA(int J2) { A(J1, J2) = f2(J1, J2, N1, N2); });
       });

   auto AH = createHostMirrorCopy(A);

   if (!hostArraysEqual(AH, RefAH)) {
      Err += Error(ErrorCode::Fail, errorMsg("parallelFor1DFor1D FAIL", N1));
   }

   return Err;
}

Error testHiparFor1DReduce1D(int N1) {
   Error Err;

   const int N2 = N1;

   HostArray1DI4 RefSumH("RefSum1H", N1);
   HostArray1DI4 RefMaxH("RefMax1H", N1);
   for (int J1 = 0; J1 < N1; ++J1) {
      I4 Sum = 0;
      I4 Max = std::numeric_limits<I4>::min();
      for (int J2 = 0; J2 < J1; ++J2) {
         Sum += f2(J1, J2, N1, N2);
         Max = std::max(Max, f2(J1, J2, N1, N2));
      }
      RefSumH(J1) = Sum;
      RefMaxH(J1) = Max;
   }

   Array1DI4 Sum1("Sum11", N1);
   Array1DI4 Max1("Max11", N1);
   parallelForOuter(
       {N1}, KOKKOS_LAMBDA(int J1, const TeamMember &Team) {
          I4 Sum;
          parallelReduceInner(
              Team, J1,
              INNER_LAMBDA(int J2, I4 &Accum) { Accum += f2(J1, J2, N1, N2); },
              Sum);
          Sum1(J1) = Sum;

          I4 Max;
          parallelReduceInner(
              Team, J1,
              INNER_LAMBDA(int J2, I4 &Accum) {
                 Accum = Kokkos::max(Accum, f2(J1, J2, N1, N2));
              },
              Kokkos::Max<I4>(Max));
          Max1(J1) = Max;
       });

   auto Sum1H = createHostMirrorCopy(Sum1);
   if (!hostArraysEqual(Sum1H, RefSumH)) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelFor1DReduce1D Sum FAIL", N1));
   }

   auto Max1H = createHostMirrorCopy(Max1);
   if (!hostArraysEqual(Max1H, RefMaxH)) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelFor1DReduce1D Max FAIL", N1));
   }

   Array1DI4 Sum2("Sum21", N1);
   Array1DI4 Max2("Max21", N1);
   parallelForOuter(
       {N1}, KOKKOS_LAMBDA(int J1, const TeamMember &Team) {
          I4 Sum, Max;
          parallelReduceInner(
              Team, J1,
              INNER_LAMBDA(int J2, I4 &AccumSum, I4 &AccumMax) {
                 AccumSum += f2(J1, J2, N1, N2);
                 AccumMax = Kokkos::max(AccumMax, f2(J1, J2, N1, N2));
              },
              Sum, Kokkos::Max<I4>(Max));
          Sum2(J1) = Sum;
          Max2(J1) = Max;
       });

   auto Sum2H = createHostMirrorCopy(Sum2);
   auto Max2H = createHostMirrorCopy(Max2);
   if (!hostArraysEqual(Sum2H, RefSumH) || !hostArraysEqual(Max2H, RefMaxH)) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelFor1DReduce1D Sum+Max FAIL", N1));
   }

   return Err;
}

Error testHiparFor1DScan1D(int N1) {
   Error Err;

   const int N2 = N1;

   HostArray2DI4 RefRSumH("RefRSum2H", N1, N2);
   for (int J1 = 0; J1 < N1; ++J1) {
      I4 RSum = 0;
      for (int J2 = 0; J2 < J1; ++J2) {
         RefRSumH(J1, J2) = RSum;
         RSum += f2(J1, J2, N1, N2);
      }
   }

   Array2DI4 RSum("RSum2", N1, N2);
   parallelForOuter(
       {N1}, KOKKOS_LAMBDA(int J1, const TeamMember &Team) {
          parallelScanInner(
              Team, J1, INNER_LAMBDA(int J2, I4 &Accum, bool IsFinal) {
                 if (IsFinal) {
                    RSum(J1, J2) = Accum;
                 }
                 Accum += f2(J1, J2, N1, N2);
              });
       });

   auto RSumH = createHostMirrorCopy(RSum);
   if (!hostArraysEqual(RSumH, RefRSumH)) {
      Err +=
          Error(ErrorCode::Fail, errorMsg("parallelFor1DScan1D Sum FAIL", N1));
   }

   return Err;
}

Error testHiparReduce1DReduce1D(int N1) {
   Error Err;

   const int N2 = N1;

   I4 RefSum = 0;
   I4 RefMax = std::numeric_limits<I4>::min();
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < J1; ++J2) {
         RefSum += f2(J1, J2, N1, N2);
         RefMax = std::max(RefMax, f2(J1, J2, N1, N2));
      }
   }

   I4 Sum1;
   parallelReduceOuter(
       {N1},
       KOKKOS_LAMBDA(int J1, const TeamMember &Team, I4 &AccumOuter) {
          I4 SumInner;
          parallelReduceInner(
              Team, J1,
              INNER_LAMBDA(int J2, I4 &AccumInner) {
                 AccumInner += f2(J1, J2, N1, N2);
              },
              SumInner);

          Kokkos::single(PerTeam(Team), [&]() { AccumOuter += SumInner; });
       },
       Sum1);

   if (Sum1 != RefSum) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelReduce1DReduce1D Sum FAIL", N1));
   }

   I4 Max1;
   parallelReduceOuter(
       {N1},
       KOKKOS_LAMBDA(int J1, const TeamMember &Team, I4 &AccumOuter) {
          I4 MaxInner;
          parallelReduceInner(
              Team, J1,
              INNER_LAMBDA(int J2, I4 &AccumInner) {
                 AccumInner = Kokkos::max(AccumInner, f2(J1, J2, N1, N2));
              },
              Kokkos::Max<I4>(MaxInner));
          AccumOuter = Kokkos::max(AccumOuter, MaxInner);
       },
       Kokkos::Max<I4>(Max1));

   if (Max1 != RefMax) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelReduce1DReduce1D Max FAIL", N1));
   }

   I4 Sum2, Max2;
   parallelReduceOuter(
       {N1},
       KOKKOS_LAMBDA(int J1, const TeamMember &Team, I4 &AccumSumOuter,
                     I4 &AccumMaxOuter) {
          I4 SumInner, MaxInner;
          parallelReduceInner(
              Team, J1,
              INNER_LAMBDA(int J2, I4 &AccumSumInner, I4 &AccumMaxInner) {
                 AccumSumInner += f2(J1, J2, N1, N2);
                 AccumMaxInner = Kokkos::max(AccumMaxInner, f2(J1, J2, N1, N2));
              },
              SumInner, Kokkos::Max<I4>(MaxInner));

          Kokkos::single(PerTeam(Team), [&]() { AccumSumOuter += SumInner; });
          AccumMaxOuter = Kokkos::max(AccumMaxOuter, MaxInner);
       },
       Sum2, Kokkos::Max<I4>(Max2));

   if (Sum2 != RefSum || Max2 != RefMax) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelReduce1DReduce1D Sum+Max FAIL", N1));
   }

   return Err;
}

Error testHiparFor1DMultiple1D(int N1, int N2) {
   Error Err;

   const int N2P1 = N2 + 1;

   HostArray2DI4 RefAH("RefA2H", N1, N2P1);
   HostArray2DI4 RefBH("RefB2H", N1, N2);
   HostArray2DI4 RefCH("RefC2H", N1, N2);
   HostArray1DI4 RefDH("RefD1H", N1);

   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2P1; ++J2) {
         RefAH(J1, J2) = f2(J1, J2, N1, N2P1);
      }

      for (int J2 = 0; J2 < N2; ++J2) {
         RefBH(J1, J2) = (RefAH(J1, J2) + RefAH(J1, J2 + 1)) / 2;
      }

      Real ScanAccum = 0;
      for (int J2 = 0; J2 < N2; ++J2) {
         ScanAccum += RefBH(J1, J2);
         RefCH(J1, J2) = ScanAccum;
      }

      if (J1 % 2 == 0) {
         for (int J2 = 0; J2 < N2; ++J2) {
            RefCH(J1, J2) += 1;
         }
      }

      Real SumAccum = 0;
      for (int J2 = 0; J2 < N2; ++J2) {
         SumAccum += RefCH(J1, J2);
      }
      RefDH(J1) = SumAccum;
   }

   Array2DI4 A("A2", N1, N2P1);
   Array2DI4 B("B2", N1, N2);
   Array2DI4 C("C2", N1, N2);
   Array1DI4 D("D1", N1);

   parallelForOuter(
       {N1}, KOKKOS_LAMBDA(int J1, const TeamMember &Team) {
          parallelForInner(
              Team, N2P1,
              INNER_LAMBDA(int J2) { A(J1, J2) = f2(J1, J2, N1, N2P1); });

          teamBarrier(Team);

          parallelForInner(
              Team, N2, INNER_LAMBDA(int J2) {
                 B(J1, J2) = (A(J1, J2) + A(J1, J2 + 1)) / 2;
              });

          teamBarrier(Team);

          parallelScanInner(
              Team, N2, INNER_LAMBDA(int J2, I4 &Accum, bool IsFinal) {
                 Accum += B(J1, J2);
                 if (IsFinal) {
                    C(J1, J2) = Accum;
                 }
              });

          if (J1 % 2 == 0) {
             parallelForInner(
                 Team, N2, INNER_LAMBDA(int J2) { C(J1, J2) += 1; });
          }

          teamBarrier(Team);

          parallelReduceInner(
              Team, N2, INNER_LAMBDA(int J2, I4 &Accum) { Accum += C(J1, J2); },
              D(J1));
       });

   auto DH = createHostMirrorCopy(D);

   if (!hostArraysEqual(DH, RefDH)) {
      Err += Error(ErrorCode::Fail, errorMsg("multiple patterns FAIL", N1, N2));
   }

   return Err;
}

Error testHiparFor2DFor1D(int N1, int N2) {
   Error Err;

   const int N3 = N1 + N2 - 1;

   HostArray3DI4 RefAH("RefA3H", N1, N2, N3);
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < J1 + J2; ++J3) {
            RefAH(J1, J2, J3) = f3(J1, J2, J3, N1, N2, N3);
         }
      }
   }

   Array3DI4 A("A3", N1, N2, N3);
   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
          parallelForInner(
              Team, J1 + J2, INNER_LAMBDA(int J3) {
                 A(J1, J2, J3) = f3(J1, J2, J3, N1, N2, N3);
              });
       });

   auto AH = createHostMirrorCopy(A);

   if (!hostArraysEqual(AH, RefAH)) {
      Err +=
          Error(ErrorCode::Fail, errorMsg("parallelFor2DFor1D FAIL", N1, N2));
   }

   return Err;
}

Error testHiparFor2DReduce1D(int N1, int N2) {
   Error Err;

   const int N3 = N1 + N2 - 1;

   HostArray2DI4 RefSumH("RefSum2H", N1, N2);
   HostArray2DI4 RefMaxH("RefMax2H", N1, N2);
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         I4 Sum = 0;
         I4 Max = std::numeric_limits<I4>::min();
         for (int J3 = 0; J3 < J1 + J2; ++J3) {
            Sum += f3(J1, J2, J3, N1, N2, N3);
            Max = std::max(Max, f3(J1, J2, J3, N1, N2, N3));
         }
         RefSumH(J1, J2) = Sum;
         RefMaxH(J1, J2) = Max;
      }
   }

   Array2DI4 Sum1("Sum12", N1, N2);
   Array2DI4 Max1("Max12", N1, N2);
   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
          I4 Sum;
          parallelReduceInner(
              Team, J1 + J2,
              INNER_LAMBDA(int J3, I4 &Accum) {
                 Accum += f3(J1, J2, J3, N1, N2, N3);
              },
              Sum);
          Sum1(J1, J2) = Sum;

          I4 Max;
          parallelReduceInner(
              Team, J1 + J2,
              INNER_LAMBDA(int J3, I4 &Accum) {
                 Accum = Kokkos::max(Accum, f3(J1, J2, J3, N1, N2, N3));
              },
              Kokkos::Max<I4>(Max));
          Max1(J1, J2) = Max;
       });

   auto Sum1H = createHostMirrorCopy(Sum1);
   if (!hostArraysEqual(Sum1H, RefSumH)) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelFor2DReduce1D Sum FAIL", N1, N2));
   }

   auto Max1H = createHostMirrorCopy(Max1);
   if (!hostArraysEqual(Max1H, RefMaxH)) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelFor2DReduce1D Max FAIL", N1, N2));
   }

   Array2DI4 Sum2("Sum22", N1, N2);
   Array2DI4 Max2("Max22", N1, N2);
   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
          I4 Sum, Max;
          parallelReduceInner(
              Team, J1 + J2,
              INNER_LAMBDA(int J3, I4 &AccumSum, I4 &AccumMax) {
                 AccumSum += f3(J1, J2, J3, N1, N2, N3);
                 AccumMax = Kokkos::max(AccumMax, f3(J1, J2, J3, N1, N2, N3));
              },
              Sum, Kokkos::Max<I4>(Max));
          Sum2(J1, J2) = Sum;
          Max2(J1, J2) = Max;
       });

   auto Sum2H = createHostMirrorCopy(Sum2);
   auto Max2H = createHostMirrorCopy(Max2);
   if (!hostArraysEqual(Sum2H, RefSumH) || !hostArraysEqual(Max2H, RefMaxH)) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelFor2DReduce1D Sum+Max FAIL", N1, N2));
   }

   return Err;
}

Error testHiparFor2DScan1D(int N1, int N2) {
   Error Err;

   const int N3 = N1 + N2 - 1;

   HostArray3DI4 RefRSumH("RefRSum3H", N1, N2, N3);
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         I4 RSum = 0;
         for (int J3 = 0; J3 < J1 + J2; ++J3) {
            RefRSumH(J1, J2, J3) = RSum;
            RSum += f3(J1, J2, J3, N1, N2, N3);
         }
      }
   }

   Array3DI4 RSum("RSum3", N1, N2, N3);
   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
          parallelScanInner(
              Team, J1 + J2, INNER_LAMBDA(int J3, I4 &Accum, bool IsFinal) {
                 if (IsFinal) {
                    RSum(J1, J2, J3) = Accum;
                 }
                 Accum += f3(J1, J2, J3, N1, N2, N3);
              });
       });

   auto RSumH = createHostMirrorCopy(RSum);
   if (!hostArraysEqual(RSumH, RefRSumH)) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelFor2DScan1D Sum FAIL", N1, N2));
   }

   return Err;
}

Error testHiparReduce2DReduce1D(int N1, int N2) {
   Error Err;

   const int N3 = N1 + N2 - 1;

   I4 RefSum = 0;
   I4 RefMax = std::numeric_limits<I4>::min();
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < J1 + J2; ++J3) {
            RefSum += f3(J1, J2, J3, N1, N2, N3);
            RefMax = std::max(RefMax, f3(J1, J2, J3, N1, N2, N3));
         }
      }
   }

   I4 Sum1;
   parallelReduceOuter(
       {N1, N2},
       KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team, I4 &AccumOuter) {
          I4 SumInner;
          parallelReduceInner(
              Team, J1 + J2,
              INNER_LAMBDA(int J3, I4 &AccumInner) {
                 AccumInner += f3(J1, J2, J3, N1, N2, N3);
              },
              SumInner);

          Kokkos::single(PerTeam(Team), [&]() { AccumOuter += SumInner; });
       },
       Sum1);

   if (Sum1 != RefSum) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelReduce2DReduce1D Sum FAIL", N1, N2));
   }

   I4 Max1;
   parallelReduceOuter(
       {N1, N2},
       KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team, I4 &AccumOuter) {
          I4 MaxInner;
          parallelReduceInner(
              Team, J1 + J2,
              INNER_LAMBDA(int J3, I4 &AccumInner) {
                 AccumInner =
                     Kokkos::max(AccumInner, f3(J1, J2, J3, N1, N2, N3));
              },
              Kokkos::Max<I4>(MaxInner));
          AccumOuter = Kokkos::max(AccumOuter, MaxInner);
       },
       Kokkos::Max<I4>(Max1));

   if (Max1 != RefMax) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelReduce2DReduce1D Max FAIL", N1, N2));
   }

   I4 Sum2, Max2;
   parallelReduceOuter(
       {N1, N2},
       KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team, I4 &AccumSumOuter,
                     I4 &AccumMaxOuter) {
          I4 SumInner, MaxInner;
          parallelReduceInner(
              Team, J1 + J2,
              INNER_LAMBDA(int J3, I4 &AccumSumInner, I4 &AccumMaxInner) {
                 AccumSumInner += f3(J1, J2, J3, N1, N2, N3);
                 AccumMaxInner =
                     Kokkos::max(AccumMaxInner, f3(J1, J2, J3, N1, N2, N3));
              },
              SumInner, Kokkos::Max<I4>(MaxInner));

          Kokkos::single(PerTeam(Team), [&]() { AccumSumOuter += SumInner; });
          AccumMaxOuter = Kokkos::max(AccumMaxOuter, MaxInner);
       },
       Sum2, Kokkos::Max<I4>(Max2));

   if (Sum2 != RefSum || Max2 != RefMax) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("parallelReduce2DReduce1D Sum+Max FAIL", N1, N2));
   }

   return Err;
}

Error testHiparFor2DMultiple1D(int N1, int N2, int N3) {
   Error Err;

   const int N3P1 = N3 + 1;

   HostArray3DI4 RefAH("RefA3H", N1, N2, N3P1);
   HostArray3DI4 RefBH("RefB3H", N1, N2, N3);
   HostArray3DI4 RefCH("RefC3H", N1, N2, N3);
   HostArray2DI4 RefDH("RefD2H", N1, N2);

   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < N3P1; ++J3) {
            RefAH(J1, J2, J3) = f3(J1, J2, J3, N1, N2, N3P1);
         }

         for (int J3 = 0; J3 < N3; ++J3) {
            RefBH(J1, J2, J3) = (RefAH(J1, J2, J3) + RefAH(J1, J2, J3 + 1)) / 2;
         }

         Real ScanAccum = 0;
         for (int J3 = 0; J3 < N3; ++J3) {
            ScanAccum += RefBH(J1, J2, J3);
            RefCH(J1, J2, J3) = ScanAccum;
         }

         if (J1 % 2 == 0 && J2 % 2 == 1) {
            for (int J3 = 0; J3 < N3; ++J3) {
               RefCH(J1, J2, J3) += 1;
            }
         }

         Real SumAccum = 0;
         for (int J3 = 0; J3 < N3; ++J3) {
            SumAccum += RefCH(J1, J2, J3);
         }
         RefDH(J1, J2) = SumAccum;
      }
   }

   Array3DI4 A("A3", N1, N2, N3P1);
   Array3DI4 B("B3", N1, N2, N3);
   Array3DI4 C("C3", N1, N2, N3);
   Array2DI4 D("D2", N1, N2);

   parallelForOuter(
       {N1, N2}, KOKKOS_LAMBDA(int J1, int J2, const TeamMember &Team) {
          parallelForInner(
              Team, N3P1, INNER_LAMBDA(int J3) {
                 A(J1, J2, J3) = f3(J1, J2, J3, N1, N2, N3P1);
              });

          teamBarrier(Team);

          parallelForInner(
              Team, N3, INNER_LAMBDA(int J3) {
                 B(J1, J2, J3) = (A(J1, J2, J3) + A(J1, J2, J3 + 1)) / 2;
              });

          teamBarrier(Team);

          parallelScanInner(
              Team, N3, INNER_LAMBDA(int J3, I4 &Accum, bool IsFinal) {
                 Accum += B(J1, J2, J3);
                 if (IsFinal) {
                    C(J1, J2, J3) = Accum;
                 }
              });

          if (J1 % 2 == 0 && J2 % 2 == 1) {
             parallelForInner(
                 Team, N3, INNER_LAMBDA(int J3) { C(J1, J2, J3) += 1; });
          }

          teamBarrier(Team);

          parallelReduceInner(
              Team, N3,
              INNER_LAMBDA(int J3, I4 &Accum) { Accum += C(J1, J2, J3); },
              D(J1, J2));
       });

   auto DH = createHostMirrorCopy(D);

   if (!hostArraysEqual(DH, RefDH)) {
      Err += Error(ErrorCode::Fail,
                   errorMsg("multiple patterns FAIL", N1, N2, N3));
   }

   return Err;
}

int main(int argc, char **argv) {
   Error Err;

   MPI_Init(&argc, &argv);
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   initLogging(DefEnv);

   try {
      Kokkos::initialize(argc, argv);
      {
         for (auto N1 : {1, 4, 9, 31, 33, 63, 65}) {
            Err += testHiparFor1DFor1D(N1);
            Err += testHiparFor1DReduce1D(N1);
            Err += testHiparFor1DScan1D(N1);
// outer reduce tests fail with SYCL and old Kokkos
// remove this check once E3SM starts using Kokkos 4.7.1
#if !defined(KOKKOS_ENABLE_SYCL) || KOKKOS_VERSION_GREATER_EQUAL(4, 7, 1)
            Err += testHiparReduce1DReduce1D(N1);
#endif

            Err += testHiparFor1DMultiple1D(1, N1);
            Err += testHiparFor1DMultiple1D((N1 + 1) / 2, N1);
            Err += testHiparFor1DMultiple1D(2 * N1, N1);
         }

         for (auto N1 : {1, 5, 10}) {
            for (auto N2 : {1, 4, 9, 31, 33, 63, 65}) {
               Err += testHiparFor2DFor1D(N1, N2);
               Err += testHiparFor2DReduce1D(N1, N2);
               Err += testHiparFor2DScan1D(N1, N2);
// outer reduce tests fail with SYCL and old Kokkos
// remove this check once E3SM starts using Kokkos 4.7.1
#if !defined(KOKKOS_ENABLE_SYCL) || KOKKOS_VERSION_GREATER_EQUAL(4, 7, 1)
               Err += testHiparReduce2DReduce1D(N1, N2);
#endif
               Err += testHiparFor2DMultiple1D(1, N1, N2);
               Err += testHiparFor2DMultiple1D((N1 + 1) / 2, N1, N2);
               Err += testHiparFor2DMultiple1D(2 * N1, N1, N2);
            }
         }
      }
      Kokkos::finalize();
   } catch (const std::exception &Ex) {
      Err += Error(ErrorCode::Fail, Ex.what() + std::string(": FAIL"));
   } catch (...) {
      Err += Error(ErrorCode::Fail, "Unknown: FAIL");
   }

   CHECK_ERROR_ABORT(Err, "Kokkos Wrappers Unit Tests FAIL");

   MPI_Finalize();

   return 0;
}
