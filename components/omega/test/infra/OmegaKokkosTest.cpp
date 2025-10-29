//===-- Test driver for OMEGA Kokkos -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Kokkos
///
/// This driver tests the Omega wrappers of Kokkos capabilities.
///
//
//===-----------------------------------------------------------------------===/

#include "OmegaKokkos.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"

#include "mpi.h"
#include <limits>

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

static KOKKOS_FUNCTION int f1(int J1, int N1) { return -N1 / 4 + J1; }

static KOKKOS_FUNCTION int f2(int J1, int J2, int N1, int N2) {
   return -(N1 * N2) / 4 + J1 + N1 * J2;
}

static KOKKOS_FUNCTION int f3(int J1, int J2, int J3, int N1, int N2, int N3) {
   return -(N1 * N2 * N3) / 4 + J1 + N1 * (J2 + N2 * J3);
}

static KOKKOS_FUNCTION int f4(int J1, int J2, int J3, int J4, int N1, int N2,
                              int N3, int N4) {
   return -(N1 * N2 * N3 * N4) / 4 + J1 + N1 * (J2 + N2 * (J3 + N3 * J4));
}

static KOKKOS_FUNCTION int f5(int J1, int J2, int J3, int J4, int J5, int N1,
                              int N2, int N3, int N4, int N5) {
   return -(N1 * N2 * N3 * N4 * N5) / 4 + J1 +
          N1 * (J2 + N2 * (J3 + N3 * (J4 + N4 * J5)));
}

Error testParallelFor1D() {
   Error Err;

   const int N1 = 127;

   HostArray1DI4 RefAH("RefA1H", N1);
   for (int J1 = 0; J1 < N1; ++J1) {
      RefAH(J1) = f1(J1, N1);
   }

   Array1DI4 A("A1", N1);
   parallelFor({N1}, KOKKOS_LAMBDA(int J1) { A(J1) = f1(J1, N1); });

   auto AH = createHostMirrorCopy(A);

   if (!hostArraysEqual(AH, RefAH)) {
      Err += Error(ErrorCode::Fail, "parallelFor 1D FAIL");
   }

   return Err;
}

Error testParallelReduce1D() {
   Error Err;

   const int N1 = 129;

   I4 RefSum = 0;
   I4 RefMax = std::numeric_limits<I4>::min();

   for (int J1 = 0; J1 < N1; ++J1) {
      RefSum += f1(J1, N1);
      RefMax = std::max(RefMax, f1(J1, N1));
   }

   // Test simple sum reduction
   I4 Sum1;
   parallelReduce(
       {N1}, KOKKOS_LAMBDA(int J1, I4 &Accum) { Accum += f1(J1, N1); }, Sum1);

   if (Sum1 != RefSum) {
      Err += Error(ErrorCode::Fail, "parallelReduce 1D sum FAIL");
   }

   // Test simple max reduction
   I4 Max1;
   parallelReduce(
       {N1},
       KOKKOS_LAMBDA(int J1, I4 &Accum) {
          Accum = Kokkos::max(Accum, f1(J1, N1));
       },
       Kokkos::Max<I4>(Max1));

   if (Max1 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 1D max FAIL");
   }

   // Test doing both reductions together
   I4 Sum2, Max2;
   parallelReduce(
       {N1},
       KOKKOS_LAMBDA(int J1, I4 &SumAccum, I4 &MaxAccum) {
          SumAccum += f1(J1, N1);
          MaxAccum = Kokkos::max(MaxAccum, f1(J1, N1));
       },
       Sum2, Kokkos::Max<I4>(Max2));

   if (Sum2 != RefSum || Max2 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 1D sum+max FAIL");
   }

   return Err;
}

Error testParallelFor2D() {
   Error Err;

   const int N1 = 63;
   const int N2 = 31;

   HostArray2DI4 RefAH("RefA2H", N1, N2);
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         RefAH(J1, J2) = f2(J1, J2, N1, N2);
      }
   }

   Array2DI4 A("A2", N1, N2);
   parallelFor(
       {N1, N2},
       KOKKOS_LAMBDA(int J1, int J2) { A(J1, J2) = f2(J1, J2, N1, N2); });

   auto AH = createHostMirrorCopy(A);

   if (!hostArraysEqual(AH, RefAH)) {
      Err += Error(ErrorCode::Fail, "parallelFor 2D FAIL");
   }

   return Err;
}

Error testParallelReduce2D() {
   Error Err;

   const int N1 = 33;
   const int N2 = 65;

   I4 RefSum = 0;
   I4 RefMax = std::numeric_limits<I4>::min();

   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         RefSum += f2(J1, J2, N1, N2);
         RefMax = std::max(RefMax, f2(J1, J2, N1, N2));
      }
   }

   // Test simple sum reduction
   I4 Sum1;
   parallelReduce(
       {N1, N2},
       KOKKOS_LAMBDA(int J1, int J2, I4 &Accum) {
          Accum += f2(J1, J2, N1, N2);
       },
       Sum1);

   if (Sum1 != RefSum) {
      Err += Error(ErrorCode::Fail, "parallelReduce 2D sum FAIL");
   }

   // Test simple max reduction
   I4 Max1;
   parallelReduce(
       {N1, N2},
       KOKKOS_LAMBDA(int J1, int J2, I4 &Accum) {
          Accum = Kokkos::max(Accum, f2(J1, J2, N1, N2));
       },
       Kokkos::Max<I4>(Max1));

   if (Max1 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 2D max FAIL");
   }

   // Test doing both reductions together
   I4 Sum2, Max2;
   parallelReduce(
       {N1, N2},
       KOKKOS_LAMBDA(int J1, int J2, I4 &SumAccum, I4 &MaxAccum) {
          SumAccum += f2(J1, J2, N1, N2);
          MaxAccum = Kokkos::max(MaxAccum, f2(J1, J2, N1, N2));
       },
       Sum2, Kokkos::Max<I4>(Max2));

   if (Sum2 != RefSum || Max2 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 2D sum+max FAIL");
   }

   return Err;
}

Error testParallelFor3D() {
   Error Err;

   const int N1 = 2;
   const int N2 = 100;
   const int N3 = 60;

   HostArray3DI4 RefAH("RefA3H", N1, N2, N3);
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < N3; ++J3) {
            RefAH(J1, J2, J3) = f3(J1, J2, J3, N1, N2, N3);
         }
      }
   }

   Array3DI4 A("A3", N1, N2, N3);
   parallelFor(
       {N1, N2, N3}, KOKKOS_LAMBDA(int J1, int J2, int J3) {
          A(J1, J2, J3) = f3(J1, J2, J3, N1, N2, N3);
       });

   auto AH = createHostMirrorCopy(A);

   if (!hostArraysEqual(AH, RefAH)) {
      Err += Error(ErrorCode::Fail, "parallelFor 3D FAIL");
   }

   return Err;
}

Error testParallelReduce3D() {
   Error Err;

   const int N1 = 10;
   const int N2 = 1;
   const int N3 = 60;

   I4 RefSum = 0;
   I4 RefMax = std::numeric_limits<I4>::min();

   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < N3; ++J3) {
            RefSum += f3(J1, J2, J3, N1, N2, N3);
            RefMax = std::max(RefMax, f3(J1, J2, J3, N1, N2, N3));
         }
      }
   }

   // Test simple sum reduction
   I4 Sum1;
   parallelReduce(
       {N1, N2, N3},
       KOKKOS_LAMBDA(int J1, int J2, int J3, I4 &Accum) {
          Accum += f3(J1, J2, J3, N1, N2, N3);
       },
       Sum1);

   if (Sum1 != RefSum) {
      Err += Error(ErrorCode::Fail, "parallelReduce 3D sum FAIL");
   }

   // Test simple max reduction
   I4 Max1;
   parallelReduce(
       {N1, N2, N3},
       KOKKOS_LAMBDA(int J1, int J2, int J3, I4 &Accum) {
          Accum = Kokkos::max(Accum, f3(J1, J2, J3, N1, N2, N3));
       },
       Kokkos::Max<I4>(Max1));

   if (Max1 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 3D max FAIL");
   }

   // Test doing both reductions together
   I4 Sum2, Max2;
   parallelReduce(
       {N1, N2, N3},
       KOKKOS_LAMBDA(int J1, int J2, int J3, I4 &SumAccum, I4 &MaxAccum) {
          SumAccum += f3(J1, J2, J3, N1, N2, N3);
          MaxAccum = Kokkos::max(MaxAccum, f3(J1, J2, J3, N1, N2, N3));
       },
       Sum2, Kokkos::Max<I4>(Max2));

   if (Sum2 != RefSum || Max2 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 3D sum+max FAIL");
   }

   return Err;
}

Error testParallelFor4D() {
   Error Err;

   const int N1 = 2;
   const int N2 = 4;
   const int N3 = 8;
   const int N4 = 16;

   HostArray4DI4 RefAH("RefA4H", N1, N2, N3, N4);
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < N3; ++J3) {
            for (int J4 = 0; J4 < N4; ++J4) {
               RefAH(J1, J2, J3, J4) = f4(J1, J2, J3, J4, N1, N2, N3, N4);
            }
         }
      }
   }

   Array4DI4 A("A4", N1, N2, N3, N4);
   parallelFor(
       {N1, N2, N3, N4}, KOKKOS_LAMBDA(int J1, int J2, int J3, int J4) {
          A(J1, J2, J3, J4) = f4(J1, J2, J3, J4, N1, N2, N3, N4);
       });

   auto AH = createHostMirrorCopy(A);

   if (!hostArraysEqual(AH, RefAH)) {
      Err += Error(ErrorCode::Fail, "parallelFor 4D FAIL");
   }

   return Err;
}

Error testParallelReduce4D() {
   Error Err;

   const int N1 = 3;
   const int N2 = 5;
   const int N3 = 7;
   const int N4 = 11;

   I4 RefSum = 0;
   I4 RefMax = std::numeric_limits<I4>::min();

   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < N3; ++J3) {
            for (int J4 = 0; J4 < N4; ++J4) {
               RefSum += f4(J1, J2, J3, J4, N1, N2, N3, N4);
               RefMax = std::max(RefMax, f4(J1, J2, J3, J4, N1, N2, N3, N4));
            }
         }
      }
   }

   // Test simple sum reduction
   I4 Sum1;
   parallelReduce(
       {N1, N2, N3, N4},
       KOKKOS_LAMBDA(int J1, int J2, int J3, int J4, I4 &Accum) {
          Accum += f4(J1, J2, J3, J4, N1, N2, N3, N4);
       },
       Sum1);

   if (Sum1 != RefSum) {
      Err += Error(ErrorCode::Fail, "parallelReduce 4D sum FAIL");
   }

   // Test simple max reduction
   I4 Max1;
   parallelReduce(
       {N1, N2, N3, N4},
       KOKKOS_LAMBDA(int J1, int J2, int J3, int J4, I4 &Accum) {
          Accum = Kokkos::max(Accum, f4(J1, J2, J3, J4, N1, N2, N3, N4));
       },
       Kokkos::Max<I4>(Max1));

   if (Max1 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 4D max FAIL");
   }

   // Test doing both reductions together
   I4 Sum2, Max2;
   parallelReduce(
       {N1, N2, N3, N4},
       KOKKOS_LAMBDA(int J1, int J2, int J3, int J4, I4 &SumAccum,
                     I4 &MaxAccum) {
          SumAccum += f4(J1, J2, J3, J4, N1, N2, N3, N4);
          MaxAccum = Kokkos::max(MaxAccum, f4(J1, J2, J3, J4, N1, N2, N3, N4));
       },
       Sum2, Kokkos::Max<I4>(Max2));

   if (Sum2 != RefSum || Max2 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 4D sum+max FAIL");
   }

   return Err;
}

Error testParallelFor5D() {
   Error Err;

   const int N1 = 33;
   const int N2 = 1;
   const int N3 = 3;
   const int N4 = 2;
   const int N5 = 129;

   HostArray5DI4 RefAH("RefA5H", N1, N2, N3, N4, N5);
   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < N3; ++J3) {
            for (int J4 = 0; J4 < N4; ++J4) {
               for (int J5 = 0; J5 < N5; ++J5) {
                  RefAH(J1, J2, J3, J4, J5) =
                      f5(J1, J2, J3, J4, J5, N1, N2, N3, N4, N5);
               }
            }
         }
      }
   }

   Array5DI4 A("A5", N1, N2, N3, N4, N5);
   parallelFor(
       {N1, N2, N3, N4, N5},
       KOKKOS_LAMBDA(int J1, int J2, int J3, int J4, int J5) {
          A(J1, J2, J3, J4, J5) = f5(J1, J2, J3, J4, J5, N1, N2, N3, N4, N5);
       });

   auto AH = createHostMirrorCopy(A);

   if (!hostArraysEqual(AH, RefAH)) {
      Err += Error(ErrorCode::Fail, "parallelFor 5D FAIL");
   }

   return Err;
}

Error testParallelReduce5D() {
   Error Err;

   const int N1 = 2;
   const int N2 = 31;
   const int N3 = 1;
   const int N4 = 127;
   const int N5 = 3;

   I4 RefSum = 0;
   I4 RefMax = std::numeric_limits<I4>::min();

   for (int J1 = 0; J1 < N1; ++J1) {
      for (int J2 = 0; J2 < N2; ++J2) {
         for (int J3 = 0; J3 < N3; ++J3) {
            for (int J4 = 0; J4 < N4; ++J4) {
               for (int J5 = 0; J5 < N5; ++J5) {
                  RefSum += f5(J1, J2, J3, J4, J5, N1, N2, N3, N4, N5);
                  RefMax = std::max(RefMax,
                                    f5(J1, J2, J3, J4, J5, N1, N2, N3, N4, N5));
               }
            }
         }
      }
   }

   // Test simple sum reduction
   I4 Sum1;
   parallelReduce(
       {N1, N2, N3, N4, N5},
       KOKKOS_LAMBDA(int J1, int J2, int J3, int J4, int J5, I4 &Accum) {
          Accum += f5(J1, J2, J3, J4, J5, N1, N2, N3, N4, N5);
       },
       Sum1);

   if (Sum1 != RefSum) {
      Err += Error(ErrorCode::Fail, "parallelReduce 5D sum FAIL");
   }

   // Test simple max reduction
   I4 Max1;
   parallelReduce(
       {N1, N2, N3, N4, N5},
       KOKKOS_LAMBDA(int J1, int J2, int J3, int J4, int J5, I4 &Accum) {
          Accum =
              Kokkos::max(Accum, f5(J1, J2, J3, J4, J5, N1, N2, N3, N4, N5));
       },
       Kokkos::Max<I4>(Max1));

   if (Max1 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 5D max FAIL");
   }

   // Test doing both reductions together
   I4 Sum2, Max2;
   parallelReduce(
       {N1, N2, N3, N4, N5},
       KOKKOS_LAMBDA(int J1, int J2, int J3, int J4, int J5, I4 &SumAccum,
                     I4 &MaxAccum) {
          SumAccum += f5(J1, J2, J3, J4, J5, N1, N2, N3, N4, N5);
          MaxAccum =
              Kokkos::max(MaxAccum, f5(J1, J2, J3, J4, J5, N1, N2, N3, N4, N5));
       },
       Sum2, Kokkos::Max<I4>(Max2));

   if (Sum2 != RefSum || Max2 != RefMax) {
      Err += Error(ErrorCode::Fail, "parallelReduce 5D sum+max FAIL");
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
         Err += testParallelFor1D();
         Err += testParallelReduce1D();

         Err += testParallelFor2D();
         Err += testParallelReduce2D();

         Err += testParallelFor3D();
         Err += testParallelReduce3D();

         Err += testParallelFor4D();
         Err += testParallelReduce4D();

         Err += testParallelFor5D();
         Err += testParallelReduce5D();
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
