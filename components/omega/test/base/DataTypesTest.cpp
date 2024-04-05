//===-- Test driver for OMEGA data types -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA data types
///
/// This driver tests the definition of various file types for the OMEGA
/// model. In particular, it tests the length of various data types and
/// outputs a PASS if all of them are the expected length. It also tests
/// that a build with SINGLE_PRECISION converts the defaul real type to
/// single precision (4-byte) floating point.
///
//
//===-----------------------------------------------------------------------===/

#include <iostream>

#include "DataTypes.h"
#include "OmegaKokkos.h"
#include "mpi.h"

using namespace OMEGA;

int main(int argc, char *argv[]) {

   // initialize environments
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {

      // declare variables of each supported type
      I4 MyInt4   = 1;
      I8 MyInt8   = 2;
      R4 MyR4     = 3.0;
      R8 MyR8     = 4.0000000000001;
      Real MyReal = 5.000001;
      // using operator""_Real;
      auto MyRealLiteral = 1._Real;
      int SizeTmp        = 0;

      // Check expected size (in bytes) for data types
      SizeTmp = sizeof(MyInt4);
      if (SizeTmp == 4)
         std::cout << "Size of I4: PASS" << std::endl;
      else
         std::cout << "Size of I4: FAIL " << SizeTmp << std::endl;

      SizeTmp = sizeof(MyInt8);
      if (SizeTmp == 8)
         std::cout << "Size of I8: PASS" << std::endl;
      else
         std::cout << "Size of I8: FAIL " << SizeTmp << std::endl;

      SizeTmp = sizeof(MyR4);
      if (SizeTmp == 4)
         std::cout << "Size of R4: PASS" << std::endl;
      else
         std::cout << "Size of R4: FAIL " << SizeTmp << std::endl;

      SizeTmp = sizeof(MyR8);
      if (SizeTmp == 8)
         std::cout << "Size of R8: PASS" << std::endl;
      else
         std::cout << "Size of R8: FAIL " << SizeTmp << std::endl;

      SizeTmp = sizeof(MyReal);
#ifdef SINGLE_PRECISION
      if (SizeTmp == 4)
         std::cout << "Size of Real is 4: PASS" << std::endl;
      else
         std::cout << "Size of Real is 4: FAIL " << SizeTmp << std::endl;
#else
      if (SizeTmp == 8)
         std::cout << "Size of Real is 8: PASS" << std::endl;
      else
         std::cout << "Size of Real is 8: FAIL " << SizeTmp << std::endl;
#endif

      SizeTmp = sizeof(MyRealLiteral);
      if (SizeTmp == sizeof(Real))
         std::cout << "Size of Real literal: PASS" << std::endl;
      else
         std::cout << "Size of Real literal: FAIL " << SizeTmp << std::endl;

      // Test creation of device arrays and copying to/from host
      // by initializing on the device, copying to host and comparing with
      // a reference host array.

      int NumCells    = 100;
      int NumVertLvls = 100;
      int NumTracers  = 4;
      int NumTimeLvls = 2;
      int NumExtra    = 2;

      // Test for 1DI4
      Array1DI4 TstArr1DI4("TstArr1DI4", NumCells);
      HostArray1DI4 RefArr1DI4("RefArr1DI4", NumCells);

      for (int i = 0; i < NumCells; ++i) {
         RefArr1DI4(i) = i;
      }

      parallelFor({NumCells}, KOKKOS_LAMBDA(int i) { TstArr1DI4(i) = i; });

      Kokkos::fence();

      auto TstHost1DI4 = createHostMirrorCopy(TstArr1DI4);

      int icount = 0;
      for (int i = 0; i < NumCells; ++i) {
         if (TstHost1DI4(i) != RefArr1DI4(i))
            ++icount;
      }

      if (icount == 0)
         std::cout << "Kokkos 1DI4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 1DI4 test: FAIL" << std::endl;

      // Test for 2DI4
      Array2DI4 TstArr2DI4("TstArr2DI4", NumCells, NumVertLvls);
      HostArray2DI4 RefArr2DI4("RefArr2DI4", NumCells, NumVertLvls);

      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr2DI4(j, i) = i + j;
         }
      }

      parallelFor(
          {NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int j, int i) { TstArr2DI4(j, i) = i + j; });

      Kokkos::fence();

      auto TstHost2DI4 = createHostMirrorCopy(TstArr2DI4);

      icount = 0;
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost2DI4(j, i) != RefArr2DI4(j, i))
               ++icount;
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 2DI4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 2DI4 test: FAIL" << std::endl;

      // Test for 3DI4
      Array3DI4 TstArr3DI4("TstArr3DI4", NumTracers, NumCells, NumVertLvls);
      HostArray3DI4 RefArr3DI4("RefArr3DI4", NumTracers, NumCells, NumVertLvls);

      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               RefArr3DI4(k, j, i) = i + j + k;
            }
         }
      }

      parallelFor(
          {NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int k, int j, int i) {
             TstArr3DI4(k, j, i) = i + j + k;
          });

      Kokkos::fence();

      auto TstHost3DI4 = createHostMirrorCopy(TstArr3DI4);

      icount = 0;
      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               if (TstHost3DI4(k, j, i) != RefArr3DI4(k, j, i))
                  ++icount;
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 3DI4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 3DI4 test: FAIL" << std::endl;

      // Test for 4DI4
      Array4DI4 TstArr4DI4("TstArr4DI4", NumTimeLvls, NumTracers, NumCells,
                           NumVertLvls);
      HostArray4DI4 RefArr4DI4("RefArr4DI4", NumTimeLvls, NumTracers, NumCells,
                               NumVertLvls);

      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  RefArr4DI4(m, k, j, i) = i + j + k + m;
               }
            }
         }
      }

      parallelFor(
          {NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
             TstArr4DI4(m, k, j, i) = i + j + k + m;
          });

      Kokkos::fence();

      auto TstHost4DI4 = createHostMirrorCopy(TstArr4DI4);

      icount = 0;
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost4DI4(m, k, j, i) != RefArr4DI4(m, k, j, i))
                     ++icount;
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 4DI4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 4DI4 test: FAIL" << std::endl;

      // Test for 5DI4
      Array5DI4 TstArr5DI4("TstArr5DI4", NumExtra, NumTimeLvls, NumTracers,
                           NumCells, NumVertLvls);
      HostArray5DI4 RefArr5DI4("RefArr5DI4", NumExtra, NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);

      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     RefArr5DI4(n, m, k, j, i) = i + j + k + m + n;
                  }
               }
            }
         }
      }

      parallelFor(
          {NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int n, int m, int k, int j, int i) {
             TstArr5DI4(n, m, k, j, i) = i + j + k + m + n;
          });

      Kokkos::fence();

      auto TstHost5DI4 = createHostMirrorCopy(TstArr5DI4);

      icount = 0;
      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     if (TstHost5DI4(n, m, k, j, i) !=
                         RefArr5DI4(n, m, k, j, i))
                        ++icount;
                  }
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 5DI4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 5DI4 test: FAIL" << std::endl;

      // Test for 1DI8
      Array1DI8 TstArr1DI8("TstArr1DI8", NumCells);
      HostArray1DI8 RefArr1DI8("RefArr1DI8", NumCells);

      for (int i = 0; i < NumCells; ++i) {
         RefArr1DI8(i) = i;
      }

      parallelFor({NumCells}, KOKKOS_LAMBDA(int i) { TstArr1DI8(i) = i; });

      Kokkos::fence();

      auto TstHost1DI8 = createHostMirrorCopy(TstArr1DI8);

      icount = 0;
      for (int i = 0; i < NumCells; ++i) {
         if (TstHost1DI8(i) != RefArr1DI8(i))
            ++icount;
      }

      if (icount == 0)
         std::cout << "Kokkos 1DI8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 1DI8 test: FAIL" << std::endl;

      // Test for 2DI8
      Array2DI8 TstArr2DI8("TstArr2DI8", NumCells, NumVertLvls);
      HostArray2DI8 RefArr2DI8("RefArr2DI8", NumCells, NumVertLvls);

      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr2DI8(j, i) = i + j;
         }
      }

      parallelFor(
          {NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int j, int i) { TstArr2DI8(j, i) = i + j; });

      Kokkos::fence();

      auto TstHost2DI8 = createHostMirrorCopy(TstArr2DI8);

      icount = 0;
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost2DI8(j, i) != RefArr2DI8(j, i))
               ++icount;
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 2DI8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 2DI8 test: FAIL" << std::endl;

      // Test for 3DI8
      Array3DI8 TstArr3DI8("TstArr3DI8", NumTracers, NumCells, NumVertLvls);
      HostArray3DI8 RefArr3DI8("RefArr3DI8", NumTracers, NumCells, NumVertLvls);

      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               RefArr3DI8(k, j, i) = i + j + k;
            }
         }
      }

      parallelFor(
          {NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int k, int j, int i) {
             TstArr3DI8(k, j, i) = i + j + k;
          });

      Kokkos::fence();

      auto TstHost3DI8 = createHostMirrorCopy(TstArr3DI8);

      icount = 0;
      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               if (TstHost3DI8(k, j, i) != RefArr3DI8(k, j, i))
                  ++icount;
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 3DI8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 3DI8 test: FAIL" << std::endl;

      // Test for 4DI8
      Array4DI8 TstArr4DI8("TstArr4DI8", NumTimeLvls, NumTracers, NumCells,
                           NumVertLvls);
      HostArray4DI8 RefArr4DI8("RefArr4DI8", NumTimeLvls, NumTracers, NumCells,
                               NumVertLvls);

      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  RefArr4DI8(m, k, j, i) = i + j + k + m;
               }
            }
         }
      }

      parallelFor(
          {NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
             TstArr4DI8(m, k, j, i) = i + j + k + m;
          });

      Kokkos::fence();

      auto TstHost4DI8 = createHostMirrorCopy(TstArr4DI8);

      icount = 0;
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost4DI8(m, k, j, i) != RefArr4DI8(m, k, j, i))
                     ++icount;
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 4DI8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 4DI8 test: FAIL" << std::endl;

      // Test for 5DI8
      Array5DI8 TstArr5DI8("TstArr5DI8", NumExtra, NumTimeLvls, NumTracers,
                           NumCells, NumVertLvls);
      HostArray5DI8 RefArr5DI8("RefArr5DI8", NumExtra, NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);

      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     RefArr5DI8(n, m, k, j, i) = i + j + k + m + n;
                  }
               }
            }
         }
      }

      parallelFor(
          {NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int n, int m, int k, int j, int i) {
             TstArr5DI8(n, m, k, j, i) = i + j + k + m + n;
          });

      Kokkos::fence();

      auto TstHost5DI8 = createHostMirrorCopy(TstArr5DI8);

      icount = 0;
      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     if (TstHost5DI8(n, m, k, j, i) !=
                         RefArr5DI8(n, m, k, j, i))
                        ++icount;
                  }
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 5DI8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 5DI8 test: FAIL" << std::endl;

      // Test for 1DR4
      Array1DR4 TstArr1DR4("TstArr1DR4", NumCells);
      HostArray1DR4 RefArr1DR4("RefArr1DR4", NumCells);

      for (int i = 0; i < NumCells; ++i) {
         RefArr1DR4(i) = i;
      }

      parallelFor({NumCells}, KOKKOS_LAMBDA(int i) { TstArr1DR4(i) = i; });

      Kokkos::fence();

      auto TstHost1DR4 = createHostMirrorCopy(TstArr1DR4);

      icount = 0;
      for (int i = 0; i < NumCells; ++i) {
         if (TstHost1DR4(i) != RefArr1DR4(i))
            ++icount;
      }

      if (icount == 0)
         std::cout << "Kokkos 1DR4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 1DR4 test: FAIL" << std::endl;

      // Test for 2DR4
      Array2DR4 TstArr2DR4("TstArr2DR4", NumCells, NumVertLvls);
      HostArray2DR4 RefArr2DR4("RefArr2DR4", NumCells, NumVertLvls);

      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr2DR4(j, i) = i + j;
         }
      }

      parallelFor(
          {NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int j, int i) { TstArr2DR4(j, i) = i + j; });

      Kokkos::fence();

      auto TstHost2DR4 = createHostMirrorCopy(TstArr2DR4);

      icount = 0;
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost2DR4(j, i) != RefArr2DR4(j, i))
               ++icount;
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 2DR4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 2DR4 test: FAIL" << std::endl;

      // Test for 3DR4
      Array3DR4 TstArr3DR4("TstArr3DR4", NumTracers, NumCells, NumVertLvls);
      HostArray3DR4 RefArr3DR4("RefArr3DR4", NumTracers, NumCells, NumVertLvls);

      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               RefArr3DR4(k, j, i) = i + j + k;
            }
         }
      }

      parallelFor(
          {NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int k, int j, int i) {
             TstArr3DR4(k, j, i) = i + j + k;
          });

      Kokkos::fence();

      auto TstHost3DR4 = createHostMirrorCopy(TstArr3DR4);

      icount = 0;
      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               if (TstHost3DR4(k, j, i) != RefArr3DR4(k, j, i))
                  ++icount;
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 3DR4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 3DR4 test: FAIL" << std::endl;

      // Test for 4DR4
      Array4DR4 TstArr4DR4("TstArr4DR4", NumTimeLvls, NumTracers, NumCells,
                           NumVertLvls);
      HostArray4DR4 RefArr4DR4("RefArr4DR4", NumTimeLvls, NumTracers, NumCells,
                               NumVertLvls);

      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  RefArr4DR4(m, k, j, i) = i + j + k + m;
               }
            }
         }
      }

      parallelFor(
          {NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
             TstArr4DR4(m, k, j, i) = i + j + k + m;
          });

      Kokkos::fence();

      auto TstHost4DR4 = createHostMirrorCopy(TstArr4DR4);

      icount = 0;
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost4DR4(m, k, j, i) != RefArr4DR4(m, k, j, i))
                     ++icount;
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 4DR4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 4DR4 test: FAIL" << std::endl;

      // Test for 5DR4
      Array5DR4 TstArr5DR4("TstArr5DR4", NumExtra, NumTimeLvls, NumTracers,
                           NumCells, NumVertLvls);
      HostArray5DR4 RefArr5DR4("RefArr5DR4", NumExtra, NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);

      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     RefArr5DR4(n, m, k, j, i) = i + j + k + m + n;
                  }
               }
            }
         }
      }

      parallelFor(
          {NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int n, int m, int k, int j, int i) {
             TstArr5DR4(n, m, k, j, i) = i + j + k + m + n;
          });

      Kokkos::fence();

      auto TstHost5DR4 = createHostMirrorCopy(TstArr5DR4);

      icount = 0;
      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     if (TstHost5DR4(n, m, k, j, i) !=
                         RefArr5DR4(n, m, k, j, i))
                        ++icount;
                  }
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 5DR4 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 5DR4 test: FAIL" << std::endl;

      // Test for 1DR8
      Array1DR8 TstArr1DR8("TstArr1DR8", NumCells);
      HostArray1DR8 RefArr1DR8("RefArr1DR8", NumCells);

      for (int i = 0; i < NumCells; ++i) {
         RefArr1DR8(i) = i;
      }

      parallelFor({NumCells}, KOKKOS_LAMBDA(int i) { TstArr1DR8(i) = i; });

      Kokkos::fence();

      auto TstHost1DR8 = createHostMirrorCopy(TstArr1DR8);

      icount = 0;
      for (int i = 0; i < NumCells; ++i) {
         if (TstHost1DR8(i) != RefArr1DR8(i))
            ++icount;
      }

      if (icount == 0)
         std::cout << "Kokkos 1DR8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 1DR8 test: FAIL" << std::endl;

      // Test for 2DR8
      Array2DR8 TstArr2DR8("TstArr2DR8", NumCells, NumVertLvls);
      HostArray2DR8 RefArr2DR8("RefArr2DR8", NumCells, NumVertLvls);

      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr2DR8(j, i) = i + j;
         }
      }

      parallelFor(
          {NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int j, int i) { TstArr2DR8(j, i) = i + j; });

      Kokkos::fence();

      auto TstHost2DR8 = createHostMirrorCopy(TstArr2DR8);

      icount = 0;
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost2DR8(j, i) != RefArr2DR8(j, i))
               ++icount;
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 2DR8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 2DR8 test: FAIL" << std::endl;

      // Test for 3DR8
      Array3DR8 TstArr3DR8("TstArr3DR8", NumTracers, NumCells, NumVertLvls);
      HostArray3DR8 RefArr3DR8("RefArr3DR8", NumTracers, NumCells, NumVertLvls);

      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               RefArr3DR8(k, j, i) = i + j + k;
            }
         }
      }

      parallelFor(
          {NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int k, int j, int i) {
             TstArr3DR8(k, j, i) = i + j + k;
          });

      Kokkos::fence();

      auto TstHost3DR8 = createHostMirrorCopy(TstArr3DR8);

      icount = 0;
      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               if (TstHost3DR8(k, j, i) != RefArr3DR8(k, j, i))
                  ++icount;
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 3DR8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 3DR8 test: FAIL" << std::endl;

      // Test for 4DR8
      Array4DR8 TstArr4DR8("TstArr4DR8", NumTimeLvls, NumTracers, NumCells,
                           NumVertLvls);
      HostArray4DR8 RefArr4DR8("RefArr4DR8", NumTimeLvls, NumTracers, NumCells,
                               NumVertLvls);

      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  RefArr4DR8(m, k, j, i) = i + j + k + m;
               }
            }
         }
      }

      parallelFor(
          {NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
             TstArr4DR8(m, k, j, i) = i + j + k + m;
          });

      Kokkos::fence();

      auto TstHost4DR8 = createHostMirrorCopy(TstArr4DR8);

      icount = 0;
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost4DR8(m, k, j, i) != RefArr4DR8(m, k, j, i))
                     ++icount;
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 4DR8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 4DR8 test: FAIL" << std::endl;

      // Test for 5DR8
      Array5DR8 TstArr5DR8("TstArr5DR8", NumExtra, NumTimeLvls, NumTracers,
                           NumCells, NumVertLvls);
      HostArray5DR8 RefArr5DR8("RefArr5DR8", NumExtra, NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);

      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     RefArr5DR8(n, m, k, j, i) = i + j + k + m + n;
                  }
               }
            }
         }
      }

      parallelFor(
          {NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int n, int m, int k, int j, int i) {
             TstArr5DR8(n, m, k, j, i) = i + j + k + m + n;
          });

      Kokkos::fence();

      auto TstHost5DR8 = createHostMirrorCopy(TstArr5DR8);

      icount = 0;
      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     if (TstHost5DR8(n, m, k, j, i) !=
                         RefArr5DR8(n, m, k, j, i))
                        ++icount;
                  }
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 5DR8 test: PASS" << std::endl;
      else
         std::cout << "Kokkos 5DR8 test: FAIL" << std::endl;

      // Test for 1DReal
      Array1DReal TstArr1DReal("TstArr1DReal", NumCells);
      HostArray1DReal RefArr1DReal("RefArr1DReal", NumCells);

      for (int i = 0; i < NumCells; ++i) {
         RefArr1DReal(i) = i;
      }

      parallelFor({NumCells}, KOKKOS_LAMBDA(int i) { TstArr1DReal(i) = i; });

      Kokkos::fence();

      auto TstHost1DReal = createHostMirrorCopy(TstArr1DReal);

      icount = 0;
      for (int i = 0; i < NumCells; ++i) {
         if (TstHost1DReal(i) != RefArr1DReal(i))
            ++icount;
      }

      if (icount == 0)
         std::cout << "Kokkos 1DReal test: PASS" << std::endl;
      else
         std::cout << "Kokkos 1DReal test: FAIL" << std::endl;

      // Test for 2DReal
      Array2DReal TstArr2DReal("TstArr2DReal", NumCells, NumVertLvls);
      HostArray2DReal RefArr2DReal("RefArr2DReal", NumCells, NumVertLvls);

      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr2DReal(j, i) = i + j;
         }
      }

      parallelFor(
          {NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int j, int i) { TstArr2DReal(j, i) = i + j; });

      Kokkos::fence();

      auto TstHost2DReal = createHostMirrorCopy(TstArr2DReal);

      icount = 0;
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost2DReal(j, i) != RefArr2DReal(j, i))
               ++icount;
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 2DReal test: PASS" << std::endl;
      else
         std::cout << "Kokkos 2DReal test: FAIL" << std::endl;

      // Test for 3DReal
      Array3DReal TstArr3DReal("TstArr3DReal", NumTracers, NumCells,
                               NumVertLvls);
      HostArray3DReal RefArr3DReal("RefArr3DReal", NumTracers, NumCells,
                                   NumVertLvls);

      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               RefArr3DReal(k, j, i) = i + j + k;
            }
         }
      }

      parallelFor(
          {NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int k, int j, int i) {
             TstArr3DReal(k, j, i) = i + j + k;
          });

      Kokkos::fence();

      auto TstHost3DReal = createHostMirrorCopy(TstArr3DReal);

      icount = 0;
      for (int k = 0; k < NumTracers; ++k) {
         for (int j = 0; j < NumCells; ++j) {
            for (int i = 0; i < NumVertLvls; ++i) {
               if (TstHost3DReal(k, j, i) != RefArr3DReal(k, j, i))
                  ++icount;
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 3DReal test: PASS" << std::endl;
      else
         std::cout << "Kokkos 3DReal test: FAIL" << std::endl;

      // Test for 4DReal
      Array4DReal TstArr4DReal("TstArr4DReal", NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);
      HostArray4DReal RefArr4DReal("RefArr4DReal", NumTimeLvls, NumTracers,
                                   NumCells, NumVertLvls);

      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  RefArr4DReal(m, k, j, i) = i + j + k + m;
               }
            }
         }
      }

      parallelFor(
          {NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int m, int k, int j, int i) {
             TstArr4DReal(m, k, j, i) = i + j + k + m;
          });

      Kokkos::fence();

      auto TstHost4DReal = createHostMirrorCopy(TstArr4DReal);

      icount = 0;
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost4DReal(m, k, j, i) != RefArr4DReal(m, k, j, i))
                     ++icount;
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 4DReal test: PASS" << std::endl;
      else
         std::cout << "Kokkos 4DReal test: FAIL" << std::endl;

      // Test for 5DReal
      Array5DReal TstArr5DReal("TstArr5DReal", NumExtra, NumTimeLvls,
                               NumTracers, NumCells, NumVertLvls);
      HostArray5DReal RefArr5DReal("RefArr5DReal", NumExtra, NumTimeLvls,
                                   NumTracers, NumCells, NumVertLvls);

      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     RefArr5DReal(n, m, k, j, i) = i + j + k + m + n;
                  }
               }
            }
         }
      }

      parallelFor(
          {NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int n, int m, int k, int j, int i) {
             TstArr5DReal(n, m, k, j, i) = i + j + k + m + n;
          });

      Kokkos::fence();

      auto TstHost5DReal = createHostMirrorCopy(TstArr5DReal);

      icount = 0;
      for (int n = 0; n < NumExtra; ++n) {
         for (int m = 0; m < NumTimeLvls; ++m) {
            for (int k = 0; k < NumTracers; ++k) {
               for (int j = 0; j < NumCells; ++j) {
                  for (int i = 0; i < NumVertLvls; ++i) {
                     if (TstHost5DReal(n, m, k, j, i) !=
                         RefArr5DReal(n, m, k, j, i))
                        ++icount;
                  }
               }
            }
         }
      }

      if (icount == 0)
         std::cout << "Kokkos 5DReal test: PASS" << std::endl;
      else
         std::cout << "Kokkos 5DReal test: FAIL" << std::endl;

      // finalize environments
      // MPI_Status status;
   }
   Kokkos::finalize();
   MPI_Finalize();

} // end of main
//===-----------------------------------------------------------------------===/
