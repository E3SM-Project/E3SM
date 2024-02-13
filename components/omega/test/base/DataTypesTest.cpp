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
#include "mpi.h"

int main(int argc, char *argv[]) {

   // initialize environments
   MPI_Init(&argc, &argv);
   yakl::init();

   // declare variables of each supported type
   OMEGA::I4 MyInt4   = 1;
   OMEGA::I8 MyInt8   = 2;
   OMEGA::R4 MyR4     = 3.0;
   OMEGA::R8 MyR8     = 4.0000000000001;
   OMEGA::Real MyReal = 5.000001;
   using OMEGA::operator""_Real;
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
   if (SizeTmp == sizeof(OMEGA::Real))
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

   using yakl::c::Bounds;
   using yakl::c::parallel_for;

   // Test for 1DI4
   OMEGA::Array1DI4 TstArr1DI4("TstArr", NumCells);
   OMEGA::ArrayHost1DI4 RefArr1DI4("RefArr", NumCells);

   for (int i = 0; i < NumCells; ++i) {
      RefArr1DI4(i) = i;
   }

   parallel_for(
       Bounds<1>(NumCells), YAKL_LAMBDA(int i) { TstArr1DI4(i) = i; });

   yakl::fence();
   auto TstHost1DI4 = TstArr1DI4.createHostCopy();

   int icount = 0;
   for (int i = 0; i < NumCells; ++i) {
      if (TstHost1DI4(i) != RefArr1DI4(i))
         ++icount;
   }
   TstHost1DI4.deallocate();
   RefArr1DI4.deallocate();
   TstArr1DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DI4 test: FAIL" << std::endl;

   // Test for 2DI4
   OMEGA::Array2DI4 TstArr2DI4("TstArr", NumCells, NumVertLvls);
   OMEGA::ArrayHost2DI4 RefArr2DI4("RefArr", NumCells, NumVertLvls);

   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         RefArr2DI4(j, i) = i + j;
      }
   }

   parallel_for(
       Bounds<2>(NumCells, NumVertLvls),
       YAKL_LAMBDA(int j, int i) { TstArr2DI4(j, i) = i + j; });

   yakl::fence();
   auto TstHost2DI4 = TstArr2DI4.createHostCopy();

   icount = 0;
   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         if (TstHost2DI4(j, i) != RefArr2DI4(j, i))
            ++icount;
      }
   }
   TstHost2DI4.deallocate();
   RefArr2DI4.deallocate();
   TstArr2DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DI4 test: FAIL" << std::endl;

   // Test for 3DI4
   OMEGA::Array3DI4 TstArr3DI4("TstArr", NumTracers, NumCells, NumVertLvls);
   OMEGA::ArrayHost3DI4 RefArr3DI4("RefArr", NumTracers, NumCells, NumVertLvls);

   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr3DI4(k, j, i) = i + j + k;
         }
      }
   }

   parallel_for(
       Bounds<3>(NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int k, int j, int i) { TstArr3DI4(k, j, i) = i + j + k; });

   yakl::fence();
   auto TstHost3DI4 = TstArr3DI4.createHostCopy();

   icount = 0;
   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost3DI4(k, j, i) != RefArr3DI4(k, j, i))
               ++icount;
         }
      }
   }
   TstHost3DI4.deallocate();
   RefArr3DI4.deallocate();
   TstArr3DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DI4 test: FAIL" << std::endl;

   // Test for 4DI4
   OMEGA::Array4DI4 TstArr4DI4("TstArr", NumTimeLvls, NumTracers, NumCells,
                               NumVertLvls);
   OMEGA::ArrayHost4DI4 RefArr4DI4("RefArr", NumTimeLvls, NumTracers, NumCells,
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

   parallel_for(
       Bounds<4>(NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int m, int k, int j, int i) {
          TstArr4DI4(m, k, j, i) = i + j + k + m;
       });

   yakl::fence();
   auto TstHost4DI4 = TstArr4DI4.createHostCopy();

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
   TstHost4DI4.deallocate();
   RefArr4DI4.deallocate();
   TstArr4DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DI4 test: FAIL" << std::endl;

   // Test for 5DI4
   OMEGA::Array5DI4 TstArr5DI4("TstArr", NumExtra, NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);
   OMEGA::ArrayHost5DI4 RefArr5DI4("RefArr", NumExtra, NumTimeLvls, NumTracers,
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

   parallel_for(
       Bounds<5>(NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int n, int m, int k, int j, int i) {
          TstArr5DI4(n, m, k, j, i) = i + j + k + m + n;
       });

   yakl::fence();
   auto TstHost5DI4 = TstArr5DI4.createHostCopy();

   icount = 0;
   for (int n = 0; n < NumExtra; ++n) {
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost5DI4(n, m, k, j, i) != RefArr5DI4(n, m, k, j, i))
                     ++icount;
               }
            }
         }
      }
   }
   TstHost5DI4.deallocate();
   RefArr5DI4.deallocate();
   TstArr5DI4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DI4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DI4 test: FAIL" << std::endl;

   // Test for 1DI8
   OMEGA::Array1DI8 TstArr1DI8("TstArr", NumCells);
   OMEGA::ArrayHost1DI8 RefArr1DI8("RefArr", NumCells);

   for (int i = 0; i < NumCells; ++i) {
      RefArr1DI8(i) = i;
   }

   parallel_for(
       Bounds<1>(NumCells), YAKL_LAMBDA(int i) { TstArr1DI8(i) = i; });

   yakl::fence();
   auto TstHost1DI8 = TstArr1DI8.createHostCopy();

   icount = 0;
   for (int i = 0; i < NumCells; ++i) {
      if (TstHost1DI8(i) != RefArr1DI8(i))
         ++icount;
   }
   TstHost1DI8.deallocate();
   RefArr1DI8.deallocate();
   TstArr1DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DI8 test: FAIL" << std::endl;

   // Test for 2DI8
   OMEGA::Array2DI8 TstArr2DI8("TstArr", NumCells, NumVertLvls);
   OMEGA::ArrayHost2DI8 RefArr2DI8("RefArr", NumCells, NumVertLvls);

   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         RefArr2DI8(j, i) = i + j;
      }
   }

   parallel_for(
       Bounds<2>(NumCells, NumVertLvls),
       YAKL_LAMBDA(int j, int i) { TstArr2DI8(j, i) = i + j; });

   yakl::fence();
   auto TstHost2DI8 = TstArr2DI8.createHostCopy();

   icount = 0;
   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         if (TstHost2DI8(j, i) != RefArr2DI8(j, i))
            ++icount;
      }
   }
   TstHost2DI8.deallocate();
   RefArr2DI8.deallocate();
   TstArr2DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DI8 test: FAIL" << std::endl;

   // Test for 3DI8
   OMEGA::Array3DI8 TstArr3DI8("TstArr", NumTracers, NumCells, NumVertLvls);
   OMEGA::ArrayHost3DI8 RefArr3DI8("RefArr", NumTracers, NumCells, NumVertLvls);

   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr3DI8(k, j, i) = i + j + k;
         }
      }
   }

   parallel_for(
       Bounds<3>(NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int k, int j, int i) { TstArr3DI8(k, j, i) = i + j + k; });

   yakl::fence();
   auto TstHost3DI8 = TstArr3DI8.createHostCopy();

   icount = 0;
   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost3DI8(k, j, i) != RefArr3DI8(k, j, i))
               ++icount;
         }
      }
   }
   TstHost3DI8.deallocate();
   RefArr3DI8.deallocate();
   TstArr3DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DI8 test: FAIL" << std::endl;

   // Test for 4DI8
   OMEGA::Array4DI8 TstArr4DI8("TstArr", NumTimeLvls, NumTracers, NumCells,
                               NumVertLvls);
   OMEGA::ArrayHost4DI8 RefArr4DI8("RefArr", NumTimeLvls, NumTracers, NumCells,
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

   parallel_for(
       Bounds<4>(NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int m, int k, int j, int i) {
          TstArr4DI8(m, k, j, i) = i + j + k + m;
       });

   yakl::fence();
   auto TstHost4DI8 = TstArr4DI8.createHostCopy();

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
   TstHost4DI8.deallocate();
   RefArr4DI8.deallocate();
   TstArr4DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DI8 test: FAIL" << std::endl;

   // Test for 5DI8
   OMEGA::Array5DI8 TstArr5DI8("TstArr", NumExtra, NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);
   OMEGA::ArrayHost5DI8 RefArr5DI8("RefArr", NumExtra, NumTimeLvls, NumTracers,
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

   parallel_for(
       Bounds<5>(NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int n, int m, int k, int j, int i) {
          TstArr5DI8(n, m, k, j, i) = i + j + k + m + n;
       });

   yakl::fence();
   auto TstHost5DI8 = TstArr5DI8.createHostCopy();

   icount = 0;
   for (int n = 0; n < NumExtra; ++n) {
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost5DI8(n, m, k, j, i) != RefArr5DI8(n, m, k, j, i))
                     ++icount;
               }
            }
         }
      }
   }
   TstHost5DI8.deallocate();
   RefArr5DI8.deallocate();
   TstArr5DI8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DI8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DI8 test: FAIL" << std::endl;

   // Test for 1DR4
   OMEGA::Array1DR4 TstArr1DR4("TstArr", NumCells);
   OMEGA::ArrayHost1DR4 RefArr1DR4("RefArr", NumCells);

   for (int i = 0; i < NumCells; ++i) {
      RefArr1DR4(i) = i;
   }

   parallel_for(
       Bounds<1>(NumCells), YAKL_LAMBDA(int i) { TstArr1DR4(i) = i; });

   yakl::fence();
   auto TstHost1DR4 = TstArr1DR4.createHostCopy();

   icount = 0;
   for (int i = 0; i < NumCells; ++i) {
      if (TstHost1DR4(i) != RefArr1DR4(i))
         ++icount;
   }
   TstHost1DR4.deallocate();
   RefArr1DR4.deallocate();
   TstArr1DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DR4 test: FAIL" << std::endl;

   // Test for 2DR4
   OMEGA::Array2DR4 TstArr2DR4("TstArr", NumCells, NumVertLvls);
   OMEGA::ArrayHost2DR4 RefArr2DR4("RefArr", NumCells, NumVertLvls);

   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         RefArr2DR4(j, i) = i + j;
      }
   }

   parallel_for(
       Bounds<2>(NumCells, NumVertLvls),
       YAKL_LAMBDA(int j, int i) { TstArr2DR4(j, i) = i + j; });

   yakl::fence();
   auto TstHost2DR4 = TstArr2DR4.createHostCopy();

   icount = 0;
   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         if (TstHost2DR4(j, i) != RefArr2DR4(j, i))
            ++icount;
      }
   }
   TstHost2DR4.deallocate();
   RefArr2DR4.deallocate();
   TstArr2DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DR4 test: FAIL" << std::endl;

   // Test for 3DR4
   OMEGA::Array3DR4 TstArr3DR4("TstArr", NumTracers, NumCells, NumVertLvls);
   OMEGA::ArrayHost3DR4 RefArr3DR4("RefArr", NumTracers, NumCells, NumVertLvls);

   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr3DR4(k, j, i) = i + j + k;
         }
      }
   }

   parallel_for(
       Bounds<3>(NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int k, int j, int i) { TstArr3DR4(k, j, i) = i + j + k; });

   yakl::fence();
   auto TstHost3DR4 = TstArr3DR4.createHostCopy();

   icount = 0;
   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost3DR4(k, j, i) != RefArr3DR4(k, j, i))
               ++icount;
         }
      }
   }
   TstHost3DR4.deallocate();
   RefArr3DR4.deallocate();
   TstArr3DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DR4 test: FAIL" << std::endl;

   // Test for 4DR4
   OMEGA::Array4DR4 TstArr4DR4("TstArr", NumTimeLvls, NumTracers, NumCells,
                               NumVertLvls);
   OMEGA::ArrayHost4DR4 RefArr4DR4("RefArr", NumTimeLvls, NumTracers, NumCells,
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

   parallel_for(
       Bounds<4>(NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int m, int k, int j, int i) {
          TstArr4DR4(m, k, j, i) = i + j + k + m;
       });

   yakl::fence();
   auto TstHost4DR4 = TstArr4DR4.createHostCopy();

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
   TstHost4DR4.deallocate();
   RefArr4DR4.deallocate();
   TstArr4DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DR4 test: FAIL" << std::endl;

   // Test for 5DR4
   OMEGA::Array5DR4 TstArr5DR4("TstArr", NumExtra, NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);
   OMEGA::ArrayHost5DR4 RefArr5DR4("RefArr", NumExtra, NumTimeLvls, NumTracers,
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

   parallel_for(
       Bounds<5>(NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int n, int m, int k, int j, int i) {
          TstArr5DR4(n, m, k, j, i) = i + j + k + m + n;
       });

   yakl::fence();
   auto TstHost5DR4 = TstArr5DR4.createHostCopy();

   icount = 0;
   for (int n = 0; n < NumExtra; ++n) {
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost5DR4(n, m, k, j, i) != RefArr5DR4(n, m, k, j, i))
                     ++icount;
               }
            }
         }
      }
   }
   TstHost5DR4.deallocate();
   RefArr5DR4.deallocate();
   TstArr5DR4.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DR4 test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DR4 test: FAIL" << std::endl;

   // Test for 1DR8
   OMEGA::Array1DR8 TstArr1DR8("TstArr", NumCells);
   OMEGA::ArrayHost1DR8 RefArr1DR8("RefArr", NumCells);

   for (int i = 0; i < NumCells; ++i) {
      RefArr1DR8(i) = i;
   }

   parallel_for(
       Bounds<1>(NumCells), YAKL_LAMBDA(int i) { TstArr1DR8(i) = i; });

   yakl::fence();
   auto TstHost1DR8 = TstArr1DR8.createHostCopy();

   icount = 0;
   for (int i = 0; i < NumCells; ++i) {
      if (TstHost1DR8(i) != RefArr1DR8(i))
         ++icount;
   }
   TstHost1DR8.deallocate();
   RefArr1DR8.deallocate();
   TstArr1DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DR8 test: FAIL" << std::endl;

   // Test for 2DR8
   OMEGA::Array2DR8 TstArr2DR8("TstArr", NumCells, NumVertLvls);
   OMEGA::ArrayHost2DR8 RefArr2DR8("RefArr", NumCells, NumVertLvls);

   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         RefArr2DR8(j, i) = i + j;
      }
   }

   parallel_for(
       Bounds<2>(NumCells, NumVertLvls),
       YAKL_LAMBDA(int j, int i) { TstArr2DR8(j, i) = i + j; });

   yakl::fence();
   auto TstHost2DR8 = TstArr2DR8.createHostCopy();

   icount = 0;
   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         if (TstHost2DR8(j, i) != RefArr2DR8(j, i))
            ++icount;
      }
   }
   TstHost2DR8.deallocate();
   RefArr2DR8.deallocate();
   TstArr2DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DR8 test: FAIL" << std::endl;

   // Test for 3DR8
   OMEGA::Array3DR8 TstArr3DR8("TstArr", NumTracers, NumCells, NumVertLvls);
   OMEGA::ArrayHost3DR8 RefArr3DR8("RefArr", NumTracers, NumCells, NumVertLvls);

   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr3DR8(k, j, i) = i + j + k;
         }
      }
   }

   parallel_for(
       Bounds<3>(NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int k, int j, int i) { TstArr3DR8(k, j, i) = i + j + k; });

   yakl::fence();
   auto TstHost3DR8 = TstArr3DR8.createHostCopy();

   icount = 0;
   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost3DR8(k, j, i) != RefArr3DR8(k, j, i))
               ++icount;
         }
      }
   }
   TstHost3DR8.deallocate();
   RefArr3DR8.deallocate();
   TstArr3DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DR8 test: FAIL" << std::endl;

   // Test for 4DR8
   OMEGA::Array4DR8 TstArr4DR8("TstArr", NumTimeLvls, NumTracers, NumCells,
                               NumVertLvls);
   OMEGA::ArrayHost4DR8 RefArr4DR8("RefArr", NumTimeLvls, NumTracers, NumCells,
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

   parallel_for(
       Bounds<4>(NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int m, int k, int j, int i) {
          TstArr4DR8(m, k, j, i) = i + j + k + m;
       });

   yakl::fence();
   auto TstHost4DR8 = TstArr4DR8.createHostCopy();

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
   TstHost4DR8.deallocate();
   RefArr4DR8.deallocate();
   TstArr4DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DR8 test: FAIL" << std::endl;

   // Test for 5DR8
   OMEGA::Array5DR8 TstArr5DR8("TstArr", NumExtra, NumTimeLvls, NumTracers,
                               NumCells, NumVertLvls);
   OMEGA::ArrayHost5DR8 RefArr5DR8("RefArr", NumExtra, NumTimeLvls, NumTracers,
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

   parallel_for(
       Bounds<5>(NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int n, int m, int k, int j, int i) {
          TstArr5DR8(n, m, k, j, i) = i + j + k + m + n;
       });

   yakl::fence();
   auto TstHost5DR8 = TstArr5DR8.createHostCopy();

   icount = 0;
   for (int n = 0; n < NumExtra; ++n) {
      for (int m = 0; m < NumTimeLvls; ++m) {
         for (int k = 0; k < NumTracers; ++k) {
            for (int j = 0; j < NumCells; ++j) {
               for (int i = 0; i < NumVertLvls; ++i) {
                  if (TstHost5DR8(n, m, k, j, i) != RefArr5DR8(n, m, k, j, i))
                     ++icount;
               }
            }
         }
      }
   }
   TstHost5DR8.deallocate();
   RefArr5DR8.deallocate();
   TstArr5DR8.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DR8 test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DR8 test: FAIL" << std::endl;

   // Test for 1DReal
   OMEGA::Array1DReal TstArr1DReal("TstArr", NumCells);
   OMEGA::ArrayHost1DReal RefArr1DReal("RefArr", NumCells);

   for (int i = 0; i < NumCells; ++i) {
      RefArr1DReal(i) = i;
   }

   parallel_for(
       Bounds<1>(NumCells), YAKL_LAMBDA(int i) { TstArr1DReal(i) = i; });

   yakl::fence();
   auto TstHost1DReal = TstArr1DReal.createHostCopy();

   icount = 0;
   for (int i = 0; i < NumCells; ++i) {
      if (TstHost1DReal(i) != RefArr1DReal(i))
         ++icount;
   }
   TstHost1DReal.deallocate();
   RefArr1DReal.deallocate();
   TstArr1DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 1DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 1DReal test: FAIL" << std::endl;

   // Test for 2DReal
   OMEGA::Array2DReal TstArr2DReal("TstArr", NumCells, NumVertLvls);
   OMEGA::ArrayHost2DReal RefArr2DReal("RefArr", NumCells, NumVertLvls);

   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         RefArr2DReal(j, i) = i + j;
      }
   }

   parallel_for(
       Bounds<2>(NumCells, NumVertLvls),
       YAKL_LAMBDA(int j, int i) { TstArr2DReal(j, i) = i + j; });

   yakl::fence();
   auto TstHost2DReal = TstArr2DReal.createHostCopy();

   icount = 0;
   for (int j = 0; j < NumCells; ++j) {
      for (int i = 0; i < NumVertLvls; ++i) {
         if (TstHost2DReal(j, i) != RefArr2DReal(j, i))
            ++icount;
      }
   }
   TstHost2DReal.deallocate();
   RefArr2DReal.deallocate();
   TstArr2DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 2DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 2DReal test: FAIL" << std::endl;

   // Test for 3DReal
   OMEGA::Array3DReal TstArr3DReal("TstArr", NumTracers, NumCells, NumVertLvls);
   OMEGA::ArrayHost3DReal RefArr3DReal("RefArr", NumTracers, NumCells,
                                       NumVertLvls);

   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            RefArr3DReal(k, j, i) = i + j + k;
         }
      }
   }

   parallel_for(
       Bounds<3>(NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int k, int j, int i) { TstArr3DReal(k, j, i) = i + j + k; });

   yakl::fence();
   auto TstHost3DReal = TstArr3DReal.createHostCopy();

   icount = 0;
   for (int k = 0; k < NumTracers; ++k) {
      for (int j = 0; j < NumCells; ++j) {
         for (int i = 0; i < NumVertLvls; ++i) {
            if (TstHost3DReal(k, j, i) != RefArr3DReal(k, j, i))
               ++icount;
         }
      }
   }
   TstHost3DReal.deallocate();
   RefArr3DReal.deallocate();
   TstArr3DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 3DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 3DReal test: FAIL" << std::endl;

   // Test for 4DReal
   OMEGA::Array4DReal TstArr4DReal("TstArr", NumTimeLvls, NumTracers, NumCells,
                                   NumVertLvls);
   OMEGA::ArrayHost4DReal RefArr4DReal("RefArr", NumTimeLvls, NumTracers,
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

   parallel_for(
       Bounds<4>(NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int m, int k, int j, int i) {
          TstArr4DReal(m, k, j, i) = i + j + k + m;
       });

   yakl::fence();
   auto TstHost4DReal = TstArr4DReal.createHostCopy();

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
   TstHost4DReal.deallocate();
   RefArr4DReal.deallocate();
   TstArr4DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 4DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 4DReal test: FAIL" << std::endl;

   // Test for 5DReal
   OMEGA::Array5DReal TstArr5DReal("TstArr", NumExtra, NumTimeLvls, NumTracers,
                                   NumCells, NumVertLvls);
   OMEGA::ArrayHost5DReal RefArr5DReal("RefArr", NumExtra, NumTimeLvls,
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

   parallel_for(
       Bounds<5>(NumExtra, NumTimeLvls, NumTracers, NumCells, NumVertLvls),
       YAKL_LAMBDA(int n, int m, int k, int j, int i) {
          TstArr5DReal(n, m, k, j, i) = i + j + k + m + n;
       });

   yakl::fence();
   auto TstHost5DReal = TstArr5DReal.createHostCopy();

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
   TstHost5DReal.deallocate();
   RefArr5DReal.deallocate();
   TstArr5DReal.deallocate();

   if (icount == 0)
      std::cout << "YAKL 5DReal test: PASS" << std::endl;
   else
      std::cout << "YAKL 5DReal test: FAIL" << std::endl;

   // finalize environments
   // MPI_Status status;
   yakl::finalize();
   MPI_Finalize();

} // end of main
//===-----------------------------------------------------------------------===/
