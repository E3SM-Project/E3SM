//===-- Test driver for OMEGA Reductions -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Reductions
///
//
//===-----------------------------------------------------------------------===/

#include <string>

#include <mpi.h>

#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Reductions.h"

using namespace OMEGA;

int main(int argc, char *argv[]) {

   int RetVal = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {

      // Create reference values based on MPI_COMM_WORLD
      MPI_Comm Comm;
      int MyTask, MySize, err;
      MachEnv::init(MPI_COMM_WORLD);
      MachEnv *DefEnv = MachEnv::getDefault();
      Comm            = DefEnv->getComm();
      MyTask          = DefEnv->getMyTask();
      MySize          = DefEnv->getNumTasks();

      // Initialize the Logging system
      OMEGA::initLogging(DefEnv);

      I4 MyInt4 = 1, MyResI4 = 0;
      I8 MyInt8 = 2, MyResI8 = 0;
      R4 MyR4 = 3.000001, MyResR4 = 0.0;
      R8 MyR8 = 4.0000000000001, MyResR8 = 0.0;
      Real MyReal = 5.000001, MyResReal = 0.0;

      // test SUM for scalars
      err             = globalSum(&MyInt4, Comm, &MyResI4);
      I4 expI4        = MyInt4 * MySize;
      char const *res = "FAIL";
      if (err == 0 && MyResI4 == expI4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum I4:    %s (exp,act=%d,%d)\n", res, expI4, MyResI4);

      err      = globalSum(&MyInt8, Comm, &MyResI8);
      I8 expI8 = MyInt8 * MySize;
      res      = "FAIL";
      if (err == 0 && MyResI8 == expI8)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum I8:    %s (exp,act=%ld,%ld)\n", res, expI8, MyResI8);

      err      = globalSum(&MyR4, Comm, &MyResR4);
      R4 expR4 = MyR4 * MySize;
      res      = "FAIL";
      if (err == 0 && MyResR4 == expR4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum R4:    %s (exp,act=%f,%f)\n", res, expR4, MyResR4);

      err      = globalSum(&MyR8, Comm, &MyResR8);
      R8 expR8 = MyR8 * MySize;
      res      = "FAIL";
      if (err == 0 && MyResR8 == expR8)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum R8:    %s (exp,act=%.15lf,%.15lf)\n", res, expR8,
             MyResR8);

      err          = globalSum(&MyReal, Comm, &MyResReal);
      Real expReal = MyReal * MySize;
      res          = "FAIL";
      if (err == 0 && MyResReal == expReal)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum real:  %s (exp,act=%.15lf,%.15lf)\n", res, expReal,
             MyResReal);

      // test SUM of I4 arrays
      int i, j, k;
      I4 NumCells = 10, NumVertLvls = 10, c = 0;
      HostArray1DI4 HostArr1DI4("HostArrD1", NumCells);
      HostArray2DI4 HostArr2DI4("HostArrD2", NumCells, NumVertLvls);
      I4 Sum1DI4 = 0, Sum2DI4 = 0;
      for (i = 0; i < NumCells; i++) {
         HostArr1DI4(i) = i;
         Sum1DI4 += i;
         for (j = 0; j < NumVertLvls; j++) {
            HostArr2DI4(i, j) = c;
            Sum2DI4 += c;
            c++;
         }
      }
      err   = globalSum(HostArr1DI4, Comm, &MyResI4);
      expI4 = Sum1DI4 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResI4 == expI4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum A1DI4: %s (exp,act=%d,%d)\n", res, expI4, MyResI4);

      err   = globalSum(HostArr2DI4, Comm, &MyResI4);
      expI4 = Sum2DI4 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResI4 == expI4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum A2DI4: %s (exp,act=%d,%d)\n", res, expI4, MyResI4);

      // test SUM of I8 arrays
      HostArray1DI8 HostArr1DI8("HostArrD1I8", NumCells);
      HostArray2DI8 HostArr2DI8("HostArrD2I8", NumCells, NumVertLvls);
      I8 Sum1DI8 = 0, Sum2DI8 = 0, c8 = 0;
      for (i = 0; i < NumCells; i++) {
         HostArr1DI8(i) = i;
         Sum1DI8 += i;
         for (j = 0; j < NumVertLvls; j++) {
            HostArr2DI8(i, j) = c8;
            Sum2DI8 += c8;
            c8++;
         }
      }
      err   = globalSum(HostArr1DI8, Comm, &MyResI8);
      expI8 = Sum1DI8 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResI8 == expI8)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum A1DI8: %s (exp,act=%ld,%ld)\n", res, expI8, MyResI8);

      err   = globalSum(HostArr2DI8, Comm, &MyResI8);
      expI8 = Sum2DI8 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResI8 == expI8)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum A2DI8: %s (exp,act=%ld,%ld)\n", res, expI8, MyResI8);

      // test SUM of R4 arrays
      HostArray1DR4 HostArr1DR4("HostArrD1R4", NumCells);
      HostArray2DR4 HostArr2DR4("HostArrD2R4", NumCells, NumVertLvls);
      R4 Sum1DR4 = 0.0, Sum2DR4 = 0.0, f = 0.0;
      for (i = 0; i < NumCells; i++) {
         HostArr1DR4(i) = i + 0.00001;
         Sum1DR4 += HostArr1DR4(i);
         for (j = 0; j < NumVertLvls; j++) {
            HostArr2DR4(i, j) = f;
            Sum2DR4 += HostArr2DR4(i, j);
            f += 1.00001;
         }
      }
      err   = globalSum(HostArr1DR4, Comm, &MyResR4);
      expR4 = Sum1DR4 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResR4 == expR4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum A1DR4: %s (exp,act=%.10f,%.10f)\n", res, expR4,
             MyResR4);

      err   = globalSum(HostArr2DR4, Comm, &MyResR4);
      expR4 = Sum2DR4 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResR4 == expR4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum A2DR4: %s (exp,act=%f,%f)\n", res, expR4, MyResR4);

      // test SUM of R8 arrays
      HostArray1DR8 HostArr1DR8("HostArrD1R8", NumCells);
      HostArray2DR8 HostArr2DR8("HostArrD2R8", NumCells, NumVertLvls);
      R8 Sum1DR8 = 0.0, Sum2DR8 = 0.0, d = 0.0;
      complex<double> LocalSum1D(0.0, 0.0), LocalSum2D(0.0, 0.0);
      double e, t1, t2, ai;
      for (i = 0; i < NumCells; i++) {
         HostArr1DR8(i) = i + 0.0000000000001;
         // local ddsum
         ai = HostArr1DR8(i);
         t1 = ai + real(LocalSum1D);
         e  = t1 - ai;
         t2 = ((real(LocalSum1D) - e) + (ai - (t1 - e))) + imag(LocalSum1D);
         LocalSum1D = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
         for (j = 0; j < NumVertLvls; j++) {
            HostArr2DR8(i, j) = d;
            d += 1.0000000000001;
            ai = HostArr2DR8(i, j);
            t1 = ai + real(LocalSum2D);
            e  = t1 - ai;
            t2 = ((real(LocalSum2D) - e) + (ai - (t1 - e))) + imag(LocalSum2D);
            LocalSum2D = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
         }
      }
      Sum1DR8 = real(LocalSum1D);
      Sum2DR8 = real(LocalSum2D);
      MyResR8 = 0.0;
      err     = globalSum(HostArr1DR8, Comm, &MyResR8);
      // perform serial sum across all MPI tasks
      complex<double> SerialSum(0.0, 0.0);
      for (i = 0; i < MySize; i++) {
         // ddsum across tasks
         t1 = Sum1DR8 + real(SerialSum);
         e  = t1 - Sum1DR8;
         t2 = ((real(SerialSum) - e) + (Sum1DR8 - (t1 - e))) + imag(SerialSum);
         SerialSum = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
      }
      expR8 = real(SerialSum);
      res   = "FAIL";
      if (err == 0 && MyResR8 == expR8)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum A1DR8: %s (exp,act=%.13lf,%.13lf)\n", res, expR8,
             MyResR8);

      err   = globalSum(HostArr2DR8, Comm, &MyResR8);
      expR8 = Sum2DR8 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResR8 == expR8)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum A2DR8: %s (exp,act=%.13lf,%.13lf)\n", res, expR8,
             MyResR8);

      //==========================================================================
      // test MIN, MAX of scalars
      MyInt4 = MyTask;
      err    = globalMin(&MyInt4, &MyResI4, Comm);
      res    = "FAIL";
      if (err == 0 && MyResI4 == 0)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global min I4:    %s (exp,act=0,%d)\n", res, MyResI4);

      MyInt8 = MyTask;
      err    = globalMax(&MyInt8, &MyResI8, Comm);
      res    = "FAIL";
      if (err == 0 && MyResI8 == MySize - 1)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global max I4:    %s (exp,act=%d,%ld)\n", res, MySize - 1,
             MyResI8);

      R8 MyR8Tmp = MyTask + MyR8;
      err        = globalMin(&MyR8Tmp, &MyResR8, Comm);
      res        = "FAIL";
      if (err == 0 && MyResR8 == MyR8)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global min R8:    %s (exp,act=%.13lf,%.13lf)\n", res, MyR8,
             MyResR8);

      err = globalMax(&MyR8Tmp, &MyResR8, Comm);
      res = "FAIL";
      if (err == 0 && MyResR8 == (MySize - 1 + MyR8))
         res = "PASS";
      else
         RetVal += 1;
      printf("Global max R8:    %s (exp,act=%.13lf,%.13lf)\n", res,
             MySize - 1 + MyR8, MyResR8);

      //==========================================================================
      // test MIN, MAX of arrays
      HostArray1DI4 HostA1DI4Work("HostA1DI4Work", NumCells * MySize);
      HostArray1DI4 HostA1DI4Min("HostA1DI4Min", NumCells * MySize);
      for (j = 0; j < NumCells; j++) {
         for (i = 0; i < MySize; i++) {
            k = i * NumCells + j;
            if (MyTask == i) {
               HostA1DI4Work(k) = (i + 1) * k; // processor-specific work array
            } else {
               HostA1DI4Work(k) = k; // default
            }
         }
      }
      err = globalMin(HostA1DI4Work, HostA1DI4Min, Comm);
      res = "PASS";
      for (i = 0; i < NumCells * MySize; i++) {
         if (HostA1DI4Min(i) != i) {
            res = "FAIL";
            RetVal += 1;
         }
      }
      printf("Global min A1DI4: %s\n", res);

      HostArray1DI4 HostA1DI4Max("HostA1DI4Max", NumCells * MySize);
      err = globalMax(HostA1DI4Work, HostA1DI4Max, Comm);
      res = "PASS";
      for (i = 0; i < MySize; i++) {
         for (j = 0; j < NumCells; j++) {
            k = i * NumCells + j;
            if (HostA1DI4Max(k) != (i + 1) * k) {
               res = "FAIL";
               RetVal += 1;
            }
         }
      }
      printf("Global max A1DI4: %s\n", res);

      // test SUM of I4 arrays on device
      Array1DI4 DevArr1DI4("DevArr1DI4", NumCells);
      Array2DI4 DevArr2DI4("DevArr2DI4", NumCells, NumVertLvls);

      parallelFor({NumCells}, KOKKOS_LAMBDA(int i) { DevArr1DI4(i) = i; });
      parallelFor(
          {NumCells, NumVertLvls},
          KOKKOS_LAMBDA(int i, int j) { DevArr2DI4(i, j) = i * 10 + j; });
      Kokkos::fence();

      err   = globalSum(DevArr1DI4, Comm, &MyResI4);
      expI4 = Sum1DI4 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResI4 == expI4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum device A1DI4: %s (exp,act=%d,%d)\n", res, expI4,
             MyResI4);

      err   = globalSum(DevArr2DI4, Comm, &MyResI4);
      expI4 = Sum2DI4 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResI4 == expI4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum device A2DI4: %s (exp,act=%d,%d)\n", res, expI4,
             MyResI4);

      // test SUM of R4 arrays on device
      Array1DR4 DevArr1DR4("DevArrD1R4", NumCells);
      parallelFor(
          {NumCells}, KOKKOS_LAMBDA(int i) { DevArr1DR4(i) = i + 0.00001; });
      Kokkos::fence();

      err   = globalSum(DevArr1DR4, Comm, &MyResR4);
      expR4 = Sum1DR4 * MySize;
      res   = "FAIL";
      if (err == 0 && MyResR4 == expR4)
         res = "PASS";
      else
         RetVal += 1;
      printf("Global sum device A1DR4: %s (exp,act=%.10f,%.10f)\n", res, expR4,
             MyResR4);
   }
   Kokkos::finalize();
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
