//===-- Test driver for OMEGA Reductions -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Reductions
///
/// This driver tests the OMEGA model Reductions module that sets up various
///
//
//===-----------------------------------------------------------------------===/

#include <string>

#include <mpi.h>

#include "MachEnv.h"
#include "Reductions.h"

//------------------------------------------------------------------------------
// The test driver for Reductions. This tests
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();

   // Create reference values based on MPI_COMM_WORLD
   MPI_Comm Comm;
   int MyTask, MySize, err;
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   Comm                   = DefEnv->getComm();
   MyTask                 = DefEnv->getMyTask();
   MySize                 = DefEnv->getNumTasks();

   OMEGA::I4 MyInt4 = 1, MyResI4 = 0;
   OMEGA::I8 MyInt8 = 2, MyResI8 = 0;
   OMEGA::R4 MyR4 = 3.000001, MyResR4 = 0.0;
   OMEGA::R8 MyR8 = 4.0000000000001, MyResR8 = 0.0, MyResReal = 0.0;
   OMEGA::Real MyReal = 5.000001;
   char *res          = "FAIL";

   // test SUM for scalars
   err = OMEGA::globalSum(&MyInt4, Comm, &MyResI4);
   if (err == 0 && MyResI4 == MyInt4 * MySize)
      res = "PASS";
   printf("Global sum I4:    %s (exp,act=%d,%d)\n", res, MyInt4 * MySize,
          MyResI4);

   err = OMEGA::globalSum(&MyInt8, Comm, &MyResI8);
   res = "FAIL";
   if (MyResI8 == MyInt8 * MySize)
      res = "PASS";
   printf("Global sum I8:    %s (exp,act=%lld,%lld)\n", res, MyInt8 * MySize,
          MyResI8);

   err = OMEGA::globalSum(&MyR4, Comm, &MyResR4);
   res = "FAIL";
   if (MyResR4 == MyR4 * MySize)
      res = "PASS";
   printf("Global sum R4:    %s (exp,act=%f,%f)\n", res, MyR4 * MySize,
          MyResR4);

   err = OMEGA::globalSum(&MyR8, Comm, &MyResR8);
   res = "FAIL";
   if (MyResR8 == MyR8 * MySize)
      res = "PASS";
   printf("Global sum R8:    %s (exp,act=%.15lf,%.15lf)\n", res, MyR8 * MySize,
          MyResR8);

   err = OMEGA::globalSum(&MyReal, Comm, &MyResReal);
   res = "FAIL";
   if (MyResReal == MyReal * MySize)
      res = "PASS";
   printf("Global sum real:  %s (exp,act=%.15lf,%.15lf)\n", res,
          MyReal * MySize, MyResReal);

   // test SUM of I4 arrays
   int i, j, k;
   OMEGA::I4 NumCells = 10, NumVertLvls = 10, c = 0;
   OMEGA::HostArray1DI4 HostArr1DI4("HostArrD1", NumCells);
   OMEGA::HostArray2DI4 HostArr2DI4("HostArrD2", NumCells, NumVertLvls);
   OMEGA::I8 Sum1D = 0, Sum2D = 0;
   for (i = 0; i < NumCells; i++) {
      HostArr1DI4(i) = i;
      Sum1D += i;
      for (j = 0; j < NumVertLvls; j++) {
         HostArr2DI4(i, j) = c;
         Sum2D += c;
         c++;
      }
   }
   err = OMEGA::globalSum(HostArr1DI4, Comm, &MyResI4);
   res = "FAIL";
   if (MyResI4 == Sum1D * MySize)
      res = "PASS";
   printf("Global sum A1DI4: %s (exp,act=%d,%d)\n", res, Sum1D * MySize,
          MyResI4);

   err = OMEGA::globalSum(HostArr2DI4, Comm, &MyResI4);
   res = "FAIL";
   if (MyResI4 == Sum2D * MySize)
      res = "PASS";
   printf("Global sum A2DI4: %s (exp,act=%d,%d)\n", res, Sum2D * MySize,
          MyResI4);
   // HostArr1DI4.deallocate();
   // HostArr2DI4.deallocate();

   // test SUM of I8 arrays
   OMEGA::HostArray1DI8 HostArr1DI8("HostArrD1I8", NumCells);
   OMEGA::HostArray2DI8 HostArr2DI8("HostArrD2I8", NumCells, NumVertLvls);
   Sum1D = 0, Sum2D = 0, c = 0;
   for (i = 0; i < NumCells; i++) {
      HostArr1DI8(i) = i;
      Sum1D += i;
      for (j = 0; j < NumVertLvls; j++) {
         HostArr2DI8(i, j) = c;
         Sum2D += c;
         c++;
      }
   }
   err = OMEGA::globalSum(HostArr1DI8, Comm, &MyResI8);
   res = "FAIL";
   if (MyResI8 == Sum1D * MySize)
      res = "PASS";
   printf("Global sum A1DI8: %s (exp,act=%d,%d)\n", res, Sum1D * MySize,
          MyResI8);

   err = OMEGA::globalSum(HostArr2DI8, Comm, &MyResI8);
   res = "FAIL";
   if (MyResI8 == Sum2D * MySize)
      res = "PASS";
   printf("Global sum A2DI8: %s (exp,act=%d,%d)\n", res, Sum2D * MySize,
          MyResI8);
   // HostArr1DI8.deallocate();
   // HostArr2DI8.deallocate();

   // test SUM of R4 arrays
   OMEGA::HostArray1DR4 HostArr1DR4("HostArrD1R4", NumCells);
   OMEGA::HostArray2DR4 HostArr2DR4("HostArrD2R4", NumCells, NumVertLvls);
   OMEGA::R4 Sum1DR4 = 0.0, Sum2DR4 = 0.0, f = 0.0;
   for (i = 0; i < NumCells; i++) {
      HostArr1DR4(i) = i + 0.00001;
      Sum1DR4 += HostArr1DR4(i);
      for (j = 0; j < NumVertLvls; j++) {
         HostArr2DR4(i, j) = f;
         Sum2DR4 += f;
         f += 1.00001;
      }
   }
   err = OMEGA::globalSum(HostArr1DR4, Comm, &MyResR4);
   res = "FAIL";
   if (MyResR4 == Sum1DR4 * MySize)
      res = "PASS";
   printf("Global sum A1DR4: %s (exp,act=%.10f,%.10f)\n", res, Sum1DR4 * MySize,
          MyResR4);

   err = OMEGA::globalSum(HostArr2DR4, Comm, &MyResR4);
   res = "FAIL";
   if (MyResR4 == Sum2DR4 * MySize)
      res = "PASS";
   printf("Global sum A2DR4: %s (exp,act=%f,%f)\n", res, Sum2DR4 * MySize,
          MyResR4);
   // HostArr1DR4.deallocate();
   // HostArr2DR4.deallocate();

   // test SUM of R8 arrays
   OMEGA::HostArray1DR8 HostArr1DR8("HostArrD1R8", NumCells);
   OMEGA::HostArray2DR8 HostArr2DR8("HostArrD2R8", NumCells, NumVertLvls);
   OMEGA::R8 Sum1DR8 = 0.0, Sum2DR8 = 0.0, d = 0.0;
   double _Complex LocalSum1D = CMPLX(0.0, 0.0);
   double _Complex LocalSum2D = CMPLX(0.0, 0.0);
   double e, t1, t2, ai;
   for (i = 0; i < NumCells; i++) {
      HostArr1DR8(i) = i + 0.0000000000001;
      // local ddsum
      ai = HostArr1DR8(i);
      t1 = ai + creal(LocalSum1D);
      e  = t1 - ai;
      t2 = ((creal(LocalSum1D) - e) + (ai - (t1 - e))) + cimag(LocalSum1D);
      LocalSum1D = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
      for (j = 0; j < NumVertLvls; j++) {
         HostArr2DR8(i, j) = d;
         d += 1.0000000000001;
         ai = HostArr2DR8(i, j);
         t1 = ai + creal(LocalSum2D);
         e  = t1 - ai;
         t2 = ((creal(LocalSum2D) - e) + (ai - (t1 - e))) + cimag(LocalSum2D);
         LocalSum2D = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
      }
   }
   Sum1DR8   = creal(LocalSum1D);
   Sum2DR8   = creal(LocalSum2D);
   MyResReal = 0.0;
   err       = OMEGA::globalSum(HostArr1DR8, Comm, &MyResReal);
   if (err != 0) {
      printf("Global sum A1DR8: globalSum call FAIL");
   }
   // perform serial sum across all MPI tasks
   double _Complex SerialSum = CMPLX(0.0, 0.0);
   for (i = 0; i < MySize; i++) {
      // ddsum across tasks
      t1 = Sum1DR8 + creal(SerialSum);
      e  = t1 - Sum1DR8;
      t2 = ((creal(SerialSum) - e) + (Sum1DR8 - (t1 - e))) + cimag(SerialSum);
      SerialSum = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
   }
   res              = "FAIL";
   OMEGA::R8 expSum = creal(SerialSum);
   if (MyResReal == expSum)
      res = "PASS";
   printf("Global sum A1DR8: %s (exp,act=%.13lf,%.13lf)\n", res, expSum,
          MyResReal);

   err = OMEGA::globalSum(HostArr2DR8, Comm, &MyResReal);
   if (err != 0) {
      printf("Global sum A2DR8: globalSum call FAIL");
   }
   res = "FAIL";
   if (MyResReal == Sum2DR8 * MySize)
      res = "PASS";
   printf("Global sum A2DR8: %s (exp,act=%.13lf,%.13lf)\n", res,
          Sum2DR8 * MySize, MyResReal);
   // HostArr1DR8.deallocate();
   // HostArr2DR8.deallocate();

   //==========================================================================
   // test MIN, MAX of scalars
   MyInt4 = MyTask;
   err    = OMEGA::globalMin(&MyInt4, &MyResI4, Comm);
   res    = "FAIL";
   if (MyResI4 == 0)
      res = "PASS";
   printf("Global min I4:    %s (exp,act=0,%d)\n", res, MyResI4);

   MyInt8 = MyTask;
   err    = OMEGA::globalMax(&MyInt8, &MyResI8, Comm);
   res    = "FAIL";
   if (MyResI8 == MySize - 1)
      res = "PASS";
   printf("Global max I4:    %s (exp,act=%d,%d)\n", res, MySize - 1, MyResI8);

   OMEGA::R8 MyR8Tmp = MyTask + MyR8;
   err               = OMEGA::globalMin(&MyR8Tmp, &MyResR8, Comm);
   res               = "FAIL";
   if (MyResR8 == MyR8)
      res = "PASS";
   printf("Global min R8:    %s (exp,act=%.13lf,%.13lf)\n", res, MyR8, MyResR8);

   err = OMEGA::globalMax(&MyR8Tmp, &MyResR8, Comm);
   res = "FAIL";
   if (MyResR8 == (MySize - 1 + MyR8))
      res = "PASS";
   printf("Global max R8:    %s (exp,act=%.13lf,%.13lf)\n", res,
          MySize - 1 + MyR8, MyResR8);

   //==========================================================================
   // test MIN, MAX of arrays
   OMEGA::HostArray1DI4 HostA1DI4Work("HostA1DI4Work", NumCells * MySize);
   OMEGA::HostArray1DI4 HostA1DI4Min("HostA1DI4Min", NumCells * MySize);
   for (i = 0; i < MySize; i++) {
      for (j = 0; j < NumCells; j++) {
         k = i * NumCells + j;
         if (MyTask == i) {
            HostA1DI4Work(k) = (i + 1) * k; // processor-specific work array
         } else {
            HostA1DI4Work(k) = k; // default
         }
      }
   }
   err = OMEGA::globalMin(HostA1DI4Work, HostA1DI4Min, Comm);
   res = "PASS";
   for (i = 0; i < NumCells * MySize; i++) {
      // printf("ReductionsTest::HostA1DI4Min(%2d)=%2d\n",i,HostA1DI4Min(i));
      if (HostA1DI4Min(i) != i)
         res = "FAIL";
   }
   printf("Global min A1DI4: %s\n", res);
   // HostA1DI4Min.deallocate();

   OMEGA::HostArray1DI4 HostA1DI4Max("HostA1DI4Max", NumCells * MySize);
   err = OMEGA::globalMax(HostA1DI4Work, HostA1DI4Max, Comm);
   res = "PASS";
   for (i = 0; i < MySize; i++) {
      for (j = 0; j < NumCells; j++) {
         k = i * NumCells + j;
         // printf("ReductionsTest::HostA1DI4Max(%2d)=%2d\n",k,HostA1DI4Max(k));
         if (HostA1DI4Max(k) != (i + 1) * k)
            res = "FAIL";
      }
   }
   printf("Global max A1DI4: %s\n", res);
   // HostA1DI4Max.deallocate();
   // HostA1DI4Work.deallocate();

   // Kokkos::finalize();
   MPI_Finalize();
} // end of main
//===-----------------------------------------------------------------------===/
