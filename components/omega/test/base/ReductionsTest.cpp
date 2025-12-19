//===-- Test driver for OMEGA Reductions -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Reductions
///
//
//===-----------------------------------------------------------------------===/

#include <algorithm>
#include <cmath>
#include <limits>
#include <string>

#include <mpi.h>

#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "Reductions.h"

using namespace OMEGA;

//------------------------------------------------------------------------------
// utility function to compute one iteration of the DD algorithm for
// reproducible double precision sums
void sumDDTest(complex<double> &ddb, double &dda) {
   double t1 = dda + real(ddb);
   double e  = t1 - dda;
   double t2 = ((real(ddb) - e) + (dda - (t1 - e))) + imag(ddb);
   ddb       = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
}

//------------------------------------------------------------------------------
// Scalar sum test function
void testScalarSums() {

   // Get some MPI values based on default environment
   MachEnv *DefEnv = MachEnv::getDefault();
   MPI_Comm Comm   = DefEnv->getComm();
   int MyTask      = DefEnv->getMyTask();
   int NTasks      = DefEnv->getNumTasks();

   // Get similar values based on subset environment for reproducibility tests
   MachEnv *SubEnv  = MachEnv::get("Subset");
   MPI_Comm CommSub = SubEnv->getComm();
   bool IsSubMember = SubEnv->isMember();

   // For reproducibility tests, define a full range of values within
   // each data type.  We want the smallest (abs) value > 0
   // Further restrict the max values so that we do not inadvertantly
   // exceed limits.
   I4 MaxEntries = 8;
   I4 MaxI4      = std::numeric_limits<I4>::max() / (10 * MaxEntries);
   I8 MaxI8      = std::numeric_limits<I8>::max() / (10 * MaxEntries);
   R4 MaxR4      = std::numeric_limits<R4>::max() / (10.0 * MaxEntries);
   R8 MaxR8      = std::numeric_limits<R8>::max() / (10.0 * MaxEntries);
   R4 EpsR4      = 0.0001;
   R8 EpsR8      = 0.0000000000001;
   // Compute max safe exponents based on adjusted max above
   I4 ExpI4 = std::log2(MaxI4);
   I4 ExpI8 = std::log2(MaxI8);
   I4 ExpR4 = std::log10(MaxR4);
   I4 ExpR8 = std::log10(MaxR8);
   // To cover a large range of values, compute a factor for test vals
   // on each rank
   std::vector<I4> FacI4(NTasks);
   std::vector<I8> FacI8(NTasks);
   std::vector<R4> FacR4(NTasks);
   std::vector<R8> FacR8(NTasks);
   for (int I = 0; I < NTasks; ++I) {
      FacI4[I] = pow(2, std::min(I, ExpI4));
      FacI8[I] = pow(2, std::min(I, ExpI8));
      FacR4[I] = pow(10.0, std::min(I, ExpR4));
      FacR8[I] = pow(10.0, std::min(I, ExpR8));
   }

   // Initialize test and reference values
   I4 TstI4 = MyTask * FacI4[MyTask];
   I4 SumI4 = 0;
   I4 RefI4 = 0;
   I4 SubI4 = 0;
   I8 TstI8 = MyTask * FacI8[MyTask];
   I8 SumI8 = 0;
   I8 RefI8 = 0;
   I8 SubI8 = 0;
   R4 TstR4 = (MyTask + EpsR4) * FacR4[MyTask];
   R4 SumR4 = 0.0;
   R4 RefR4 = 0.0;
   R4 SubR4 = 0.0;
   R8 TstR8 = (MyTask + EpsR8) * FacR8[MyTask];
   R8 SumR8 = 0.0;
   R8 RefR8 = 0.0;
   R8 SubR8 = 0.0;
   R8 TmpR1 = 0.0; // for reproducible R4 sums
   R8 TmpR2 = 0.0; // for reproducible R4 sums
   // For reference values, compute a serial sum of all test values on default
   // and subset domains. Use special reproducible sum function for R8.
   double DDValRef;
   double DDValSub;
   complex<double> DDSumRef(0.0, 0.0);
   complex<double> DDSumSub(0.0, 0.0);
   for (int Task = 0; Task < NTasks; ++Task) {
      RefI4 += Task * FacI4[Task];
      RefI8 += Task * FacI8[Task];
      TmpR1 += (Task + EpsR4) * FacR4[Task];
      DDValRef = (Task + EpsR8) * FacR8[Task];
      sumDDTest(DDSumRef, DDValRef);
      if (Task < 4) {
         SubI4 += Task * FacI4[Task];
         SubI8 += Task * FacI8[Task];
         TmpR2 += (Task + EpsR4) * FacR4[Task];
         DDValSub = (Task + EpsR8) * FacR8[Task];
         sumDDTest(DDSumSub, DDValSub);
      }
   }
   RefR4 = TmpR1;
   RefR8 = real(DDSumRef);
   if (MyTask < 4) {
      SubR4 = TmpR2;
      SubR8 = real(DDSumSub);
   }

   // Scalar sum sanity checks
   // Test global scalar sums against ref values
   SumI4 = globalSum(TstI4, Comm);
   SumI8 = globalSum(TstI8, Comm);
   SumR4 = globalSum(TstR4, Comm);
   SumR8 = globalSum(TstR8, Comm);

   if (SumI4 != RefI4)
      ABORT_ERROR("ReductionsTest: FAIL globalSum (I4 scalar)"
                  "Expected = {} Actual = {}",
                  RefI4, SumI4);
   if (SumI8 != RefI8)
      ABORT_ERROR("ReductionsTest: FAIL globalSum (I8 scalar)"
                  "Expected = {} Actual = {}",
                  RefI8, SumI8);
   if (SumR4 != RefR4)
      ABORT_ERROR("ReductionsTest: FAIL globalSum (R4 scalar)"
                  "Expected = {} Actual = {}",
                  RefR4, SumR4);
   if (SumR8 != RefR8)
      ABORT_ERROR("ReductionsTest: FAIL globalSum (R8 scalar)"
                  "Expected = {} Actual = {}",
                  RefR8, SumR8);

   // Test scalar sums with subset communicator
   // Also tests reproducibility using different processor count
   if (IsSubMember) {
      SumI4 = globalSum(TstI4, CommSub);
      SumI8 = globalSum(TstI8, CommSub);
      SumR4 = globalSum(TstR4, CommSub);
      SumR8 = globalSum(TstR8, CommSub);

      if (SumI4 != SubI4)
         ABORT_ERROR("ReductionsTest: FAIL globalSum (I4 scalar) subset"
                     "Expected = {} Actual = {}",
                     SubI4, SumI4);
      if (SumI8 != SubI8)
         ABORT_ERROR("ReductionsTest: FAIL globalSum (I8 scalar) subset"
                     "Expected = {} Actual = {}",
                     SubI8, SumI8);
      if (SumR4 != SubR4)
         ABORT_ERROR("ReductionsTest: FAIL globalSum (R4 scalar) subset"
                     "Expected = {} Actual = {}",
                     SubR4, SumR4);
      if (SumR8 != SubR8)
         ABORT_ERROR("ReductionsTest: FAIL globalSum (R8 scalar) subset"
                     "Expected = {} Actual = {}",
                     SubR8, SumR8);
   } // end if subset member

   // Scalar sum multi-field checks - checks both value and reproducibility
   // For reproducibility, we set the second field to have values in the
   // reverse order of the original
   int NFields = 2;

   std::vector<I4> LocVecI4(NFields);
   std::vector<I8> LocVecI8(NFields);
   std::vector<R4> LocVecR4(NFields);
   std::vector<R8> LocVecR8(NFields);
   std::vector<I4> SumVecI4(NFields);
   std::vector<I8> SumVecI8(NFields);
   std::vector<R4> SumVecR4(NFields);
   std::vector<R8> SumVecR8(NFields);
   I4 TaskRev  = NTasks - 1 - MyTask;
   LocVecI4[0] = TstI4;
   LocVecI8[0] = TstI8;
   LocVecR4[0] = TstR4;
   LocVecR8[0] = TstR8;
   LocVecI4[1] = TaskRev * FacI4[TaskRev];
   LocVecI8[1] = TaskRev * FacI8[TaskRev];
   LocVecR4[1] = (TaskRev + EpsR4) * FacR4[TaskRev];
   LocVecR8[1] = (TaskRev + EpsR8) * FacR8[TaskRev];

   // Test multi-field sums
   SumVecI4 = globalSum(LocVecI4, Comm);
   SumVecI8 = globalSum(LocVecI8, Comm);
   SumVecR4 = globalSum(LocVecR4, Comm);
   SumVecR8 = globalSum(LocVecR8, Comm);
   // Test sum for each field
   for (int I = 0; I < NFields; ++I) {
      if (SumVecI4[I] != RefI4)
         ABORT_ERROR("ReductionsTest: FAIL globalSum (I4 scalar multifield)"
                     "Field {} Expected = {} Actual = {}",
                     I, RefI4, SumVecI4[I]);
      if (SumVecI8[I] != RefI8)
         ABORT_ERROR("ReductionsTest: FAIL globalSum (I8 scalar multifield)"
                     "Field {} Expected = {} Actual = {}",
                     I, RefI8, SumVecI8[I]);
      if (SumVecR4[I] != RefR4)
         ABORT_ERROR("ReductionsTest: FAIL globalSum (R4 scalar multifield)"
                     "Field {} Expected = {} Actual = {}",
                     I, RefR4, SumVecR4[I]);
      if (SumVecR8[I] != RefR8)
         ABORT_ERROR("ReductionsTest: FAIL globalSum (R8 scalar multifield)"
                     "Field {} Expected = {} Actual = {}",
                     I, RefR8, SumVecR8[I]);
   }

} // End testScalarSums

//------------------------------------------------------------------------------
// Individual array test functions
template <typename AT, typename T>
void testArray(const std::string &TestLabel, ///< [in] label for test
               const AT Array,               ///< [in] array to be summed
               const T RefResult,            ///< [in] reference result
               const std::vector<I4> *IndxRange = nullptr ///< [in] index range
) {
   MachEnv *DefEnv = MachEnv::getDefault();
   MPI_Comm Comm   = DefEnv->getComm();
   T TestResult    = 0;

   if (IndxRange == nullptr) {
      TestResult = globalSum(Array, Comm);
   } else {
      TestResult = globalSum(Array, Comm, IndxRange);
   }

   // temporarily document result
   LOG_INFO(" GlobalSum Test {}: Expected {} Actual {}", TestLabel, RefResult,
            TestResult);
   // Exit on error
   if (TestResult != RefResult)
      ABORT_ERROR("GlobalSum {} Test FAIL: Expected {} Actual {}", TestLabel,
                  RefResult, TestResult);
}
//------------------------------------------------------------------------------
// Individual array product test functions
template <typename AT, typename T>
void testArrayProd(
    const std::string &TestLabel,              ///< [in] label for test
    const AT Arr1,                             ///< [in] first array in product
    const AT Arr2,                             ///< [in] second array in product
    const T RefResult,                         ///< [in] reference result
    const std::vector<I4> *IndxRange = nullptr ///< [in] index range
) {
   MachEnv *DefEnv = MachEnv::getDefault();
   MPI_Comm Comm   = DefEnv->getComm();
   T TestResult    = 0;

   if (IndxRange == nullptr) {
      TestResult = globalSum(Arr1, Arr2, Comm);
   } else {
      TestResult = globalSum(Arr1, Arr2, Comm, IndxRange);
   }

   // temporarily document result
   LOG_INFO(" GlobalSum Test {}: Expected {} Actual {}", TestLabel, RefResult,
            TestResult);
   // Exit on error
   if (TestResult != RefResult)
      ABORT_ERROR("GlobalSum {} Test FAIL: Expected {} Actual {}", TestLabel,
                  RefResult, TestResult);
}
//------------------------------------------------------------------------------
// Array sum test function
void testArraySums() {

   // Get some MPI values based on default environment
   MachEnv *DefEnv = MachEnv::getDefault();
   int MyTask      = DefEnv->getMyTask();
   int NTasks      = DefEnv->getNumTasks();

   // Set model size for array tests
   int Nx        = 5;
   int Ny        = 5;
   int Nz        = 3;
   int Nm        = 3;
   int Nn        = 3;
   int NxGlob    = Nx * NTasks;
   I8 MaxEntries = NxGlob * Ny * Nz * Nm * Nn;

   // For reproducibility tests, define a full range of values within
   // each data type.  We want the smallest (abs) value > 0
   // Further restrict the max values so that we do not inadvertantly
   // exceed limits.
   I4 MaxI4 = std::numeric_limits<I4>::max() / (10 * MaxEntries);
   I8 MaxI8 = std::numeric_limits<I8>::max() / (10 * MaxEntries);
   R4 MaxR4 = std::numeric_limits<R4>::max() / (10.0 * MaxEntries);
   R8 MaxR8 = std::numeric_limits<R8>::max() / (10.0 * MaxEntries);
   R4 EpsR4 = 0.0001;
   R8 EpsR8 = 0.0000000000001;
   // Compute max safe exponents based on adjusted max above
   I4 ExpI4 = std::log2(MaxI4);
   I4 ExpI8 = std::log2(MaxI8);
   I4 ExpR4 = std::log10(MaxR4);
   I4 ExpR8 = std::log10(MaxR8);
   // To cover a large range of values, compute a factor for test vals
   // on each rank
   std::vector<I4> FacI4(NTasks);
   std::vector<I8> FacI8(NTasks);
   std::vector<R4> FacR4(NTasks);
   std::vector<R8> FacR8(NTasks);
   for (int I = 0; I < NTasks; ++I) {
      FacI4[I] = pow(2, std::min(I, ExpI4));
      FacI8[I] = pow(2, std::min(I, ExpI8));
      FacR4[I] = pow(10.0, std::min(I, ExpR4));
      FacR8[I] = pow(10.0, std::min(I, ExpR8));
   }

   // Allocate test arrays
   HostArray1DI4 TestHost1DI4("Host1DI4", Nx);
   HostArray2DI4 TestHost2DI4("Host2DI4", Nx, Ny);
   HostArray3DI4 TestHost3DI4("Host3DI4", Nx, Ny, Nz);
   HostArray4DI4 TestHost4DI4("Host4DI4", Nx, Ny, Nz, Nm);
   HostArray5DI4 TestHost5DI4("Host5DI4", Nx, Ny, Nz, Nm, Nn);
   HostArray1DI8 TestHost1DI8("Host1DI8", Nx);
   HostArray2DI8 TestHost2DI8("Host2DI8", Nx, Ny);
   HostArray3DI8 TestHost3DI8("Host3DI8", Nx, Ny, Nz);
   HostArray4DI8 TestHost4DI8("Host4DI8", Nx, Ny, Nz, Nm);
   HostArray5DI8 TestHost5DI8("Host5DI8", Nx, Ny, Nz, Nm, Nn);
   HostArray1DR4 TestHost1DR4("Host1DR4", Nx);
   HostArray2DR4 TestHost2DR4("Host2DR4", Nx, Ny);
   HostArray3DR4 TestHost3DR4("Host3DR4", Nx, Ny, Nz);
   HostArray4DR4 TestHost4DR4("Host4DR4", Nx, Ny, Nz, Nm);
   HostArray5DR4 TestHost5DR4("Host5DR4", Nx, Ny, Nz, Nm, Nn);
   HostArray1DR8 TestHost1DR8("Host1DR8", Nx);
   HostArray2DR8 TestHost2DR8("Host2DR8", Nx, Ny);
   HostArray3DR8 TestHost3DR8("Host3DR8", Nx, Ny, Nz);
   HostArray4DR8 TestHost4DR8("Host4DR8", Nx, Ny, Nz, Nm);
   HostArray5DR8 TestHost5DR8("Host5DR8", Nx, Ny, Nz, Nm, Nn);
   Array1DI4 Test1DI4("Test1DI4", Nx);
   Array2DI4 Test2DI4("Test2DI4", Nx, Ny);
   Array3DI4 Test3DI4("Test3DI4", Nx, Ny, Nz);
   Array4DI4 Test4DI4("Test4DI4", Nx, Ny, Nz, Nm);
   Array5DI4 Test5DI4("Test5DI4", Nx, Ny, Nz, Nm, Nn);
   Array1DI8 Test1DI8("Test1DI8", Nx);
   Array2DI8 Test2DI8("Test2DI8", Nx, Ny);
   Array3DI8 Test3DI8("Test3DI8", Nx, Ny, Nz);
   Array4DI8 Test4DI8("Test4DI8", Nx, Ny, Nz, Nm);
   Array5DI8 Test5DI8("Test5DI8", Nx, Ny, Nz, Nm, Nn);
   Array1DR4 Test1DR4("Test1DR4", Nx);
   Array2DR4 Test2DR4("Test2DR4", Nx, Ny);
   Array3DR4 Test3DR4("Test3DR4", Nx, Ny, Nz);
   Array4DR4 Test4DR4("Test4DR4", Nx, Ny, Nz, Nm);
   Array5DR4 Test5DR4("Test5DR4", Nx, Ny, Nz, Nm, Nn);
   Array1DR8 Test1DR8("Test1DR8", Nx);
   Array2DR8 Test2DR8("Test2DR8", Nx, Ny);
   Array3DR8 Test3DR8("Test3DR8", Nx, Ny, Nz);
   Array4DR8 Test4DR8("Test4DR8", Nx, Ny, Nz, Nm);
   Array5DR8 Test5DR8("Test5DR8", Nx, Ny, Nz, Nm, Nn);

   // Compute reference values using a serial reproducible sum
   I4 Ref1DI4 = 0;
   I8 Ref1DI8 = 0;
   R4 Ref1DR4 = 0.0;
   R8 Ref1DR8 = 0.0;
   I4 Ref2DI4 = 0;
   I8 Ref2DI8 = 0;
   R4 Ref2DR4 = 0.0;
   R8 Ref2DR8 = 0.0;
   I4 Ref3DI4 = 0;
   I8 Ref3DI8 = 0;
   R4 Ref3DR4 = 0.0;
   R8 Ref3DR8 = 0.0;
   I4 Ref4DI4 = 0;
   I8 Ref4DI8 = 0;
   R4 Ref4DR4 = 0.0;
   R8 Ref4DR8 = 0.0;
   I4 Ref5DI4 = 0;
   I8 Ref5DI8 = 0;
   R4 Ref5DR4 = 0.0;
   R8 Ref5DR8 = 0.0;
   R8 Tmp1DR4 = 0.0;
   R8 Tmp2DR4 = 0.0;
   R8 Tmp3DR4 = 0.0;
   R8 Tmp4DR4 = 0.0;
   R8 Tmp5DR4 = 0.0;
   R8 DDVal1D = 0.0;
   R8 DDVal2D = 0.0;
   R8 DDVal3D = 0.0;
   R8 DDVal4D = 0.0;
   R8 DDVal5D = 0.0;
   complex<double> DDSumRef1D(0.0, 0.0);
   complex<double> DDSumRef2D(0.0, 0.0);
   complex<double> DDSumRef3D(0.0, 0.0);
   complex<double> DDSumRef4D(0.0, 0.0);
   complex<double> DDSumRef5D(0.0, 0.0);
   // Compute reference sums
   for (int Task = 0; Task < NTasks; ++Task) {
      for (int I = 0; I < Nx; ++I) {
         int IGlob = Task * Nx + I;
         Ref1DI4 += IGlob * FacI4[Task];
         Ref1DI8 += IGlob * FacI8[Task];
         Tmp1DR4 += (IGlob + EpsR4) * FacR4[Task];
         DDVal1D = (IGlob + EpsR8) * FacR8[Task];
         sumDDTest(DDSumRef1D, DDVal1D); // local repro sum
         for (int J = 0; J < Ny; ++J) {
            int Jindx = IGlob + J;
            Ref2DI4 += Jindx * FacI4[Task];
            Ref2DI8 += Jindx * FacI8[Task];
            Tmp2DR4 += (Jindx + EpsR4) * FacR4[Task];
            DDVal2D = (Jindx + EpsR8) * FacR8[Task];
            sumDDTest(DDSumRef2D, DDVal2D); // local repro sum
            for (int K = 0; K < Nz; ++K) {
               int Kindx = IGlob + J + K;
               Ref3DI4 += Kindx * FacI4[Task];
               Ref3DI8 += Kindx * FacI8[Task];
               Tmp3DR4 += (Kindx + EpsR4) * FacR4[Task];
               DDVal3D = (Kindx + EpsR8) * FacR8[Task];
               sumDDTest(DDSumRef3D, DDVal3D); // local repro sum
               for (int M = 0; M < Nm; ++M) {
                  int Mindx = IGlob + J + K + M;
                  Ref4DI4 += Mindx * FacI4[Task];
                  Ref4DI8 += Mindx * FacI8[Task];
                  Tmp4DR4 += (Mindx + EpsR4) * FacR4[Task];
                  DDVal4D = (Mindx + EpsR8) * FacR8[Task];
                  sumDDTest(DDSumRef4D, DDVal4D); // local repro sum
                  for (int N = 0; N < Nn; ++N) {
                     int Nindx = IGlob + J + K + M + N;
                     Ref5DI4 += Nindx * FacI4[Task];
                     Ref5DI8 += Nindx * FacI8[Task];
                     Tmp5DR4 += (Nindx + EpsR4) * FacR4[Task];
                     DDVal5D = (Nindx + EpsR8) * FacR8[Task];
                     sumDDTest(DDSumRef5D, DDVal5D); // local repro sum
                  }
               }
            }
         }
      }
   }
   Ref1DR4 = Tmp1DR4;
   Ref2DR4 = Tmp2DR4;
   Ref3DR4 = Tmp3DR4;
   Ref4DR4 = Tmp4DR4;
   Ref5DR4 = Tmp5DR4;
   Ref1DR8 = real(DDSumRef1D);
   Ref2DR8 = real(DDSumRef2D);
   Ref3DR8 = real(DDSumRef3D);
   Ref4DR8 = real(DDSumRef4D);
   Ref5DR8 = real(DDSumRef5D);

   // Fill host arrays
   for (int I = 0; I < Nx; ++I) {
      int IGlobal     = MyTask * Nx + I;
      TestHost1DI4(I) = IGlobal * FacI4[MyTask];
      TestHost1DI8(I) = IGlobal * FacI8[MyTask];
      TestHost1DR4(I) = (IGlobal + EpsR4) * FacR4[MyTask];
      TestHost1DR8(I) = (IGlobal + EpsR8) * FacR8[MyTask];
      for (int J = 0; J < Ny; ++J) {
         int JIndx          = IGlobal + J;
         TestHost2DI4(I, J) = JIndx * FacI4[MyTask];
         TestHost2DI8(I, J) = JIndx * FacI8[MyTask];
         TestHost2DR4(I, J) = (JIndx + EpsR4) * FacR4[MyTask];
         TestHost2DR8(I, J) = (JIndx + EpsR8) * FacR8[MyTask];
         for (int K = 0; K < Nz; ++K) {
            int KIndx             = IGlobal + J + K;
            TestHost3DI4(I, J, K) = KIndx * FacI4[MyTask];
            TestHost3DI8(I, J, K) = KIndx * FacI8[MyTask];
            TestHost3DR4(I, J, K) = (KIndx + EpsR4) * FacR4[MyTask];
            TestHost3DR8(I, J, K) = (KIndx + EpsR8) * FacR8[MyTask];
            for (int M = 0; M < Nm; ++M) {
               int MIndx                = IGlobal + J + K + M;
               TestHost4DI4(I, J, K, M) = MIndx * FacI4[MyTask];
               TestHost4DI8(I, J, K, M) = MIndx * FacI8[MyTask];
               TestHost4DR4(I, J, K, M) = (MIndx + EpsR4) * FacR4[MyTask];
               TestHost4DR8(I, J, K, M) = (MIndx + EpsR8) * FacR8[MyTask];
               for (int N = 0; N < Nn; ++N) {
                  int NIndx                   = IGlobal + J + K + M + N;
                  TestHost5DI4(I, J, K, M, N) = NIndx * FacI4[MyTask];
                  TestHost5DI8(I, J, K, M, N) = NIndx * FacI8[MyTask];
                  TestHost5DR4(I, J, K, M, N) = (NIndx + EpsR4) * FacR4[MyTask];
                  TestHost5DR8(I, J, K, M, N) = (NIndx + EpsR8) * FacR8[MyTask];
               }
            }
         }
      }
   }

   // Copy values to device test arrays
   deepCopy(Test1DI4, TestHost1DI4);
   deepCopy(Test2DI4, TestHost2DI4);
   deepCopy(Test3DI4, TestHost3DI4);
   deepCopy(Test4DI4, TestHost4DI4);
   deepCopy(Test5DI4, TestHost5DI4);
   deepCopy(Test1DI8, TestHost1DI8);
   deepCopy(Test2DI8, TestHost2DI8);
   deepCopy(Test3DI8, TestHost3DI8);
   deepCopy(Test4DI8, TestHost4DI8);
   deepCopy(Test5DI8, TestHost5DI8);
   deepCopy(Test1DR4, TestHost1DR4);
   deepCopy(Test2DR4, TestHost2DR4);
   deepCopy(Test3DR4, TestHost3DR4);
   deepCopy(Test4DR4, TestHost4DR4);
   deepCopy(Test5DR4, TestHost5DR4);
   deepCopy(Test1DR8, TestHost1DR8);
   deepCopy(Test2DR8, TestHost2DR8);
   deepCopy(Test3DR8, TestHost3DR8);
   deepCopy(Test4DR8, TestHost4DR8);
   deepCopy(Test5DR8, TestHost5DR8);

   //-----
   // Test results for full arrays
   testArray("1DI4Host Full", TestHost1DI4, Ref1DI4);
   testArray("1DI8Host Full", TestHost1DI8, Ref1DI8);
   testArray("1DR4Host Full", TestHost1DR4, Ref1DR4);
   testArray("1DR8Host Full", TestHost1DR8, Ref1DR8);
   testArray("2DI4Host Full", TestHost2DI4, Ref2DI4);
   testArray("2DI8Host Full", TestHost2DI8, Ref2DI8);
   testArray("2DR4Host Full", TestHost2DR4, Ref2DR4);
   testArray("2DR8Host Full", TestHost2DR8, Ref2DR8);
   testArray("3DI4Host Full", TestHost3DI4, Ref3DI4);
   testArray("3DI8Host Full", TestHost3DI8, Ref3DI8);
   testArray("3DR4Host Full", TestHost3DR4, Ref3DR4);
   testArray("3DR8Host Full", TestHost3DR8, Ref3DR8);
   testArray("4DI4Host Full", TestHost4DI4, Ref4DI4);
   testArray("4DI8Host Full", TestHost4DI8, Ref4DI8);
   testArray("4DR4Host Full", TestHost4DR4, Ref4DR4);
   testArray("4DR8Host Full", TestHost4DR8, Ref4DR8);
   testArray("5DI4Host Full", TestHost5DI4, Ref5DI4);
   testArray("5DI8Host Full", TestHost5DI8, Ref5DI8);
   testArray("5DR4Host Full", TestHost5DR4, Ref5DR4);
   testArray("5DR8Host Full", TestHost5DR8, Ref5DR8);

   testArray("1DI4Dev Full", Test1DI4, Ref1DI4);
   testArray("1DI8Dev Full", Test1DI8, Ref1DI8);
   testArray("1DR4Dev Full", Test1DR4, Ref1DR4);
   testArray("1DR8Dev Full", Test1DR8, Ref1DR8);
   testArray("2DI4Dev Full", Test2DI4, Ref2DI4);
   testArray("2DI8Dev Full", Test2DI8, Ref2DI8);
   testArray("2DR4Dev Full", Test2DR4, Ref2DR4);
   testArray("2DR8Dev Full", Test2DR8, Ref2DR8);
   testArray("3DI4Dev Full", Test3DI4, Ref3DI4);
   testArray("3DI8Dev Full", Test3DI8, Ref3DI8);
   testArray("3DR4Dev Full", Test3DR4, Ref3DR4);
   testArray("3DR8Dev Full", Test3DR8, Ref3DR8);
   testArray("4DI4Dev Full", Test4DI4, Ref4DI4);
   testArray("4DI8Dev Full", Test4DI8, Ref4DI8);
   testArray("4DR4Dev Full", Test4DR4, Ref4DR4);
   testArray("4DR8Dev Full", Test4DR8, Ref4DR8);
   testArray("5DI4Dev Full", Test5DI4, Ref5DI4);
   testArray("5DI8Dev Full", Test5DI8, Ref5DI8);
   testArray("5DR4Dev Full", Test5DR4, Ref5DR4);
   testArray("5DR8Dev Full", Test5DR8, Ref5DR8);

   //-----
   // Test arrays with index range and restrictions and test sums with product
   // by creating a mask that is only 1 where index range is valid and compare.
   //-----

   // Initialize mask arrays
   HostArray1DI4 MaskHost1DI4("MaskH1DI4", Nx);
   HostArray1DI8 MaskHost1DI8("MaskH1DI8", Nx);
   HostArray1DR4 MaskHost1DR4("MaskH1DR4", Nx);
   HostArray1DR8 MaskHost1DR8("MaskH1DR8", Nx);
   HostArray2DI4 MaskHost2DI4("MaskH2DI4", Nx, Ny);
   HostArray2DI8 MaskHost2DI8("MaskH2DI8", Nx, Ny);
   HostArray2DR4 MaskHost2DR4("MaskH2DR4", Nx, Ny);
   HostArray2DR8 MaskHost2DR8("MaskH2DR8", Nx, Ny);
   HostArray3DI4 MaskHost3DI4("MaskH3DI4", Nx, Ny, Nz);
   HostArray3DI8 MaskHost3DI8("MaskH3DI8", Nx, Ny, Nz);
   HostArray3DR4 MaskHost3DR4("MaskH3DR4", Nx, Ny, Nz);
   HostArray3DR8 MaskHost3DR8("MaskH3DR8", Nx, Ny, Nz);
   HostArray4DI4 MaskHost4DI4("MaskH4DI4", Nx, Ny, Nz, Nm);
   HostArray4DI8 MaskHost4DI8("MaskH4DI8", Nx, Ny, Nz, Nm);
   HostArray4DR4 MaskHost4DR4("MaskH4DR4", Nx, Ny, Nz, Nm);
   HostArray4DR8 MaskHost4DR8("MaskH4DR8", Nx, Ny, Nz, Nm);
   HostArray5DI4 MaskHost5DI4("MaskH5DI4", Nx, Ny, Nz, Nm, Nn);
   HostArray5DI8 MaskHost5DI8("MaskH5DI8", Nx, Ny, Nz, Nm, Nn);
   HostArray5DR4 MaskHost5DR4("MaskH5DR4", Nx, Ny, Nz, Nm, Nn);
   HostArray5DR8 MaskHost5DR8("MaskH5DR8", Nx, Ny, Nz, Nm, Nn);
   Array1DI4 Mask1DI4("Mask1DI4", Nx);
   Array1DI8 Mask1DI8("Mask1DI8", Nx);
   Array1DR4 Mask1DR4("Mask1DR4", Nx);
   Array1DR8 Mask1DR8("Mask1DR8", Nx);
   Array2DI4 Mask2DI4("Mask2DI4", Nx, Ny);
   Array2DI8 Mask2DI8("Mask2DI8", Nx, Ny);
   Array2DR4 Mask2DR4("Mask2DR4", Nx, Ny);
   Array2DR8 Mask2DR8("Mask2DR8", Nx, Ny);
   Array3DI4 Mask3DI4("Mask3DI4", Nx, Ny, Nz);
   Array3DI8 Mask3DI8("Mask3DI8", Nx, Ny, Nz);
   Array3DR4 Mask3DR4("Mask3DR4", Nx, Ny, Nz);
   Array3DR8 Mask3DR8("Mask3DR8", Nx, Ny, Nz);
   Array4DI4 Mask4DI4("Mask4DI4", Nx, Ny, Nz, Nm);
   Array4DI8 Mask4DI8("Mask4DI8", Nx, Ny, Nz, Nm);
   Array4DR4 Mask4DR4("Mask4DR4", Nx, Ny, Nz, Nm);
   Array4DR8 Mask4DR8("Mask4DR8", Nx, Ny, Nz, Nm);
   Array5DI4 Mask5DI4("Mask5DI4", Nx, Ny, Nz, Nm, Nn);
   Array5DI8 Mask5DI8("Mask5DI8", Nx, Ny, Nz, Nm, Nn);
   Array5DR4 Mask5DR4("Mask5DR4", Nx, Ny, Nz, Nm, Nn);
   Array5DR8 Mask5DR8("Mask5DR8", Nx, Ny, Nz, Nm, Nn);

   // Set restricted index range
   int IMin = 2;
   int IMax = Nx - 2;
   int JMin = 1;
   int JMax = Ny - 1;
   int KMin = 1;
   int KMax = 1;
   int MMin = 1;
   int MMax = 1;
   int NMin = 1;
   int NMax = 1;
   std::vector<int> AddRange(10);
   AddRange[0] = IMin;
   AddRange[1] = IMax;
   AddRange[2] = JMin;
   AddRange[3] = JMax;
   AddRange[4] = KMin;
   AddRange[5] = KMax;
   AddRange[6] = MMin;
   AddRange[7] = MMax;
   AddRange[8] = NMin;
   AddRange[9] = NMax;
   // Reset various sums and compute new reference values
   Ref1DI4    = 0;
   Ref1DI8    = 0;
   Ref1DR4    = 0.0;
   Ref1DR8    = 0.0;
   Tmp1DR4    = 0.0;
   DDVal1D    = 0.0;
   Ref2DI4    = 0;
   Ref2DI8    = 0;
   Ref2DR4    = 0.0;
   Ref2DR8    = 0.0;
   Tmp2DR4    = 0.0;
   DDVal2D    = 0.0;
   Ref3DI4    = 0;
   Ref3DI8    = 0;
   Ref3DR4    = 0.0;
   Ref3DR8    = 0.0;
   Tmp3DR4    = 0.0;
   DDVal3D    = 0.0;
   Ref4DI4    = 0;
   Ref4DI8    = 0;
   Ref4DR4    = 0.0;
   Ref4DR8    = 0.0;
   Tmp4DR4    = 0.0;
   DDVal4D    = 0.0;
   Ref5DI4    = 0;
   Ref5DI8    = 0;
   Ref5DR4    = 0.0;
   Ref5DR8    = 0.0;
   Tmp5DR4    = 0.0;
   DDVal5D    = 0.0;
   DDSumRef1D = complex<double>(0.0, 0.0);
   DDSumRef2D = complex<double>(0.0, 0.0);
   DDSumRef3D = complex<double>(0.0, 0.0);
   DDSumRef4D = complex<double>(0.0, 0.0);
   DDSumRef5D = complex<double>(0.0, 0.0);
   // Compute new reference sums for restricted range
   for (int Task = 0; Task < NTasks; ++Task) {
      for (int I = IMin; I <= IMax; ++I) {
         int IGlob = Task * Nx + I;
         Ref1DI4 += IGlob * FacI4[Task];
         Ref1DI8 += IGlob * FacI8[Task];
         Tmp1DR4 += (IGlob + EpsR4) * FacR4[Task];
         DDVal1D = (IGlob + EpsR8) * FacR8[Task];
         sumDDTest(DDSumRef1D, DDVal1D); // local repro sum
         for (int J = JMin; J <= JMax; ++J) {
            int Jindx = IGlob + J;
            Ref2DI4 += Jindx * FacI4[Task];
            Ref2DI8 += Jindx * FacI8[Task];
            Tmp2DR4 += (Jindx + EpsR4) * FacR4[Task];
            DDVal2D = (Jindx + EpsR8) * FacR8[Task];
            sumDDTest(DDSumRef2D, DDVal2D); // local repro sum
            for (int K = KMin; K <= KMax; ++K) {
               int Kindx = IGlob + J + K;
               Ref3DI4 += Kindx * FacI4[Task];
               Ref3DI8 += Kindx * FacI8[Task];
               Tmp3DR4 += (Kindx + EpsR4) * FacR4[Task];
               DDVal3D = (Kindx + EpsR8) * FacR8[Task];
               sumDDTest(DDSumRef3D, DDVal3D); // local repro sum
               for (int M = MMin; M <= MMax; ++M) {
                  int Mindx = IGlob + J + K + M;
                  Ref4DI4 += Mindx * FacI4[Task];
                  Ref4DI8 += Mindx * FacI8[Task];
                  Tmp4DR4 += (Mindx + EpsR4) * FacR4[Task];
                  DDVal4D = (Mindx + EpsR8) * FacR8[Task];
                  sumDDTest(DDSumRef4D, DDVal4D); // local repro sum
                  for (int N = NMin; N <= NMax; ++N) {
                     int Nindx = IGlob + J + K + M + N;
                     Ref5DI4 += Nindx * FacI4[Task];
                     Ref5DI8 += Nindx * FacI8[Task];
                     Tmp5DR4 += (Nindx + EpsR4) * FacR4[Task];
                     DDVal5D = (Nindx + EpsR8) * FacR8[Task];
                     sumDDTest(DDSumRef5D, DDVal5D); // local repro sum
                  }
               }
            }
         }
      }
   }
   Ref1DR4 = Tmp1DR4;
   Ref2DR4 = Tmp2DR4;
   Ref3DR4 = Tmp3DR4;
   Ref4DR4 = Tmp4DR4;
   Ref5DR4 = Tmp5DR4;
   Ref1DR8 = real(DDSumRef1D);
   Ref2DR8 = real(DDSumRef2D);
   Ref3DR8 = real(DDSumRef3D);
   Ref4DR8 = real(DDSumRef4D);
   Ref5DR8 = real(DDSumRef5D);

   // Fill mask arrays on host
   for (int I = IMin; I <= IMax; ++I) {
      MaskHost1DI4(I) = 1;
      MaskHost1DI8(I) = 1;
      MaskHost1DR4(I) = 1.0;
      MaskHost1DR8(I) = 1.0;
      for (int J = JMin; J <= JMax; ++J) {
         MaskHost2DI4(I, J) = 1;
         MaskHost2DI8(I, J) = 1;
         MaskHost2DR4(I, J) = 1.0;
         MaskHost2DR8(I, J) = 1.0;
         for (int K = KMin; K <= KMax; ++K) {
            MaskHost3DI4(I, J, K) = 1;
            MaskHost3DI8(I, J, K) = 1;
            MaskHost3DR4(I, J, K) = 1.0;
            MaskHost3DR8(I, J, K) = 1.0;
            for (int M = MMin; M <= MMax; ++M) {
               MaskHost4DI4(I, J, K, M) = 1;
               MaskHost4DI8(I, J, K, M) = 1;
               MaskHost4DR4(I, J, K, M) = 1.0;
               MaskHost4DR8(I, J, K, M) = 1.0;
               for (int N = NMin; N <= NMax; ++N) {
                  MaskHost5DI4(I, J, K, M, N) = 1;
                  MaskHost5DI8(I, J, K, M, N) = 1;
                  MaskHost5DR4(I, J, K, M, N) = 1.0;
                  MaskHost5DR8(I, J, K, M, N) = 1.0;
               }
            }
         }
      }
   }

   // create equivalent device arrays
   deepCopy(Mask1DI4, MaskHost1DI4);
   deepCopy(Mask1DI8, MaskHost1DI8);
   deepCopy(Mask1DR4, MaskHost1DR4);
   deepCopy(Mask1DR8, MaskHost1DR8);
   deepCopy(Mask2DI4, MaskHost2DI4);
   deepCopy(Mask2DI8, MaskHost2DI8);
   deepCopy(Mask2DR4, MaskHost2DR4);
   deepCopy(Mask2DR8, MaskHost2DR8);
   deepCopy(Mask3DI4, MaskHost3DI4);
   deepCopy(Mask3DI8, MaskHost3DI8);
   deepCopy(Mask3DR4, MaskHost3DR4);
   deepCopy(Mask3DR8, MaskHost3DR8);
   deepCopy(Mask4DI4, MaskHost4DI4);
   deepCopy(Mask4DI8, MaskHost4DI8);
   deepCopy(Mask4DR4, MaskHost4DR4);
   deepCopy(Mask4DR8, MaskHost4DR8);
   deepCopy(Mask5DI4, MaskHost5DI4);
   deepCopy(Mask5DI8, MaskHost5DI8);
   deepCopy(Mask5DR4, MaskHost5DR4);
   deepCopy(Mask5DR8, MaskHost5DR8);

   //-----
   // Test sums with range limits
   testArray("1DI4Host Add Range", TestHost1DI4, Ref1DI4, &AddRange);
   testArray("1DI8Host Add Range", TestHost1DI8, Ref1DI8, &AddRange);
   testArray("1DR4Host Add Range", TestHost1DR4, Ref1DR4, &AddRange);
   testArray("1DR8Host Add Range", TestHost1DR8, Ref1DR8, &AddRange);
   testArray("2DI4Host Add Range", TestHost2DI4, Ref2DI4, &AddRange);
   testArray("2DI8Host Add Range", TestHost2DI8, Ref2DI8, &AddRange);
   testArray("2DR4Host Add Range", TestHost2DR4, Ref2DR4, &AddRange);
   testArray("2DR8Host Add Range", TestHost2DR8, Ref2DR8, &AddRange);
   testArray("3DI4Host Add Range", TestHost3DI4, Ref3DI4, &AddRange);
   testArray("3DI8Host Add Range", TestHost3DI8, Ref3DI8, &AddRange);
   testArray("3DR4Host Add Range", TestHost3DR4, Ref3DR4, &AddRange);
   testArray("3DR8Host Add Range", TestHost3DR8, Ref3DR8, &AddRange);
   testArray("4DI4Host Add Range", TestHost4DI4, Ref4DI4, &AddRange);
   testArray("4DI8Host Add Range", TestHost4DI8, Ref4DI8, &AddRange);
   testArray("4DR4Host Add Range", TestHost4DR4, Ref4DR4, &AddRange);
   testArray("4DR8Host Add Range", TestHost4DR8, Ref4DR8, &AddRange);
   testArray("5DI4Host Add Range", TestHost5DI4, Ref5DI4, &AddRange);
   testArray("5DI8Host Add Range", TestHost5DI8, Ref5DI8, &AddRange);
   testArray("5DR4Host Add Range", TestHost5DR4, Ref5DR4, &AddRange);
   testArray("5DR8Host Add Range", TestHost5DR8, Ref5DR8, &AddRange);

   testArray("1DI4Dev Add Range", Test1DI4, Ref1DI4, &AddRange);
   testArray("1DI8Dev Add Range", Test1DI8, Ref1DI8, &AddRange);
   testArray("1DR4Dev Add Range", Test1DR4, Ref1DR4, &AddRange);
   testArray("1DR8Dev Add Range", Test1DR8, Ref1DR8, &AddRange);
   testArray("2DI4Dev Add Range", Test2DI4, Ref2DI4, &AddRange);
   testArray("2DI8Dev Add Range", Test2DI8, Ref2DI8, &AddRange);
   testArray("2DR4Dev Add Range", Test2DR4, Ref2DR4, &AddRange);
   testArray("2DR8Dev Add Range", Test2DR8, Ref2DR8, &AddRange);
   testArray("3DI4Dev Add Range", Test3DI4, Ref3DI4, &AddRange);
   testArray("3DI8Dev Add Range", Test3DI8, Ref3DI8, &AddRange);
   testArray("3DR4Dev Add Range", Test3DR4, Ref3DR4, &AddRange);
   testArray("3DR8Dev Add Range", Test3DR8, Ref3DR8, &AddRange);
   testArray("4DI4Dev Add Range", Test4DI4, Ref4DI4, &AddRange);
   testArray("4DI8Dev Add Range", Test4DI8, Ref4DI8, &AddRange);
   testArray("4DR4Dev Add Range", Test4DR4, Ref4DR4, &AddRange);
   testArray("4DR8Dev Add Range", Test4DR8, Ref4DR8, &AddRange);
   testArray("5DI4Dev Add Range", Test5DI4, Ref5DI4, &AddRange);
   testArray("5DI8Dev Add Range", Test5DI8, Ref5DI8, &AddRange);
   testArray("5DR4Dev Add Range", Test5DR4, Ref5DR4, &AddRange);
   testArray("5DR8Dev Add Range", Test5DR8, Ref5DR8, &AddRange);

   //-----
   // Test sums with product - full arrays
   testArrayProd("1DI4Host Product Full", TestHost1DI4, MaskHost1DI4, Ref1DI4);
   testArrayProd("1DI8Host Product Full", TestHost1DI8, MaskHost1DI8, Ref1DI8);
   testArrayProd("1DR4Host Product Full", TestHost1DR4, MaskHost1DR4, Ref1DR4);
   testArrayProd("1DR8Host Product Full", TestHost1DR8, MaskHost1DR8, Ref1DR8);
   testArrayProd("2DI4Host Product Full", TestHost2DI4, MaskHost2DI4, Ref2DI4);
   testArrayProd("2DI8Host Product Full", TestHost2DI8, MaskHost2DI8, Ref2DI8);
   testArrayProd("2DR4Host Product Full", TestHost2DR4, MaskHost2DR4, Ref2DR4);
   testArrayProd("2DR8Host Product Full", TestHost2DR8, MaskHost2DR8, Ref2DR8);
   testArrayProd("3DI4Host Product Full", TestHost3DI4, MaskHost3DI4, Ref3DI4);
   testArrayProd("3DI8Host Product Full", TestHost3DI8, MaskHost3DI8, Ref3DI8);
   testArrayProd("3DR4Host Product Full", TestHost3DR4, MaskHost3DR4, Ref3DR4);
   testArrayProd("3DR8Host Product Full", TestHost3DR8, MaskHost3DR8, Ref3DR8);
   testArrayProd("4DI4Host Product Full", TestHost4DI4, MaskHost4DI4, Ref4DI4);
   testArrayProd("4DI8Host Product Full", TestHost4DI8, MaskHost4DI8, Ref4DI8);
   testArrayProd("4DR4Host Product Full", TestHost4DR4, MaskHost4DR4, Ref4DR4);
   testArrayProd("4DR8Host Product Full", TestHost4DR8, MaskHost4DR8, Ref4DR8);
   testArrayProd("5DI4Host Product Full", TestHost5DI4, MaskHost5DI4, Ref5DI4);
   testArrayProd("5DI8Host Product Full", TestHost5DI8, MaskHost5DI8, Ref5DI8);
   testArrayProd("5DR4Host Product Full", TestHost5DR4, MaskHost5DR4, Ref5DR4);
   testArrayProd("5DR8Host Product Full", TestHost5DR8, MaskHost5DR8, Ref5DR8);

   testArrayProd("1DI4Dev Product Full", Test1DI4, Mask1DI4, Ref1DI4);
   testArrayProd("1DI8Dev Product Full", Test1DI8, Mask1DI8, Ref1DI8);
   testArrayProd("1DR4Dev Product Full", Test1DR4, Mask1DR4, Ref1DR4);
   testArrayProd("1DR8Dev Product Full", Test1DR8, Mask1DR8, Ref1DR8);
   testArrayProd("2DI4Dev Product Full", Test2DI4, Mask2DI4, Ref2DI4);
   testArrayProd("2DI8Dev Product Full", Test2DI8, Mask2DI8, Ref2DI8);
   testArrayProd("2DR4Dev Product Full", Test2DR4, Mask2DR4, Ref2DR4);
   testArrayProd("2DR8Dev Product Full", Test2DR8, Mask2DR8, Ref2DR8);
   testArrayProd("3DI4Dev Product Full", Test3DI4, Mask3DI4, Ref3DI4);
   testArrayProd("3DI8Dev Product Full", Test3DI8, Mask3DI8, Ref3DI8);
   testArrayProd("3DR4Dev Product Full", Test3DR4, Mask3DR4, Ref3DR4);
   testArrayProd("3DR8Dev Product Full", Test3DR8, Mask3DR8, Ref3DR8);
   testArrayProd("4DI4Dev Product Full", Test4DI4, Mask4DI4, Ref4DI4);
   testArrayProd("4DI8Dev Product Full", Test4DI8, Mask4DI8, Ref4DI8);
   testArrayProd("4DR4Dev Product Full", Test4DR4, Mask4DR4, Ref4DR4);
   testArrayProd("4DR8Dev Product Full", Test4DR8, Mask4DR8, Ref4DR8);
   testArrayProd("5DI4Dev Product Full", Test5DI4, Mask5DI4, Ref5DI4);
   testArrayProd("5DI8Dev Product Full", Test5DI8, Mask5DI8, Ref5DI8);
   testArrayProd("5DR4Dev Product Full", Test5DR4, Mask5DR4, Ref5DR4);
   testArrayProd("5DR8Dev Product Full", Test5DR8, Mask5DR8, Ref5DR8);

   //------
   // Test sums with product and subset range
   testArrayProd("1DI4Host Prod Range", TestHost1DI4, MaskHost1DI4, Ref1DI4,
                 &AddRange);
   testArrayProd("1DI8Host Prod Range", TestHost1DI8, MaskHost1DI8, Ref1DI8,
                 &AddRange);
   testArrayProd("1DR4Host Prod Range", TestHost1DR4, MaskHost1DR4, Ref1DR4,
                 &AddRange);
   testArrayProd("1DR8Host Prod Range", TestHost1DR8, MaskHost1DR8, Ref1DR8,
                 &AddRange);
   testArrayProd("2DI4Host Prod Range", TestHost2DI4, MaskHost2DI4, Ref2DI4,
                 &AddRange);
   testArrayProd("2DI8Host Prod Range", TestHost2DI8, MaskHost2DI8, Ref2DI8,
                 &AddRange);
   testArrayProd("2DR4Host Prod Range", TestHost2DR4, MaskHost2DR4, Ref2DR4,
                 &AddRange);
   testArrayProd("2DR8Host Prod Range", TestHost2DR8, MaskHost2DR8, Ref2DR8,
                 &AddRange);
   testArrayProd("3DI4Host Prod Range", TestHost3DI4, MaskHost3DI4, Ref3DI4,
                 &AddRange);
   testArrayProd("3DI8Host Prod Range", TestHost3DI8, MaskHost3DI8, Ref3DI8,
                 &AddRange);
   testArrayProd("3DR4Host Prod Range", TestHost3DR4, MaskHost3DR4, Ref3DR4,
                 &AddRange);
   testArrayProd("3DR8Host Prod Range", TestHost3DR8, MaskHost3DR8, Ref3DR8,
                 &AddRange);
   testArrayProd("4DI4Host Prod Range", TestHost4DI4, MaskHost4DI4, Ref4DI4,
                 &AddRange);
   testArrayProd("4DI8Host Prod Range", TestHost4DI8, MaskHost4DI8, Ref4DI8,
                 &AddRange);
   testArrayProd("4DR4Host Prod Range", TestHost4DR4, MaskHost4DR4, Ref4DR4,
                 &AddRange);
   testArrayProd("4DR8Host Prod Range", TestHost4DR8, MaskHost4DR8, Ref4DR8,
                 &AddRange);
   testArrayProd("5DI4Host Prod Range", TestHost5DI4, MaskHost5DI4, Ref5DI4,
                 &AddRange);
   testArrayProd("5DI8Host Prod Range", TestHost5DI8, MaskHost5DI8, Ref5DI8,
                 &AddRange);
   testArrayProd("5DR4Host Prod Range", TestHost5DR4, MaskHost5DR4, Ref5DR4,
                 &AddRange);
   testArrayProd("5DR8Host Prod Range", TestHost5DR8, MaskHost5DR8, Ref5DR8,
                 &AddRange);

   testArrayProd("1DI4Dev Prod Range", Test1DI4, Mask1DI4, Ref1DI4, &AddRange);
   testArrayProd("1DI8Dev Prod Range", Test1DI8, Mask1DI8, Ref1DI8, &AddRange);
   testArrayProd("1DR4Dev Prod Range", Test1DR4, Mask1DR4, Ref1DR4, &AddRange);
   testArrayProd("1DR8Dev Prod Range", Test1DR8, Mask1DR8, Ref1DR8, &AddRange);
   testArrayProd("2DI4Dev Prod Range", Test2DI4, Mask2DI4, Ref2DI4, &AddRange);
   testArrayProd("2DI8Dev Prod Range", Test2DI8, Mask2DI8, Ref2DI8, &AddRange);
   testArrayProd("2DR4Dev Prod Range", Test2DR4, Mask2DR4, Ref2DR4, &AddRange);
   testArrayProd("2DR8Dev Prod Range", Test2DR8, Mask2DR8, Ref2DR8, &AddRange);
   testArrayProd("3DI4Dev Prod Range", Test3DI4, Mask3DI4, Ref3DI4, &AddRange);
   testArrayProd("3DI8Dev Prod Range", Test3DI8, Mask3DI8, Ref3DI8, &AddRange);
   testArrayProd("3DR4Dev Prod Range", Test3DR4, Mask3DR4, Ref3DR4, &AddRange);
   testArrayProd("3DR8Dev Prod Range", Test3DR8, Mask3DR8, Ref3DR8, &AddRange);
   testArrayProd("4DI4Dev Prod Range", Test4DI4, Mask4DI4, Ref4DI4, &AddRange);
   testArrayProd("4DI8Dev Prod Range", Test4DI8, Mask4DI8, Ref4DI8, &AddRange);
   testArrayProd("4DR4Dev Prod Range", Test4DR4, Mask4DR4, Ref4DR4, &AddRange);
   testArrayProd("4DR8Dev Prod Range", Test4DR8, Mask4DR8, Ref4DR8, &AddRange);
   testArrayProd("5DI4Dev Prod Range", Test5DI4, Mask5DI4, Ref5DI4, &AddRange);
   testArrayProd("5DI8Dev Prod Range", Test5DI8, Mask5DI8, Ref5DI8, &AddRange);
   testArrayProd("5DR4Dev Prod Range", Test5DR4, Mask5DR4, Ref5DR4, &AddRange);
   testArrayProd("5DR8Dev Prod Range", Test5DR8, Mask5DR8, Ref5DR8, &AddRange);

} // End testArraySums

//------------------------------------------------------------------------------
// Main test driver
int main(int argc, char *argv[]) {

   // Initialize various environments and utilities
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   OMEGA::initLogging(DefEnv);
   LOG_INFO("------ Global Reductions Unit Tests ------");

   // For reproducibility across tasks, create a smaller sub-environment
   int NTasks = DefEnv->getNumTasks();
   if (NTasks < 8)
      ABORT_ERROR(
          "ReductionsTest: FAIL must run unit test with at least 8 tasks");
   MachEnv::create("Subset", DefEnv, 4); // contiguous subset environment

   // Call individual test routines
   testScalarSums();
   testArraySums();
   // testMinMax();

   //------------------------------------------------------------------------
   // Multifield sums
   // TODO

   //------------------------------------------------------------------------
   // Multifield sums with subset index range
   // TODO

   //------------------------------------------------------------------------
   // Multifield sums with product
   // TODO

   //------------------------------------------------------------------------
   // Multifield sums with product and subset index range
   // TODO

   //------------------------------------------------------------------------
   // Reproducibility test Array sums
   // TODO

   //------------------------------------------------------------------------
   // Reproducibility test Arrays sums with product
   // TODO

   LOG_INFO("------ Global Reductions Unit Tests Sucessful ------");

   // clean up
   MachEnv::removeAll();
   Kokkos::finalize();
   MPI_Finalize();

   return 0; // if we made it here, return success

} // end of main
//===-----------------------------------------------------------------------===/
