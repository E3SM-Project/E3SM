#ifndef OMEGA_REDUCTIONS_H
#define OMEGA_REDUCTIONS_H
//===-- base/Reductions.h - MPI reduction definitions -----------*- C++ -*-===//
//
/// \file
/// \brief Defines
///
//
//===----------------------------------------------------------------------===//

#include <complex>
using std::complex;

#include "DataTypes.h"
#include "OmegaKokkos.h"

namespace OMEGA {

static int R8SumInitialized = 0;

static MPI_Op MPI_SUMDD; // special MPI operator for reproducible R8 sum

void ddSum(void *InBuffer, void *OutBuffer, int *Len, MPI_Datatype *DataType) {
   complex<double> *dda = (complex<double> *)InBuffer;
   complex<double> *ddb = (complex<double> *)OutBuffer;
   double e, t1, t2;
   for (int i = 0; i < *Len; i++) {
      t1 = real(dda[i]) + real(ddb[i]);
      e  = t1 - real(dda[i]);
      t2 = ((real(ddb[i]) - e) + (real(dda[i]) - (t1 - e))) + imag(dda[i]) +
           imag(ddb[i]);
      ddb[i] = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
   }
}

int globalSumInit() {
   int ierr = MPI_Op_create(&ddSum, 1, &MPI_SUMDD);
   if (ierr == 0)
      R8SumInitialized = 1;
   return ierr;
}

///-----------------------------------------------------------------------------
/// Sum given values across all Comm's MPI processors
///-----------------------------------------------------------------------------

//////////
// Global sum scalars
//////////
// I4
int globalSum(const I4 *Val, const MPI_Comm Comm, I4 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_SUM, Comm);
}

// I8
int globalSum(const I8 *Val, const MPI_Comm Comm, I8 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_SUM, Comm);
}

// R4
int globalSum(const R4 *Val, const MPI_Comm Comm, R4 *Res) {
   R8 LocalTmp, GlobalTmp;
   LocalTmp = *Val;
   int ierr =
       MPI_Allreduce(&LocalTmp, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
   *Res = GlobalTmp;
   return ierr;
}

// R8
int globalSum(const R8 *Val, const MPI_Comm Comm, R8 *Res) {
   // initialize reproducible MPI_SUMDD operator
   if (!R8SumInitialized) {
      globalSumInit();
   }

   complex<double> LocalTmp(*Val, 0.0);
   complex<double> GlobalTmp(0.0, 0.0);

   int ierr = MPI_Allreduce(&LocalTmp, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                            MPI_SUMDD, Comm);
   *Res     = real(GlobalTmp);
   return ierr;
}

//////////
// Global sum arrays
//////////
// I4 or I8 array
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_integral_v<typename Kokkos::View<T>::value_type>, int>
globalSum(const Kokkos::View<T, ML, MS> arr, const MPI_Comm Comm, IT *GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   IT LocalSum = 0;
   int dim     = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   if (Kokkos::SpaceAccessibility<MS,
                                  Kokkos::HostSpace>::accessible) { // on host
      for (i = imin; i < imax; i++) {
         LocalSum += arr.data()[i];
      }
   } else { // on device
      parallelReduce(
          {imax}, KOKKOS_LAMBDA(int i, IT &Accum) { Accum += arr.data()[i]; },
          LocalSum);
   }
   if (typeid(IT) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalSum, GlobalSum, 1, MPI_INT32_T, MPI_SUM, Comm);
   } else {
      ierr = MPI_Allreduce(&LocalSum, GlobalSum, 1, MPI_INT64_T, MPI_SUM, Comm);
   }
   return ierr;
}

// R4 array
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, int>
globalSum(const Kokkos::View<T, ML, MS> arr, const MPI_Comm Comm, R4 *GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   R8 GlobalTmp = 0.0, LocalSum = 0.0;
   int dim = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   if (Kokkos::SpaceAccessibility<MS, Kokkos::HostSpace>::accessible) {
      for (i = imin; i < imax; i++) {
         LocalSum += arr.data()[i];
      }
   } else {
      parallelReduce(
          {imax}, KOKKOS_LAMBDA(int i, R8 &Accum) { Accum += arr.data()[i]; },
          LocalSum);
   }
   ierr = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
   *GlobalSum = GlobalTmp;
   return ierr;
}

// R8 array
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, int>
globalSum(const Kokkos::View<T, ML, MS> arr, const MPI_Comm Comm, R8 *GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int dim = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }

   if (Kokkos::SpaceAccessibility<MS, Kokkos::HostSpace>::accessible) {
      // Accumulate the local sum using Knuth's algorithm
      complex<double> LocalSum(0.0, 0.0), GlobalTmp(0.0, 0.0);
      double e, t1, t2, ai;
      for (i = imin; i < imax; i++) {
         ai = arr.data()[i];
         t1 = ai + real(LocalSum);
         e  = t1 - ai;
         t2 = ((real(LocalSum) - e) + (ai - (t1 - e))) + imag(LocalSum);
         // The result is t1 + t2, after normalization.
         LocalSum = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
      }
      ierr       = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                                 MPI_SUMDD, Comm);
      *GlobalSum = real(GlobalTmp);
   } else {
      R8 LocalSum = 0.0, GlobalTmp = 0.0;
      parallelReduce(
          {imax}, KOKKOS_LAMBDA(int i, R8 &Accum) { Accum += arr.data()[i]; },
          LocalSum);
      ierr = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
      *GlobalSum = GlobalTmp;
   }
   return ierr;
}

//////////
// Global sum with product
//////////
// I4 or I8 array
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_integral_v<typename Kokkos::View<T>::value_type>, int>
globalSum(const Kokkos::View<T, ML, MS> arr, const Kokkos::View<T, ML, MS> arr2,
          const MPI_Comm Comm, IT *GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   T LocalSum = 0;
   int dim    = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   if (Kokkos::SpaceAccessibility<MS,
                                  Kokkos::HostSpace>::accessible) { // on host
      for (i = imin; i < imax; i++) {
         LocalSum += arr.data()[i] * arr2.data()[i];
      }
   } else { // on device
      parallelReduce(
          {imax},
          KOKKOS_LAMBDA(int i, IT &Accum) {
             Accum += arr.data()[i] * arr2.data()[i];
          },
          LocalSum);
   }
   if (typeid(IT) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalSum, GlobalSum, 1, MPI_INT32_T, MPI_SUM, Comm);
   } else {
      ierr = MPI_Allreduce(&LocalSum, GlobalSum, 1, MPI_INT64_T, MPI_SUM, Comm);
   }
   return ierr;
}

// R4 array
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, int>
globalSum(const Kokkos::View<T, ML, MS> arr, const Kokkos::View<T, ML, MS> arr2,
          const MPI_Comm Comm, R4 *GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   R8 GlobalTmp = 0.0, LocalSum = 0.0;
   int dim = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   if (Kokkos::SpaceAccessibility<MS, Kokkos::HostSpace>::accessible) {
      for (i = imin; i < imax; i++) {
         LocalSum += arr.data()[i] * arr2.data()[i];
      }
   } else {
      parallelReduce(
          {imax},
          KOKKOS_LAMBDA(int i, R8 &Accum) {
             Accum += arr.data()[i] * arr2.data()[i];
          },
          LocalSum);
   }
   ierr = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
   *GlobalSum = GlobalTmp;
   return ierr;
}

// R8 array
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, int>
globalSum(const Kokkos::View<T, ML, MS> arr, const Kokkos::View<T, ML, MS> arr2,
          const MPI_Comm Comm, R8 *GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int dim = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }

   if (Kokkos::SpaceAccessibility<MS, Kokkos::HostSpace>::accessible) {
      // Accumulate the local sum using Knuth's algorithm
      complex<double> LocalSum(0.0, 0.0), GlobalTmp(0.0, 0.0);
      double e, t1, t2, ai;
      for (i = imin; i < imax; i++) {
         ai = arr.data()[i] * arr2.data()[i];
         t1 = ai + real(LocalSum);
         e  = t1 - ai;
         t2 = ((real(LocalSum) - e) + (ai - (t1 - e))) + imag(LocalSum);
         // The result is t1 + t2, after normalization.
         LocalSum = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
      }
      ierr       = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                                 MPI_SUMDD, Comm);
      *GlobalSum = real(GlobalTmp);
   } else {
      R8 LocalSum = 0.0, GlobalTmp = 0.0;
      parallelReduce(
          {imax},
          KOKKOS_LAMBDA(int i, R8 &Accum) {
             Accum += arr.data()[i] * arr2.data()[i];
          },
          LocalSum);
      ierr = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
      *GlobalSum = GlobalTmp;
   }
   return ierr;
}

//////////
// Global sum multi-field
//////////
// I4 scalars
int globalSum(const std::vector<I4> scalars, const MPI_Comm Comm,
              std::vector<I4> GlobalSum) {
   int nFlds = scalars.size();
   return MPI_Allreduce(&scalars[0], &GlobalSum[0], nFlds, MPI_INT32_T, MPI_SUM,
                        Comm);
}

// I8 scalars
int globalSum(const std::vector<I8> scalars, const MPI_Comm Comm,
              std::vector<I8> GlobalSum) {
   int nFlds = scalars.size();
   return MPI_Allreduce(&scalars[0], &GlobalSum[0], nFlds, MPI_INT64_T, MPI_SUM,
                        Comm);
}

// R4 scalars
int globalSum(const std::vector<R4> scalars, const MPI_Comm Comm,
              std::vector<R4> GlobalSum) {
   int nFlds = scalars.size();
   R8 LocalTmp[nFlds], GlobalTmp[nFlds];
   int i, ierr;
   for (i = 0; i < nFlds; i++) {
      LocalTmp[i] = scalars[i]; // R8<-R4
   }
   ierr = MPI_Allreduce(LocalTmp, GlobalTmp, nFlds, MPI_DOUBLE, MPI_SUM, Comm);
   for (i = 0; i < nFlds; i++) {
      GlobalSum[i] = GlobalTmp[i]; // R4<-R8
   }
   return ierr;
}

// R8 scalars
int globalSum(const std::vector<R8> scalars, const MPI_Comm Comm,
              std::vector<R8> GlobalSum) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int nFlds = scalars.size();
   complex<double> LocalTmp[nFlds], GlobalTmp[nFlds];
   int i, ierr;
   for (i = 0; i < nFlds; i++) {
      LocalTmp[i]  = complex<double>(scalars[i], 0.0);
      GlobalTmp[i] = complex<double>(0.0, 0.0);
   }
   ierr = MPI_Allreduce(LocalTmp, GlobalTmp, nFlds, MPI_C_DOUBLE_COMPLEX,
                        MPI_SUMDD, Comm);
   for (i = 0; i < nFlds; i++) {
      GlobalSum[i] = real(GlobalTmp[i]);
   }
   return ierr;
}

// I4 arrays
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I4, typename Kokkos::View<T>::value_type>, int>
globalSum(const std::vector<Kokkos::View<T, ML, MS>> arrays,
          const MPI_Comm Comm, std::vector<I4> GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   int dim   = arrays[0].rank;
   I4 LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = 0;
      for (i = imin; i < imax; i++) {
         LocalSum[ifld] += arrays[ifld].data()[i];
      }
   }
   return MPI_Allreduce(LocalSum, &GlobalSum[0], nFlds, MPI_INT32_T, MPI_SUM,
                        Comm);
}

// I8 arrays
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I8, typename Kokkos::View<T>::value_type>, int>
globalSum(const std::vector<Kokkos::View<T, ML, MS>> arrays,
          const MPI_Comm Comm, std::vector<I8> GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   int dim   = arrays[0].rank;
   I8 LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = 0;
      for (i = imin; i < imax; i++) {
         LocalSum[ifld] += arrays[ifld].data()[i];
      }
   }
   return MPI_Allreduce(LocalSum, &GlobalSum[0], nFlds, MPI_INT64_T, MPI_SUM,
                        Comm);
}

// R4 arrays
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, int>
globalSum(const std::vector<Kokkos::View<T, ML, MS>> arrays,
          const MPI_Comm Comm, std::vector<R4> GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld, ierr;
   int nFlds = arrays.size();
   int dim   = arrays[0].rank;
   R8 GlobalTmp[nFlds], LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = 0.0;
      for (i = imin; i < imax; i++) {
         LocalSum[ifld] += arrays[ifld].data()[i];
      }
   }
   ierr = MPI_Allreduce(LocalSum, GlobalTmp, nFlds, MPI_DOUBLE, MPI_SUM, Comm);
   for (ifld = 0; ifld < nFlds; ifld++) {
      GlobalSum[ifld] = GlobalTmp[ifld]; // R4<-R8
   }
   return ierr;
}

// R8 arrays
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, int>
globalSum(const std::vector<Kokkos::View<T, ML, MS>> arrays,
          const MPI_Comm Comm, std::vector<R8> GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int i, imin, imax, ifld, ierr;
   int nFlds = arrays.size();
   int dim   = arrays[0].rank;
   complex<double> GlobalTmp[nFlds], LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   double e, t1, t2, ai;
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = complex<double>(0.0, 0.0);
      for (i = imin; i < imax; i++) {
         ai = arrays[ifld].data()[i];
         t1 = ai + real(LocalSum[ifld]);
         e  = t1 - ai;
         t2 = ((real(LocalSum[ifld]) - e) + (ai - (t1 - e))) +
              imag(LocalSum[ifld]);
         // The result is t1 + t2, after normalization.
         LocalSum[ifld] = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
      }
      GlobalTmp[ifld] = complex<double>(0.0, 0.0);
   }
   ierr = MPI_Allreduce(LocalSum, GlobalTmp, nFlds, MPI_C_DOUBLE_COMPLEX,
                        MPI_SUMDD, Comm);
   for (ifld = 0; ifld < nFlds; ifld++) {
      GlobalSum[ifld] = real(GlobalTmp[ifld]);
   }
   return ierr;
}

//////////
// Global sum multi-field with product
//////////
// I4 arrays
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I4, typename Kokkos::View<T>::value_type>, int>
globalSum(const std::vector<Kokkos::View<T, ML, MS>> arrays,
          const std::vector<Kokkos::View<T, ML, MS>> arrays2,
          const MPI_Comm Comm, std::vector<I4> GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   int dim   = arrays[0].rank;
   I4 LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = 0;
      for (i = imin; i < imax; i++) {
         LocalSum[ifld] += arrays[ifld].data()[i] * arrays2[ifld].data()[i];
      }
   }
   return MPI_Allreduce(LocalSum, &GlobalSum[0], nFlds, MPI_INT32_T, MPI_SUM,
                        Comm);
}

// I8 arrays
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I8, typename Kokkos::View<T>::value_type>, int>
globalSum(const std::vector<Kokkos::View<T, ML, MS>> arrays,
          const std::vector<Kokkos::View<T, ML, MS>> arrays2,
          const MPI_Comm Comm, std::vector<I8> GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   int dim   = arrays[0].rank;
   I8 LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = 0;
      for (i = imin; i < imax; i++) {
         LocalSum[ifld] += arrays[ifld].data()[i] * arrays2[ifld].data()[i];
      }
   }
   return MPI_Allreduce(LocalSum, &GlobalSum[0], nFlds, MPI_INT64_T, MPI_SUM,
                        Comm);
}

// R4 arrays
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, int>
globalSum(const std::vector<Kokkos::View<T, ML, MS>> arrays,
          const std::vector<Kokkos::View<T, ML, MS>> arrays2,
          const MPI_Comm Comm, std::vector<R4> GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld, ierr;
   int nFlds = arrays.size();
   int dim   = arrays[0].rank;
   R8 GlobalTmp[nFlds], LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = 0.0;
      for (i = imin; i < imax; i++) {
         LocalSum[ifld] += arrays[ifld].data()[i] * arrays2[ifld].data()[i];
      }
   }
   ierr = MPI_Allreduce(LocalSum, GlobalTmp, nFlds, MPI_DOUBLE, MPI_SUM, Comm);
   for (ifld = 0; ifld < nFlds; ifld++) {
      GlobalSum[ifld] = GlobalTmp[ifld]; // R4<-R8
   }
   return ierr;
}

// R8 arrays
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, int>
globalSum(const std::vector<Kokkos::View<T, ML, MS>> arrays,
          const std::vector<Kokkos::View<T, ML, MS>> arrays2,
          const MPI_Comm Comm, std::vector<R8> GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int i, imin, imax, ifld, ierr;
   int nFlds = arrays.size();
   int dim   = arrays[0].rank;
   complex<double> GlobalTmp[nFlds], LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   double e, t1, t2, ai;
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = complex<double>(0.0, 0.0);
      for (i = imin; i < imax; i++) {
         ai = arrays[ifld].data()[i] * arrays2[ifld].data()[i];
         t1 = ai + real(LocalSum[ifld]);
         e  = t1 - ai;
         t2 = ((real(LocalSum[ifld]) - e) + (ai - (t1 - e))) +
              imag(LocalSum[ifld]);
         // The result is t1 + t2, after normalization.
         LocalSum[ifld] = complex<double>(t1 + t2, t2 - ((t1 + t2) - t1));
      }
      GlobalTmp[ifld] = complex<double>(0.0, 0.0);
   }
   ierr = MPI_Allreduce(LocalSum, GlobalTmp, nFlds, MPI_C_DOUBLE_COMPLEX,
                        MPI_SUMDD, Comm);
   for (ifld = 0; ifld < nFlds; ifld++) {
      GlobalSum[ifld] = real(GlobalTmp[ifld]);
   }
   return ierr;
}

//////////
// Global minval
//////////
// Array
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, int>
globalMinVal(const Kokkos::View<T, ML, MS> arr, const MPI_Comm Comm,
             IT GlobalMinVal, const std::vector<I4> *IndxRange = nullptr) {
   int dim = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   IT LocalMinVal = arr.data()[imin];
   for (i = imin + 1; i < imax; i++) {
      if (LocalMinVal > arr.data()[i]) {
         LocalMinVal = arr.data()[i];
      }
   }

   if (typeid(IT) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_INT32_T, MPI_MIN,
                           Comm);
   } else if (typeid(IT) == typeid(I8)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_INT64_T, MPI_MIN,
                           Comm);
   } else if (typeid(IT) == typeid(R4)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_FLOAT, MPI_MIN,
                           Comm);
   } else if (typeid(IT) == typeid(R8)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_DOUBLE, MPI_MIN,
                           Comm);
   }
   return ierr;
}

// Array with mask
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, int>
globalMinVal(const Kokkos::View<T, ML, MS> arr,
             const Kokkos::View<T, ML, MS> arr2, const MPI_Comm Comm,
             IT GlobalMinVal, const std::vector<I4> *IndxRange = nullptr) {
   int dim = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   IT tmp, LocalMinVal = arr.data()[imin] * arr2.data()[imin];
   for (i = imin + 1; i < imax; i++) {
      tmp = arr.data()[i] * arr2.data()[i];
      if (LocalMinVal > tmp) {
         LocalMinVal = tmp;
      }
   }

   if (typeid(IT) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_INT32_T, MPI_MIN,
                           Comm);
   } else if (typeid(IT) == typeid(I8)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_INT64_T, MPI_MIN,
                           Comm);
   } else if (typeid(IT) == typeid(R4)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_FLOAT, MPI_MIN,
                           Comm);
   } else if (typeid(IT) == typeid(R8)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_DOUBLE, MPI_MIN,
                           Comm);
   }
   return ierr;
}

// Array multi-field
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, int>
globalMinVal(const std::vector<Kokkos::View<T, ML, MS>> arrays,
             const MPI_Comm Comm, std::vector<IT> GlobalMinVal,
             const std::vector<I4> *IndxRange = nullptr) {
   int dim = arrays[0].rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   int ifld, nFlds = arrays.size();
   IT LocalMinVal[nFlds];
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalMinVal[ifld] = arrays[ifld].data()[imin];
      for (i = imin + 1; i < imax; i++) {
         if (LocalMinVal[ifld] > arrays[ifld].data()[i]) {
            LocalMinVal[ifld] = arrays[ifld].data()[i];
         }
      }
   }
   if (typeid(IT) == typeid(I4)) {
      ierr = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], nFlds, MPI_INT32_T,
                           MPI_MIN, Comm);
   } else if (typeid(IT) == typeid(I8)) {
      ierr = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], nFlds, MPI_INT64_T,
                           MPI_MIN, Comm);
   } else if (typeid(IT) == typeid(R4)) {
      ierr = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], nFlds, MPI_FLOAT,
                           MPI_MIN, Comm);
   } else if (typeid(IT) == typeid(R8)) {
      ierr = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], nFlds, MPI_DOUBLE,
                           MPI_MIN, Comm);
   }
   return ierr;
}

//////////
// Global maxval
//////////
// Array
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, int>
globalMaxVal(const Kokkos::View<T, ML, MS> arr, const MPI_Comm Comm,
             IT GlobalMaxVal, const std::vector<I4> *IndxRange = nullptr) {
   int dim = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   IT LocalMaxVal = arr.data()[imin];
   for (i = imin + 1; i < imax; i++) {
      if (LocalMaxVal < arr.data()[i]) {
         LocalMaxVal = arr.data()[i];
      }
   }

   if (typeid(IT) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_INT32_T, MPI_MAX,
                           Comm);
   } else if (typeid(IT) == typeid(I8)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_INT64_T, MPI_MAX,
                           Comm);
   } else if (typeid(IT) == typeid(R4)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_FLOAT, MPI_MAX,
                           Comm);
   } else if (typeid(IT) == typeid(R8)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_DOUBLE, MPI_MAX,
                           Comm);
   }
   return ierr;
}

// Array with mask
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, int>
globalMaxVal(const Kokkos::View<T, ML, MS> arr,
             const Kokkos::View<T, ML, MS> arr2, const MPI_Comm Comm,
             IT GlobalMaxVal, const std::vector<I4> *IndxRange = nullptr) {
   int dim = arr.rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   IT tmp, LocalMaxVal = arr.data()[imin] * arr2.data()[imin];
   for (i = imin + 1; i < imax; i++) {
      tmp = arr.data()[i] * arr2.data()[i];
      if (LocalMaxVal < tmp) {
         LocalMaxVal = tmp;
      }
   }

   if (typeid(IT) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_INT32_T, MPI_MAX,
                           Comm);
   } else if (typeid(IT) == typeid(I8)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_INT64_T, MPI_MAX,
                           Comm);
   } else if (typeid(IT) == typeid(R4)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_FLOAT, MPI_MAX,
                           Comm);
   } else if (typeid(IT) == typeid(R8)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_DOUBLE, MPI_MAX,
                           Comm);
   }
   return ierr;
}

// Array multi-field
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, int>
globalMaxVal(const std::vector<Kokkos::View<T, ML, MS>> arrays,
             const MPI_Comm Comm, std::vector<IT> GlobalMaxVal,
             const std::vector<I4> *IndxRange = nullptr) {
   int dim = arrays[0].rank;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   int ifld, nFlds = arrays.size();
   IT LocalMaxVal[nFlds];
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalMaxVal[ifld] = arrays[ifld].data()[imin];
      for (i = imin + 1; i < imax; i++) {
         if (LocalMaxVal[ifld] < arrays[ifld].data()[i]) {
            LocalMaxVal[ifld] = arrays[ifld].data()[i];
         }
      }
   }
   if (typeid(IT) == typeid(I4)) {
      ierr = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], nFlds, MPI_INT32_T,
                           MPI_MAX, Comm);
   } else if (typeid(IT) == typeid(I8)) {
      ierr = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], nFlds, MPI_INT64_T,
                           MPI_MAX, Comm);
   } else if (typeid(IT) == typeid(R4)) {
      ierr = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], nFlds, MPI_FLOAT,
                           MPI_MAX, Comm);
   } else if (typeid(IT) == typeid(R8)) {
      ierr = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], nFlds, MPI_DOUBLE,
                           MPI_MAX, Comm);
   }
   return ierr;
}

///-----------------------------------------------------------------------------
/// Get MIN-value across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------
int globalMin(const I4 *Val, I4 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_MIN, Comm);
}

int globalMin(const I8 *Val, I8 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_MIN, Comm);
}

int globalMin(const R4 *Val, R4 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_MIN, Comm);
}

int globalMin(const R8 *Val, R8 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_MIN, Comm);
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I4, typename Kokkos::View<T>::value_type>, int>
globalMin(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT32_T, MPI_MIN,
                        Comm);
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I8, typename Kokkos::View<T>::value_type>, int>
globalMin(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT64_T, MPI_MIN,
                        Comm);
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, int>
globalMin(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_FLOAT, MPI_MIN,
                        Comm);
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, int>
globalMin(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_DOUBLE, MPI_MIN,
                        Comm);
}

///-----------------------------------------------------------------------------
/// Get MAX-value across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------
int globalMax(const I4 *Val, I4 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_MAX, Comm);
}

int globalMax(const I8 *Val, I8 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_MAX, Comm);
}

int globalMax(const R4 *Val, R4 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_MAX, Comm);
}

int globalMax(const R8 *Val, R8 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_MAX, Comm);
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I4, typename Kokkos::View<T>::value_type>, int>
globalMax(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT32_T, MPI_MAX,
                        Comm);
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I8, typename Kokkos::View<T>::value_type>, int>
globalMax(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT64_T, MPI_MAX,
                        Comm);
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, int>
globalMax(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_FLOAT, MPI_MAX,
                        Comm);
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, int>
globalMax(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_DOUBLE, MPI_MAX,
                        Comm);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_REDUCTIONS_H
