#ifndef OMEGA_REDUCTIONS_H
#define OMEGA_REDUCTIONS_H
//===-- base/Reductions.h - MPI reduction definitions -----------*- C++ -*-===//
//
/// \file
/// \brief Defines
///
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include <complex.h>

namespace OMEGA {

static int R8SumInitialized = 0;

static MPI_Op MPI_SUMDD; // special MPI operator for reproducible R8 sum

void ddSum(void *InBuffer, void *OutBuffer, int *Len, MPI_Datatype *DataType) {
   double _Complex *dda = (double _Complex *)InBuffer;
   double _Complex *ddb = (double _Complex *)OutBuffer;
   double e, t1, t2;
   for (int i = 0; i < *Len; i++) {
      t1 = creal(dda[i]) + creal(ddb[i]);
      e  = t1 - creal(dda[i]);
      t2 = ((creal(ddb[i]) - e) + (creal(dda[i]) - (t1 - e))) + cimag(dda[i]) +
           cimag(ddb[i]);
      ddb[i] = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
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

   double _Complex LocalTmp  = CMPLX(*Val, 0.0);
   double _Complex GlobalTmp = CMPLX(0.0, 0.0);

   int ierr = MPI_Allreduce(&LocalTmp, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                            MPI_SUMDD, Comm);
   *Res     = creal(GlobalTmp);
   return ierr;
}

//////////
// Global sum arrays
//////////
// I4 or I8 array
template <class T, int dim>
typename std::enable_if<std::is_integral<T>::value, // value is true if T is an
                                                    // integer type
                        T>::type                    // return type is T
globalSum(yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr,
          const MPI_Comm Comm, T *GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   T LocalSum = 0;
   if (IndxRange == nullptr) {
      LocalSum = yakl::intrinsics::sum(arr);
   } else {
      for (int i = (*IndxRange)[0]; i < (*IndxRange)[dim * 2 - 1]; i++) {
         LocalSum += arr.data()[i];
      }
   }
   return MPI_Allreduce(&LocalSum, GlobalSum, 1, MPI_INT64_T, MPI_SUM, Comm);
}

// R4 array
template <int dim>
int globalSum(yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const arr,
              const MPI_Comm Comm, R4 *GlobalSum,
              const std::vector<I4> *IndxRange = nullptr) {
   R8 GlobalTmp = 0.0, LocalSum = 0.0;
   if (IndxRange == nullptr) {
      LocalSum = yakl::intrinsics::sum(arr);
   } else {
      for (int i = (*IndxRange)[0]; i < (*IndxRange)[dim * 2 - 1]; i++) {
         LocalSum += arr.data()[i];
      }
   }
   int ierr =
       MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
   *GlobalSum = GlobalTmp;
   return ierr;
}

// R8 array
template <int dim>
int globalSum(yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const arr,
              const MPI_Comm Comm, R8 *GlobalSum,
              const std::vector<I4> *IndxRange = nullptr) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }

   // Accumulate the local sum using Knuth's algorithm
   double _Complex LocalSum = CMPLX(0.0, 0.0), GlobalTmp = CMPLX(0.0, 0.0);
   double e, t1, t2, ai;
   for (i = imin; i < imax; i++) {
      ai = arr.data()[i];
      t1 = ai + creal(LocalSum);
      e  = t1 - ai;
      t2 = ((creal(LocalSum) - e) + (ai - (t1 - e))) + cimag(LocalSum);
      // The result is t1 + t2, after normalization.
      LocalSum = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
   }
   ierr       = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                              MPI_SUMDD, Comm);
   *GlobalSum = creal(GlobalTmp);
   return ierr;
}

//////////
// Global sum with product
//////////
// I4 or I8 array
template <class T, int dim>
typename std::enable_if<std::is_integral<T>::value, T>::type
globalSum(yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr,
          yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr2,
          const MPI_Comm Comm, T *GlobalSum,
          const std::vector<I4> *IndxRange = nullptr) {
   T LocalSum = 0;
   int i, imin, imax;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   for (i = imin; i < imax; i++) {
      LocalSum += arr.data()[i] * arr2.data()[i];
   }
   return MPI_Allreduce(&LocalSum, GlobalSum, 1, MPI_INT64_T, MPI_SUM, Comm);
}

// R4 array
template <int dim>
int globalSum(yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const arr,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const arr2,
              const MPI_Comm Comm, R4 *GlobalSum,
              const std::vector<I4> *IndxRange = nullptr) {
   R8 GlobalTmp = 0.0, LocalSum = 0.0;
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   for (i = imin; i < imax; i++) {
      LocalSum += arr.data()[i] * arr2.data()[i];
   }
   ierr = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
   *GlobalSum = GlobalTmp;
   return ierr;
}

// R8 array
template <int dim>
int globalSum(yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const arr,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const arr2,
              const MPI_Comm Comm, R8 *GlobalSum,
              const std::vector<I4> *IndxRange = nullptr) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }

   // Accumulate the local sum using Knuth's algorithm
   double _Complex LocalSum = CMPLX(0.0, 0.0), GlobalTmp = CMPLX(0.0, 0.0);
   double e, t1, t2, ai;
   for (i = imin; i < imax; i++) {
      ai = arr.data()[i] * arr2.data()[i];
      t1 = ai + creal(LocalSum);
      e  = t1 - ai;
      t2 = ((creal(LocalSum) - e) + (ai - (t1 - e))) + cimag(LocalSum);
      // The result is t1 + t2, after normalization.
      LocalSum = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
   }
   ierr       = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                              MPI_SUMDD, Comm);
   *GlobalSum = creal(GlobalTmp);
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
   double _Complex LocalTmp[nFlds], GlobalTmp[nFlds];
   int i, ierr;
   for (i = 0; i < nFlds; i++) {
      LocalTmp[i]  = CMPLX(scalars[i], 0.0);
      GlobalTmp[i] = CMPLX(0.0, 0.0);
   }
   ierr = MPI_Allreduce(LocalTmp, GlobalTmp, nFlds, MPI_C_DOUBLE_COMPLEX,
                        MPI_SUMDD, Comm);
   for (i = 0; i < nFlds; i++) {
      GlobalSum[i] = creal(GlobalTmp[i]);
   }
   return ierr;
}

// I4 arrays
template <int dim>
int globalSum(
    const std::vector<yakl::Array<I4, dim, yakl::memHost, yakl::styleC>> arrays,
    const MPI_Comm Comm, std::vector<I4> GlobalSum,
    const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   I4 LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
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
template <int dim>
int globalSum(
    const std::vector<yakl::Array<I8, dim, yakl::memHost, yakl::styleC>> arrays,
    const MPI_Comm Comm, std::vector<I8> GlobalSum,
    const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   I8 LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
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
template <int dim>
int globalSum(
    const std::vector<yakl::Array<R4, dim, yakl::memHost, yakl::styleC>> arrays,
    const MPI_Comm Comm, std::vector<R4> GlobalSum,
    const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   R8 GlobalTmp[nFlds], LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
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
template <int dim>
int globalSum(
    const std::vector<yakl::Array<R8, dim, yakl::memHost, yakl::styleC>> arrays,
    const MPI_Comm Comm, std::vector<R8> GlobalSum,
    const std::vector<I4> *IndxRange = nullptr) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   double _Complex GlobalTmp[nFlds], LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   double e, t1, t2, ai;
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = CMPLX(0.0, 0.0);
      for (i = imin; i < imax; i++) {
         ai = arrays[ifld].data()[i];
         t1 = ai + creal(LocalSum[ifld]);
         e  = t1 - ai;
         t2 = ((creal(LocalSum[ifld]) - e) + (ai - (t1 - e))) +
              cimag(LocalSum[ifld]);
         // The result is t1 + t2, after normalization.
         LocalSum[ifld] = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
      }
      GlobalTmp[ifld] = CMPLX(0.0, 0.0);
   }
   ierr = MPI_Allreduce(LocalSum, GlobalTmp, nFlds, MPI_C_DOUBLE_COMPLEX,
                        MPI_SUMDD, Comm);
   for (ifld = 0; ifld < nFlds; ifld++) {
      GlobalSum[ifld] = creal(GlobalTmp[ifld]);
   }
   return ierr;
}

//////////
// Global sum multi-field with product
//////////
// I4 arrays
template <int dim>
int globalSum(
    const std::vector<yakl::Array<I4, dim, yakl::memHost, yakl::styleC>> arrays,
    const std::vector<yakl::Array<I4, dim, yakl::memHost, yakl::styleC>>
        arrays2,
    const MPI_Comm Comm, std::vector<I4> GlobalSum,
    const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   I4 LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
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
template <int dim>
int globalSum(
    const std::vector<yakl::Array<I8, dim, yakl::memHost, yakl::styleC>> arrays,
    const std::vector<yakl::Array<I8, dim, yakl::memHost, yakl::styleC>>
        arrays2,
    const MPI_Comm Comm, std::vector<I8> GlobalSum,
    const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   I8 LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
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
template <int dim>
int globalSum(
    const std::vector<yakl::Array<R4, dim, yakl::memHost, yakl::styleC>> arrays,
    const std::vector<yakl::Array<R4, dim, yakl::memHost, yakl::styleC>>
        arrays2,
    const MPI_Comm Comm, std::vector<R4> GlobalSum,
    const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   R8 GlobalTmp[nFlds], LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
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
template <int dim>
int globalSum(
    const std::vector<yakl::Array<R8, dim, yakl::memHost, yakl::styleC>> arrays,
    const std::vector<yakl::Array<R8, dim, yakl::memHost, yakl::styleC>>
        arrays2,
    const MPI_Comm Comm, std::vector<R8> GlobalSum,
    const std::vector<I4> *IndxRange = nullptr) {
   if (!R8SumInitialized) {
      globalSumInit();
   }
   int i, imin, imax, ifld;
   int nFlds = arrays.size();
   double _Complex GlobalTmp[nFlds], LocalSum[nFlds];
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   double e, t1, t2, ai;
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalSum[ifld] = CMPLX(0.0, 0.0);
      for (i = imin; i < imax; i++) {
         ai = arrays[ifld].data()[i] * arrays2[ifld].data()[i];
         t1 = ai + creal(LocalSum[ifld]);
         e  = t1 - ai;
         t2 = ((creal(LocalSum[ifld]) - e) + (ai - (t1 - e))) +
              cimag(LocalSum[ifld]);
         // The result is t1 + t2, after normalization.
         LocalSum[ifld] = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
      }
      GlobalTmp[ifld] = CMPLX(0.0, 0.0);
   }
   ierr = MPI_Allreduce(LocalSum, GlobalTmp, nFlds, MPI_C_DOUBLE_COMPLEX,
                        MPI_SUMDD, Comm);
   for (ifld = 0; ifld < nFlds; ifld++) {
      GlobalSum[ifld] = creal(GlobalTmp[ifld]);
   }
   return ierr;
}

//////////
// Global minval
//////////
// Array
template <class T, int dim>
int globalMinVal(yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr,
                 const MPI_Comm Comm, T GlobalMinVal,
                 const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   T LocalMinVal = arr.data()[imin];
   for (i = imin + 1; i < imax; i++) {
      if (LocalMinVal > arr.data()[i]) {
         LocalMinVal = arr.data()[i];
      }
   }

   if (typeid(T) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_INT32_T, MPI_MIN,
                           Comm);
   } else if (typeid(T) == typeid(I8)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_INT64_T, MPI_MIN,
                           Comm);
   } else if (typeid(T) == typeid(R4)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_FLOAT, MPI_MIN,
                           Comm);
   } else if (typeid(T) == typeid(R8)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_DOUBLE, MPI_MIN,
                           Comm);
   }
   return ierr;
}

// Array with mask
template <class T, int dim>
int globalMinVal(yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr,
                 yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr2,
                 const MPI_Comm Comm, T GlobalMinVal,
                 const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   T tmp, LocalMinVal = arr.data()[imin] * arr2.data()[imin];
   for (i = imin + 1; i < imax; i++) {
      tmp = arr.data()[i] * arr2.data()[i];
      if (LocalMinVal > tmp) {
         LocalMinVal = tmp;
      }
   }

   if (typeid(T) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_INT32_T, MPI_MIN,
                           Comm);
   } else if (typeid(T) == typeid(I8)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_INT64_T, MPI_MIN,
                           Comm);
   } else if (typeid(T) == typeid(R4)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_FLOAT, MPI_MIN,
                           Comm);
   } else if (typeid(T) == typeid(R8)) {
      ierr = MPI_Allreduce(&LocalMinVal, &GlobalMinVal, 1, MPI_DOUBLE, MPI_MIN,
                           Comm);
   }
   return ierr;
}

// Array multi-field
template <class T, int dim>
int globalMinVal(
    const std::vector<yakl::Array<T, dim, yakl::memHost, yakl::styleC>> arrays,
    const MPI_Comm Comm, std::vector<T> GlobalMinVal,
    const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   int ifld, nFlds = arrays.size();
   T LocalMinVal[nFlds];
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalMinVal[ifld] = arrays[ifld].data()[imin];
      for (i = imin + 1; i < imax; i++) {
         if (LocalMinVal[ifld] > arrays[ifld].data()[i]) {
            LocalMinVal[ifld] = arrays[ifld].data()[i];
         }
      }
   }
   if (typeid(T) == typeid(I4)) {
      ierr = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], nFlds, MPI_INT32_T,
                           MPI_MIN, Comm);
   } else if (typeid(T) == typeid(I8)) {
      ierr = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], nFlds, MPI_INT64_T,
                           MPI_MIN, Comm);
   } else if (typeid(T) == typeid(R4)) {
      ierr = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], nFlds, MPI_FLOAT,
                           MPI_MIN, Comm);
   } else if (typeid(T) == typeid(R8)) {
      ierr = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], nFlds, MPI_DOUBLE,
                           MPI_MIN, Comm);
   }
   return ierr;
}

//////////
// Global maxval
//////////
// Array
template <class T, int dim>
int globalMaxVal(yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr,
                 const MPI_Comm Comm, T GlobalMaxVal,
                 const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   T LocalMaxVal = arr.data()[imin];
   for (i = imin + 1; i < imax; i++) {
      if (LocalMaxVal < arr.data()[i]) {
         LocalMaxVal = arr.data()[i];
      }
   }

   if (typeid(T) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_INT32_T, MPI_MAX,
                           Comm);
   } else if (typeid(T) == typeid(I8)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_INT64_T, MPI_MAX,
                           Comm);
   } else if (typeid(T) == typeid(R4)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_FLOAT, MPI_MAX,
                           Comm);
   } else if (typeid(T) == typeid(R8)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_DOUBLE, MPI_MAX,
                           Comm);
   }
   return ierr;
}

// Array with mask
template <class T, int dim>
int globalMaxVal(yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr,
                 yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr2,
                 const MPI_Comm Comm, T GlobalMaxVal,
                 const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arr.totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   T tmp, LocalMaxVal = arr.data()[imin] * arr2.data()[imin];
   for (i = imin + 1; i < imax; i++) {
      tmp = arr.data()[i] * arr2.data()[i];
      if (LocalMaxVal < tmp) {
         LocalMaxVal = tmp;
      }
   }

   if (typeid(T) == typeid(I4)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_INT32_T, MPI_MAX,
                           Comm);
   } else if (typeid(T) == typeid(I8)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_INT64_T, MPI_MAX,
                           Comm);
   } else if (typeid(T) == typeid(R4)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_FLOAT, MPI_MAX,
                           Comm);
   } else if (typeid(T) == typeid(R8)) {
      ierr = MPI_Allreduce(&LocalMaxVal, &GlobalMaxVal, 1, MPI_DOUBLE, MPI_MAX,
                           Comm);
   }
   return ierr;
}

// Array multi-field
template <class T, int dim>
int globalMaxVal(
    const std::vector<yakl::Array<T, dim, yakl::memHost, yakl::styleC>> arrays,
    const MPI_Comm Comm, std::vector<T> GlobalMaxVal,
    const std::vector<I4> *IndxRange = nullptr) {
   int i, imin, imax, ierr;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].totElems();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   int ifld, nFlds = arrays.size();
   T LocalMaxVal[nFlds];
   for (ifld = 0; ifld < nFlds; ifld++) {
      LocalMaxVal[ifld] = arrays[ifld].data()[imin];
      for (i = imin + 1; i < imax; i++) {
         if (LocalMaxVal[ifld] < arrays[ifld].data()[i]) {
            LocalMaxVal[ifld] = arrays[ifld].data()[i];
         }
      }
   }
   if (typeid(T) == typeid(I4)) {
      ierr = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], nFlds, MPI_INT32_T,
                           MPI_MAX, Comm);
   } else if (typeid(T) == typeid(I8)) {
      ierr = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], nFlds, MPI_INT64_T,
                           MPI_MAX, Comm);
   } else if (typeid(T) == typeid(R4)) {
      ierr = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], nFlds, MPI_FLOAT,
                           MPI_MAX, Comm);
   } else if (typeid(T) == typeid(R8)) {
      ierr = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], nFlds, MPI_DOUBLE,
                           MPI_MAX, Comm);
   }
   return ierr;
}

///-----------------------------------------------------------------------------
/// Get MIN-value across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------
int GlobalMin(const I4 *Val, I4 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_MIN, Comm);
}

int GlobalMin(const I8 *Val, I8 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_MIN, Comm);
}

int GlobalMin(const R4 *Val, R4 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_MIN, Comm);
}

int GlobalMin(const R8 *Val, R8 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_MIN, Comm);
}

template <int dim>
int GlobalMin(yakl::Array<I4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I4, dim, yakl::memHost, yakl::styleC> out,
              const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT32_T, MPI_MIN,
                        Comm);
}

template <int dim>
int GlobalMin(yakl::Array<I8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC> out,
              const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT64_T, MPI_MIN,
                        Comm);
}

template <int dim>
int GlobalMin(yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> out,
              const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_FLOAT, MPI_MIN,
                        Comm);
}

template <int dim>
int GlobalMin(yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> out,
              const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_DOUBLE, MPI_MIN,
                        Comm);
}

///-----------------------------------------------------------------------------
/// Get MAX-value across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------
int GlobalMax(const I4 *Val, I4 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_MAX, Comm);
}

int GlobalMax(const I8 *Val, I8 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_MAX, Comm);
}

int GlobalMax(const R4 *Val, R4 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_MAX, Comm);
}

int GlobalMax(const R8 *Val, R8 *Res, const MPI_Comm Comm) {
   return MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_MAX, Comm);
}

template <int dim>
int GlobalMax(yakl::Array<I4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I4, dim, yakl::memHost, yakl::styleC> out,
              const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT32_T, MPI_MAX,
                        Comm);
}

template <int dim>
int GlobalMax(yakl::Array<I8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC> out,
              const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT64_T, MPI_MAX,
                        Comm);
}

template <int dim>
int GlobalMax(yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> out,
              const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_FLOAT, MPI_MAX,
                        Comm);
}

template <int dim>
int GlobalMax(yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> out,
              const MPI_Comm Comm) {
   return MPI_Allreduce(in.data(), out.data(), in.size(), MPI_DOUBLE, MPI_MAX,
                        Comm);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_REDUCTIONS_H
