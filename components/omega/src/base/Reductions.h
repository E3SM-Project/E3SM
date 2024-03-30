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
// Scalars
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
// Arrays: sum each array locally, then sum scalars across all MPI tasks
//////////
// I4 or I8 array wout indxRange
template <class T, int dim>
typename std::enable_if<std::is_integral<T>::value, // value is true if T is an
                                                    // integer type
                        T>::type                    // return type is T
globalSum(yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr,
          const MPI_Comm Comm, T *val) {
   T LocalSum = yakl::intrinsics::sum(arr);
   return MPI_Allreduce(&LocalSum, val, 1, MPI_INT64_T, MPI_SUM, Comm);
}

// I4 or I8 array with indxRange
template <class T, int dim>
typename std::enable_if<std::is_integral<T>::value, // value is true if T is an
                                                    // integer type
                        T>::type                    // return type is T
globalSum(yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr,
          const MPI_Comm Comm, const std::vector<I4> IndxRange, T *val) {
   T LocalSum = 0;
   for (int i = IndxRange[0]; i < IndxRange[1]; i++) {
      LocalSum += arr.data()[i];
   }
   return MPI_Allreduce(&LocalSum, val, 1, MPI_INT64_T, MPI_SUM, Comm);
}

// R4 array wout indxRange
template <int dim>
int globalSum(yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const arr,
              const MPI_Comm Comm, R4 *val) {
   R8 GlobalSum = 0.0, LocalSum = yakl::intrinsics::sum(arr);
   int ierr =
       MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPI_DOUBLE, MPI_SUM, Comm);
   *val = GlobalSum;
   return ierr;
}

// R4 array with indxRange
template <int dim>
int globalSum(yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const arr,
              const MPI_Comm Comm, const std::vector<I4> IndxRange, R4 *val) {
   R8 GlobalSum = 0.0, LocalSum = 0.0;
   for (int i = IndxRange[0]; i < IndxRange[1]; i++) {
      LocalSum += arr.data()[i];
   }
   int ierr =
       MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPI_DOUBLE, MPI_SUM, Comm);
   *val = GlobalSum;
   return ierr;
}

// R8 array wout indxRange
template <int dim>
int globalSum(yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const arr,
              const MPI_Comm Comm, R8 *val) {
   // initialize reproducible MPI_SUMDD operator
   if (!R8SumInitialized) {
      globalSumInit();
   }

   // Accumulate the local sum using Knuth's algorithm
   double _Complex LocalSum = CMPLX(0.0, 0.0);
   double e, t1, t2, ai;
   for (int i = 0; i < arr.totElems(); i++) {
      ai = arr.data()[i];
      t1 = ai + creal(LocalSum);
      e  = t1 - ai;
      t2 = ((creal(LocalSum) - e) + (ai - (t1 - e))) + cimag(LocalSum);
      // The result is t1 + t2, after normalization.
      LocalSum = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
   }
   double _Complex GlobalSum = CMPLX(0.0, 0.0);
   int ierr = MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPI_C_DOUBLE_COMPLEX,
                            MPI_SUMDD, Comm);
   *val     = creal(GlobalSum);
   return ierr;
}

// R8 array with indxRange
template <int dim>
int globalSum(yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const arr,
              const MPI_Comm Comm, const std::vector<I4> IndxRange, R8 *val) {
   // initialize reproducible MPI_SUMDD operator
   if (!R8SumInitialized) {
      globalSumInit();
   }

   // Accumulate the local sum using Knuth's algorithm
   double _Complex LocalSum = CMPLX(0.0, 0.0);
   double e, t1, t2, ai;
   for (int i = IndxRange[0]; i < IndxRange[1]; i++) {
      ai = arr.data()[i];
      t1 = ai + creal(LocalSum);
      e  = t1 - ai;
      t2 = ((creal(LocalSum) - e) + (ai - (t1 - e))) + cimag(LocalSum);
      // The result is t1 + t2, after normalization.
      LocalSum = CMPLX(t1 + t2, t2 - ((t1 + t2) - t1));
   }
   double _Complex GlobalSum = CMPLX(0.0, 0.0);
   int ierr = MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPI_C_DOUBLE_COMPLEX,
                            MPI_SUMDD, Comm);
   *val     = creal(GlobalSum);
   return ierr;
}

///-----------------------------------------------------------------------------
/// Get MIN-value across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------
int GlobalMin(const MachEnv *InEnv, const I4 *Val, I4 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_MIN, InEnv->getComm());
}

int GlobalMin(const MachEnv *InEnv, const I8 *Val, I8 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_MIN, InEnv->getComm());
}

int GlobalMin(const MachEnv *InEnv, const R4 *Val, R4 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_MIN, InEnv->getComm());
}

int GlobalMin(const MachEnv *InEnv, const R8 *Val, R8 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_MIN, InEnv->getComm());
}

template <int dim>
int GlobalMin(const MachEnv *InEnv,
              yakl::Array<I4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I4, dim, yakl::memHost, yakl::styleC> out) {
   int i, err, sz = in.size();
   I4 snd[sz], rcv[sz];
   for (i = 0; i < sz; i++)
      snd[i] = in(i);
   err = MPI_Allreduce(&snd, &rcv, sz, MPI_INT32_T, MPI_MIN, InEnv->getComm());
   for (i = 0; i < sz; i++)
      out(i) = rcv[i];
   return err;
}

template <int dim>
int GlobalMin(const MachEnv *InEnv,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC> out) {
   int i, err, sz = in.size();
   I8 snd[sz], rcv[sz];
   for (i = 0; i < sz; i++)
      snd[i] = in(i);
   err = MPI_Allreduce(&snd, &rcv, sz, MPI_INT64_T, MPI_MIN, InEnv->getComm());
   for (i = 0; i < sz; i++)
      out(i) = rcv[i];
   return err;
}

template <int dim>
int GlobalMin(const MachEnv *InEnv,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> out) {
   int i, err, sz = in.size();
   R4 snd[sz], rcv[sz];
   for (i = 0; i < sz; i++)
      snd[i] = in(i);
   err = MPI_Allreduce(&snd, &rcv, sz, MPI_FLOAT, MPI_MIN, InEnv->getComm());
   for (i = 0; i < sz; i++)
      out(i) = rcv[i];
   return err;
}

template <int dim>
int GlobalMin(const MachEnv *InEnv,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> out) {
   int i, err, sz = in.size();
   R4 snd[sz], rcv[sz];
   for (i = 0; i < sz; i++)
      snd[i] = in(i);
   err = MPI_Allreduce(&snd, &rcv, sz, MPI_DOUBLE, MPI_MIN, InEnv->getComm());
   for (i = 0; i < sz; i++)
      out(i) = rcv[i];
   return err;
}

///-----------------------------------------------------------------------------
/// Get MAX-value across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------
int GlobalMax(const MachEnv *InEnv, const I4 *Val, I4 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_MAX, InEnv->getComm());
}

int GlobalMax(const MachEnv *InEnv, const I8 *Val, I8 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_MAX, InEnv->getComm());
}

int GlobalMax(const MachEnv *InEnv, const R4 *Val, R4 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_MAX, InEnv->getComm());
}

int GlobalMax(const MachEnv *InEnv, const R8 *Val, R8 *Res) {
   return MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_MAX, InEnv->getComm());
}

template <int dim>
int GlobalMax(const MachEnv *InEnv,
              yakl::Array<I4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I4, dim, yakl::memHost, yakl::styleC> out) {
   int i, err, sz = in.size();
   I4 snd[sz], rcv[sz];
   for (i = 0; i < sz; i++)
      snd[i] = in(i);
   err = MPI_Allreduce(&snd, &rcv, sz, MPI_INT32_T, MPI_MAX, InEnv->getComm());
   for (i = 0; i < sz; i++)
      out(i) = rcv[i];
   return err;
}

template <int dim>
int GlobalMax(const MachEnv *InEnv,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC> out) {
   int i, err, sz = in.size();
   I8 snd[sz], rcv[sz];
   for (i = 0; i < sz; i++)
      snd[i] = in(i);
   err = MPI_Allreduce(&snd, &rcv, sz, MPI_INT64_T, MPI_MAX, InEnv->getComm());
   for (i = 0; i < sz; i++)
      out(i) = rcv[i];
   return err;
}

template <int dim>
int GlobalMax(const MachEnv *InEnv,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> out) {
   int i, err, sz = in.size();
   R4 snd[sz], rcv[sz];
   for (i = 0; i < sz; i++)
      snd[i] = in(i);
   err = MPI_Allreduce(&snd, &rcv, sz, MPI_FLOAT, MPI_MAX, InEnv->getComm());
   for (i = 0; i < sz; i++)
      out(i) = rcv[i];
   return err;
}

template <int dim>
int GlobalMax(const MachEnv *InEnv,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> out) {
   int i, err, sz = in.size();
   R8 snd[sz], rcv[sz];
   for (i = 0; i < sz; i++)
      snd[i] = in(i);
   err = MPI_Allreduce(&snd, &rcv, sz, MPI_DOUBLE, MPI_MAX, InEnv->getComm());
   for (i = 0; i < sz; i++)
      out(i) = rcv[i];
   return err;
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_REDUCTIONS_H
