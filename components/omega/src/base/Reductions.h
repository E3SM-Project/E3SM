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

namespace OMEGA {

///-----------------------------------------------------------------------------
/// Sum given values across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------

// Scalars
int GlobalSum(const MachEnv *InEnv, const I4 *Val, I4 *Res) {
  return MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_SUM, InEnv->getComm());
}

int GlobalSum(const MachEnv *InEnv, const I8 *Val, I8 *Res) {
  return MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_SUM, InEnv->getComm());
}

int GlobalSum(const MachEnv *InEnv, const R4 *Val, R4 *Res) {
  return MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_SUM, InEnv->getComm());
}

int GlobalSum(const MachEnv *InEnv, const R8 *Val, R8 *Res) {
  return MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_SUM, InEnv->getComm());
}

// Arrays: sum each array locally, then sum scalars across all MPI tasks
template <class T, int dim>
typename std::enable_if<
  std::is_integral<T>::value, // value is true if T is an integer type
  T>::type                    // return type is T
GlobalSum(const MachEnv *InEnv, yakl::Array<T, dim, yakl::memHost, yakl::styleC> const arr, T *val) {
  T LocalSum = yakl::intrinsics::sum(arr);
  return MPI_Allreduce(&LocalSum, val, 1, MPI_INT64_T, MPI_SUM, InEnv->getComm());
}

template <int dim>
int GlobalSum(const MachEnv *InEnv, yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const arr, R4 *val) {
  R4 LocalSum = yakl::intrinsics::sum(arr);
  return MPI_Allreduce(&LocalSum, val, 1, MPI_FLOAT, MPI_SUM, InEnv->getComm());
}

template <int dim>
int GlobalSum(const MachEnv *InEnv, yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const arr, R8 *val) {
  R8 LocalSum = yakl::intrinsics::sum(arr);
  return MPI_Allreduce(&LocalSum, val, 1, MPI_DOUBLE, MPI_SUM, InEnv->getComm());
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
              yakl::Array<I4, dim, yakl::memHost, yakl::styleC>      out) {
  int i, err, sz = in.size();
  I4 snd[sz], rcv[sz];
  for (i = 0; i < sz; i++) snd[i] = in(i);
  err = MPI_Allreduce(&snd, &rcv, sz, MPI_INT32_T, MPI_MIN, InEnv->getComm());
  for (i = 0; i < sz; i++) out(i) = rcv[i];
  return err;
}

template <int dim>
int GlobalMin(const MachEnv *InEnv,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC>      out) {
  int i, err, sz = in.size();
  I8 snd[sz], rcv[sz];
  for (i = 0; i < sz; i++) snd[i] = in(i);
  err = MPI_Allreduce(&snd, &rcv, sz, MPI_INT64_T, MPI_MIN, InEnv->getComm());
  for (i = 0; i < sz; i++) out(i) = rcv[i];
  return err;
}

template <int dim>
int GlobalMin(const MachEnv *InEnv,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC>      out) {
  int i, err, sz = in.size();
  R4 snd[sz], rcv[sz];
  for (i = 0; i < sz; i++) snd[i] = in(i);
  err = MPI_Allreduce(&snd, &rcv, sz, MPI_FLOAT, MPI_MIN, InEnv->getComm());
  for (i = 0; i < sz; i++) out(i) = rcv[i];
  return err;
}

template <int dim>
int GlobalMin(const MachEnv *InEnv,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC>      out) {
  int i, err, sz = in.size();
  R4 snd[sz], rcv[sz];
  for (i = 0; i < sz; i++) snd[i] = in(i);
  err = MPI_Allreduce(&snd, &rcv, sz, MPI_DOUBLE, MPI_MIN, InEnv->getComm());
  for (i = 0; i < sz; i++) out(i) = rcv[i];
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
              yakl::Array<I4, dim, yakl::memHost, yakl::styleC>      out) {
  int i, err, sz = in.size();
  I4 snd[sz], rcv[sz];
  for (i = 0; i < sz; i++) snd[i] = in(i);
  err = MPI_Allreduce(&snd, &rcv, sz, MPI_INT32_T, MPI_MAX, InEnv->getComm());
  for (i = 0; i < sz; i++) out(i) = rcv[i];
  return err;
}

template <int dim>
int GlobalMax(const MachEnv *InEnv,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<I8, dim, yakl::memHost, yakl::styleC>      out) {
  int i, err, sz = in.size();
  I8 snd[sz], rcv[sz];
  for (i = 0; i < sz; i++) snd[i] = in(i);
  err = MPI_Allreduce(&snd, &rcv, sz, MPI_INT64_T, MPI_MAX, InEnv->getComm());
  for (i = 0; i < sz; i++) out(i) = rcv[i];
  return err;
}

template <int dim>
int GlobalMax(const MachEnv *InEnv,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R4, dim, yakl::memHost, yakl::styleC>      out) {
  int i, err, sz = in.size();
  R4 snd[sz], rcv[sz];
  for (i = 0; i < sz; i++) snd[i] = in(i);
  err = MPI_Allreduce(&snd, &rcv, sz, MPI_FLOAT, MPI_MAX, InEnv->getComm());
  for (i = 0; i < sz; i++) out(i) = rcv[i];
  return err;
}

template <int dim>
int GlobalMax(const MachEnv *InEnv,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC> const in,
              yakl::Array<R8, dim, yakl::memHost, yakl::styleC>      out) {
  int i, err, sz = in.size();
  R8 snd[sz], rcv[sz];
  for (i = 0; i < sz; i++) snd[i] = in(i);
  err = MPI_Allreduce(&snd, &rcv, sz, MPI_DOUBLE, MPI_MAX, InEnv->getComm());
  for (i = 0; i < sz; i++) out(i) = rcv[i];
  return err;
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_REDUCTIONS_H
