#ifndef OMEGA_REDUCTIONS_H
#define OMEGA_REDUCTIONS_H
//===-- base/Reductions.h - MPI reduction definitions -----------*- C++ -*-===//
//
/// \file
/// \brief Functions for computing global reductions of distributed values
///
/// Global reductions compute global sums, minimum value or maximum value
/// of distributed data across MPI tasks. All sums are reproducible with
/// respect to the order in which the values are computed (eg due to different
/// numbers of tasks or different decompositions of data). For real (single
/// precision) values, reproducibility is obtained by computing the sums in
/// double precision and truncating the result. For double precision values,
/// the sums are computed using the double-double algorithms of Knuth and
/// Bailey following the implementation of He and Ding. See:
///   Knuth, DE., 1969, The Art of Computer Programming, chapter 4, vol 2.
///   Yozo Hida and S Xiaoye and David H Bailey, 2008, Library for double-double
///      and quad-double arithmetic,
///      https://www.davidhbailey.com/dhbpapers/qd.pdf
///   He, Yun and Chris Ding, 2001, Using Accurate Arithmetics to Improve
///      Numerical Reproducibility and Stability in Parallel Applications,
///      J. Supercomputing, vol 18, pp 259-277
///
/// Array functions are templated on memory local (host/device) because these
/// are the same on CPU-only machines. However, the sums are calculated
/// differently for each data type and we cannot assume linear storage for
/// Kokkos arrays, so there are separate function definitions for type and
/// array dimension.
//
//===----------------------------------------------------------------------===//

#include <complex>
using std::complex;

#include "DataTypes.h"
#include "Error.h"
#include "OmegaKokkos.h"

namespace OMEGA {

static bool R8SumNotInitialized = true;

static MPI_Op MPI_SUMDD; // special MPI operator for reproducible R8 sum

//------------------------------------------------------------------------------
// Special sum function to use in custom double precision MPI reduction.
// Implements the DD algorithm for reproducible sums as implemented in
//
void sumDD(void *InBuffer, void *OutBuffer, int *Len, MPI_Datatype *DataType) {
   complex<double> *DDa = (complex<double> *)InBuffer;
   complex<double> *DDb = (complex<double> *)OutBuffer;
   double E, T1, T2;
   for (int I = 0; I < *Len; I++) {
      T1 = real(DDa[I]) + real(DDb[I]);
      E  = T1 - real(DDa[I]);
      T2 = ((real(DDb[I]) - E) + (real(DDa[I]) - (T1 - E))) + imag(DDa[I]) +
           imag(DDb[I]);
      DDb[I] = complex<double>(T1 + T2, T2 - ((T1 + T2) - T1));
   }
}

// Local iteration of a DD sum for use in local sums
void sumDDLocal(complex<double> &DDb, //< [inout] local sum and residual
                const double &DDa     //< [in] local double to add to sum
) {
   double T1 = DDa + real(DDb);
   double E  = T1 - DDa;
   double T2 = ((real(DDb) - E) + (DDa - (T1 - E))) + imag(DDb);
   DDb       = complex<double>(T1 + T2, T2 - ((T1 + T2) - T1));
}

//------------------------------------------------------------------------------
// Special functions to extract array properties and index range
/// Determines whether a Kokkos array is on host
template <typename T, typename ML, typename MS>
bool isReduceArrayOnHost(
    const Kokkos::View<T, ML, MS> Array ///< [in] array to extract info
) {
   // Determine where array is located
   return Kokkos::SpaceAccessibility<MS, Kokkos::HostSpace>::accessible;
}

/// Extracts range and strides from Kokkos array
template <typename T, typename ML, typename MS>
void getReduceArrayInfo(
    std::vector<I4> &IRange,              ///< [out] index range for reduction
    std::vector<I8> &Strides,             ///< [out] array stride for each dim
    const Kokkos::View<T, ML, MS> &Array, ///< [in] array to extract info
    const std::vector<I4> *IndxRange      ///< [in] input index range
) {
   // Get array and layout information
   int ArrDim = Array.rank;
   OMEGA_REQUIRE((ArrDim > 0) and (ArrDim < 6),
                 "Reductions: Array with {} dimensions not supported", ArrDim);
   // Kokkos may pad array dims so extract the actual stride used for storage
   // For some reason, need to use plain array for the stride function rather
   // than the data pointer from std::vector and copy to vector later...
   size_t TmpStrides[ArrDim + 1];
   Array.stride(TmpStrides);
   // Strides returns non-zero stride for last dims so only set up to Dim - 1
   for (int IDim = 0; IDim < ArrDim; ++IDim)
      Strides[IDim] = TmpStrides[IDim];
   // Set index range
   std::fill(IRange.begin(), IRange.end(), 0); // init to zero
   if (IndxRange == nullptr) {
      if (ArrDim > 0)
         IRange[1] = Array.extent(0) - 1;
      if (ArrDim > 1)
         IRange[3] = Array.extent(1) - 1;
      if (ArrDim > 2)
         IRange[5] = Array.extent(2) - 1;
      if (ArrDim > 3)
         IRange[7] = Array.extent(3) - 1;
      if (ArrDim > 4)
         IRange[9] = Array.extent(4) - 1;
   } else {
      if (ArrDim > 0) {
         IRange[0] = (*IndxRange)[0];
         IRange[1] = (*IndxRange)[1];
      }
      if (ArrDim > 1) {
         IRange[2] = (*IndxRange)[2];
         IRange[3] = (*IndxRange)[3];
      }
      if (ArrDim > 2) {
         IRange[4] = (*IndxRange)[4];
         IRange[5] = (*IndxRange)[5];
      }
      if (ArrDim > 3) {
         IRange[6] = (*IndxRange)[6];
         IRange[7] = (*IndxRange)[7];
      }
      if (ArrDim > 4) {
         IRange[8] = (*IndxRange)[8];
         IRange[9] = (*IndxRange)[9];
      }
   }
}

/// Copies some array index info to device
void copyReduceInfoToDevice(
    Array1DI4 &DevRange,     ///< [out] index range for reduction
    Array1DI8 &DevStrides,   ///< [out] array stride for each dim
    std::vector<I4> &IRange, ///< [in] index range for reduction
    std::vector<I8> &Strides ///< [in] array stride for each dim
) {
   HostArray1DI4 HostRange("IRange", 10);
   HostArray1DI8 HostStrides("Strides", 5);
   for (int I = 0; I < 5; ++I) {
      HostStrides(I)       = Strides[I];
      HostRange(2 * I)     = IRange[2 * I];
      HostRange(2 * I + 1) = IRange[2 * I + 1];
   }
   deepCopy(DevRange, HostRange);
   deepCopy(DevStrides, HostStrides);
}

//------------------------------------------------------------------------------
// Initialize the special DD sum operator for double precision reproducible sums
void globalSumInit() {
   int Err = MPI_Op_create(&sumDD, 1, &MPI_SUMDD);
   if (Err == 0)
      R8SumNotInitialized = false;
}

///-----------------------------------------------------------------------------
/// Scalar sum - sums local scalar values across all MPI processors associated
/// with the communicator.  Return value is sum across all tasks and sums are
/// bit reproducible.
/// I4 specific interface
I4 globalSum(const I4 &Val,      ///< [in] local scalar value to be summed
             const MPI_Comm Comm ///< [in] MPI communicator
) {
   I4 Result;
   int Err = MPI_Allreduce(&Val, &Result, 1, MPI_INT32_T, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (I4 scalar): Error in MPI_Allreduce");
   return Result;
}

/// I8 specific interface
I8 globalSum(const I8 &Val,      ///< [in] local scalar value to be summed
             const MPI_Comm Comm ///< [in] MPI communicator
) {
   I8 Result;
   int Err = MPI_Allreduce(&Val, &Result, 1, MPI_INT64_T, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (I8 scalar): Error in MPI_Allreduce");
   return Result;
}

/// R4 specific interface
R4 globalSum(const R4 &Val,      ///< [in] local scalar value to be summed
             const MPI_Comm Comm ///< [in] MPI communicator
) {
   R8 LocalTmp, GlobalTmp;
   LocalTmp = Val; // convert to double for reproducibility
   int Err = MPI_Allreduce(&LocalTmp, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (R4 scalar): Error in MPI_Allreduce");
   R4 Result = GlobalTmp; // convert back to R4
   return Result;
}

/// R8 specific interface
R8 globalSum(const R8 &Val,      ///< [in] local scalar value to be summed
             const MPI_Comm Comm ///< [in] MPI communicator
) {
   // initialize reproducible MPI_SUMDD operator
   if (R8SumNotInitialized)
      globalSumInit();

   // For reproducibility, store sum and residual in a complex number pair
   // and use special MPI sum operator defined at init
   complex<double> LocalTmp(Val, 0.0);
   complex<double> GlobalTmp(0.0, 0.0);

   int Err = MPI_Allreduce(&LocalTmp, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                           MPI_SUMDD, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (R8 scalar): Error in MPI_Allreduce");
   // Final sum is the real part of the complex pair
   R8 Result = real(GlobalTmp);
   return Result;
}

///-----------------------------------------------------------------------------
/// Multifield scalar sums - sums a set of scalar values across MPI tasks
/// Return value is a vector containg the sum of each scalar across all tasks
/// I4 specific interface
std::vector<I4> globalSum(
    const std::vector<I4> &Scalars, ///< [in] vector of scalars to be summed
    const MPI_Comm Comm             ///< [in] MPI communicator
) {
   int NFields = Scalars.size();
   std::vector<I4> Results(NFields);
   int Err = MPI_Allreduce(&Scalars[0], &Results[0], NFields, MPI_INT32_T,
                           MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (I4 scalar multifield): Error in MPI_Allreduce");
   return Results;
}

/// I8 specific interface
std::vector<I8> globalSum(
    const std::vector<I8> &Scalars, ///< [in] vector of scalars to be summed
    const MPI_Comm Comm             ///< [in] MPI communicator
) {
   int NFields = Scalars.size();
   std::vector<I8> Results(NFields);
   int Err = MPI_Allreduce(&Scalars[0], &Results[0], NFields, MPI_INT64_T,
                           MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (I8 scalar multifield): Error in MPI_Allreduce");
   return Results;
}

/// R4 specific interface
std::vector<R4> globalSum(
    const std::vector<R4> &Scalars, ///< [in] vector of scalars to be summed
    const MPI_Comm Comm             ///< [in] MPI communicator
) {
   int NFields = Scalars.size();
   // For reproducibility, perform sum in double precision
   std::vector<R8> LocalTmp(NFields);
   std::vector<R8> GlobalTmp(NFields);
   for (int i = 0; i < NFields; i++) {
      LocalTmp[i] = Scalars[i]; // R8<-R4
   }
   int Err = MPI_Allreduce(&LocalTmp[0], &GlobalTmp[0], NFields, MPI_DOUBLE,
                           MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (R4 scalar multifield): Error in MPI_Allreduce");
   std::vector<R4> Results(NFields);
   for (int i = 0; i < NFields; i++) {
      Results[i] = GlobalTmp[i]; // R4<-R8
   }
   return Results;
}

/// R8 specific interface
std::vector<R8> globalSum(
    const std::vector<R8> &Scalars, ///< [in] vector of scalars to be summed
    const MPI_Comm Comm             ///< [in] MPI communicator
) {
   if (R8SumNotInitialized)
      globalSumInit();

   int NFields = Scalars.size();
   // The DD algorithm uses a complex number to store sums and residuals
   complex<double> LocalTmp[NFields], GlobalTmp[NFields];
   for (int i = 0; i < NFields; i++) {
      LocalTmp[i]  = complex<double>(Scalars[i], 0.0);
      GlobalTmp[i] = complex<double>(0.0, 0.0);
   }
   int Err = MPI_Allreduce(LocalTmp, GlobalTmp, NFields, MPI_C_DOUBLE_COMPLEX,
                           MPI_SUMDD, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (R8 scalar multifield): Error in MPI_Allreduce");
   std::vector<R8> Results(NFields);
   for (int i = 0; i < NFields; i++) {
      Results[i] = real(GlobalTmp[i]); // result is in real part of complex
   }
   return Results;
}

///-----------------------------------------------------------------------------
/// Computes sums of distributed arrays over the local indices and across MPI
/// tasks. The optional IndxRange argument is a vector that supplies the
/// min and max local indices along each dimension over which the sum is
/// computed. If not provided, the sum is computed over the full array.
/// The final scalar sum is returned as the result. The sums are bit-by-bit
/// reproducible.
/// I4 specific interface
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I4, typename Kokkos::View<T>::value_type>, I4>
globalSum(const Kokkos::View<T, ML, MS> Array, ///< [in] array to be summed
          const MPI_Comm Comm,                 ///< [in] MPI Communicator
          const std::vector<I4> *IndxRange = nullptr ///< [in] index range
) {
   // Get array and layout information
   bool IsHost = isReduceArrayOnHost(Array);
   std::vector<I8> Strides(5, 0);
   std::vector<I4> IRange(10, 0);
   getReduceArrayInfo(IRange, Strides, Array, IndxRange);

   // Compute local sum on host or device
   I4 LocalSum = 0;
   if (IsHost) {
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     LocalSum += Array.data()[LinearAdd];
                  }
               }
            }
         }
      }
   } else { // on device
      Array1DI4 DevRange("IRange", 10);
      Array1DI8 DevStrides("Strides", 5);
      copyReduceInfoToDevice(DevRange, DevStrides, IRange, Strides);
      OMEGA_SCOPE(LocStrides, DevStrides);
      OMEGA_SCOPE(LocRange, DevRange);
      OMEGA_SCOPE(LocArray, Array);
      parallelReduce(
          {IRange[1] + 1, IRange[3] + 1, IRange[5] + 1, IRange[7] + 1,
           IRange[9] + 1},
          KOKKOS_LAMBDA(int I, int J, int K, int L, int M, I4 &Accum) {
             if (I >= LocRange(0) and J >= LocRange(2) and K >= LocRange(4) and
                 L >= LocRange(6) and M >= LocRange(8)) {
                size_t LinearAdd = I * LocStrides(0) + J * LocStrides(1) +
                                   K * LocStrides(2) + L * LocStrides(3) +
                                   M * LocStrides(4);
                Accum += LocArray.data()[LinearAdd];
             }
          },
          LocalSum);
   } // end if onHost
   // Compute final sum by adding local sums from each MPI task
   I4 GlobalSum;
   int Err =
       MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPI_INT32_T, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (I4 Array): Error in MPI_Allreduce");
   return GlobalSum;
}
/// I8 specific interface
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I8, typename Kokkos::View<T>::value_type>, I8>
globalSum(const Kokkos::View<T, ML, MS> Array, ///< [in] array to be summed
          const MPI_Comm Comm,                 ///< [in] MPI Communicator
          const std::vector<I4> *IndxRange = nullptr ///< [in] index range
) {
   // Get array and layout information
   bool IsHost = isReduceArrayOnHost(Array);
   std::vector<I8> Strides(5, 0);
   std::vector<I4> IRange(10, 0);
   getReduceArrayInfo(IRange, Strides, Array,
                      IndxRange); ///< [out] true if host array

   // Compute local sum on host or device
   I8 LocalSum = 0;
   if (IsHost) {
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     LocalSum += Array.data()[LinearAdd];
                  }
               }
            }
         }
      }
   } else { // on device
      // Transfer some info on device
      Array1DI4 DevRange("IRange", 10);
      Array1DI8 DevStrides("Strides", 5);
      copyReduceInfoToDevice(DevRange, DevStrides, IRange, Strides);
      OMEGA_SCOPE(LocStrides, DevStrides);
      OMEGA_SCOPE(LocRange, DevRange);
      OMEGA_SCOPE(LocArray, Array);
      parallelReduce(
          {IRange[1] + 1, IRange[3] + 1, IRange[5] + 1, IRange[7] + 1,
           IRange[9] + 1},
          KOKKOS_LAMBDA(int I, int J, int K, int L, int M, I8 &Accum) {
             if (I >= LocRange(0) and J >= LocRange(2) and K >= LocRange(4) and
                 L >= LocRange(6) and M >= LocRange(8)) {
                size_t LinearAdd = I * LocStrides(0) + J * LocStrides(1) +
                                   K * LocStrides(2) + L * LocStrides(3) +
                                   M * LocStrides(4);
                Accum += LocArray.data()[LinearAdd];
             }
          },
          LocalSum);
   } // end if onHost
   // Compute final sum by adding local sums from each MPI task
   I8 GlobalSum;
   int Err =
       MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPI_INT64_T, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (I8 Array): Error in MPI_Allreduce");
   return GlobalSum;
}
/// R4 specific interface
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, R4>
globalSum(const Kokkos::View<T, ML, MS> Array, ///< [in] array to be summed
          const MPI_Comm Comm,                 ///< [in] MPI Communicator
          const std::vector<I4> *IndxRange = nullptr ///< [in] index range
) {
   // Get array and layout information
   bool IsHost = isReduceArrayOnHost(Array);
   std::vector<I8> Strides(5, 0);
   std::vector<I4> IRange(10, 0);
   getReduceArrayInfo(IRange, Strides, Array,
                      IndxRange); ///< [out] true if host array

   // Compute local sum on host or device
   R8 LocalSum = 0;
   if (IsHost) {
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     LocalSum += Array.data()[LinearAdd];
                  }
               }
            }
         }
      }
   } else { // on device
      // Transfer some info on device
      Array1DI4 DevRange("IRange", 10);
      Array1DI8 DevStrides("Strides", 5);
      copyReduceInfoToDevice(DevRange, DevStrides, IRange, Strides);
      OMEGA_SCOPE(LocStrides, DevStrides);
      OMEGA_SCOPE(LocRange, DevRange);
      OMEGA_SCOPE(LocArray, Array);
      parallelReduce(
          {IRange[1] + 1, IRange[3] + 1, IRange[5] + 1, IRange[7] + 1,
           IRange[9] + 1},
          KOKKOS_LAMBDA(int I, int J, int K, int L, int M, R8 &Accum) {
             if (I >= LocRange(0) and J >= LocRange(2) and K >= LocRange(4) and
                 L >= LocRange(6) and M >= LocRange(8)) {
                size_t LinearAdd = I * LocStrides[0] + J * LocStrides[1] +
                                   K * LocStrides[2] + L * LocStrides[3] +
                                   M * LocStrides[4];
                Accum += LocArray.data()[LinearAdd];
             }
          },
          LocalSum);
   } // end if onHost
   // Compute final sum by adding local sums from each MPI task
   R8 GlobalTmp;
   int Err = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (R4 Array): Error in MPI_Allreduce");
   R4 GlobalSum = GlobalTmp; // convert back to R4
   return GlobalSum;
}
/// R8 specific interface
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, R8>
globalSum(const Kokkos::View<T, ML, MS> Array, ///< [in] array to be summed
          const MPI_Comm Comm,                 ///< [in] MPI Communicator
          const std::vector<I4> *IndxRange = nullptr ///< [in] index range
) {
   // Get array and layout information
   bool IsHost = isReduceArrayOnHost(Array);
   std::vector<I8> Strides(5, 0);
   std::vector<I4> IRange(10, 0);
   getReduceArrayInfo(IRange, Strides, Array,
                      IndxRange); ///< [out] true if host array

   // Compute local sum on host or device
   complex<double> DDTmp; // Sum and residual for DD algorithm
   if (IsHost) {
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     R8 DataTmp = Array.data()[LinearAdd];
                     sumDDLocal(DDTmp, DataTmp);
                  }
               }
            }
         }
      }
   } else { // on device
      // Kokkos sums not reproducible so transfer data to host and compute
      // reproducible sum there
      auto ArrayH = createHostMirrorCopy(Array);
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     R8 DataTmp = ArrayH.data()[LinearAdd];
                     sumDDLocal(DDTmp, DataTmp);
                  }
               }
            }
         }
      }
      // Original Kokkos sum code
      //R8 LocalSum = 0;
      // Transfer some info on device
      //Array1DI4 DevRange("IRange", 10);
      //Array1DI8 DevStrides("Strides", 5);
      //copyReduceInfoToDevice(DevRange, DevStrides, IRange, Strides);
      //OMEGA_SCOPE(LocStrides, DevStrides);
      //OMEGA_SCOPE(LocRange, DevRange);
      //OMEGA_SCOPE(LocArray, Array);
      //parallelReduce(
      //    {IRange[1] + 1, IRange[3] + 1, IRange[5] + 1, IRange[7] + 1,
      //     IRange[9] + 1},
      //    KOKKOS_LAMBDA(int I, int J, int K, int L, int M, R8 &Accum) {
      //       if (I >= LocRange(0) and J >= LocRange(2) and K >= LocRange(4) and
      //           L >= LocRange(6) and M >= LocRange(8)) {
      //          size_t LinearAdd = I * LocStrides[0] + J * LocStrides[1] +
      //                             K * LocStrides[2] + L * LocStrides[3] +
      //                             M * LocStrides[4];
      //          Accum += LocArray.data()[LinearAdd];
      //       }
      //    },
      //    LocalSum);
      //DDTmp = complex<double>(LocalSum, 0.0);
   } // end if onHost
   // Compute final sum by adding local sums from each MPI task
   complex<double> GlobalTmp(0.0, 0.0);
   int Err = MPI_Allreduce(&DDTmp, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                           MPI_SUMDD, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (R8 Array): Error in MPI_Allreduce");
   R8 GlobalSum = real(GlobalTmp); // final sum is in real part
   return GlobalSum;
}
///-----------------------------------------------------------------------------
/// Sums a product of two arrays (eg a mask), computing the global sum of a
/// product of two distributed arrays over the local index dimensions and
/// across MPI tasks. The two arrays must be the same size and type. An optional
/// index range can be specified through a vector of min, max indices for each
/// dimension of the arrays.
/// I4 specific interface
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I4, typename Kokkos::View<T>::value_type>, I4>
globalSum(const Kokkos::View<T, ML, MS> Arr1, ///< [in] first  array in product
          const Kokkos::View<T, ML, MS> Arr2, ///< [in] second array in product
          const MPI_Comm Comm,                ///< [in] MPI Communicator
          const std::vector<I4> *IndxRange = nullptr ///< [in] opt index range
) {
   // Some error checks
   OMEGA_REQUIRE(Arr1.rank == Arr2.rank,
                 "globalSum (I4 Array product): Arrays must have same rank");
   OMEGA_REQUIRE(Arr1.size() == Arr2.size(),
                 "globalSum (I4 Array product): Arrays must have same size");
   // Get array and layout information
   bool IsHost = isReduceArrayOnHost(Arr1);
   std::vector<I8> Strides(5, 0);
   std::vector<I4> IRange(10, 0);
   getReduceArrayInfo(IRange, Strides, Arr1,
                      IndxRange); ///< [out] true if host array

   // Compute local sum on host or device
   I4 LocalSum = 0;
   if (IsHost) {
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     LocalSum +=
                         Arr1.data()[LinearAdd] * Arr2.data()[LinearAdd];
                  }
               }
            }
         }
      }
   } else { // on device
      // Transfer some info on device
      Array1DI4 DevRange("IRange", 10);
      Array1DI8 DevStrides("Strides", 5);
      copyReduceInfoToDevice(DevRange, DevStrides, IRange, Strides);
      OMEGA_SCOPE(LocStrides, DevStrides);
      OMEGA_SCOPE(LocRange, DevRange);
      OMEGA_SCOPE(LocArr1, Arr1);
      OMEGA_SCOPE(LocArr2, Arr2);
      parallelReduce(
          {IRange[1] + 1, IRange[3] + 1, IRange[5] + 1, IRange[7] + 1,
           IRange[9] + 1},
          KOKKOS_LAMBDA(int I, int J, int K, int L, int M, I4 &Accum) {
             if (I >= LocRange(0) and J >= LocRange(2) and K >= LocRange(4) and
                 L >= LocRange(6) and M >= LocRange(8)) {
                size_t LinearAdd = I * LocStrides[0] + J * LocStrides[1] +
                                   K * LocStrides[2] + L * LocStrides[3] +
                                   M * LocStrides[4];
                Accum += LocArr1.data()[LinearAdd] * LocArr2.data()[LinearAdd];
             }
          },
          LocalSum);
   } // end if onHost
   // Compute final sum by adding local sums from each MPI task
   I4 GlobalSum;
   int Err =
       MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPI_INT32_T, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (I4 Array product): Error in MPI_Allreduce");
   return GlobalSum;
}
/// I8 specific interface
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I8, typename Kokkos::View<T>::value_type>, I8>
globalSum(const Kokkos::View<T, ML, MS> Arr1, // [in] first  array in product
          const Kokkos::View<T, ML, MS> Arr2, // [in] second array in product
          const MPI_Comm Comm,                // [in] MPI Communicator
          const std::vector<I4> *IndxRange = nullptr ///< [in] opt index range
) {
   // Some error checks
   OMEGA_REQUIRE(Arr1.rank == Arr2.rank,
                 "globalSum (I8 Array product): Arrays must have same rank");
   OMEGA_REQUIRE(Arr1.size() == Arr2.size(),
                 "globalSum (I8 Array product): Arrays must have same size");
   // Get array and layout information
   bool IsHost = isReduceArrayOnHost(Arr1);
   std::vector<I8> Strides(5, 0);
   std::vector<I4> IRange(10, 0);
   getReduceArrayInfo(IRange, Strides, Arr1,
                      IndxRange); ///< [out] true if host array

   // Compute local sum on host or device
   I8 LocalSum = 0;
   if (IsHost) {
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     LocalSum +=
                         Arr1.data()[LinearAdd] * Arr2.data()[LinearAdd];
                  }
               }
            }
         }
      }
   } else { // on device
      // Transfer some info on device
      Array1DI4 DevRange("IRange", 10);
      Array1DI8 DevStrides("Strides", 5);
      copyReduceInfoToDevice(DevRange, DevStrides, IRange, Strides);
      OMEGA_SCOPE(LocStrides, DevStrides);
      OMEGA_SCOPE(LocRange, DevRange);
      OMEGA_SCOPE(LocArr1, Arr1);
      OMEGA_SCOPE(LocArr2, Arr2);
      parallelReduce(
          {IRange[1] + 1, IRange[3] + 1, IRange[5] + 1, IRange[7] + 1,
           IRange[9] + 1},
          KOKKOS_LAMBDA(int I, int J, int K, int L, int M, I8 &Accum) {
             if (I >= LocRange(0) and J >= LocRange(2) and K >= LocRange(4) and
                 L >= LocRange(6) and M >= LocRange(8)) {
                size_t LinearAdd = I * LocStrides[0] + J * LocStrides[1] +
                                   K * LocStrides[2] + L * LocStrides[3] +
                                   M * LocStrides[4];
                Accum += LocArr1.data()[LinearAdd] * LocArr2.data()[LinearAdd];
             }
          },
          LocalSum);
   } // end if onHost
   // Compute final sum by adding local sums from each MPI task
   I8 GlobalSum;
   int Err =
       MPI_Allreduce(&LocalSum, &GlobalSum, 1, MPI_INT64_T, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (I8 array product): Error in MPI_Allreduce");
   return GlobalSum;
}
/// R4 specific interface
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, R4>
globalSum(const Kokkos::View<T, ML, MS> Arr1, // [in] first  array in product
          const Kokkos::View<T, ML, MS> Arr2, // [in] second array in product
          const MPI_Comm Comm,                // [in] MPI Communicator
          const std::vector<I4> *IndxRange = nullptr ///< [in] opt index range
) {
   // Some error checks
   OMEGA_REQUIRE(Arr1.rank == Arr2.rank,
                 "globalSum (R4 Array product): Arrays must have same rank");
   OMEGA_REQUIRE(Arr1.size() == Arr2.size(),
                 "globalSum (R4 Array product): Arrays must have same size");
   // Get array and layout information
   bool IsHost = isReduceArrayOnHost(Arr1);
   std::vector<I8> Strides(5, 0);
   std::vector<I4> IRange(10, 0);
   getReduceArrayInfo(IRange, Strides, Arr1,
                      IndxRange); ///< [out] true if host array

   // Compute local sum on host or device
   R8 LocalSum = 0;
   if (IsHost) {
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     // convert each to R8 to be sure prod is computed in R8
                     R8 DTmp1 = Arr1.data()[LinearAdd];
                     R8 DTmp2 = Arr2.data()[LinearAdd];
                     LocalSum += DTmp1 * DTmp2;
                  }
               }
            }
         }
      }
   } else { // on device
      // Transfer some info on device
      Array1DI4 DevRange("IRange", 10);
      Array1DI8 DevStrides("Strides", 5);
      copyReduceInfoToDevice(DevRange, DevStrides, IRange, Strides);
      OMEGA_SCOPE(LocStrides, DevStrides);
      OMEGA_SCOPE(LocRange, DevRange);
      OMEGA_SCOPE(LocArr1, Arr1);
      OMEGA_SCOPE(LocArr2, Arr2);
      parallelReduce(
          {IRange[1] + 1, IRange[3] + 1, IRange[5] + 1, IRange[7] + 1,
           IRange[9] + 1},
          KOKKOS_LAMBDA(int I, int J, int K, int L, int M, R8 &Accum) {
             if (I >= LocRange(0) and J >= LocRange(2) and K >= LocRange(4) and
                 L >= LocRange(6) and M >= LocRange(8)) {
                size_t LinearAdd = I * LocStrides[0] + J * LocStrides[1] +
                                   K * LocStrides[2] + L * LocStrides[3] +
                                   M * LocStrides[4];
                // convert each to R8 to be sure prod is computed in R8
                R8 DTmp1 = LocArr1.data()[LinearAdd];
                R8 DTmp2 = LocArr2.data()[LinearAdd];
                Accum += DTmp1 * DTmp2;
             }
          },
          LocalSum);
   } // end if onHost
   // Compute final sum by adding local sums from each MPI task
   R8 GlobalTmp;
   int Err = MPI_Allreduce(&LocalSum, &GlobalTmp, 1, MPI_DOUBLE, MPI_SUM, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (R4 array product): Error in MPI_Allreduce");
   R4 GlobalSum = GlobalTmp;
   return GlobalSum;
}
//*** R8
template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, R8>
globalSum(const Kokkos::View<T, ML, MS> Arr1, // [in] first  array in product
          const Kokkos::View<T, ML, MS> Arr2, // [in] second array in product
          const MPI_Comm Comm,                // [in] MPI Communicator
          const std::vector<I4> *IndxRange = nullptr ///< [in] opt index range
) {
   // Some error checks
   OMEGA_REQUIRE(Arr1.rank == Arr2.rank,
                 "globalSum (R8 Array product): Arrays must have same rank");
   OMEGA_REQUIRE(Arr1.size() == Arr2.size(),
                 "globalSum (R8 Array product): Arrays must have same size");
   // Get array and layout information
   bool IsHost = isReduceArrayOnHost(Arr1);
   std::vector<I8> Strides(5, 0);
   std::vector<I4> IRange(10, 0);
   getReduceArrayInfo(IRange, Strides, Arr1,
                      IndxRange); ///< [out] true if host array

   // Compute local sum on host or device
   complex<double> DDTmp; // Sum and residual for DD algorithm
   if (IsHost) {
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     R8 ProdTmp =
                         Arr1.data()[LinearAdd] * Arr2.data()[LinearAdd];
                     sumDDLocal(DDTmp, ProdTmp);
                  }
               }
            }
         }
      }
   } else { // on device
      // Kokkos sums not reproducible so transfer data to host and compute
      // reproducible sum there
      auto Arr1H = createHostMirrorCopy(Arr1);
      auto Arr2H = createHostMirrorCopy(Arr2);
      for (I4 I = IRange[0]; I <= IRange[1]; ++I) {
         for (I4 J = IRange[2]; J <= IRange[3]; ++J) {
            for (I4 K = IRange[4]; K <= IRange[5]; ++K) {
               for (I4 L = IRange[6]; L <= IRange[7]; ++L) {
                  for (I4 M = IRange[8]; M <= IRange[9]; ++M) {
                     size_t LinearAdd = I * Strides[0] + J * Strides[1] +
                                        K * Strides[2] + L * Strides[3] +
                                        M * Strides[4];
                     R8 ProdTmp =
                         Arr1H.data()[LinearAdd] * Arr2H.data()[LinearAdd];
                     sumDDLocal(DDTmp, ProdTmp);
                  }
               }
            }
         }
      }
      //R8 LocalSum = 0;
      // Transfer some info on device
      //Array1DI4 DevRange("IRange", 10);
      //Array1DI8 DevStrides("Strides", 5);
      //copyReduceInfoToDevice(DevRange, DevStrides, IRange, Strides);
      //OMEGA_SCOPE(LocStrides, DevStrides);
      //OMEGA_SCOPE(LocRange, DevRange);
      //OMEGA_SCOPE(LocArr1, Arr1);
      //OMEGA_SCOPE(LocArr2, Arr2);
      //parallelReduce(
      //    {IRange[1] + 1, IRange[3] + 1, IRange[5] + 1, IRange[7] + 1,
      //     IRange[9] + 1},
      //    KOKKOS_LAMBDA(int I, int J, int K, int L, int M, R8 &Accum) {
      //       if (I >= LocRange(0) and J >= LocRange(2) and K >= LocRange(4) and
      //           L >= LocRange(6) and M >= LocRange(8)) {
      //          size_t LinearAdd = I * LocStrides[0] + J * LocStrides[1] +
      //                             K * LocStrides[2] + L * LocStrides[3] +
      //                             M * LocStrides[4];
      //          Accum += LocArr1.data()[LinearAdd] * LocArr2.data()[LinearAdd];
      //       }
      //    },
      //    LocalSum);
      //DDTmp = complex<double>(LocalSum, 0.0);
   } // end if onHost
   // Compute final sum by adding local sums from each MPI task
   complex<double> GlobalTmp(0.0, 0.0);
   int Err = MPI_Allreduce(&DDTmp, &GlobalTmp, 1, MPI_C_DOUBLE_COMPLEX,
                           MPI_SUMDD, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalSum (R8 array product): Error in MPI_Allreduce");
   R8 GlobalSum = real(GlobalTmp);
   return GlobalSum;
}
///-----------------------------------------------------------------------------
/*
//////////
// Global minval
//////////
// Array
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, void>
globalMinVal(const Kokkos::View<T, ML, MS> arr, const MPI_Comm Comm,
             IT *GlobalMinVal, const std::vector<I4> *IndxRange = nullptr) {
   int dim = arr.rank;
   int i, imin, imax;
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

   int Err;
   if (typeid(IT) == typeid(I4)) {
      Err = MPI_Allreduce(&LocalMinVal, GlobalMinVal, 1, MPI_INT32_T, MPI_MIN,
                          Comm);
   } else if (typeid(IT) == typeid(I8)) {
      Err = MPI_Allreduce(&LocalMinVal, GlobalMinVal, 1, MPI_INT64_T, MPI_MIN,
                          Comm);
   } else if (typeid(IT) == typeid(R4)) {
      Err = MPI_Allreduce(&LocalMinVal, GlobalMinVal, 1, MPI_FLOAT, MPI_MIN,
                          Comm);
   } else if (typeid(IT) == typeid(R8)) {
      Err = MPI_Allreduce(&LocalMinVal, GlobalMinVal, 1, MPI_DOUBLE, MPI_MIN,
                          Comm);
   }
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMinval: Error in MPI_Allreduce");
}

// Array with mask
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, void>
globalMinVal(const Kokkos::View<T, ML, MS> arr,
             const Kokkos::View<T, ML, MS> arr2, const MPI_Comm Comm,
             IT *GlobalMinVal, const std::vector<I4> *IndxRange = nullptr) {
   int dim = arr.rank;
   int i, imin, imax;
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

   int Err;
   if (typeid(IT) == typeid(I4)) {
      Err = MPI_Allreduce(&LocalMinVal, GlobalMinVal, 1, MPI_INT32_T, MPI_MIN,
                          Comm);
   } else if (typeid(IT) == typeid(I8)) {
      Err = MPI_Allreduce(&LocalMinVal, GlobalMinVal, 1, MPI_INT64_T, MPI_MIN,
                          Comm);
   } else if (typeid(IT) == typeid(R4)) {
      Err = MPI_Allreduce(&LocalMinVal, GlobalMinVal, 1, MPI_FLOAT, MPI_MIN,
                          Comm);
   } else if (typeid(IT) == typeid(R8)) {
      Err = MPI_Allreduce(&LocalMinVal, GlobalMinVal, 1, MPI_DOUBLE, MPI_MIN,
                          Comm);
   }
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMinval with mask: Error in MPI_Allreduce");
}

// Array multi-field
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, void>
globalMinVal(const std::vector<Kokkos::View<T, ML, MS>> arrays,
             const MPI_Comm Comm, std::vector<IT> GlobalMinVal,
             const std::vector<I4> *IndxRange = nullptr) {
   int dim = arrays[0].rank;
   int i, imin, imax;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   int ifld, NFields = arrays.size();
   IT LocalMinVal[NFields];
   for (ifld = 0; ifld < NFields; ifld++) {
      LocalMinVal[ifld] = arrays[ifld].data()[imin];
      for (i = imin + 1; i < imax; i++) {
         if (LocalMinVal[ifld] > arrays[ifld].data()[i]) {
            LocalMinVal[ifld] = arrays[ifld].data()[i];
         }
      }
   }
   int Err;
   if (typeid(IT) == typeid(I4)) {
      Err = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], NFields, MPI_INT32_T,
                          MPI_MIN, Comm);
   } else if (typeid(IT) == typeid(I8)) {
      Err = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], NFields, MPI_INT64_T,
                          MPI_MIN, Comm);
   } else if (typeid(IT) == typeid(R4)) {
      Err = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], NFields, MPI_FLOAT,
                          MPI_MIN, Comm);
   } else if (typeid(IT) == typeid(R8)) {
      Err = MPI_Allreduce(LocalMinVal, &GlobalMinVal[0], NFields, MPI_DOUBLE,
                          MPI_MIN, Comm);
   }
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMinval multifield: Error in MPI_Allreduce");
}

//////////
// Global maxval
//////////
// Array
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, void>
globalMaxVal(const Kokkos::View<T, ML, MS> arr, const MPI_Comm Comm,
             IT *GlobalMaxVal, const std::vector<I4> *IndxRange = nullptr) {
   int dim = arr.rank;
   int i, imin, imax;
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

   int Err;
   if (typeid(IT) == typeid(I4)) {
      Err = MPI_Allreduce(&LocalMaxVal, GlobalMaxVal, 1, MPI_INT32_T, MPI_MAX,
                          Comm);
   } else if (typeid(IT) == typeid(I8)) {
      Err = MPI_Allreduce(&LocalMaxVal, GlobalMaxVal, 1, MPI_INT64_T, MPI_MAX,
                          Comm);
   } else if (typeid(IT) == typeid(R4)) {
      Err = MPI_Allreduce(&LocalMaxVal, GlobalMaxVal, 1, MPI_FLOAT, MPI_MAX,
                          Comm);
   } else if (typeid(IT) == typeid(R8)) {
      Err = MPI_Allreduce(&LocalMaxVal, GlobalMaxVal, 1, MPI_DOUBLE, MPI_MAX,
                          Comm);
   }
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMaxval: Error in MPI_Allreduce");
}

// Array with mask
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, void>
globalMaxVal(const Kokkos::View<T, ML, MS> arr,
             const Kokkos::View<T, ML, MS> arr2, const MPI_Comm Comm,
             IT *GlobalMaxVal, const std::vector<I4> *IndxRange = nullptr) {
   int dim = arr.rank;
   int i, imin, imax;
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

   int Err;
   if (typeid(IT) == typeid(I4)) {
      Err = MPI_Allreduce(&LocalMaxVal, GlobalMaxVal, 1, MPI_INT32_T, MPI_MAX,
                          Comm);
   } else if (typeid(IT) == typeid(I8)) {
      Err = MPI_Allreduce(&LocalMaxVal, GlobalMaxVal, 1, MPI_INT64_T, MPI_MAX,
                          Comm);
   } else if (typeid(IT) == typeid(R4)) {
      Err = MPI_Allreduce(&LocalMaxVal, GlobalMaxVal, 1, MPI_FLOAT, MPI_MAX,
                          Comm);
   } else if (typeid(IT) == typeid(R8)) {
      Err = MPI_Allreduce(&LocalMaxVal, GlobalMaxVal, 1, MPI_DOUBLE, MPI_MAX,
                          Comm);
   }
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMaxval with mask: Error in MPI_Allreduce");
}

// Array multi-field
template <typename T, typename IT, typename ML, typename MS>
std::enable_if_t<std::is_same_v<IT, typename Kokkos::View<T>::value_type>, void>
globalMaxVal(const std::vector<Kokkos::View<T, ML, MS>> arrays,
             const MPI_Comm Comm, std::vector<IT> GlobalMaxVal,
             const std::vector<I4> *IndxRange = nullptr) {
   int dim = arrays[0].rank;
   int i, imin, imax;
   if (IndxRange == nullptr) {
      imin = 0;
      imax = arrays[0].size();
   } else {
      imin = (*IndxRange)[0];
      imax = (*IndxRange)[dim * 2 - 1];
   }
   int ifld, NFields = arrays.size();
   IT LocalMaxVal[NFields];
   for (ifld = 0; ifld < NFields; ifld++) {
      LocalMaxVal[ifld] = arrays[ifld].data()[imin];
      for (i = imin + 1; i < imax; i++) {
         if (LocalMaxVal[ifld] < arrays[ifld].data()[i]) {
            LocalMaxVal[ifld] = arrays[ifld].data()[i];
         }
      }
   }
   int Err;
   if (typeid(IT) == typeid(I4)) {
      Err = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], NFields, MPI_INT32_T,
                          MPI_MAX, Comm);
   } else if (typeid(IT) == typeid(I8)) {
      Err = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], NFields, MPI_INT64_T,
                          MPI_MAX, Comm);
   } else if (typeid(IT) == typeid(R4)) {
      Err = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], NFields, MPI_FLOAT,
                          MPI_MAX, Comm);
   } else if (typeid(IT) == typeid(R8)) {
      Err = MPI_Allreduce(LocalMaxVal, &GlobalMaxVal[0], NFields, MPI_DOUBLE,
                          MPI_MAX, Comm);
   }
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMaxVal multifield: Error in MPI_Allreduce");
}

///-----------------------------------------------------------------------------
/// Get MIN-value across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------
void globalMin(const I4 *Val, I4 *Res, const MPI_Comm Comm) {
   int Err = MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_MIN, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMin (I4 scalar): Error in MPI_Allreduce");
}

void globalMin(const I8 *Val, I8 *Res, const MPI_Comm Comm) {
   int Err = MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_MIN, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMin (I8 scalar): Error in MPI_Allreduce");
}

void globalMin(const R4 *Val, R4 *Res, const MPI_Comm Comm) {
   int Err = MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_MIN, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMin (R4 scalar): Error in MPI_Allreduce");
}

void globalMin(const R8 *Val, R8 *Res, const MPI_Comm Comm) {
   int Err = MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_MIN, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMin (R8 scalar): Error in MPI_Allreduce");
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I4, typename Kokkos::View<T>::value_type>, void>
globalMin(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   int Err = MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT32_T,
                           MPI_MIN, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMin (I4 array): Error in MPI_Allreduce");
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I8, typename Kokkos::View<T>::value_type>, void>
globalMin(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   int Err = MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT64_T,
                           MPI_MIN, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMin (I8 array): Error in MPI_Allreduce");
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, void>
globalMin(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   int Err = MPI_Allreduce(in.data(), out.data(), in.size(), MPI_FLOAT, MPI_MIN,
                           Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMin (R4 array): Error in MPI_Allreduce");
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, void>
globalMin(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   int Err = MPI_Allreduce(in.data(), out.data(), in.size(), MPI_DOUBLE,
                           MPI_MIN, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMin (R8 array): Error in MPI_Allreduce");
}

///-----------------------------------------------------------------------------
/// Get MAX-value across all MPI processors in the MachEnv
///-----------------------------------------------------------------------------
void globalMax(const I4 *Val, I4 *Res, const MPI_Comm Comm) {
   int Err = MPI_Allreduce(Val, Res, 1, MPI_INT32_T, MPI_MAX, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMax (I4 scalar): Error in MPI_Allreduce");
}

void globalMax(const I8 *Val, I8 *Res, const MPI_Comm Comm) {
   int Err = MPI_Allreduce(Val, Res, 1, MPI_INT64_T, MPI_MAX, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMax (I8 scalar): Error in MPI_Allreduce");
}

void globalMax(const R4 *Val, R4 *Res, const MPI_Comm Comm) {
   int Err = MPI_Allreduce(Val, Res, 1, MPI_FLOAT, MPI_MAX, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMax (R4 scalar): Error in MPI_Allreduce");
}

void globalMax(const R8 *Val, R8 *Res, const MPI_Comm Comm) {
   int Err = MPI_Allreduce(Val, Res, 1, MPI_DOUBLE, MPI_MAX, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMax (R8 scalar): Error in MPI_Allreduce");
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I4, typename Kokkos::View<T>::value_type>, void>
globalMax(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   int Err = MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT32_T,
                           MPI_MAX, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMax (I4 array): Error in MPI_Allreduce");
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<I8, typename Kokkos::View<T>::value_type>, void>
globalMax(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   int Err = MPI_Allreduce(in.data(), out.data(), in.size(), MPI_INT64_T,
                           MPI_MAX, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMax (I8 array): Error in MPI_Allreduce");
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R4, typename Kokkos::View<T>::value_type>, void>
globalMax(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   int Err = MPI_Allreduce(in.data(), out.data(), in.size(), MPI_FLOAT, MPI_MAX,
                           Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMax (R4 array): Error in MPI_Allreduce");
}

template <typename T, typename ML, typename MS>
std::enable_if_t<std::is_same_v<R8, typename Kokkos::View<T>::value_type>, void>
globalMax(Kokkos::View<T, ML, MS> const in, Kokkos::View<T, ML, MS> out,
          const MPI_Comm Comm) {
   int Err = MPI_Allreduce(in.data(), out.data(), in.size(), MPI_DOUBLE,
                           MPI_MAX, Comm);
   if (Err != MPI_SUCCESS)
      ABORT_ERROR("globalMax (R8 array): Error in MPI_Allreduce");
}
*/
} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // defined OMEGA_REDUCTIONS_H
