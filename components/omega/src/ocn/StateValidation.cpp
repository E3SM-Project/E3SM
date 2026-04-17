//===-- ocn/StateValidation.cpp - ocean state validation --------*- C++ -*-===//
//
// Validates ocean state fields by checking for NaN values and
// out-of-bounds conditions. Any failure triggers a critical error log with
// backtrace and MPI_Abort on the local communicator.
//
//===----------------------------------------------------------------------===//

#include "StateValidation.h"

#include "AuxiliaryState.h"
#include "DataTypes.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Tracers.h"
#include "mpi.h"

#include <cmath>
#include <cpptrace/cpptrace.hpp>
#include <string>
#include <utility>

namespace OMEGA {

//------------------------------------------------------------------------------
// Helper: abort on the local Omega communicator with a message and backtrace
static void abortWithMessage(const std::string &Msg) {
   LOG_CRITICAL("{}", Msg);
   cpptrace::generate_trace().print();
   MPI_Comm Comm = MachEnv::getDefault()->getComm();
   MPI_Abort(Comm, static_cast<int>(ErrorCode::Critical));
}

//------------------------------------------------------------------------------
// Helper: count NaN entries and out-of-range entries in a 2-D Real device
// array over the first NCells/NEdges rows and NVert columns.
// Returns {NaNCount, OutOfRangeCount}.
static std::pair<I4, I4> checkArray2D(const Array2DReal &Arr, I4 NRows,
                                      I4 NCols, Real MinVal, Real MaxVal,
                                      bool CheckMin) {
   I4 NaNCount        = 0;
   I4 OutOfRangeCount = 0;

   parallelReduce(
       "CheckNaN", {NRows, NCols},
       KOKKOS_LAMBDA(int Row, int Col, int &Accum) {
          Real Val = Arr(Row, Col);
          if (Kokkos::isnan(Val)) {
             ++Accum;
          }
       },
       NaNCount);

   parallelReduce(
       "CheckBounds", {NRows, NCols},
       KOKKOS_LAMBDA(int Row, int Col, int &Accum) {
          Real Val = Arr(Row, Col);
          if (!Kokkos::isnan(Val)) {
             if (Val > MaxVal) {
                ++Accum;
             } else if (CheckMin && Val < MinVal) {
                ++Accum;
             }
          }
       },
       OutOfRangeCount);

   return {NaNCount, OutOfRangeCount};
}

//------------------------------------------------------------------------------
// Helper: count NaN and out-of-range entries for a single tracer (row = cell,
// col = vert) extracted from the 3-D tracer array at the given tracer index.
static std::pair<I4, I4> checkTracerArray(const Array3DReal &Tracers3D,
                                          I4 TracerIdx, I4 NCells, I4 NVert,
                                          Real MinVal, Real MaxVal) {
   I4 NaNCount        = 0;
   I4 OutOfRangeCount = 0;

   parallelReduce(
       "CheckTracerNaN", {NCells, NVert},
       KOKKOS_LAMBDA(int Cell, int K, int &Accum) {
          Real Val = Tracers3D(TracerIdx, Cell, K);
          if (Kokkos::isnan(Val)) {
             ++Accum;
          }
       },
       NaNCount);

   parallelReduce(
       "CheckTracerBounds", {NCells, NVert},
       KOKKOS_LAMBDA(int Cell, int K, int &Accum) {
          Real Val = Tracers3D(TracerIdx, Cell, K);
          if (!Kokkos::isnan(Val)) {
             if (Val < MinVal || Val > MaxVal) {
                ++Accum;
             }
          }
       },
       OutOfRangeCount);

   return {NaNCount, OutOfRangeCount};
}

//------------------------------------------------------------------------------
/// Validate ocean state fields for NaN and out-of-bounds conditions.
/// Aborts via MPI_Abort on failure.
void validateOceanState(const OceanState *State, const AuxiliaryState *AuxState,
                        I4 TimeLevel) {

   bool AnyFailure = false;

   // -------------------------------------------------------------------------
   // LayerThickness: valid range [1e-10, 1000]
   // -------------------------------------------------------------------------
   {
      Array2DReal LayerThick = State->getLayerThickness(TimeLevel);
      auto [NaNs, OOB] =
          checkArray2D(LayerThick, State->NCellsOwned, State->NVertLayers,
                       static_cast<Real>(1e-10), static_cast<Real>(1000.0),
                       /*CheckMin=*/true);

      if (NaNs > 0) {
         LOG_CRITICAL(
             "StateValidation: LayerThickness contains {} NaN value(s)", NaNs);
         AnyFailure = true;
      }
      if (OOB > 0) {
         LOG_CRITICAL("StateValidation: LayerThickness has {} value(s) outside "
                      "valid range [1e-10, 1000]",
                      OOB);
         AnyFailure = true;
      }
   }

   // -------------------------------------------------------------------------
   // KineticEnergyCell: valid range [0, 10]
   // -------------------------------------------------------------------------
   {
      const Array2DReal &KE = AuxState->KineticAux.KineticEnergyCell;
      auto [NaNs, OOB] =
          checkArray2D(KE, State->NCellsOwned, State->NVertLayers,
                       static_cast<Real>(0.0), static_cast<Real>(10.0),
                       /*CheckMin=*/true);

      if (NaNs > 0) {
         LOG_CRITICAL(
             "StateValidation: KineticEnergyCell contains {} NaN value(s)",
             NaNs);
         AnyFailure = true;
      }
      if (OOB > 0) {
         LOG_CRITICAL(
             "StateValidation: KineticEnergyCell has {} value(s) outside "
             "valid range [0, 10]",
             OOB);
         AnyFailure = true;
      }
   }

   // -------------------------------------------------------------------------
   // Temperature tracer: valid range [-10, 50]
   // -------------------------------------------------------------------------
   if (Tracers::IndxTemp != -1) {
      Array3DReal AllTracers = Tracers::getAll(TimeLevel);
      auto [NaNs, OOB]       = checkTracerArray(
          AllTracers, Tracers::IndxTemp, State->NCellsOwned, State->NVertLayers,
          static_cast<Real>(-10.0), static_cast<Real>(50.0));

      if (NaNs > 0) {
         LOG_CRITICAL("StateValidation: Temperature contains {} NaN value(s)",
                      NaNs);
         AnyFailure = true;
      }
      if (OOB > 0) {
         LOG_CRITICAL("StateValidation: Temperature has {} value(s) outside "
                      "valid range [-10, 50]",
                      OOB);
         AnyFailure = true;
      }
   }

   // -------------------------------------------------------------------------
   // Salinity tracer: valid range [-2, 60]
   // -------------------------------------------------------------------------
   if (Tracers::IndxSalt != -1) {
      Array3DReal AllTracers = Tracers::getAll(TimeLevel);
      auto [NaNs, OOB]       = checkTracerArray(
          AllTracers, Tracers::IndxSalt, State->NCellsOwned, State->NVertLayers,
          static_cast<Real>(-2.0), static_cast<Real>(60.0));

      if (NaNs > 0) {
         LOG_CRITICAL("StateValidation: Salinity contains {} NaN value(s)",
                      NaNs);
         AnyFailure = true;
      }
      if (OOB > 0) {
         LOG_CRITICAL("StateValidation: Salinity has {} value(s) outside "
                      "valid range [-2, 60]",
                      OOB);
         AnyFailure = true;
      }
   }

   // -------------------------------------------------------------------------
   // Abort if any check failed
   // -------------------------------------------------------------------------
   if (AnyFailure) {
      abortWithMessage("StateValidation: Ocean state validation failed. "
                       "See critical messages above for details.");
   }
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
