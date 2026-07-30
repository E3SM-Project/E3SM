//===-- ocn/StateValidation.cpp - ocean state validation --------*- C++ -*-===//
//
// Validates ocean state fields by checking for NaN values and
// out-of-bounds conditions. Abort behaviour is controlled at runtime via the
// YAML options Omega::State::Validation::AbortOnNan (default true) and
// Omega::State::Validation::AbortOnOOB (default false). When an abort flag is
// set a critical error log with backtrace is emitted and MPI_Abort is called;
// otherwise a CRITICAL log message is written and execution continues.
//
//===----------------------------------------------------------------------===//

#include "StateValidation.h"

#include "AuxiliaryState.h"
#include "Config.h"
#include "DataTypes.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OceanState.h"
#include "OmegaKokkos.h"
#include "Tracers.h"
#include "VertCoord.h"
#include "mpi.h"

#include <cmath>
#include <cpptrace/cpptrace.hpp>
#include <string>
#include <utility>

namespace OMEGA {

struct ValidationAbortConfig {
   bool AbortOnOOB = false;
   bool AbortOnNan = true;
};

static ValidationAbortConfig getValidationAbortConfig() {
   ValidationAbortConfig AbortConfig;

   Config *OmegaConfig = Config::getOmegaConfig();
   if (!OmegaConfig) {
      return AbortConfig;
   }

   Config StateConfig("State");
   if (!OmegaConfig->get(StateConfig).isSuccess()) {
      return AbortConfig;
   }

   Config ValidationConfig("Validation");
   if (!StateConfig.get(ValidationConfig).isSuccess()) {
      return AbortConfig;
   }

   ValidationConfig.get("AbortOnOOB", AbortConfig.AbortOnOOB);
   ValidationConfig.get("AbortOnNan", AbortConfig.AbortOnNan);

   return AbortConfig;
}

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
// array over the first NCells/NEdges rows and NVert columns, restricted to
// active cells (CellMask > 0).
// Returns {NaNCount, OutOfRangeCount}.
static std::pair<I4, I4> checkArray2D(const Array2DReal &Arr, I4 NRows,
                                      I4 NCols, Real MinVal, Real MaxVal,
                                      bool CheckMin,
                                      const Array2DReal &CellMask) {
   I4 NaNCount        = 0;
   I4 OutOfRangeCount = 0;

   parallelReduce(
       "CheckNaN", {NRows, NCols},
       KOKKOS_LAMBDA(int Row, int Col, int &Accum) {
          if (CellMask(Row, Col) == 0) {
             return;
          }
          Real Val = Arr(Row, Col);
          if (Kokkos::isnan(Val)) {
             ++Accum;
          }
       },
       NaNCount);

   parallelReduce(
       "CheckBounds", {NRows, NCols},
       KOKKOS_LAMBDA(int Row, int Col, int &Accum) {
          if (CellMask(Row, Col) == 0) {
             return;
          }
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
// col = vert) extracted from the 3-D tracer array at the given tracer index,
// restricted to active cells (CellMask > 0).
static std::pair<I4, I4> checkTracerArray(const Array3DReal &Tracers3D,
                                          I4 TracerIdx, I4 NCells, I4 NVert,
                                          Real MinVal, Real MaxVal,
                                          const Array2DReal &CellMask) {
   I4 NaNCount        = 0;
   I4 OutOfRangeCount = 0;

   parallelReduce(
       "CheckTracerNaN", {NCells, NVert},
       KOKKOS_LAMBDA(int Cell, int K, int &Accum) {
          if (CellMask(Cell, K) == 0) {
             return;
          }
          Real Val = Tracers3D(TracerIdx, Cell, K);
          if (Kokkos::isnan(Val)) {
             ++Accum;
          }
       },
       NaNCount);

   parallelReduce(
       "CheckTracerBounds", {NCells, NVert},
       KOKKOS_LAMBDA(int Cell, int K, int &Accum) {
          if (CellMask(Cell, K) == 0) {
             return;
          }
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
/// Check ocean state fields for NaN and out-of-bounds conditions.
/// Only active cells (where CellMask > 0) are checked.
/// Returns the total count of errors found; does not abort.
std::pair<I4, I4> checkOceanState(const OceanState *State,
                                  const AuxiliaryState *AuxState,
                                  const VertCoord *VCoord, I4 TimeLevel) {

   I4 TotalNaNs = 0;
   I4 TotalOOBs = 0;

   const Array2DReal &CellMask = VCoord->CellMask;

   // -------------------------------------------------------------------------
   // PseudoThickness: valid range [1e-10, 1000]
   // -------------------------------------------------------------------------
   {
      Array2DReal PseudoThick = State->getPseudoThickness(TimeLevel);
      auto [NaNs, OOB] =
          checkArray2D(PseudoThick, State->NCellsOwned, State->NVertLayers,
                       static_cast<Real>(1e-10), static_cast<Real>(1000.0),
                       /*CheckMin=*/true, CellMask);

      if (NaNs > 0) {
         LOG_CRITICAL(
             "StateValidation: PseudoThickness contains {} NaN value(s)", NaNs);
         TotalNaNs += NaNs;
      }
      if (OOB > 0) {
         LOG_CRITICAL(
             "StateValidation: PseudoThickness has {} value(s) outside "
             "valid range [1e-10, 1000]",
             OOB);
         TotalOOBs += OOB;
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
                       /*CheckMin=*/true, CellMask);

      if (NaNs > 0) {
         LOG_CRITICAL(
             "StateValidation: KineticEnergyCell contains {} NaN value(s)",
             NaNs);
         TotalNaNs += NaNs;
      }
      if (OOB > 0) {
         LOG_CRITICAL(
             "StateValidation: KineticEnergyCell has {} value(s) outside "
             "valid range [0, 10]",
             OOB);
         TotalOOBs += OOB;
      }
   }

   // -------------------------------------------------------------------------
   // Temperature tracer: valid range [-10, 50]
   // -------------------------------------------------------------------------
   if (Tracers::IndxTemp != -1) {
      Array3DReal AllTracers = Tracers::getAll(TimeLevel);
      auto [NaNs, OOB]       = checkTracerArray(
          AllTracers, Tracers::IndxTemp, State->NCellsOwned, State->NVertLayers,
          static_cast<Real>(-10.0), static_cast<Real>(50.0), CellMask);

      if (NaNs > 0) {
         LOG_CRITICAL("StateValidation: Temperature contains {} NaN value(s)",
                      NaNs);
         TotalNaNs += NaNs;
      }
      if (OOB > 0) {
         LOG_CRITICAL("StateValidation: Temperature has {} value(s) outside "
                      "valid range [-10, 50]",
                      OOB);
         TotalOOBs += OOB;
      }
   }

   // -------------------------------------------------------------------------
   // Salinity tracer: valid range [-2, 60]
   // -------------------------------------------------------------------------
   if (Tracers::IndxSalt != -1) {
      Array3DReal AllTracers = Tracers::getAll(TimeLevel);
      auto [NaNs, OOB]       = checkTracerArray(
          AllTracers, Tracers::IndxSalt, State->NCellsOwned, State->NVertLayers,
          static_cast<Real>(-2.0), static_cast<Real>(60.0), CellMask);

      if (NaNs > 0) {
         LOG_CRITICAL("StateValidation: Salinity contains {} NaN value(s)",
                      NaNs);
         TotalNaNs += NaNs;
      }
      if (OOB > 0) {
         LOG_CRITICAL("StateValidation: Salinity has {} value(s) outside "
                      "valid range [-2, 60]",
                      OOB);
         TotalOOBs += OOB;
      }
   }

   return {TotalNaNs, TotalOOBs};
}

//------------------------------------------------------------------------------
/// Validate ocean state fields for NaN and out-of-bounds conditions.
/// Only active cells (where CellMask > 0) are checked.
/// Aborts via MPI_Abort on failure.
void validateOceanState(const OceanState *State, const AuxiliaryState *AuxState,
                        const VertCoord *VCoord, I4 TimeLevel) {

   auto [NaNs, OOBs]      = checkOceanState(State, AuxState, VCoord, TimeLevel);
   const auto AbortConfig = getValidationAbortConfig();
   bool abort             = false;

   if (NaNs > 0 && AbortConfig.AbortOnNan) {
      abort = true;
   }
   if (OOBs > 0 && AbortConfig.AbortOnOOB) {
      abort = true;
   }
   if (abort) {
      abortWithMessage(
          "StateValidation: Ocean state validation aborting with " +
          std::to_string(NaNs) + " NaNs and " + std::to_string(OOBs) +
          " OOBs." + "See critical messages above for details.");
   } else if ((NaNs + OOBs) > 0) {
      // Only report the flags that are actually suppressing an abort, i.e.
      // those corresponding to the issues that were detected
      std::string Flags;
      if (NaNs > 0) {
         Flags = "AbortOnNan is false";
      }
      if (OOBs > 0) {
         if (!Flags.empty())
            Flags += " and ";
         Flags += "AbortOnOOB is false";
      }
      LOG_CRITICAL("StateValidation: Ocean state validation failed with {} "
                   "NaNs and {} OOBs. {}; continuing after logging.",
                   NaNs, OOBs, Flags);
   }
}

} // namespace OMEGA

//===----------------------------------------------------------------------===//
