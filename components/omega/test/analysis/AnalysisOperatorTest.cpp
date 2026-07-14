//===-- Test 5.1: Multi-type operator correctness tests --------*- C++ -*-===//
//
// Comprehensive unit tests for Analysis operators across all supported
// array types (ranks 1D/2D/3D and scalar types I4/I8/R4/R8)
//
//===-----------------------------------------------------------------------===//

#include "Analysis.h"
#include "AnalysisOpFactory.h"
#include "Decomp.h"
#include "Field.h"
#include "Forcing.h"
#include "Halo.h"
#include "HorzMesh.h"
#include "IO.h"
#include "IOStream.h"
#include "Logging.h"
#include "TimeStepper.h"
#include "VertCoord.h"
#include "VertMix.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <type_traits>
#include <vector>

using namespace OMEGA;

// Test result tracking
int NumTests  = 0;
int NumPassed = 0;
int NumFailed = 0;

//===----------------------------------------------------------------------===//
// Generic Helper Template Struct
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Template struct consolidating all test helper functions
template <typename ArrayType> struct TestHelper {
   using ScalarT             = typename ArrayType::non_const_value_type;
   static constexpr int Rank = ArrayType::rank;

   // Type-aware tolerance for floating point comparisons
   static ScalarT getTolerance() {
      if constexpr (std::is_integral_v<ScalarT>) {
         return 0; // Exact equality for integers
      } else if constexpr (std::is_same_v<ScalarT, float>) {
         return 1.0e-4f; // Single precision tolerance
      } else {
         return 1.0e-8; // Double precision tolerance
      }
   }

   // Get dimensions based on rank
   static std::vector<I4> getDims(const HorzMesh *Mesh,
                                  const VertCoord *VCoord) {
      if constexpr (Rank == 1) {
         return {Mesh->NCellsSize}; // 1D horizontal array over cells
      } else if constexpr (Rank == 2) {
         return {Mesh->NCellsSize, VCoord->NVertLayers};
      } else if constexpr (Rank == 3) {
         return {Tracers::getNumTracers(), Mesh->NCellsSize,
                 VCoord->NVertLayers};
      }
      return {};
   }

   // Get dimension names
   static std::vector<std::string> getDimNames() {
      if constexpr (Rank == 1) {
         return {"NCells"};
      } else if constexpr (Rank == 2) {
         return {"NCells", "NVertLayers"};
      } else if constexpr (Rank == 3) {
         return {"NTracers", "NCells", "NVertLayers"};
      }
      return {};
   }

   // Create test field for 1D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 1, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4)> ValueFunc) {

      auto DimNames = getDimNames();
      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", -1.0e30, 1.0e30, 1, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0]);
      TestField->attachData<ArrayType>(TestData);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         TestDataHost(i) = ValueFunc(i);
      }
      Kokkos::deep_copy(TestData, TestDataHost);
   }

   // Create test field for 2D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 2, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4, I4)> ValueFunc) {

      auto DimNames = getDimNames();
      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", -1.0e30, 1.0e30, 2, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1]);
      TestField->attachData<ArrayType>(TestData);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Dims[1]; ++j) {
            TestDataHost(i, j) = ValueFunc(i, j);
         }
      }
      Kokkos::deep_copy(TestData, TestDataHost);
   }

   // Create test field for 3D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 3, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4, I4, I4)> ValueFunc) {

      auto DimNames = getDimNames();
      auto TestField =
          Field::create(FieldName, "Test field for multi-type validation", "m",
                        "", -1.0e30, 1.0e30, 3, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1], Dims[2]);
      TestField->attachData<ArrayType>(TestData);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Dims[1]; ++j) {
            for (I4 k = 0; k < Dims[2]; ++k) {
               TestDataHost(i, j, k) = ValueFunc(i, j, k);
            }
         }
      }
      Kokkos::deep_copy(TestData, TestDataHost);
   }
};

//------------------------------------------------------------------------------
// Helper function to report test results
void reportTest(const std::string &TestName, bool Passed) {
   NumTests++;
   if (Passed) {
      NumPassed++;
   } else {
      NumFailed++;
      LOG_ERROR("FAIL: {}", TestName);
   }
}

//===----------------------------------------------------------------------===//
// Operator Test Templates
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Template for testing SpatialMaxOp with any array type
template <typename ArrayType>
void testSpatialMaxOpType(const std::string &TypeName, const MachEnv *Env,
                          const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldMax_" + TypeName;

   // Create test field with values based on global cell IDs for MPI correctness
   // Get global cell IDs to ensure unique values across all ranks
   auto Decomp  = Decomp::getDefault();
   auto CellIDH = Decomp->CellIDH; // Global cell IDs (1-based)

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [CellIDH](I4 i) -> ScalarT {
         return static_cast<ScalarT>(CellIDH(i) - 1); // Convert to 0-based
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims,
                          [CellIDH, VCoord](I4 i, I4 j) -> ScalarT {
                             // Unique value: cellID * NVertLayers + layerIndex
                             return static_cast<ScalarT>(
                                 (CellIDH(i) - 1) * VCoord->NVertLayers + j);
                          });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          FieldName, Dims,
          [CellIDH, VCoord, Mesh](I4 i, I4 j, I4 k) -> ScalarT {
             // Unique value: tracerIdx * (NCellsGlobal * NVertLayers) + cellID
             // * NVertLayers + layerIdx
             return static_cast<ScalarT>(
                 i * (Mesh->NCellsGlobal * VCoord->NVertLayers) +
                 (CellIDH(j) - 1) * VCoord->NVertLayers + k);
          });
   }

   // Compute expected max. The operator performs a global reduction across all
   // ranks via MPI_Allreduce, so the expected maximum is based on global mesh
   // properties. The maximum value corresponds to the largest indices in each
   // dimension.
   ScalarT ExpectedMax = 0;

   if constexpr (Rank == 1) {
      // 1D: max is just the highest cell ID (0-based)
      ExpectedMax = static_cast<ScalarT>(Mesh->NCellsGlobal - 1);
   } else if constexpr (Rank == 2) {
      // 2D: max = (NCellsGlobal - 1) * NVertLayers + (NVertLayers - 1)
      ExpectedMax =
          static_cast<ScalarT>((Mesh->NCellsGlobal - 1) * VCoord->NVertLayers +
                               (VCoord->NVertLayers - 1));
   } else if constexpr (Rank == 3) {
      // 3D: max = (NTracers - 1) * (NCellsGlobal * NVertLayers) +
      //           (NCellsGlobal - 1) * NVertLayers + (NVertLayers - 1)
      I4 NTracers = Tracers::getNumTracers();
      ExpectedMax = static_cast<ScalarT>(
          (NTracers - 1) * (Mesh->NCellsGlobal * VCoord->NVertLayers) +
          (Mesh->NCellsGlobal - 1) * VCoord->NVertLayers +
          (VCoord->NVertLayers - 1));
   }

   // Create and compute operator
   Config EmptyConfig;
   auto MaxOp =
       AnalysisOpFactory::createOp("SpatialMax", {FieldName}, EmptyConfig);
   MaxOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MaxOp->compute(TestTime);

   // Get result. The operator attaches output as Array1D<ScalarT>, so retrieve
   // with the matching type to avoid reinterpreting bits via
   // static_pointer_cast.
   auto ResultField = Field::get(FieldName + "_SpatialMax");
   auto ResultData =
       ResultField->getDataArray<typename Array1D<ScalarT>::type>();
   auto ResultHost = Kokkos::create_mirror_view(ResultData);
   Kokkos::deep_copy(ResultHost, ResultData);

   Real ComputedMax     = static_cast<Real>(ResultHost(0));
   Real ExpectedMaxReal = static_cast<Real>(ExpectedMax);

   // Verify
   bool Passed = (std::abs(ComputedMax - ExpectedMaxReal) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialMaxOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Got: {}", ExpectedMaxReal, ComputedMax);
   }
}

//------------------------------------------------------------------------------
// Template for testing SpatialMinOp with any array type
template <typename ArrayType>
void testSpatialMinOpType(const std::string &TypeName, const MachEnv *Env,
                          const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldMin_" + TypeName;

   // Create test field with values based on global cell IDs for MPI correctness
   auto Decomp  = Decomp::getDefault();
   auto CellIDH = Decomp->CellIDH; // Global cell IDs (1-based)

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [CellIDH](I4 i) -> ScalarT {
         return static_cast<ScalarT>(CellIDH(i) - 1 +
                                     100); // Convert to 0-based, offset by 100
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(
          FieldName, Dims, [CellIDH, VCoord](I4 i, I4 j) -> ScalarT {
             // Unique value: cellID * NVertLayers + layerIndex + offset
             return static_cast<ScalarT>(
                 (CellIDH(i) - 1) * VCoord->NVertLayers + j + 100);
          });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          FieldName, Dims,
          [CellIDH, VCoord, Mesh](I4 i, I4 j, I4 k) -> ScalarT {
             // Unique value: tracerIdx * (NCellsGlobal * NVertLayers) + cellID
             // * NVertLayers + layerIdx + offset
             return static_cast<ScalarT>(
                 i * (Mesh->NCellsGlobal * VCoord->NVertLayers) +
                 (CellIDH(j) - 1) * VCoord->NVertLayers + k + 100);
          });
   }

   // Compute expected min. The operator performs a global reduction across all
   // ranks via MPI_Allreduce, so the expected minimum is based on global mesh
   // properties. The minimum value is always at i=0, j=0 (cell with global ID
   // 0), k=0, plus offset 100.
   ScalarT ExpectedMin = static_cast<ScalarT>(100);

   // Create and compute operator
   Config EmptyConfig;
   auto MinOp =
       AnalysisOpFactory::createOp("SpatialMin", {FieldName}, EmptyConfig);
   MinOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MinOp->compute(TestTime);

   // Get result. The operator attaches output as Array1D<ScalarT>, so retrieve
   // with the matching type to avoid reinterpreting bits via
   // static_pointer_cast.
   auto ResultField = Field::get(FieldName + "_SpatialMin");
   auto ResultData =
       ResultField->getDataArray<typename Array1D<ScalarT>::type>();
   auto ResultHost = Kokkos::create_mirror_view(ResultData);
   Kokkos::deep_copy(ResultHost, ResultData);

   Real ComputedMin     = static_cast<Real>(ResultHost(0));
   Real ExpectedMinReal = static_cast<Real>(ExpectedMin);

   // Verify
   bool Passed = (std::abs(ComputedMin - ExpectedMinReal) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialMinOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Got: {}", ExpectedMinReal, ComputedMin);
   }
}

//------------------------------------------------------------------------------
// Template for testing SpatialMeanOp with any array type
template <typename ArrayType>
void testSpatialMeanOpType(const std::string &TypeName, const MachEnv *Env,
                           const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldMean_" + TypeName;

   // Create test field with alternating values to properly test mean
   // calculation
   ScalarT Value1 = static_cast<ScalarT>(10);
   ScalarT Value2 = static_cast<ScalarT>(20);

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [Value1, Value2](I4 i) -> ScalarT {
         return ((i % 2) == 0) ? Value1 : Value2;
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims,
                          [Value1, Value2](I4 i, I4 j) -> ScalarT {
                             return (((i + j) % 2) == 0) ? Value1 : Value2;
                          });
   } else if constexpr (Rank == 3) {
      Helper::createField(FieldName, Dims,
                          [Value1, Value2](I4 i, I4 j, I4 k) -> ScalarT {
                             return (((i + j + k) % 2) == 0) ? Value1 : Value2;
                          });
   }

   // Calculate expected mean by counting actual Value1 and Value2 elements
   // in the active region (accounting for masked layers)
   I8 LocalCount1 = 0, LocalCount2 = 0;

   if constexpr (Rank == 1) {
      // 1D: count based on horizontal index pattern
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         if ((i % 2) == 0)
            LocalCount1++;
         else
            LocalCount2++;
      }
   } else if constexpr (Rank == 2) {
      // 2D: count based on (i+j) pattern within active vertical layers
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
              ++j) {
            if (((i + j) % 2) == 0)
               LocalCount1++;
            else
               LocalCount2++;
         }
      }
   } else if constexpr (Rank == 3) {
      // 3D: count based on (t+i+j) pattern within active layers
      I4 NTracers = Tracers::getNumTracers();
      for (I4 t = 0; t < NTracers; ++t) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
            for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
                 ++j) {
               if (((t + i + j) % 2) == 0)
                  LocalCount1++;
               else
                  LocalCount2++;
            }
         }
      }
   }

   // Global sum across MPI ranks
   I8 Count1         = globalSum(LocalCount1, Env->getComm());
   I8 Count2         = globalSum(LocalCount2, Env->getComm());
   I8 TotalElements  = Count1 + Count2;
   Real ExpectedMean = (static_cast<Real>(Value1) * static_cast<Real>(Count1) +
                        static_cast<Real>(Value2) * static_cast<Real>(Count2)) /
                       static_cast<Real>(TotalElements);

   // Create and compute operator
   Config EmptyConfig;
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, EmptyConfig);
   MeanOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MeanOp->compute(TestTime);

   // Get result. The operator attaches output as Array1DReal (always Real type
   // regardless of input type).
   auto ResultField = Field::get(FieldName + "_SpatialMean");
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = Kokkos::create_mirror_view(ResultData);
   Kokkos::deep_copy(ResultHost, ResultData);

   Real ComputedMean = ResultHost(0);

   // Verify
   bool Passed = (std::abs(ComputedMean - ExpectedMean) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialMeanOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Got: {}", ExpectedMean, ComputedMean);
   }
}

//------------------------------------------------------------------------------
// Template for testing SpatialStdDevOp with any array type
template <typename ArrayType>
void testSpatialStdDevOpType(const std::string &TypeName, const MachEnv *Env,
                             const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldStdDev_" + TypeName;

   // Create test field with alternating values to properly test std dev
   // calculation
   ScalarT Value1 = static_cast<ScalarT>(10);
   ScalarT Value2 = static_cast<ScalarT>(20);

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [Value1, Value2](I4 i) -> ScalarT {
         return ((i % 2) == 0) ? Value1 : Value2;
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims,
                          [Value1, Value2](I4 i, I4 j) -> ScalarT {
                             return (((i + j) % 2) == 0) ? Value1 : Value2;
                          });
   } else if constexpr (Rank == 3) {
      Helper::createField(FieldName, Dims,
                          [Value1, Value2](I4 i, I4 j, I4 k) -> ScalarT {
                             return (((i + j + k) % 2) == 0) ? Value1 : Value2;
                          });
   }

   // Calculate expected standard deviation by counting actual Value1 and Value2
   // elements in the active region (accounting for masked layers)
   I8 LocalCount1 = 0, LocalCount2 = 0;

   if constexpr (Rank == 1) {
      // 1D: count based on horizontal index pattern
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         if ((i % 2) == 0)
            LocalCount1++;
         else
            LocalCount2++;
      }
   } else if constexpr (Rank == 2) {
      // 2D: count based on (i+j) pattern within active vertical layers
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
              ++j) {
            if (((i + j) % 2) == 0)
               LocalCount1++;
            else
               LocalCount2++;
         }
      }
   } else if constexpr (Rank == 3) {
      // 3D: count based on (t+i+j) pattern within active layers
      I4 NTracers = Tracers::getNumTracers();
      for (I4 t = 0; t < NTracers; ++t) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
            for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
                 ++j) {
               if (((t + i + j) % 2) == 0)
                  LocalCount1++;
               else
                  LocalCount2++;
            }
         }
      }
   }

   // Global sum across MPI ranks
   I8 Count1        = globalSum(LocalCount1, Env->getComm());
   I8 Count2        = globalSum(LocalCount2, Env->getComm());
   I8 TotalElements = Count1 + Count2;
   Real Mean        = (static_cast<Real>(Value1) * static_cast<Real>(Count1) +
                static_cast<Real>(Value2) * static_cast<Real>(Count2)) /
               static_cast<Real>(TotalElements);

   // Standard deviation: sqrt(sum((x_i - mean)^2) / N)
   Real SumSquaredDiff = static_cast<Real>(Count1) *
                             std::pow(static_cast<Real>(Value1) - Mean, 2.0) +
                         static_cast<Real>(Count2) *
                             std::pow(static_cast<Real>(Value2) - Mean, 2.0);
   Real ExpectedStdDev =
       std::sqrt(SumSquaredDiff / static_cast<Real>(TotalElements));

   // SpatialStdDevOp requires a pre-existing _SpatialMean field for the input.
   // Create and compute a SpatialMeanOp first so that field is registered.
   Config EmptyConfig;
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, EmptyConfig);
   MeanOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MeanOp->compute(TestTime);

   // Now create and compute the StdDev operator
   auto StdDevOp =
       AnalysisOpFactory::createOp("SpatialStdDev", {FieldName}, EmptyConfig);
   StdDevOp->initialize(Env, Mesh, VCoord, EmptyConfig);
   StdDevOp->compute(TestTime);

   // Get result. The operator attaches output as Array1DReal (always Real type
   // regardless of input type).
   auto ResultField = Field::get(FieldName + "_SpatialStdDev");
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = Kokkos::create_mirror_view(ResultData);
   Kokkos::deep_copy(ResultHost, ResultData);

   Real ComputedStdDev = ResultHost(0);

   // Verify
   bool Passed = (std::abs(ComputedStdDev - ExpectedStdDev) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialStdDevOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Got: {}", ExpectedStdDev, ComputedStdDev);
   }
}

//------------------------------------------------------------------------------
// Template for testing TimeMeanOp with any array type
template <typename ArrayType>
void testTimeMeanOpType(const std::string &TypeName, const MachEnv *Env,
                        const HorzMesh *Mesh, const VertCoord *VCoord,
                        Clock *ModelClock) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldTimeMean_" + TypeName;

   // Create test field with initial value
   ScalarT BaseValue = static_cast<ScalarT>(5);

   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims,
                          [BaseValue](I4 i) -> ScalarT { return BaseValue; });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims, [BaseValue](I4 i, I4 j) -> ScalarT {
         return BaseValue;
      });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          FieldName, Dims,
          [BaseValue](I4 i, I4 j, I4 k) -> ScalarT { return BaseValue; });
   }

   // Get the field and its data array for updating during time loop
   auto TestField = Field::get(FieldName);
   auto TestData  = TestField->template getDataArray<ArrayType>();

   // Set up time stepping parameters
   const int NumSteps        = 5; // Accumulate over 5 timesteps
   TimeInterval StepInterval = ModelClock->getTimeStep();

   // Calculate period interval (NumSteps * timestep)
   R8 StepSeconds;
   StepInterval.get(StepSeconds, TimeUnits::Seconds);
   R8 PeriodSeconds = StepSeconds * NumSteps;
   TimeInterval PeriodInterval(PeriodSeconds, TimeUnits::Seconds);

   // Create TimeMeanOp with a valid period string (e.g., "5seconds")
   // The Period string is just a label used in the output field name
   Config OpConfig;
   std::string PeriodLabel =
       std::to_string(static_cast<int>(PeriodSeconds)) + "seconds";
   //   OpConfig.set("Period", PeriodLabel);

   auto TimeMeanOp = AnalysisOpFactory::createOp(
       "TimeMean", {FieldName}, makeOpConfig(opParam("Period", PeriodLabel)));
   TimeMeanOp->initialize(Env, Mesh, VCoord, OpConfig);

   // Create a period alarm that rings after NumSteps
   TimeInstant StartTime = ModelClock->getCurrentTime();
   Alarm PeriodAlarm("TestPeriodAlarm_" + TypeName, PeriodInterval, StartTime);
   TimeMeanOp->setPeriodAlarm(&PeriodAlarm);

   // Time-stepping loop: update field values and compute mean at each step
   std::vector<ScalarT> ValuesAtEachStep;

   for (int step = 0; step < NumSteps; ++step) {
      // Update field values to simulate time evolution
      // Value at each step = BaseValue + step (e.g., 5, 6, 7, 8, 9)
      ScalarT CurrentValue = static_cast<ScalarT>(static_cast<Real>(BaseValue) +
                                                  static_cast<Real>(step));
      ValuesAtEachStep.push_back(CurrentValue);

      // Update the field data on device
      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      Kokkos::deep_copy(TestDataHost, TestData);

      if constexpr (Rank == 1) {
         for (I4 i = 0; i < Dims[0]; ++i) {
            TestDataHost(i) = CurrentValue;
         }
      } else if constexpr (Rank == 2) {
         for (I4 i = 0; i < Dims[0]; ++i) {
            for (I4 j = 0; j < Dims[1]; ++j) {
               TestDataHost(i, j) = CurrentValue;
            }
         }
      } else if constexpr (Rank == 3) {
         for (I4 i = 0; i < Dims[0]; ++i) {
            for (I4 j = 0; j < Dims[1]; ++j) {
               for (I4 k = 0; k < Dims[2]; ++k) {
                  TestDataHost(i, j, k) = CurrentValue;
               }
            }
         }
      }

      Kokkos::deep_copy(TestData, TestDataHost);

      // Advance clock to next timestep
      ModelClock->advance();
      TimeInstant CurrentTime = ModelClock->getCurrentTime();

      // Update alarm status based on current time
      PeriodAlarm.updateStatus(CurrentTime);

      // Compute the time mean (accumulates internally)
      TimeMeanOp->compute(CurrentTime);

      // Check if alarm is ringing (should ring after last step)
      if (PeriodAlarm.isRinging()) {
         // Mean should now be finalized
         break;
      }
   }

   // Calculate expected mean: average of [BaseValue, BaseValue+1, ...,
   // BaseValue+(NumSteps-1)] For BaseValue=5 and NumSteps=5: avg of [5, 6, 7,
   // 8, 9] = 7.0
   Real Sum = 0.0;
   for (const auto &val : ValuesAtEachStep) {
      Sum += static_cast<Real>(val);
   }
   Real ExpectedMean = Sum / static_cast<Real>(NumSteps);

   // Get result field - output field name includes the Period label
   std::string ResultFieldName = FieldName + "_TimeMean" + PeriodLabel;
   auto ResultField            = Field::get(ResultFieldName);

   // Verify a sample of values. The TimeMeanOp output field is always Real type
   // regardless of input type.
   bool Passed = true;
   if constexpr (Rank == 1) {
      auto ResultData = ResultField->getDataArray<Array1D_t<Real>>();
      auto ResultHost = Kokkos::create_mirror_view(ResultData);
      Kokkos::deep_copy(ResultHost, ResultData);

      for (I4 i = 0; i < std::min(10, Dims[0]); ++i) {
         Real ComputedValue = ResultHost(i);
         if (std::abs(ComputedValue - ExpectedMean) >
             static_cast<Real>(Helper::getTolerance())) {
            Passed = false;
            LOG_ERROR("  At index {}: Expected {}, Got {}", i, ExpectedMean,
                      ComputedValue);
            break;
         }
      }
   } else if constexpr (Rank == 2) {
      auto ResultData = ResultField->getDataArray<Array2D_t<Real>>();
      auto ResultHost = Kokkos::create_mirror_view(ResultData);
      Kokkos::deep_copy(ResultHost, ResultData);

      for (I4 i = 0; i < std::min(5, Dims[0]); ++i) {
         for (I4 j = 0; j < std::min(5, Dims[1]); ++j) {
            Real ComputedValue = ResultHost(i, j);
            if (std::abs(ComputedValue - ExpectedMean) >
                static_cast<Real>(Helper::getTolerance())) {
               Passed = false;
               LOG_ERROR("  At index ({}, {}): Expected {}, Got {}", i, j,
                         ExpectedMean, ComputedValue);
               break;
            }
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      auto ResultData = ResultField->getDataArray<Array3D_t<Real>>();
      auto ResultHost = Kokkos::create_mirror_view(ResultData);
      Kokkos::deep_copy(ResultHost, ResultData);

      for (I4 i = 0; i < std::min(3, Dims[0]); ++i) {
         for (I4 j = 0; j < std::min(3, Dims[1]); ++j) {
            for (I4 k = 0; k < std::min(3, Dims[2]); ++k) {
               Real ComputedValue = ResultHost(i, j, k);
               if (std::abs(ComputedValue - ExpectedMean) >
                   static_cast<Real>(Helper::getTolerance())) {
                  Passed = false;
                  LOG_ERROR("  At index ({}, {}, {}): Expected {}, Got {}", i,
                            j, k, ExpectedMean, ComputedValue);
                  break;
               }
            }
            if (!Passed)
               break;
         }
         if (!Passed)
            break;
      }
   }

   reportTest("TimeMeanOp: " + TypeName, Passed);

   if (!Passed) {
      LOG_ERROR("  Expected mean: {}", ExpectedMean);
      LOG_ERROR("  Period label: {}", PeriodLabel);
   }
}

//===----------------------------------------------------------------------===//
// Main Test Functions
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Test SpatialMaxOp with all array types
void testSpatialMaxOp(const MachEnv *Env, const HorzMesh *Mesh,
                      const VertCoord *VCoord) {

   // 1D arrays - 4 scalar types
   testSpatialMaxOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays - 4 scalar types
   testSpatialMaxOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays - 4 scalar types
   testSpatialMaxOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testSpatialMaxOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test SpatialMinOp with all array types
void testSpatialMinOp(const MachEnv *Env, const HorzMesh *Mesh,
                      const VertCoord *VCoord) {

   // 1D arrays
   testSpatialMinOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testSpatialMinOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays
   testSpatialMinOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testSpatialMinOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays
   testSpatialMinOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testSpatialMinOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testSpatialMinOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test SpatialMeanOp with all array types
void testSpatialMeanOp(const MachEnv *Env, const HorzMesh *Mesh,
                       const VertCoord *VCoord) {

   // 1D arrays
   testSpatialMeanOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays
   testSpatialMeanOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays
   testSpatialMeanOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testSpatialMeanOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test SpatialStdDevOp with all array types
void testSpatialStdDevOp(const MachEnv *Env, const HorzMesh *Mesh,
                         const VertCoord *VCoord) {

   // 1D arrays
   testSpatialStdDevOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays
   testSpatialStdDevOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays
   testSpatialStdDevOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testSpatialStdDevOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test TimeMeanOp with all array types
void testTimeMeanOp(const MachEnv *Env, const HorzMesh *Mesh,
                    const VertCoord *VCoord, Clock *ModelClock) {

   // 1D arrays
   testTimeMeanOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord, ModelClock);

   // 2D arrays
   testTimeMeanOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord, ModelClock);

   // 3D arrays
   testTimeMeanOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord, ModelClock);
   testTimeMeanOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord, ModelClock);
}

//===----------------------------------------------------------------------===//
// Initialization and finalization functions
//===----------------------------------------------------------------------===//

//------------------------------------------------------------------------------
// Initialize needed modules
void initAnalysisTest() {

   I4 Err;

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // First step of time stepper initialization needed for IOstream
   TimeStepper::init1();

   // Get the model clock
   TimeStepper *DefStepper = TimeStepper::getDefault();
   Clock *ModelClock       = DefStepper->getClock();

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   // Initialize streams
   IOStream::init(ModelClock);

   // Initialize fields
   Field::init(ModelClock);

   // Initialize the default halo
   Err = Halo::init();
   if (Err != 0)
      ABORT_ERROR("AnalysisOperatorTest: error initializing default halo");

   // Initialize the default mesh
   HorzMesh::init(ModelClock);

   // Initialize the default vertical coordinate
   VertCoord::init();

   // Initialize tracers
   Tracers::init();

   // Initialize auxiliary state
   AuxiliaryState::init();

   // Initialize equation of state
   Eos::init();

   // Initialize pressure gradient
   PressureGrad::init();

   // Initialize forcing
   Forcing::init();

   // Initialize vertical mixing
   VertMix::init();

   // Initialize tendencies
   Tendencies::init();

   // Initialize vertical advection
   VertAdv::init();

   // Second step of time stepper initialization
   TimeStepper::init2();

   // Initialize ocean state
   Err = OceanState::init();
   if (Err != 0)
      ABORT_ERROR("AnalysisOperatorTest: error initializing default state");

   // Register all analysis operators
   Analysis::init();
}

//------------------------------------------------------------------------------
// Clean-up modules
void finalizeAnalysisTest() {

   Analysis::finalize();
   IOStream::finalize();
   VertMix::destroyInstance();
   Forcing::clear();
   OceanState::clear();
   Tracers::clear();
   AuxiliaryState::clear();
   PressureGrad::clear();
   Tendencies::clear();
   VertAdv::clear();
   VertCoord::clear();
   TimeStepper::clear();
   HorzMesh::clear();
   Field::clear();
   Dimension::clear();
   Halo::clear();
   Decomp::clear();
   MachEnv::removeAll();
}

//===----------------------------------------------------------------------===//
// Main test driver
//===----------------------------------------------------------------------===//

int main(int argc, char *argv[]) {

   int Err = 0;

   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      initAnalysisTest();

      auto DefEnv     = MachEnv::getDefault();
      auto DefStepper = TimeStepper::getDefault();
      auto Mesh       = HorzMesh::getDefault();
      auto VCoord     = VertCoord::getDefault();
      auto ModelClock = DefStepper->getClock();

      testSpatialMaxOp(DefEnv, Mesh, VCoord);

      testSpatialMinOp(DefEnv, Mesh, VCoord);

      testSpatialMeanOp(DefEnv, Mesh, VCoord);

      testSpatialStdDevOp(DefEnv, Mesh, VCoord);

      testTimeMeanOp(DefEnv, Mesh, VCoord, ModelClock);

      if (NumFailed > 0) {
         Err = 1;
         LOG_ERROR("AnalysisOperatorTest failure");
         LOG_ERROR("  Total tests: {}", NumTests);
         LOG_ERROR("  Passed: {}", NumPassed);
         LOG_ERROR("  Failed: {}", NumFailed);
      }

      finalizeAnalysisTest();
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   return Err;
}
