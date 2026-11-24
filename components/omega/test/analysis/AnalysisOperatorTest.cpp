//===-- Test 5.1: Multi-type operator correctness tests ---------*- C++ -*-===//
//
// Comprehensive unit tests for Analysis operators across supported array types
//
//===----------------------------------------------------------------------===//

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
#include "OceanState.h"
#include "TimeStepper.h"
#include "VertAdv.h"
#include "VertCoord.h"
#include "VertMix.h"

#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <type_traits>
#include <vector>

using namespace OMEGA;

// Test result tracking
int NumTests  = 0;
int NumPassed = 0;
int NumFailed = 0;

// Units and standard name given to every test field, which the operators
// must pass through to their outputs
const std::string TestUnits   = "m";
const std::string TestStdName = "test_field_standard_name";

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

   // Relative tolerance for a weighted statistic, whose summation order
   // differs from the host-side check
   static Real getRelTolerance() {
      if constexpr (std::is_same_v<ScalarT, float>) {
         return 1.0e-5;
      } else {
         return 1.0e-12;
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

      // Set ValidMin/ValidMax based on scalar type to avoid bad_any_cast
      ScalarT ValidMin, ValidMax;
      if constexpr (std::is_integral_v<ScalarT>) {
         ValidMin = std::numeric_limits<ScalarT>::lowest();
         ValidMax = std::numeric_limits<ScalarT>::max();
      } else if constexpr (std::is_same_v<ScalarT, float>) {
         ValidMin = -1.0e30f;
         ValidMax = 1.0e30f;
      } else {
         ValidMin = -1.0e30;
         ValidMax = 1.0e30;
      }

      auto TestField = Field::create(
          FieldName, "Test field for multi-type validation", TestUnits,
          TestStdName, ValidMin, ValidMax, 1, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0]);
      TestField->template attachData<ArrayType>(TestData, false);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         TestDataHost(i) = ValueFunc(i);
      }
      deepCopy(TestData, TestDataHost);
   }

   // Create test field for 2D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 2, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4, I4)> ValueFunc) {

      auto DimNames = getDimNames();

      // Set ValidMin/ValidMax based on scalar type to avoid bad_any_cast
      ScalarT ValidMin, ValidMax;
      if constexpr (std::is_integral_v<ScalarT>) {
         ValidMin = std::numeric_limits<ScalarT>::lowest();
         ValidMax = std::numeric_limits<ScalarT>::max();
      } else if constexpr (std::is_same_v<ScalarT, float>) {
         ValidMin = -1.0e30f;
         ValidMax = 1.0e30f;
      } else {
         ValidMin = -1.0e30;
         ValidMax = 1.0e30;
      }

      auto TestField = Field::create(
          FieldName, "Test field for multi-type validation", TestUnits,
          TestStdName, ValidMin, ValidMax, 2, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1]);
      TestField->template attachData<ArrayType>(TestData, false);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Dims[1]; ++j) {
            TestDataHost(i, j) = ValueFunc(i, j);
         }
      }
      deepCopy(TestData, TestDataHost);
   }

   // Create test field for 3D arrays
   template <int R = Rank>
   static typename std::enable_if<R == 3, void>::type
   createField(const std::string &FieldName, const std::vector<I4> &Dims,
               std::function<ScalarT(I4, I4, I4)> ValueFunc) {

      auto DimNames = getDimNames();

      // Set ValidMin/ValidMax based on scalar type to avoid bad_any_cast
      ScalarT ValidMin, ValidMax;
      if constexpr (std::is_integral_v<ScalarT>) {
         ValidMin = std::numeric_limits<ScalarT>::lowest();
         ValidMax = std::numeric_limits<ScalarT>::max();
      } else if constexpr (std::is_same_v<ScalarT, float>) {
         ValidMin = -1.0e30f;
         ValidMax = 1.0e30f;
      } else {
         ValidMin = -1.0e30;
         ValidMax = 1.0e30;
      }

      auto TestField = Field::create(
          FieldName, "Test field for multi-type validation", TestUnits,
          TestStdName, ValidMin, ValidMax, 3, DimNames);

      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1], Dims[2]);
      TestField->template attachData<ArrayType>(TestData, false);

      auto TestDataHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Dims[1]; ++j) {
            for (I4 k = 0; k < Dims[2]; ++k) {
               TestDataHost(i, j, k) = ValueFunc(i, j, k);
            }
         }
      }
      deepCopy(TestData, TestDataHost);
   }
};

//------------------------------------------------------------------------------
// Helper function to report test results
void reportTest(const std::string &TestName, bool Passed) {
   NumTests++;
   if (Passed) {
      NumPassed++;
      LOG_DEBUG("PASS: {}", TestName);
   } else {
      NumFailed++;
      LOG_ERROR("FAIL: {}", TestName);
   }
}

//------------------------------------------------------------------------------
// Host-side weights for checking the spatial statistics, computed from the
// host mesh, mask and pseudo-thickness arrays independently of SpatialWeights:
// area for horizontal fields, mass (area times pseudo-thickness) for layered
// fields with the thickness interpolated to edges and vertices as the
// auxiliary variables do, and half the mass of each adjacent layer for
// interface fields. Inactive entries have zero weight.
struct HostWeights {
   const HorzMesh *Mesh;
   const VertCoord *VCoord;
   HostArray2DReal ThickH;

   HostWeights(const HorzMesh *InMesh, const VertCoord *InVCoord)
       : Mesh(InMesh), VCoord(InVCoord) {
      auto State      = OceanState::getDefault();
      auto ThickField = Field::get(State->PseudoThicknessFldName);
      ThickH = createHostMirrorCopy(ThickField->getDataArray<Array2DReal>());
   }

   bool cellActive(I4 ICell, I4 K) const {
      return VCoord->CellMaskH(ICell, K) > 0;
   }

   // Mass of one layer of one entity
   Real cellLayer(I4 ICell, I4 K) const {
      return cellActive(ICell, K) ? Mesh->AreaCellH(ICell) * ThickH(ICell, K)
                                  : 0;
   }
   Real edgeLayer(I4 IEdge, I4 K) const {
      if (VCoord->EdgeMaskH(IEdge, K) <= 0)
         return 0;
      Real Thick = 0;
      for (I4 J = 0; J < 2; ++J) {
         I4 JCell = Mesh->CellsOnEdgeH(IEdge, J);
         if (cellActive(JCell, K))
            Thick += 0.5 * ThickH(JCell, K);
      }
      return Mesh->DcEdgeH(IEdge) * Mesh->DvEdgeH(IEdge) * Thick;
   }
   Real vertexLayer(I4 IVertex, I4 K) const {
      if (VCoord->VertexMaskH(IVertex, K) <= 0)
         return 0;
      Real Mass = 0;
      for (I4 J = 0; J < Mesh->VertexDegree; ++J) {
         I4 JCell = Mesh->CellsOnVertexH(IVertex, J);
         if (cellActive(JCell, K))
            Mass += Mesh->KiteAreasOnVertexH(IVertex, J) * ThickH(JCell, K);
      }
      return Mass;
   }

   // Weight of one interface of one cell: half of each adjacent layer
   Real cellInterface(I4 ICell, I4 K) const {
      Real Above = (K > 0) ? cellLayer(ICell, K - 1) : 0;
      Real Below = (K < VCoord->NVertLayers) ? cellLayer(ICell, K) : 0;
      return 0.5 * (Above + Below);
   }

   // Area weights of horizontal fields
   Real cellArea(I4 ICell) const {
      return cellActive(ICell, 0) ? Mesh->AreaCellH(ICell) : 0;
   }
   Real edgeArea(I4 IEdge) const {
      return (VCoord->EdgeMaskH(IEdge, 0) > 0)
                 ? Mesh->DcEdgeH(IEdge) * Mesh->DvEdgeH(IEdge)
                 : 0;
   }
};

//------------------------------------------------------------------------------
// Computes the expected weighted mean and standard deviation of a field from
// a visitor that calls back with the weight and value of every local entry.
// Two passes, so the deviations are taken from the global mean.
using WeightedVisitor = std::function<void(std::function<void(Real, Real)>)>;

void expectedStats(const WeightedVisitor &Visit, MPI_Comm Comm, Real &Mean,
                   Real &StdDev) {
   Real SumW = 0, SumWX = 0;
   Visit([&](Real W, Real X) {
      if (W > 0) {
         SumW += W;
         SumWX += W * X;
      }
   });
   Real W = globalSum(SumW, Comm);
   Mean   = globalSum(SumWX, Comm) / W;

   Real SumWDev = 0;
   Visit([&](Real Wt, Real X) {
      if (Wt > 0)
         SumWDev += Wt * (X - Mean) * (X - Mean);
   });
   StdDev = std::sqrt(globalSum(SumWDev, Comm) / W);
}

//------------------------------------------------------------------------------
// Checks a computed statistic against its expected value to a relative
// tolerance
bool checkClose(const std::string &TestName, Real Computed, Real Expected,
                Real RelTol) {
   Real Scale  = std::max(std::abs(Expected), static_cast<Real>(1));
   bool Passed = std::abs(Computed - Expected) <= RelTol * Scale;
   reportTest(TestName, Passed);
   if (!Passed)
      LOG_ERROR("  Expected: {}, Got: {}", Expected, Computed);
   return Passed;
}

//------------------------------------------------------------------------------
// Returns a string attribute of a field, or an empty string if the field has
// no such attribute. Field::create may store an empty attribute or omit it;
// both mean the field has none.
std::string getAttribute(const std::string &FieldName,
                         const std::string &AttName) {
   auto FieldPtr = Field::get(FieldName);
   std::string Value;
   if (FieldPtr && FieldPtr->hasMetadata(AttName))
      FieldPtr->getMetadata(AttName, Value);
   return Value;
}

//------------------------------------------------------------------------------
// Checks that an operator's output field carries the expected units, standard
// name and cell methods
void checkOutputMetadata(const std::string &TestName,
                         const std::string &OutputName,
                         const std::string &ExpectedUnits,
                         const std::string &ExpectedStdName,
                         const std::string &ExpectedCellMethods) {

   std::string Units       = getAttribute(OutputName, "units");
   std::string StdName     = getAttribute(OutputName, "standard_name");
   std::string CellMethods = getAttribute(OutputName, "cell_methods");

   bool Passed = (Units == ExpectedUnits) && (StdName == ExpectedStdName) &&
                 (CellMethods == ExpectedCellMethods);
   reportTest(TestName + " metadata", Passed);

   if (!Passed) {
      LOG_ERROR("  units: '{}', expected '{}'", Units, ExpectedUnits);
      LOG_ERROR("  standard_name: '{}', expected '{}'", StdName,
                ExpectedStdName);
      LOG_ERROR("  cell_methods: '{}', expected '{}'", CellMethods,
                ExpectedCellMethods);
   }
}

//------------------------------------------------------------------------------
// Returns the cell method the spatial operators record for a test field of
// the given rank: 1D fields are horizontal only, higher ranks have depth
std::string expectedSpatialCellMethod(int Rank, const std::string &Method) {
   if (Rank >= 2)
      return "area: depth: " + Method;
   return "area: " + Method;
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
   auto MaxOp =
       AnalysisOpFactory::createOp("SpatialMax", {FieldName}, makeOpConfig());
   MaxOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MaxOp->compute(TestTime);

   // Get result. The operator attaches output as Array1D<ScalarT>, so retrieve
   // with the matching type to avoid reinterpreting bits via
   // static_pointer_cast.
   auto ResultField = Field::get(FieldName + "_SpatialMax");
   auto ResultData =
       ResultField->getDataArray<typename Array1D<ScalarT>::type>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   Real ComputedMax     = static_cast<Real>(ResultHost(0));
   Real ExpectedMaxReal = static_cast<Real>(ExpectedMax);

   // Verify
   bool Passed = (std::abs(ComputedMax - ExpectedMaxReal) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialMaxOp: " + TypeName, Passed);
   checkOutputMetadata("SpatialMaxOp: " + TypeName, FieldName + "_SpatialMax",
                       TestUnits, TestStdName,
                       expectedSpatialCellMethod(Rank, "maximum"));

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
   auto MinOp =
       AnalysisOpFactory::createOp("SpatialMin", {FieldName}, makeOpConfig());
   MinOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MinOp->compute(TestTime);

   // Get result. The operator attaches output as Array1D<ScalarT>, so retrieve
   // with the matching type to avoid reinterpreting bits via
   // static_pointer_cast.
   auto ResultField = Field::get(FieldName + "_SpatialMin");
   auto ResultData =
       ResultField->getDataArray<typename Array1D<ScalarT>::type>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   Real ComputedMin     = static_cast<Real>(ResultHost(0));
   Real ExpectedMinReal = static_cast<Real>(ExpectedMin);

   // Verify
   bool Passed = (std::abs(ComputedMin - ExpectedMinReal) <=
                  static_cast<Real>(Helper::getTolerance()));
   reportTest("SpatialMinOp: " + TypeName, Passed);
   checkOutputMetadata("SpatialMinOp: " + TypeName, FieldName + "_SpatialMin",
                       TestUnits, TestStdName,
                       expectedSpatialCellMethod(Rank, "minimum"));

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

   // Expected weighted mean from the host-side weights: area for the 1D
   // field, mass for the layered ones, over owned cells and active layers
   HostWeights Weights(Mesh, VCoord);
   auto Value = [Value1, Value2](I4 Parity) -> Real {
      return static_cast<Real>((Parity % 2) == 0 ? Value1 : Value2);
   };
   WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
      if constexpr (Rank == 1) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
            Add(Weights.cellArea(i), Value(i));
      } else if constexpr (Rank == 2) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
            for (I4 j = 0; j < VCoord->NVertLayers; ++j)
               Add(Weights.cellLayer(i, j), Value(i + j));
      } else if constexpr (Rank == 3) {
         I4 NTracers = Tracers::getNumTracers();
         for (I4 t = 0; t < NTracers; ++t)
            for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
               for (I4 j = 0; j < VCoord->NVertLayers; ++j)
                  Add(Weights.cellLayer(i, j), Value(t + i + j));
      }
   };
   Real ExpectedMean, ExpectedStdDev;
   expectedStats(Visit, Env->getComm(), ExpectedMean, ExpectedStdDev);

   // Create and compute operator
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, makeOpConfig());
   MeanOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MeanOp->compute(TestTime);

   // Get result. The operator attaches output as Array1DReal (always Real type
   // regardless of input type).
   auto ResultField = Field::get(FieldName + "_SpatialMean");
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   Real ComputedMean = ResultHost(0);

   // Verify
   checkClose("SpatialMeanOp: " + TypeName, ComputedMean, ExpectedMean,
              Helper::getRelTolerance());
   checkOutputMetadata("SpatialMeanOp: " + TypeName, FieldName + "_SpatialMean",
                       TestUnits, TestStdName,
                       expectedSpatialCellMethod(Rank, "mean"));
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

   // Expected weighted standard deviation from the host-side weights
   HostWeights Weights(Mesh, VCoord);
   auto Value = [Value1, Value2](I4 Parity) -> Real {
      return static_cast<Real>((Parity % 2) == 0 ? Value1 : Value2);
   };
   WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
      if constexpr (Rank == 1) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
            Add(Weights.cellArea(i), Value(i));
      } else if constexpr (Rank == 2) {
         for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
            for (I4 j = 0; j < VCoord->NVertLayers; ++j)
               Add(Weights.cellLayer(i, j), Value(i + j));
      } else if constexpr (Rank == 3) {
         I4 NTracers = Tracers::getNumTracers();
         for (I4 t = 0; t < NTracers; ++t)
            for (I4 i = 0; i < Mesh->NCellsOwned; ++i)
               for (I4 j = 0; j < VCoord->NVertLayers; ++j)
                  Add(Weights.cellLayer(i, j), Value(t + i + j));
      }
   };
   Real ExpectedMean, ExpectedStdDev;
   expectedStats(Visit, Env->getComm(), ExpectedMean, ExpectedStdDev);

   // SpatialStdDevOp requires a pre-existing _SpatialMean field for the input.
   // Create and compute a SpatialMeanOp first so that field is registered.
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, makeOpConfig());
   MeanOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MeanOp->compute(TestTime);

   // Now create and compute the StdDev operator
   auto StdDevOp = AnalysisOpFactory::createOp("SpatialStdDev", {FieldName},
                                               makeOpConfig());
   StdDevOp->initialize(Env, Mesh, VCoord, makeOpConfig());
   StdDevOp->compute(TestTime);

   // Get result. The operator attaches output as Array1DReal (always Real type
   // regardless of input type).
   auto ResultField = Field::get(FieldName + "_SpatialStdDev");
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   Real ComputedStdDev = ResultHost(0);

   // Verify
   checkClose("SpatialStdDevOp: " + TypeName, ComputedStdDev, ExpectedStdDev,
              Helper::getRelTolerance());
   checkOutputMetadata("SpatialStdDevOp: " + TypeName,
                       FieldName + "_SpatialStdDev", TestUnits, TestStdName,
                       expectedSpatialCellMethod(Rank, "standard_deviation"));
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
      auto TestDataHost = createHostMirrorCopy(TestData);

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

      deepCopy(TestData, TestDataHost);

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
      auto ResultHost = createHostMirrorCopy(ResultData);

      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
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
      auto ResultHost = createHostMirrorCopy(ResultData);

      // Active layers [MinLayerCell, MaxLayerCell] must equal the analytic
      // time-mean; inactive layers must retain FillValueReal (the refactored
      // TimeMeanOp only writes active layers, preserving fill values).
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         const I4 KMin = VCoord->MinLayerCellH(i);
         const I4 KMax = VCoord->MaxLayerCellH(i);
         for (I4 j = 0; j < VCoord->NVertLayers; ++j) {
            Real ComputedValue = ResultHost(i, j);
            bool Active        = (j >= KMin && j <= KMax);
            Real Expected      = Active ? ExpectedMean : FillValueReal;
            if (std::abs(ComputedValue - Expected) >
                static_cast<Real>(Helper::getTolerance())) {
               Passed = false;
               LOG_ERROR("  At index ({}, {}) [active={}]: Expected {}, Got {}",
                         i, j, Active, Expected, ComputedValue);
               break;
            }
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      auto ResultData = ResultField->getDataArray<Array3D_t<Real>>();
      auto ResultHost = createHostMirrorCopy(ResultData);

      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Mesh->NCellsOwned; ++j) {
            const I4 KMin = VCoord->MinLayerCellH(j);
            const I4 KMax = VCoord->MaxLayerCellH(j);
            for (I4 k = 0; k < VCoord->NVertLayers; ++k) {
               Real ComputedValue = ResultHost(i, j, k);
               bool Active        = (k >= KMin && k <= KMax);
               Real Expected      = Active ? ExpectedMean : FillValueReal;
               if (std::abs(ComputedValue - Expected) >
                   static_cast<Real>(Helper::getTolerance())) {
                  Passed = false;
                  LOG_ERROR("  At index ({}, {}, {}) [active={}]: Expected {}, "
                            "Got {}",
                            i, j, k, Active, Expected, ComputedValue);
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
   checkOutputMetadata("TimeMeanOp: " + TypeName,
                       FieldName + "_TimeMean" + PeriodLabel, TestUnits,
                       TestStdName, "time: mean");

   if (!Passed) {
      LOG_ERROR("  Expected mean: {}", ExpectedMean);
      LOG_ERROR("  Period label: {}", PeriodLabel);
   }
}

//------------------------------------------------------------------------------
// Verify that TimeMeanOp includes the bottom interface for NVertLayersP1 fields
void testTimeMeanOpInterfaceField(const MachEnv *Env, const HorzMesh *Mesh,
                                  const VertCoord *VCoord, Clock *ModelClock) {

   using ArrayType = Array2DR8;

   const std::string FieldName = "TestFieldTimeMeanInterface";
   const std::vector<I4> Dims  = {Mesh->NCellsSize, VCoord->NVertLayersP1};
   auto TestField =
       Field::create(FieldName, "Test interface field for TimeMeanOp", "", "",
                     -1.0e30, 1.0e30, 2, {"NCells", "NVertLayersP1"});
   ArrayType TestData(FieldName + "_data", Dims[0], Dims[1]);
   TestField->attachData<ArrayType>(TestData);

   const int NumSteps        = 5;
   TimeInterval StepInterval = ModelClock->getTimeStep();
   R8 StepSeconds;
   StepInterval.get(StepSeconds, TimeUnits::Seconds);
   R8 PeriodSeconds = StepSeconds * NumSteps;
   TimeInterval PeriodInterval(PeriodSeconds, TimeUnits::Seconds);
   std::string PeriodLabel =
       std::to_string(static_cast<int>(PeriodSeconds)) + "seconds";

   auto TimeMeanOp = AnalysisOpFactory::createOp(
       "TimeMean", {FieldName}, makeOpConfig(opParam("Period", PeriodLabel)));
   TimeMeanOp->initialize(Env, Mesh, VCoord, Config{});

   TimeInstant StartTime = ModelClock->getCurrentTime();
   Alarm PeriodAlarm("TestPeriodAlarmInterface", PeriodInterval, StartTime);
   TimeMeanOp->setPeriodAlarm(&PeriodAlarm);

   for (int Step = 0; Step < NumSteps; ++Step) {
      const Real CurrentValue = 5.0_Real + static_cast<Real>(Step);
      auto TestDataHost       = createHostMirrorCopy(TestData);
      for (I4 ICell = 0; ICell < Dims[0]; ++ICell) {
         for (I4 K = 0; K < Dims[1]; ++K) {
            TestDataHost(ICell, K) = CurrentValue;
         }
      }
      deepCopy(TestData, TestDataHost);

      ModelClock->advance();
      TimeInstant CurrentTime = ModelClock->getCurrentTime();
      PeriodAlarm.updateStatus(CurrentTime);
      TimeMeanOp->compute(CurrentTime);
   }

   const Real ExpectedMean = 7.0_Real;
   auto ResultField        = Field::get(FieldName + "_TimeMean" + PeriodLabel);
   auto ResultData         = ResultField->getDataArray<Array2D_t<Real>>();
   auto ResultHost         = createHostMirrorCopy(ResultData);

   bool Passed = true;
   for (I4 ICell = 0; ICell < Mesh->NCellsOwned && Passed; ++ICell) {
      const I4 KMin = VCoord->MinLayerCellH(ICell);
      const I4 KMax = VCoord->MaxLayerCellH(ICell) + 1;
      for (I4 K = 0; K < VCoord->NVertLayersP1; ++K) {
         const bool Active   = (K >= KMin && K <= KMax);
         const Real Expected = Active ? ExpectedMean : FillValueReal;
         const Real Actual   = ResultHost(ICell, K);
         if (std::abs(Actual - Expected) > 1.0e-8) {
            Passed = false;
            LOG_ERROR("  Interface field at ({}, {}) [active={}]: Expected {}, "
                      "Got {}",
                      ICell, K, Active, Expected, Actual);
            break;
         }
      }
   }

   reportTest("TimeMeanOp: NVertLayersP1 interface field", Passed);
}

//------------------------------------------------------------------------------
// Creates a 2D Real field on the given entities with the given vertical
// dimension, filled with a pattern of values on active owned entries and
// NaN everywhere else (inactive layers and halo entities), so that a
// statistic that reads an inactive or halo entry comes out NaN. Returns the
// host copy of the values.
HostArray2DReal createEntityField(const std::string &FieldName,
                                  const std::string &HorizDim,
                                  const std::string &VertDim, I4 NEntities,
                                  I4 NOwned, I4 NVert,
                                  std::function<bool(I4, I4)> IsActive,
                                  std::function<Real(I4, I4)> ValueFunc) {

   std::vector<std::string> DimNames = {HorizDim, VertDim};
   auto TestField = Field::create(FieldName, "Weighted statistics test field",
                                  "m", "", -1.0e30, 1.0e30, 2, DimNames);
   Array2DReal Data(FieldName + "_data", NEntities, NVert);
   TestField->attachData<Array2DReal>(Data);

   auto DataHost  = Kokkos::create_mirror_view(Data);
   const Real NaN = std::numeric_limits<Real>::quiet_NaN();
   for (I4 I = 0; I < NEntities; ++I) {
      for (I4 K = 0; K < NVert; ++K) {
         bool Use       = (I < NOwned) && IsActive(I, K);
         DataHost(I, K) = Use ? ValueFunc(I, K) : NaN;
      }
   }
   Kokkos::deep_copy(Data, DataHost);
   return DataHost;
}

//------------------------------------------------------------------------------
// Computes the mean and standard deviation of a field with the operators and
// checks them against the expected values
void checkMeanAndStdDev(const std::string &TestName,
                        const std::string &FieldName, const MachEnv *Env,
                        const HorzMesh *Mesh, const VertCoord *VCoord,
                        Real ExpectedMean, Real ExpectedStdDev, Real RelTol) {

   Config EmptyConfig;
   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, EmptyConfig);
   MeanOp->initialize(Env, Mesh, VCoord, EmptyConfig);
   auto StdDevOp =
       AnalysisOpFactory::createOp("SpatialStdDev", {FieldName}, EmptyConfig);
   StdDevOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   TimeInstant TestTime;
   MeanOp->compute(TestTime);
   StdDevOp->compute(TestTime);

   auto MeanHost = createHostMirrorCopy(
       Field::get(FieldName + "_SpatialMean")->getDataArray<Array1DReal>());
   auto StdDevHost = createHostMirrorCopy(
       Field::get(FieldName + "_SpatialStdDev")->getDataArray<Array1DReal>());

   checkClose(TestName + " mean", MeanHost(0), ExpectedMean, RelTol);
   checkClose(TestName + " std dev", StdDevHost(0), ExpectedStdDev, RelTol);
}

//------------------------------------------------------------------------------
// Tests the weighted mean and standard deviation on edges, vertices and cell
// interfaces, and on a horizontal edge field, against host-side weights.
// Inactive and halo entries hold NaN, so a statistic that touched them would
// fail. A constant field must give its value and a zero standard deviation
// whatever the weights.
void testWeightedStats(const MachEnv *Env, const HorzMesh *Mesh,
                       const VertCoord *VCoord) {

   HostWeights Weights(Mesh, VCoord);
   MPI_Comm Comm      = Env->getComm();
   const Real RelTol  = 1.0e-12;
   const I4 NVertLyrs = VCoord->NVertLayers;
   auto Pattern       = [](I4 I, I4 K) -> Real {
      return static_cast<Real>(10 + (7 * I + 3 * K) % 5);
   };
   auto Constant = [](I4 I, I4 K) -> Real { return static_cast<Real>(7.5); };
   Real ExpectedMean, ExpectedStdDev;

   // Layered field on edges, mass-weighted with the edge thickness
   {
      auto Active = [&](I4 E, I4 K) { return VCoord->EdgeMaskH(E, K) > 0; };
      createEntityField("TestWeightedEdge", "NEdges", "NVertLayers",
                        Mesh->NEdgesSize, Mesh->NEdgesOwned, NVertLyrs, Active,
                        Pattern);
      WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
         for (I4 E = 0; E < Mesh->NEdgesOwned; ++E)
            for (I4 K = 0; K < NVertLyrs; ++K)
               Add(Weights.edgeLayer(E, K), Pattern(E, K));
      };
      expectedStats(Visit, Comm, ExpectedMean, ExpectedStdDev);
      checkMeanAndStdDev("Weighted stats: edge layers", "TestWeightedEdge", Env,
                         Mesh, VCoord, ExpectedMean, ExpectedStdDev, RelTol);
   }

   // Layered field on vertices, mass-weighted with the kite areas
   {
      auto Active = [&](I4 V, I4 K) { return VCoord->VertexMaskH(V, K) > 0; };
      createEntityField("TestWeightedVertex", "NVertices", "NVertLayers",
                        Mesh->NVerticesSize, Mesh->NVerticesOwned, NVertLyrs,
                        Active, Pattern);
      WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
         for (I4 V = 0; V < Mesh->NVerticesOwned; ++V)
            for (I4 K = 0; K < NVertLyrs; ++K)
               Add(Weights.vertexLayer(V, K), Pattern(V, K));
      };
      expectedStats(Visit, Comm, ExpectedMean, ExpectedStdDev);
      checkMeanAndStdDev("Weighted stats: vertex layers", "TestWeightedVertex",
                         Env, Mesh, VCoord, ExpectedMean, ExpectedStdDev,
                         RelTol);
   }

   // Interface field on cells, half the mass of each adjacent layer. The
   // bottom interface of every active column takes part.
   {
      auto Active = [&](I4 C, I4 K) {
         bool Above = (K > 0) && Weights.cellActive(C, K - 1);
         bool Below = (K < NVertLyrs) && Weights.cellActive(C, K);
         return Above || Below;
      };
      createEntityField("TestWeightedInterface", "NCells", "NVertLayersP1",
                        Mesh->NCellsSize, Mesh->NCellsOwned, NVertLyrs + 1,
                        Active, Pattern);
      WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
         for (I4 C = 0; C < Mesh->NCellsOwned; ++C)
            for (I4 K = 0; K <= NVertLyrs; ++K)
               Add(Weights.cellInterface(C, K), Pattern(C, K));
      };
      expectedStats(Visit, Comm, ExpectedMean, ExpectedStdDev);
      checkMeanAndStdDev("Weighted stats: cell interfaces",
                         "TestWeightedInterface", Env, Mesh, VCoord,
                         ExpectedMean, ExpectedStdDev, RelTol);
   }

   // Horizontal field on edges, area-weighted
   {
      std::vector<std::string> DimNames = {"NEdges"};
      auto TestField =
          Field::create("TestWeightedEdge1D", "Weighted statistics test field",
                        "m", "", -1.0e30, 1.0e30, 1, DimNames);
      Array1DReal Data("TestWeightedEdge1D_data", Mesh->NEdgesSize);
      TestField->attachData<Array1DReal>(Data);
      auto DataHost  = Kokkos::create_mirror_view(Data);
      const Real NaN = std::numeric_limits<Real>::quiet_NaN();
      for (I4 E = 0; E < Mesh->NEdgesSize; ++E) {
         bool Use    = (E < Mesh->NEdgesOwned) && (Weights.edgeArea(E) > 0);
         DataHost(E) = Use ? Pattern(E, 0) : NaN;
      }
      Kokkos::deep_copy(Data, DataHost);
      WeightedVisitor Visit = [&](std::function<void(Real, Real)> Add) {
         for (I4 E = 0; E < Mesh->NEdgesOwned; ++E)
            Add(Weights.edgeArea(E), Pattern(E, 0));
      };
      expectedStats(Visit, Comm, ExpectedMean, ExpectedStdDev);
      checkMeanAndStdDev("Weighted stats: horizontal edges",
                         "TestWeightedEdge1D", Env, Mesh, VCoord, ExpectedMean,
                         ExpectedStdDev, RelTol);
   }

   // A constant layered field on cells: the mean is the constant and the
   // standard deviation zero, whatever the weights
   {
      auto Active = [&](I4 C, I4 K) { return Weights.cellActive(C, K); };
      createEntityField("TestWeightedConstant", "NCells", "NVertLayers",
                        Mesh->NCellsSize, Mesh->NCellsOwned, NVertLyrs, Active,
                        Constant);
      checkMeanAndStdDev("Weighted stats: constant field",
                         "TestWeightedConstant", Env, Mesh, VCoord, 7.5, 0.0,
                         RelTol);
   }
}

//------------------------------------------------------------------------------
// Tests that metadata is inherited along an operator chain: a time mean of a
// spatial mean keeps the units and standard name and lists both reductions in
// order in cell_methods. Also tests that a field without units or a standard
// name yields outputs without them rather than with invented ones.
void testInheritedMetadata(const MachEnv *Env, const HorzMesh *Mesh,
                           const VertCoord *VCoord, Clock *ModelClock) {

   using Helper = TestHelper<Array2DR8>;

   std::vector<I4> Dims = Helper::getDims(Mesh, VCoord);
   Config EmptyConfig;

   // Chain: SpatialMean then TimeMean of a 2D field with units "m"
   std::string FieldName = "TestFieldChain";
   Helper::createField(FieldName, Dims,
                       [](I4 i, I4 j) -> R8 { return static_cast<R8>(1); });

   auto MeanOp =
       AnalysisOpFactory::createOp("SpatialMean", {FieldName}, EmptyConfig);
   MeanOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   std::string MeanName = FieldName + "_SpatialMean";
   std::string Period   = "1day";
   auto ChainOp         = AnalysisOpFactory::createOp(
       "TimeMean", {MeanName}, makeOpConfig(opParam("Period", Period)));
   ChainOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   checkOutputMetadata("TimeMean of SpatialMean",
                       MeanName + "_TimeMean" + Period, TestUnits, TestStdName,
                       "area: depth: mean time: mean");

   // A field with no units and no standard name: the outputs have none either
   std::string BareName              = "TestFieldBare";
   std::vector<std::string> DimNames = Helper::getDimNames();
   auto BareField = Field::create(BareName, "Test field without units", "", "",
                                  -1.0e30, 1.0e30, 2, DimNames);
   Array2DR8 BareData(BareName + "_data", Dims[0], Dims[1]);
   BareField->attachData<Array2DR8>(BareData);

   auto BareMaxOp =
       AnalysisOpFactory::createOp("SpatialMax", {BareName}, EmptyConfig);
   BareMaxOp->initialize(Env, Mesh, VCoord, EmptyConfig);

   checkOutputMetadata("SpatialMax of field without units",
                       BareName + "_SpatialMax", "", "",
                       "area: depth: maximum");
}

//------------------------------------------------------------------------------
// Tests the units operators derive along a chain like the MOC's: a velocity
// times an area is a volume flux, and multiplying the flux by 1e-6 to get
// Sverdrups scales its units by 1e6. Also tests that a product with a field
// without units has no units, and that a negative scalar drops the standard
// name of the quantity it reverses.
void testDerivedUnits(const MachEnv *Env, const HorzMesh *Mesh,
                      const VertCoord *VCoord) {

   std::vector<std::string> DimNames = {"NCells"};
   auto createUnitsField             = [&](const std::string &Name,
                               const std::string &Units,
                               const std::string &StdName) {
      auto NewField = Field::create(Name, "Derived units test field", Units,
                                                StdName, -1.0e30, 1.0e30, 1, DimNames);
      Array1DReal Data(Name + "_data", Mesh->NCellsSize);
      NewField->attachData<Array1DReal>(Data);
   };
   createUnitsField("TestUnitsVelocity", "m s-1", "upward_sea_water_velocity");
   createUnitsField("TestUnitsArea", "m2", "cell_area");
   createUnitsField("TestUnitsBare", "", "");

   Config EmptyConfig;
   auto scalarConfig = [](const std::string &Scalar) {
      return makeOpConfig(opParam("Scalar", Scalar));
   };

   // A velocity times an area is a volume flux with no standard name
   std::string FluxName = "TestUnitsVelocity_BinaryMultiply(TestUnitsArea)";
   auto FluxOp          = AnalysisOpFactory::createOp(
       "BinaryMultiply", {"TestUnitsVelocity", "TestUnitsArea"}, EmptyConfig);
   checkOutputMetadata("BinaryMultiply of velocity and area", FluxName,
                       "m3 s-1", "", "");

   // Values multiplied by 1e-6 are in units of 1e6 m3 s-1
   auto SvOp = AnalysisOpFactory::createOp("ScalarMultiply", {FluxName},
                                           scalarConfig("1.0e-6"));
   checkOutputMetadata("ScalarMultiply of flux to Sverdrups",
                       FluxName + "_ScalarMultiply(1.0e-6)", "1e+06 m3 s-1", "",
                       "");

   // A scalar of 1 keeps the units and standard name; -1 keeps the units but
   // reverses the quantity, so it drops the standard name
   auto OneOp = AnalysisOpFactory::createOp(
       "ScalarMultiply", {"TestUnitsVelocity"}, scalarConfig("1"));
   checkOutputMetadata("ScalarMultiply by 1",
                       "TestUnitsVelocity_ScalarMultiply(1)", "m s-1",
                       "upward_sea_water_velocity", "");
   auto NegOp = AnalysisOpFactory::createOp(
       "ScalarMultiply", {"TestUnitsVelocity"}, scalarConfig("-1"));
   checkOutputMetadata("ScalarMultiply by -1",
                       "TestUnitsVelocity_ScalarMultiply(-1)", "m s-1", "", "");

   // A product with a field without units has no units, rather than the
   // units of the other field
   auto BareOp = AnalysisOpFactory::createOp(
       "BinaryMultiply", {"TestUnitsVelocity", "TestUnitsBare"}, EmptyConfig);
   checkOutputMetadata("BinaryMultiply with a field without units",
                       "TestUnitsVelocity_BinaryMultiply(TestUnitsBare)", "",
                       "", "");
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
// Verify the TimeMeanOp SamplingInterval gate: with "2steps" over a 4-step
// period, only every other step is accumulated, so the mean equals the simple
// average of the SAMPLED values (not all steps). Also exercises the three-form
// equivalence ("2steps", "2perPeriod", and a raw TimeInterval string all
// normalize to the same 2-timestep interval) and the every-step default
// (bit-for-bit with the ungated path).
void testTimeMeanOpSamplingInterval(const MachEnv *Env, const HorzMesh *Mesh,
                                    const VertCoord *VCoord,
                                    Clock *ModelClock) {

   using ArrayType = Array2DR8;

   // 4-step period; timestep comes from the model clock
   const int NumSteps        = 4;
   TimeInterval StepInterval = ModelClock->getTimeStep();
   R8 StepSeconds;
   StepInterval.get(StepSeconds, TimeUnits::Seconds);
   const R8 PeriodSeconds = StepSeconds * NumSteps;
   TimeInterval PeriodInterval(PeriodSeconds, TimeUnits::Seconds);
   const std::string PeriodLabel =
       std::to_string(static_cast<int>(PeriodSeconds)) + "seconds";

   // A raw TimeInterval string equal to two timesteps, for form equivalence
   const std::string TwoStepSecondsStr =
       std::to_string(static_cast<int>(2 * StepSeconds)) + "seconds";

   // Field values per step: 10, 20, 30, 40. The SampleAlarm is anchored to the
   // period start and fires on the first step at/after PeriodStart+2*dt, i.e.
   // the 2nd and 4th steps, so the gated mean is (20 + 40) / 2 = 30. The
   // ungated (every-step) mean would be (10+20+30+40)/4 = 25.
   const std::vector<Real> StepValues = {10.0_Real, 20.0_Real, 30.0_Real,
                                         40.0_Real};
   const Real GatedExpected           = 30.0_Real;
   const Real UngatedExpected         = 25.0_Real;

   // Helper that runs one 4-step period with the given SamplingInterval config
   // (empty => every-step default) and returns the resulting mean at cell 0,
   // layer 0. Each sub-run re-anchors its PeriodAlarm/SampleAlarm to the
   // clock's current time at the sub-run start, so the 2-step sampling pattern
   // is preserved regardless of the advancing clock.
   auto runPeriod = [&](const std::string &FieldSuffix,
                        const std::string &SamplingInterval) -> Real {
      const std::string FieldName = "TestFieldTMSI_" + FieldSuffix;
      const std::vector<I4> Dims  = {Mesh->NCellsSize, VCoord->NVertLayers};
      auto TestField              = Field::create(
          FieldName, "Test field for SamplingInterval", TestUnits, TestStdName,
          -1.0e30, 1.0e30, 2, {"NCells", "NVertLayers"});
      ArrayType TestData(FieldName + "_data", Dims[0], Dims[1]);
      TestField->attachData<ArrayType>(TestData);

      Config OpConfig = makeOpConfig(opParam("Period", PeriodLabel));
      if (!SamplingInterval.empty())
         OpConfig.add("SamplingInterval", SamplingInterval);

      auto Op = AnalysisOpFactory::createOp("TimeMean", {FieldName}, OpConfig);

      // The operator anchors its SampleAlarm to the absolute clock start
      // (ModelClock->getStartTime()) with a 2-timestep stride, so its sample
      // boundaries fall at start + {0, 2*dt, 4*dt, ...}. This shared clock has
      // already been advanced an uncontrolled number of steps by the preceding
      // TimeMean tests, so a sub-run that begins on an ODD step offset relative
      // to that grid would be out of phase: samples would land on the 1st and
      // 3rd steps instead of the intended 2nd and 4th. Snap each sub-run to an
      // even-step offset so the "fires on the 2nd and 4th steps" assumption
      // holds and the gated mean is (20 + 40) / 2 = 30.
      {
         TimeInterval SinceStart =
             ModelClock->getCurrentTime() - ModelClock->getStartTime();
         R8 SinceStartSeconds;
         SinceStart.get(SinceStartSeconds, TimeUnits::Seconds);
         const int StepsSinceStart =
             static_cast<int>(std::llround(SinceStartSeconds / StepSeconds));
         if (StepsSinceStart % 2 != 0)
            ModelClock->advance();
      }

      TimeInstant StartTime = ModelClock->getCurrentTime();
      Alarm PeriodAlarm("TMSIPeriodAlarm_" + FieldSuffix, PeriodInterval,
                        StartTime);
      Op->initialize(Env, Mesh, VCoord, Config{});
      Op->setPeriodAlarm(&PeriodAlarm);

      for (int Step = 0; Step < NumSteps; ++Step) {
         auto TestDataHost = createHostMirrorCopy(TestData);
         for (I4 ICell = 0; ICell < Dims[0]; ++ICell)
            for (I4 K = 0; K < Dims[1]; ++K)
               TestDataHost(ICell, K) = StepValues[Step];
         deepCopy(TestData, TestDataHost);

         ModelClock->advance();
         TimeInstant CurrentTime = ModelClock->getCurrentTime();
         PeriodAlarm.updateStatus(CurrentTime);
         Op->compute(CurrentTime);
      }

      auto ResultField = Field::get(FieldName + "_TimeMean" + PeriodLabel);
      auto ResultData  = ResultField->getDataArray<Array2D_t<Real>>();
      auto ResultHost  = createHostMirrorCopy(ResultData);
      return ResultHost(0, 0);
   };

   bool Passed = true;

   // 1. Gated mean with "2steps"
   const Real StepsMean = runPeriod("steps", "2steps");
   if (std::abs(StepsMean - GatedExpected) > 1.0e-8) {
      Passed = false;
      LOG_ERROR("  SamplingInterval '2steps': expected {}, got {}",
                GatedExpected, StepsMean);
   }

   // 2. Three-form equivalence: "2perPeriod" and a raw TimeInterval string
   //    must give the identical gated mean as "2steps".
   const Real PerPeriodMean = runPeriod("perperiod", "2perPeriod");
   if (std::abs(PerPeriodMean - GatedExpected) > 1.0e-8) {
      Passed = false;
      LOG_ERROR("  SamplingInterval '2perPeriod': expected {}, got {}",
                GatedExpected, PerPeriodMean);
   }

   const Real RawIntervalMean = runPeriod("rawinterval", TwoStepSecondsStr);
   if (std::abs(RawIntervalMean - GatedExpected) > 1.0e-8) {
      Passed = false;
      LOG_ERROR("  SamplingInterval '{}': expected {}, got {}",
                TwoStepSecondsStr, GatedExpected, RawIntervalMean);
   }

   // 3. Every-step default (empty SamplingInterval) is bit-for-bit with the
   //    ungated accumulation of all four steps.
   const Real DefaultMean = runPeriod("default", "");
   if (std::abs(DefaultMean - UngatedExpected) > 1.0e-8) {
      Passed = false;
      LOG_ERROR("  Default (every-step) mean: expected {}, got {}",
                UngatedExpected, DefaultMean);
   }

   reportTest("TimeMeanOp: SamplingInterval gated mean + form equivalence",
              Passed);
}

//------------------------------------------------------------------------------
// Runs one 1Day averaging period of a new TimeMeanOp with the given
// SamplingInterval and returns the mean at cell 0, layer 0. The field holds
// the step number (1, 2, ...) on each step, and the period starts on a 6-hour
// boundary of the clock, which is what the sample alarm is anchored to.
Real runOneDaySampledMean(const std::string &FieldSuffix,
                          const std::string &SamplingInterval,
                          const MachEnv *Env, const HorzMesh *Mesh,
                          const VertCoord *VCoord, Clock *ModelClock) {

   using ArrayType = Array2DR8;

   const std::string FieldName = "TestFieldTMDay_" + FieldSuffix;
   const std::vector<I4> Dims  = {Mesh->NCellsSize, VCoord->NVertLayers};
   auto TestField              = Field::create(
       FieldName, "Test field for a 1Day SamplingInterval", TestUnits,
       TestStdName, -1.0e30, 1.0e30, 2, {"NCells", "NVertLayers"});
   ArrayType TestData(FieldName + "_data", Dims[0], Dims[1]);
   TestField->attachData<ArrayType>(TestData);

   Config OpConfig = makeOpConfig(opParam("Period", std::string("1Day")));
   OpConfig.add("SamplingInterval", SamplingInterval);
   auto Op = AnalysisOpFactory::createOp("TimeMean", {FieldName}, OpConfig);

   R8 StepSeconds;
   ModelClock->getTimeStep().get(StepSeconds, TimeUnits::Seconds);
   const R8 SixHourSeconds = 6.0 * SECONDS_PER_HOUR;
   while (true) {
      TimeInterval SinceStart =
          ModelClock->getCurrentTime() - ModelClock->getStartTime();
      R8 SinceStartSeconds;
      SinceStart.get(SinceStartSeconds, TimeUnits::Seconds);
      if (std::fmod(SinceStartSeconds, SixHourSeconds) == 0.0)
         break;
      ModelClock->advance();
   }

   TimeInstant StartTime = ModelClock->getCurrentTime();
   Alarm PeriodAlarm("TMDayPeriodAlarm_" + FieldSuffix,
                     TimeInterval(1, TimeUnits::Days), StartTime);
   Op->initialize(Env, Mesh, VCoord, Config{});
   Op->setPeriodAlarm(&PeriodAlarm);

   const I4 NumSteps =
       static_cast<I4>(std::llround(SECONDS_PER_DAY / StepSeconds));
   for (I4 Step = 1; Step <= NumSteps; ++Step) {
      auto TestDataHost = createHostMirrorCopy(TestData);
      for (I4 ICell = 0; ICell < Dims[0]; ++ICell)
         for (I4 K = 0; K < Dims[1]; ++K)
            TestDataHost(ICell, K) = static_cast<Real>(Step);
      deepCopy(TestData, TestDataHost);

      ModelClock->advance();
      TimeInstant CurrentTime = ModelClock->getCurrentTime();
      PeriodAlarm.updateStatus(CurrentTime);
      Op->compute(CurrentTime);
   }

   auto ResultField = Field::get(FieldName + "_TimeMean1Day");
   auto ResultData  = ResultField->getDataArray<Array2D_t<Real>>();
   auto ResultHost  = createHostMirrorCopy(ResultData);
   return ResultHost(0, 0);
}

//------------------------------------------------------------------------------
// Verify that "<n>perPeriod" divides a period given in days rather than
// truncating it: "4perPeriod" of a 1Day period samples every 6 hours, giving
// the same mean as "6hour". The samples fall on the steps S, 2S, 3S and 4S,
// where S is the number of steps in 6 hours, so the mean is 2.5 S.
void testTimeMeanOpSamplingPerDay(const MachEnv *Env, const HorzMesh *Mesh,
                                  const VertCoord *VCoord, Clock *ModelClock) {

   R8 StepSeconds;
   ModelClock->getTimeStep().get(StepSeconds, TimeUnits::Seconds);
   const R8 StepsPerSixHours = 6.0 * SECONDS_PER_HOUR / StepSeconds;
   if (StepsPerSixHours != std::floor(StepsPerSixHours)) {
      LOG_INFO("  Skipping the 1Day SamplingInterval test: the time step "
               "does not divide 6 hours");
      return;
   }
   const Real Expected = 2.5_Real * static_cast<Real>(StepsPerSixHours);

   bool Passed = true;
   for (const auto &[Suffix, Interval] :
        std::vector<std::pair<std::string, std::string>>{
            {"perperiod", "4perPeriod"}, {"rawinterval", "6hour"}}) {
      const Real Mean =
          runOneDaySampledMean(Suffix, Interval, Env, Mesh, VCoord, ModelClock);
      if (std::abs(Mean - Expected) > 1.0e-8) {
         Passed = false;
         LOG_ERROR("  SamplingInterval '{}' over 1Day: expected {}, got {}",
                   Interval, Expected, Mean);
      }
   }

   reportTest("TimeMeanOp: SamplingInterval perPeriod over a 1Day period",
              Passed);
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

   testTimeMeanOpInterfaceField(Env, Mesh, VCoord, ModelClock);

   testTimeMeanOpSamplingInterval(Env, Mesh, VCoord, ModelClock);

   testTimeMeanOpSamplingPerDay(Env, Mesh, VCoord, ModelClock);
}

//------------------------------------------------------------------------------
// Template for testing ScalarMultiplyOp with any array type
template <typename ArrayType>
void testScalarMultiplyOpType(const std::string &TypeName, const MachEnv *Env,
                              const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldScalarMult_" + TypeName;

   // Create test field with known values
   ScalarT BaseValue = static_cast<ScalarT>(10);
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
   // Test scalar multiplication
   std::string ScalarStr = "2.5";
   Real Scalar           = 2.5;
   ScalarT ExpectedValue =
       static_cast<ScalarT>(static_cast<Real>(BaseValue) * Scalar);
   // Create and compute operator
   Config OpConfig;
   OpConfig.add("Scalar", Scalar);

   auto MultOp =
       AnalysisOpFactory::createOp("ScalarMultiply", {FieldName}, OpConfig);
   MultOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MultOp->compute(TestTime);

   // Get result
   auto ResultField =
       Field::get(FieldName + "_ScalarMultiply(" + ScalarStr + ")");
   auto ResultData = ResultField->getDataArray<ArrayType>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   // Verify all owned cells
   bool Passed = true;
   if constexpr (Rank == 1) {
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         if (std::abs(static_cast<Real>(ResultHost(i)) -
                      static_cast<Real>(ExpectedValue)) >
             static_cast<Real>(Helper::getTolerance())) {
            Passed = false;
            break;
         }
      }
   } else if constexpr (Rank == 2) {
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
              ++j) {
            if (std::abs(static_cast<Real>(ResultHost(i, j)) -
                         static_cast<Real>(ExpectedValue)) >
                static_cast<Real>(Helper::getTolerance())) {
               Passed = false;
               break;
            }
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Mesh->NCellsOwned; ++j) {
            for (I4 k = VCoord->MinLayerCellH(j); k <= VCoord->MaxLayerCellH(j);
                 ++k) {
               if (std::abs(static_cast<Real>(ResultHost(i, j, k)) -
                            static_cast<Real>(ExpectedValue)) >
                   static_cast<Real>(Helper::getTolerance())) {
                  Passed = false;
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

   reportTest("ScalarMultiplyOp: " + TypeName, Passed);
   checkOutputMetadata("ScalarMultiplyOp: " + TypeName,
                       FieldName + "_ScalarMultiply(" + ScalarStr + ")",
                       "0.4 m", TestStdName, "");

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Scalar: {}", static_cast<Real>(ExpectedValue),
                Scalar);
   }
}

//------------------------------------------------------------------------------
// Test ScalarMultiplyOp with all array types
void testScalarMultiplyOp(const MachEnv *Env, const HorzMesh *Mesh,
                          const VertCoord *VCoord) {

   // 1D arrays
   testScalarMultiplyOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // 2D arrays
   testScalarMultiplyOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // 3D arrays
   testScalarMultiplyOpType<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testScalarMultiplyOpType<Array3DR8>("3D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Template for testing BinaryMultiplyOp with same-rank arrays
template <typename ArrayType>
void testBinaryMultiplyOpSameRank(const std::string &TypeName,
                                  const MachEnv *Env, const HorzMesh *Mesh,
                                  const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   std::vector<I4> Dims   = Helper::getDims(Mesh, VCoord);
   std::string Field1Name = "TestField1BinaryMult_" + TypeName;
   std::string Field2Name = "TestField2BinaryMult_" + TypeName;

   // Create two test fields with known values
   ScalarT Value1 = static_cast<ScalarT>(3);
   ScalarT Value2 = static_cast<ScalarT>(4);

   if constexpr (Rank == 1) {
      Helper::createField(Field1Name, Dims,
                          [Value1](I4 i) -> ScalarT { return Value1; });
      Helper::createField(Field2Name, Dims,
                          [Value2](I4 i) -> ScalarT { return Value2; });
   } else if constexpr (Rank == 2) {
      Helper::createField(Field1Name, Dims,
                          [Value1](I4 i, I4 j) -> ScalarT { return Value1; });
      Helper::createField(Field2Name, Dims,
                          [Value2](I4 i, I4 j) -> ScalarT { return Value2; });
   } else if constexpr (Rank == 3) {
      Helper::createField(
          Field1Name, Dims,
          [Value1](I4 i, I4 j, I4 k) -> ScalarT { return Value1; });
      Helper::createField(
          Field2Name, Dims,
          [Value2](I4 i, I4 j, I4 k) -> ScalarT { return Value2; });
   }

   ScalarT ExpectedValue = static_cast<ScalarT>(static_cast<Real>(Value1) *
                                                static_cast<Real>(Value2));

   // Create and compute operator
   auto MultOp = AnalysisOpFactory::createOp(
       "BinaryMultiply", {Field1Name, Field2Name}, makeOpConfig());
   MultOp->initialize(Env, Mesh, VCoord, makeOpConfig());
   TimeInstant TestTime;
   MultOp->compute(TestTime);

   // Get result
   auto ResultField =
       Field::get(Field1Name + "_BinaryMultiply(" + Field2Name + ")");
   auto ResultData = ResultField->getDataArray<ArrayType>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   // Verify all owned cells
   bool Passed = true;
   if constexpr (Rank == 1) {
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         if (std::abs(static_cast<Real>(ResultHost(i)) -
                      static_cast<Real>(ExpectedValue)) >
             static_cast<Real>(Helper::getTolerance())) {
            Passed = false;
            break;
         }
      }
   } else if constexpr (Rank == 2) {
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
              ++j) {
            if (std::abs(static_cast<Real>(ResultHost(i, j)) -
                         static_cast<Real>(ExpectedValue)) >
                static_cast<Real>(Helper::getTolerance())) {
               Passed = false;
               break;
            }
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Mesh->NCellsOwned; ++j) {
            for (I4 k = VCoord->MinLayerCellH(j); k <= VCoord->MaxLayerCellH(j);
                 ++k) {
               if (std::abs(static_cast<Real>(ResultHost(i, j, k)) -
                            static_cast<Real>(ExpectedValue)) >
                   static_cast<Real>(Helper::getTolerance())) {
                  Passed = false;
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

   reportTest("BinaryMultiplyOp (same-rank): " + TypeName, Passed);
   checkOutputMetadata("BinaryMultiplyOp (same-rank): " + TypeName,
                       Field1Name + "_BinaryMultiply(" + Field2Name + ")", "m2",
                       "", "");

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Value1: {}, Value2: {}",
                static_cast<Real>(ExpectedValue), static_cast<Real>(Value1),
                static_cast<Real>(Value2));
   }
}

//------------------------------------------------------------------------------
// Test BinaryMultiplyOp with vertical expansion (2D x 1D)
template <typename Array2DType>
void testBinaryMultiplyOpVerticalExpansion(const std::string &TypeName,
                                           const MachEnv *Env,
                                           const HorzMesh *Mesh,
                                           const VertCoord *VCoord) {

   using Helper2D       = TestHelper<Array2DType>;
   using ScalarT        = typename Helper2D::ScalarT;
   constexpr int Rank2D = Helper2D::Rank;

   // Only test 2D arrays for vertical expansion
   if constexpr (Rank2D != 2) {
      return;
   }

   std::vector<I4> Dims2D = Helper2D::getDims(Mesh, VCoord);
   std::vector<I4> Dims1D = {Dims2D[0]}; // Just horizontal dimension

   std::string Field2DName = "TestField2DVertExp_" + TypeName;
   std::string Field1DName = "TestField1DVertExp_" + TypeName;

   // Create 2D field with known values
   ScalarT Value2D = static_cast<ScalarT>(5);
   Helper2D::createField(Field2DName, Dims2D,
                         [Value2D](I4 i, I4 j) -> ScalarT { return Value2D; });

   // Create 1D field (horizontal only)
   ScalarT Value1D = static_cast<ScalarT>(2);
   using Helper1D  = TestHelper<Array1DR8>; // Use same scalar type base
   Helper1D::createField(Field1DName, Dims1D,
                         [Value1D](I4 i) -> ScalarT { return Value1D; });

   ScalarT ExpectedValue = static_cast<ScalarT>(static_cast<Real>(Value2D) *
                                                static_cast<Real>(Value1D));

   // Create and compute operator (2D x 1D = 2D with vertical expansion)
   auto MultOp = AnalysisOpFactory::createOp(
       "BinaryMultiply", {Field2DName, Field1DName}, makeOpConfig());
   MultOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   MultOp->compute(TestTime);

   // Get result
   auto ResultField =
       Field::get(Field2DName + "_BinaryMultiply(" + Field1DName + ")");
   auto ResultData = ResultField->getDataArray<Array2DType>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   // Verify: all active vertical layers should have same value (1D value
   // replicated)
   bool Passed = true;
   for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
      for (I4 j = VCoord->MinLayerCellH(i); j <= VCoord->MaxLayerCellH(i);
           ++j) {
         if (std::abs(static_cast<Real>(ResultHost(i, j)) -
                      static_cast<Real>(ExpectedValue)) >
             static_cast<Real>(Helper2D::getTolerance())) {
            Passed = false;
            LOG_ERROR("  At ({}, {}): Expected {}, Got {}", i, j,
                      static_cast<Real>(ExpectedValue),
                      static_cast<Real>(ResultHost(i, j)));
            break;
         }
      }
      if (!Passed)
         break;
   }

   reportTest("BinaryMultiplyOp (vertical expansion): " + TypeName, Passed);
   checkOutputMetadata("BinaryMultiplyOp (vertical expansion): " + TypeName,
                       Field2DName + "_BinaryMultiply(" + Field1DName + ")",
                       "m2", "", "");

   if (!Passed) {
      LOG_ERROR("  Expected: {}, Value2D: {}, Value1D: {}",
                static_cast<Real>(ExpectedValue), static_cast<Real>(Value2D),
                static_cast<Real>(Value1D));
   }
}

//------------------------------------------------------------------------------
// Test BinaryMultiplyOp with all array types
void testBinaryMultiplyOp(const MachEnv *Env, const HorzMesh *Mesh,
                          const VertCoord *VCoord) {

   // Test same-rank multiplication for 1D arrays
   testBinaryMultiplyOpSameRank<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // Test same-rank multiplication for 2D arrays
   testBinaryMultiplyOpSameRank<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array2DR8>("2D-R8", Env, Mesh, VCoord);

   // Test same-rank multiplication for 3D arrays
   testBinaryMultiplyOpSameRank<Array3DI4>("3D-I4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array3DI8>("3D-I8", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array3DR4>("3D-R4", Env, Mesh, VCoord);
   testBinaryMultiplyOpSameRank<Array3DR8>("3D-R8", Env, Mesh, VCoord);

   // Test vertical expansion for 2D × 1D
   testBinaryMultiplyOpVerticalExpansion<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testBinaryMultiplyOpVerticalExpansion<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testBinaryMultiplyOpVerticalExpansion<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testBinaryMultiplyOpVerticalExpansion<Array2DR8>("2D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Template for testing PrefixSumOp with 1D and 2D arrays
template <typename ArrayType>
void testPrefixSumOpType(const std::string &TypeName, const MachEnv *Env,
                         const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   const I4 NVertLayers = VCoord->NVertLayers;

   if constexpr (Rank == 1) {
      // --- 1D forward scan ---
      std::string FwdName = "TestFieldPrefixSum1D_Fwd_" + TypeName;
      auto FwdField =
          Field::create(FwdName, "Test 1D field for PrefixSum forward", "units",
                        "", static_cast<ScalarT>(-1000),
                        static_cast<ScalarT>(1000), 1, {"NVertLayers"});
      ArrayType FwdData(FwdName + "_data", NVertLayers);
      auto FwdHost = Kokkos::create_mirror_view(FwdData);
      for (I4 k = 0; k < NVertLayers; ++k) {
         FwdHost(k) = static_cast<ScalarT>(k + 1); // 1, 2, 3, ...
      }
      deepCopy(FwdData, FwdHost);
      FwdField->template attachData<ArrayType>(FwdData, false);

      auto FwdOp = AnalysisOpFactory::createOp(
          "PrefixSum", {FwdName},
          makeOpConfig(opParam("Dimension", 0), opParam("Reverse", false)));
      FwdOp->initialize(Env, Mesh, VCoord, makeOpConfig());
      TimeInstant TestTime;
      FwdOp->compute(TestTime);

      auto FwdResultField = Field::get(FwdName + "_PrefixSum");
      auto FwdResultData  = FwdResultField->getDataArray<ArrayType>();
      auto FwdResultHost  = createHostMirrorCopy(FwdResultData);

      bool FwdPassed = true;
      for (I4 k = 0; k < NVertLayers; ++k) {
         // Inclusive forward sum: Output[k] = sum(1..k+1) = (k+1)(k+2)/2
         Real Expected = static_cast<Real>((k + 1) * (k + 2)) / 2.0;
         Real Computed = static_cast<Real>(FwdResultHost(k));
         if (std::abs(Computed - Expected) >
             static_cast<Real>(Helper::getTolerance())) {
            FwdPassed = false;
            LOG_ERROR("  1D forward at {}: Expected {}, Got {}", k, Expected,
                      Computed);
            break;
         }
      }
      reportTest("PrefixSumOp: " + TypeName + " forward", FwdPassed);

      // --- 1D reverse scan ---
      std::string RevName = "TestFieldPrefixSum1D_Rev_" + TypeName;
      auto RevField =
          Field::create(RevName, "Test 1D field for PrefixSum reverse", "units",
                        "", static_cast<ScalarT>(-1000),
                        static_cast<ScalarT>(1000), 1, {"NVertLayers"});
      ArrayType RevData(RevName + "_data", NVertLayers);
      auto RevHost = Kokkos::create_mirror_view(RevData);
      for (I4 k = 0; k < NVertLayers; ++k) {
         RevHost(k) = static_cast<ScalarT>(k + 1); // 1, 2, 3, ...
      }
      deepCopy(RevData, RevHost);
      RevField->template attachData<ArrayType>(RevData, false);

      auto RevOp = AnalysisOpFactory::createOp(
          "PrefixSum", {RevName},
          makeOpConfig(opParam("Dimension", 0), opParam("Reverse", true)));
      RevOp->initialize(Env, Mesh, VCoord, makeOpConfig());
      RevOp->compute(TestTime);

      auto RevResultField = Field::get(RevName + "_PrefixSum");
      auto RevResultData  = RevResultField->getDataArray<ArrayType>();
      auto RevResultHost  = createHostMirrorCopy(RevResultData);

      bool RevPassed = true;
      for (I4 k = 0; k < NVertLayers; ++k) {
         // Inclusive reverse sum: Output[k] = sum(k+1..N)
         Real Expected = 0.0;
         for (I4 j = k; j < NVertLayers; ++j) {
            Expected += static_cast<Real>(j + 1);
         }
         Real Computed = static_cast<Real>(RevResultHost(k));
         if (std::abs(Computed - Expected) >
             static_cast<Real>(Helper::getTolerance())) {
            RevPassed = false;
            LOG_ERROR("  1D reverse at {}: Expected {}, Got {}", k, Expected,
                      Computed);
            break;
         }
      }
      reportTest("PrefixSumOp: " + TypeName + " reverse", RevPassed);

   } else if constexpr (Rank == 2) {
      // Use bins (non-distributed dimension) for testing
      // Create a simple 2D field: NBins x NVertLayers
      const I4 NBins = 10; // Small number for testing

      std::string FieldName = "TestFieldPrefixSum_" + TypeName;

      if (!Dimension::exists("TestBins")) {
         Dimension::create("TestBins", NBins);
      }
      // Create field with bins × vertical levels
      auto TestField =
          Field::create(FieldName, "Test field for PrefixSum", "units", "",
                        static_cast<ScalarT>(-1000), static_cast<ScalarT>(1000),
                        2, {"TestBins", "NVertLayers"});
      // Allocate and initialize data: value = (bin_idx + 1) * (level + 1)
      ArrayType TestData(FieldName + "_data", NBins, NVertLayers);
      auto TestHost = Kokkos::create_mirror_view(TestData);
      for (I4 i = 0; i < NBins; ++i) {
         for (I4 j = 0; j < NVertLayers; ++j) {
            TestHost(i, j) = static_cast<ScalarT>((i + 1) * (j + 1));
         }
      }
      deepCopy(TestData, TestHost);
      TestField->template attachData<ArrayType>(TestData, false);
      // Test vertical integration (dimension 1, reverse = true for MOC)
      auto PrefixOp = AnalysisOpFactory::createOp(
          "PrefixSum", {FieldName},
          makeOpConfig(opParam("Dimension", 1), opParam("Reverse", true)));
      PrefixOp->initialize(Env, Mesh, VCoord, makeOpConfig());
      TimeInstant TestTime;
      PrefixOp->compute(TestTime);
      // Get result
      auto ResultField = Field::get(FieldName + "_PrefixSum");
      auto ResultData  = ResultField->getDataArray<ArrayType>();
      auto ResultHost  = createHostMirrorCopy(ResultData);
      // Verify reverse cumulative sum
      bool Passed = true;

      for (I4 i = 0; i < NBins; ++i) {
         for (I4 j = 0; j < NVertLayers; ++j) {
            // Compute expected: sum from j to NVertLayers-1
            Real Expected = 0.0;
            for (I4 k = j; k < NVertLayers; ++k) {
               Expected += static_cast<Real>((i + 1) * (k + 1));
            }

            Real Computed = static_cast<Real>(ResultHost(i, j));

            if (std::abs(Computed - Expected) >
                static_cast<Real>(Helper::getTolerance())) {
               Passed = false;
               LOG_ERROR("  At bin {} level {}: Expected {}, Got {}", i, j,
                         Expected, Computed);
               break;
            }
         }
         if (!Passed)
            break;
      }
      reportTest("PrefixSumOp: " + TypeName, Passed);
      checkOutputMetadata("PrefixSumOp: " + TypeName, FieldName + "_PrefixSum",
                          "units", "", "");
   }
}

//------------------------------------------------------------------------------
// Test PrefixSumOp with 1D and 2D array types
void testPrefixSumOp(const MachEnv *Env, const HorzMesh *Mesh,
                     const VertCoord *VCoord) {

   // Test 1D arrays (e.g., output from TransectAccumulatorOp)
   testPrefixSumOpType<Array1DI4>("1D-I4", Env, Mesh, VCoord);
   testPrefixSumOpType<Array1DI8>("1D-I8", Env, Mesh, VCoord);
   testPrefixSumOpType<Array1DR4>("1D-R4", Env, Mesh, VCoord);
   testPrefixSumOpType<Array1DR8>("1D-R8", Env, Mesh, VCoord);

   // Test 2D arrays (most relevant for MOC vertical integration)
   testPrefixSumOpType<Array2DI4>("2D-I4", Env, Mesh, VCoord);
   testPrefixSumOpType<Array2DI8>("2D-I8", Env, Mesh, VCoord);
   testPrefixSumOpType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testPrefixSumOpType<Array2DR8>("2D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Template for testing PrefixSumOp with boundary condition (BC)
template <typename ArrayType>
void testPrefixSumOpWithBCType(const std::string &TypeName, const MachEnv *Env,
                               const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   // Only test 2D arrays for PrefixSum (most relevant for MOC)
   if constexpr (Rank != 2) {
      return;
   }

   // Only test floating point types (BC feature used with real transport
   // fields)
   if constexpr (!std::is_floating_point_v<ScalarT>) {
      return;
   }

   // Use bins (non-distributed dimension) for testing horizontal scan with BC
   const I4 NBins       = 10; // Small number for testing
   const I4 NVertLayers = VCoord->NVertLayers;

   std::string InputFieldName = "TestFieldPrefixSumBC_" + TypeName;
   std::string BCFieldName    = "TestBCPrefixSum_" + TypeName;

   // Ensure TestBins dimension exists
   if (!Dimension::exists("TestBins")) {
      Dimension::create("TestBins", NBins);
   }

   // Create BC field (1D: NVertLayers)
   // BC values represent boundary transport at southern edge (for regional MOC)
   auto BCField = Field::create(BCFieldName, "Test BC for PrefixSum", "m3/s",
                                "", static_cast<ScalarT>(-1e30),
                                static_cast<ScalarT>(1e30), 1, {"NVertLayers"});

   Array1DReal BCData(BCFieldName + "_data", NVertLayers);
   BCField->template attachData<Array1DReal>(BCData, false);

   auto BCHost = Kokkos::create_mirror_view(BCData);
   for (I4 k = 0; k < NVertLayers; ++k) {
      BCHost(k) = 100.0 * (k + 1); // Simple BC values: 100, 200, 300, ...
   }
   deepCopy(BCData, BCHost);

   // Create input field (2D: NBins × NVertLayers)
   auto InputField = Field::create(
       InputFieldName, "Test field for PrefixSum with BC", "m3/s", "",
       static_cast<ScalarT>(-1e30), static_cast<ScalarT>(1e30), 2,
       {"TestBins", "NVertLayers"});

   ArrayType InputData(InputFieldName + "_data", NBins, NVertLayers);
   auto InputHost = Kokkos::create_mirror_view(InputData);
   for (I4 i = 0; i < NBins; ++i) {
      for (I4 k = 0; k < NVertLayers; ++k) {
         InputHost(i, k) =
             static_cast<ScalarT>((i + 1) * (k + 1)); // 1,2,3... at bin 0,
                                                      // 2,4,6... at bin 1, etc.
      }
   }
   deepCopy(InputData, InputHost);
   InputField->template attachData<ArrayType>(InputData, false);

   // Test horizontal integration with BC (dimension 0, forward scan)
   // This simulates regional MOC where BC provides southern boundary transport
   auto PrefixOp = AnalysisOpFactory::createOp(
       "PrefixSum", {InputFieldName, BCFieldName},
       makeOpConfig(opParam("Dimension", 0), opParam("Reverse", false)));
   PrefixOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   PrefixOp->compute(TestTime);

   // Get result
   auto ResultField =
       Field::get(InputFieldName + "_PrefixSum(BC=" + BCFieldName + ")");
   auto ResultData = ResultField->getDataArray<ArrayType>();
   auto ResultHost = createHostMirrorCopy(ResultData);

   // Verify forward cumulative sum with BC seeding
   bool Passed = true;

   for (I4 i = 0; i < NBins; ++i) {
      for (I4 k = 0; k < NVertLayers; ++k) {
         // Compute expected: BC(k) + sum from 0 to i of Input(i,k)
         // At i=0: BC(k) + Input(0,k)
         // At i=1: BC(k) + Input(0,k) + Input(1,k)
         // etc.
         Real Expected = static_cast<Real>(BCHost(k));
         for (I4 bin = 0; bin <= i; ++bin) {
            Expected += static_cast<Real>(InputHost(bin, k));
         }

         Real Computed = static_cast<Real>(ResultHost(i, k));

         if (std::abs(Computed - Expected) >
             static_cast<Real>(Helper::getTolerance())) {
            Passed = false;
            LOG_ERROR("  At bin {} level {}: Expected {}, Got {}", i, k,
                      Expected, Computed);
            break;
         }
      }
      if (!Passed)
         break;
   }

   reportTest("PrefixSumOpWithBC: " + TypeName, Passed);
   checkOutputMetadata("PrefixSumOpWithBC: " + TypeName,
                       InputFieldName + "_PrefixSum(BC=" + BCFieldName + ")",
                       "m3/s", "", "");
}

//------------------------------------------------------------------------------
// Test PrefixSumOp with boundary condition using 2D array types
void testPrefixSumOpWithBC(const MachEnv *Env, const HorzMesh *Mesh,
                           const VertCoord *VCoord) {

   // Test 2D floating-point arrays (BC feature used for MOC regional transport)
   testPrefixSumOpWithBCType<Array2DR4>("2D-R4", Env, Mesh, VCoord);
   testPrefixSumOpWithBCType<Array2DR8>("2D-R8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test CoordinateBinningOp with 1D coordinate arrays
void testCoordinateBinningOp(const MachEnv *Env, const HorzMesh *Mesh,
                             const VertCoord *VCoord) {

   int MyRank;
   MPI_Comm_rank(Env->getComm(), &MyRank);

   std::string FieldName = "TestCoordBinning";

   // Create a simple 1D coordinate field with values from -90 to 90 (like
   // latitude)
   I4 NCells = Mesh->NCellsSize;
   auto TestField =
       Field::create(FieldName, "Test coordinate field for binning", "degrees",
                     "", -180.0, 180.0, 1, {"NCells"});

   Array1DR8 TestData(FieldName + "_data", NCells);
   TestField->attachData<Array1DR8>(TestData, false);

   auto TestDataHost = Kokkos::create_mirror_view(TestData);
   auto Decomp       = Decomp::getDefault();
   auto CellIDH      = Decomp->CellIDH;

   // Create latitude-like distribution: map global cell ID to latitude
   for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
      Real Fraction = static_cast<Real>(CellIDH(i) - 1) /
                      static_cast<Real>(Mesh->NCellsGlobal - 1);
      TestDataHost(i) = -90.0 + Fraction * 180.0; // -90 to +90
   }
   deepCopy(TestData, TestDataHost);

   // Create CoordinateBinning operator
   I4 NumBins     = 18; // 10-degree bins
   auto BinningOp = AnalysisOpFactory::createOp(
       "CoordinateBinning", {FieldName},
       makeOpConfig(opParam("NumBins", NumBins), opParam("MinBin", -90.0),
                    opParam("MaxBin", 90.0)));

   BinningOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   BinningOp->compute(TestTime);

   // Get result (bin indices)
   auto ResultField = Field::get(FieldName + "_BinIndex");
   auto ResultData  = ResultField->getDataArray<Array1DI4>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Verify bin assignments
   bool Passed = true;
   for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
      Real Coord     = TestDataHost(i);
      I4 ComputedBin = ResultHost(i);
      // Replicate exact operator logic including margins
      // (operator adds margin = (MaxBin-MinBin)*1e-6 to both ends)
      Real RawMin    = -90.0;
      Real RawMax    = 90.0;
      Real Margin    = (RawMax - RawMin) * 1.0e-6;
      Real AdjMinBin = RawMin - Margin;
      Real AdjMaxBin = RawMax + Margin;
      Real BinWidth  = (AdjMaxBin - AdjMinBin) / static_cast<Real>(NumBins);

      I4 ExpectedBin = static_cast<I4>((Coord - AdjMinBin) / BinWidth);

      // Clamp to valid range [0, NumBins-1] as the operator does
      if (ExpectedBin < 0) {
         ExpectedBin = 0;
      } else if (ExpectedBin >= NumBins) {
         ExpectedBin = NumBins - 1;
      }

      if (ComputedBin != ExpectedBin) {
         Passed = false;
         LOG_ERROR("  At index {}: Coord={}, Expected bin {}, Got bin {}", i,
                   Coord, ExpectedBin, ComputedBin);
         break;
      }
   }

   reportTest("CoordinateBinningOp", Passed);
   checkOutputMetadata("CoordinateBinningOp", FieldName + "_BinIndex", "1", "",
                       "");
}

//------------------------------------------------------------------------------
// Test BinnedAccumulatorOp with 2D field and bin indices
void testBinnedAccumulatorOp(const MachEnv *Env, const HorzMesh *Mesh,
                             const VertCoord *VCoord) {

   // Create bin index field (assign cells to bins based on pattern)
   I4 NCells  = Mesh->NCellsSize;
   I4 NumBins = 5;

   std::string BinIndexFieldName = "TestBinIndices";
   auto BinField = Field::create(BinIndexFieldName, "Test bin indices", "1", "",
                                 0, NumBins, 1, {"NCells"});

   Array1DI4 BinData(BinIndexFieldName + "_data", NCells);
   BinField->attachData<Array1DI4>(BinData, false);

   auto BinDataHost = Kokkos::create_mirror_view(BinData);
   for (I4 i = 0; i < NCells; ++i) {
      BinDataHost(i) = i % NumBins; // Distribute cells across bins
   }
   deepCopy(BinData, BinDataHost);

   // Create value field to accumulate (2D: NCells × NVertLayers)
   std::string ValueFieldName = "TestAccumValues";
   I4 NVertLayers             = VCoord->NVertLayers;

   auto ValueField =
       Field::create(ValueFieldName, "Test values for accumulation", "m3/s", "",
                     -1e30, 1e30, 2, {"NCells", "NVertLayers"});

   Array2DR8 ValueData(ValueFieldName + "_data", NCells, NVertLayers);
   ValueField->attachData<Array2DR8>(ValueData, false);

   auto ValueDataHost = Kokkos::create_mirror_view(ValueData);
   for (I4 i = 0; i < NCells; ++i) {
      for (I4 j = 0; j < NVertLayers; ++j) {
         ValueDataHost(i, j) = 1.0; // Uniform value for easy verification
      }
   }
   deepCopy(ValueData, ValueDataHost);

   // Create BinnedAccumulator operator
   auto AccumOp = AnalysisOpFactory::createOp(
       "BinnedAccumulator", {ValueFieldName, BinIndexFieldName},
       makeOpConfig(opParam("NumBins", NumBins)));
   AccumOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   AccumOp->compute(TestTime);

   // Get result (should be NumBins × NVertLayers)
   auto ResultField = Field::get(ValueFieldName + "_BinnedAccumulator(" +
                                 BinIndexFieldName + ")");
   auto ResultData  = ResultField->getDataArray<Array2DR8>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Calculate expected sums per bin per layer, respecting MinLayer/MaxLayer
   // (only cells where layer K is active contribute to that layer's bin sum)
   // Local counts per bin per layer
   std::vector<std::vector<I8>> LocalCounts(NumBins,
                                            std::vector<I8>(NVertLayers, 0));
   for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
      I4 BinIdx = BinDataHost(i);
      I4 KMin   = VCoord->MinLayerCellH(i);
      I4 KMax   = VCoord->MaxLayerCellH(i);
      for (I4 k = KMin; k <= KMax; ++k) {
         LocalCounts[BinIdx][k]++;
      }
   }

   // Global sum across MPI ranks
   std::vector<std::vector<I8>> GlobalCounts(NumBins,
                                             std::vector<I8>(NVertLayers, 0));
   for (I4 bin = 0; bin < NumBins; ++bin) {
      for (I4 k = 0; k < NVertLayers; ++k) {
         GlobalCounts[bin][k] = globalSum(LocalCounts[bin][k], Env->getComm());
      }
   }

   // Verify: each bin/layer sum should equal count * 1.0 (uniform value field)
   bool Passed = true;
   for (I4 bin = 0; bin < NumBins; ++bin) {
      for (I4 k = 0; k < NVertLayers; ++k) {
         Real Expected = static_cast<Real>(GlobalCounts[bin][k]) * 1.0;
         Real Computed = ResultHost(bin, k);

         if (std::abs(Computed - Expected) > 1.0e-8) {
            Passed = false;
            LOG_ERROR("  Bin {}, Layer {}: Expected {}, Got {}", bin, k,
                      Expected, Computed);
            break;
         }
      }
      if (!Passed)
         break;
   }

   reportTest("BinnedAccumulatorOp", Passed);
   checkOutputMetadata("BinnedAccumulatorOp",
                       ValueFieldName + "_BinnedAccumulator(" +
                           BinIndexFieldName + ")",
                       "m3/s", "", "");
}

//------------------------------------------------------------------------------
// Template for testing PseudoToGeometricOp with any array type
template <typename ArrayType>
void testPseudoToGeometricOpType(const std::string &TypeName,
                                 const MachEnv *Env, const HorzMesh *Mesh,
                                 const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   // Only test Real types (R4, R8) since SpecVol multiplication requires reals
   if constexpr (!std::is_floating_point_v<ScalarT>) {
      return;
   }

   std::vector<I4> Dims  = Helper::getDims(Mesh, VCoord);
   std::string FieldName = "TestFieldPseudoToGeom_" + TypeName;

   // Create and populate SpecVol Field with test data
   auto SpecVolField = Field::get("SpecVol");
   if (!SpecVolField) {
      // Create SpecVol field if it doesn't exist
      SpecVolField = Field::create("SpecVol", "Test specific volume", "m^3/kg",
                                   "", 0.0, 10.0, 2, {"NCells", "NVertLayers"});

      Array2DReal SpecVolDataNew("SpecVol_data", Mesh->NCellsSize,
                                 VCoord->NVertLayers);
      SpecVolField->attachData<Array2DReal>(SpecVolDataNew, false);
   }

   auto SpecVolData = SpecVolField->getDataArray<Array2DReal>();
   auto SpecVolHost = Kokkos::create_mirror_view(SpecVolData);

   // Fill with test values: SpecVol(i,k) = 1.0 + 0.01*i + 0.001*k
   for (I4 i = 0; i < Mesh->NCellsSize; ++i) {
      for (I4 k = 0; k < VCoord->NVertLayers; ++k) {
         SpecVolHost(i, k) = 1.0 + 0.01 * i + 0.001 * k;
      }
   }
   deepCopy(SpecVolData, SpecVolHost);

   // Create pseudo-coordinate test field with known values
   if constexpr (Rank == 1) {
      Helper::createField(FieldName, Dims, [](I4 i) -> ScalarT {
         return static_cast<ScalarT>(i + 1); // Simple linear values
      });
   } else if constexpr (Rank == 2) {
      Helper::createField(FieldName, Dims, [](I4 i, I4 k) -> ScalarT {
         return static_cast<ScalarT>((i + 1) * (k + 1)); // Product of indices
      });
   } else if constexpr (Rank == 3) {
      Helper::createField(FieldName, Dims, [](I4 i, I4 j, I4 k) -> ScalarT {
         return static_cast<ScalarT>((i + 1) * (j + 1) * (k + 1));
      });
   }

   // Create and compute operator
   auto PseudoToGeomOp = AnalysisOpFactory::createOp(
       "PseudoToGeometric", {FieldName}, makeOpConfig());
   PseudoToGeomOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   PseudoToGeomOp->compute(TestTime);

   // Get result
   auto ResultField = Field::get(FieldName + "_PseudoToGeometric");
   auto ResultData  = ResultField->getDataArray<ArrayType>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Get input field for verification
   auto InputField = Field::get(FieldName);
   auto InputData  = InputField->getDataArray<ArrayType>();
   auto InputHost  = createHostMirrorCopy(InputData);

   // Verify all elements
   bool Passed     = true;
   Real Tolerance  = static_cast<Real>(Helper::getTolerance());
   int NumVerified = 0;

   if constexpr (Rank == 1) {
      // 1D: check all horizontal points
      for (I4 i = 0; i < Mesh->NCellsOwned; i++) {
         // Compute expected using same precision path as operator
         Real PseudoVal  = static_cast<Real>(InputHost(i));
         Real SpecVolVal = SpecVolHost(i, 0);
         ScalarT Expected =
             static_cast<ScalarT>(RhoSw * SpecVolVal * PseudoVal);
         ScalarT Computed = ResultHost(i);

         if (std::abs(static_cast<Real>(Computed) -
                      static_cast<Real>(Expected)) /
                 static_cast<Real>(Expected) >
             Tolerance) {
            Passed = false;
            break;
         }
         NumVerified++;
      }
   } else if constexpr (Rank == 2) {
      // 2D: check all (i,k) points in active layers
      for (I4 i = 0; i < Mesh->NCellsOwned; i++) {
         // Respect active layer bounds for this cell, clamped to array bounds
         I4 KMin = VCoord->MinLayerCellH(i);
         I4 KMax = VCoord->MaxLayerCellH(i);
         for (I4 k = KMin; k <= KMax; k++) {
            // Compute expected using same precision path as operator
            Real PseudoVal  = static_cast<Real>(InputHost(i, k));
            Real SpecVolVal = SpecVolHost(i, k);
            ScalarT Expected =
                static_cast<ScalarT>(RhoSw * SpecVolVal * PseudoVal);
            ScalarT Computed = ResultHost(i, k);

            if (std::abs(static_cast<Real>(Computed) -
                         static_cast<Real>(Expected)) /
                    static_cast<Real>(Expected) >
                Tolerance) {
               Passed = false;
               break;
            }
            NumVerified++;
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      // 3D: check all (i,j,k) points in active layers
      for (I4 i = 0; i < Dims[0]; i++) {
         for (I4 j = 0; j < Mesh->NCellsOwned; j++) {
            // Respect active layer bounds for this cell, clamped to array
            // bounds
            I4 KMin = VCoord->MinLayerCellH(j);
            I4 KMax = VCoord->MaxLayerCellH(j);
            for (I4 k = KMin; k <= KMax; k++) {
               // Compute expected using same precision path as operator
               Real PseudoVal  = static_cast<Real>(InputHost(i, j, k));
               Real SpecVolVal = SpecVolHost(j, k);
               ScalarT Expected =
                   static_cast<ScalarT>(RhoSw * SpecVolVal * PseudoVal);
               ScalarT Computed = ResultHost(i, j, k);

               if (std::abs(static_cast<Real>(Computed) -
                            static_cast<Real>(Expected)) /
                       static_cast<Real>(Expected) >
                   Tolerance) {
                  Passed = false;
                  break;
               }
               NumVerified++;
            }
            if (!Passed)
               break;
         }
         if (!Passed)
            break;
      }
   }

   reportTest("PseudoToGeometricOp: " + TypeName, Passed);
   checkOutputMetadata("PseudoToGeometricOp: " + TypeName,
                       FieldName + "_PseudoToGeometric", TestUnits, TestStdName,
                       "");

   if (!Passed) {
      LOG_ERROR("PseudoToGeometricOp {} test failed after verifying {} points",
                TypeName, NumVerified);
   }
}

//------------------------------------------------------------------------------
// Test PseudoToGeometricOp with all supported array types
void testPseudoToGeometricOp(const MachEnv *Env, const HorzMesh *Mesh,
                             const VertCoord *VCoord) {

   // Test 1D arrays
   testPseudoToGeometricOpType<Array1DR4>("Array1DR4", Env, Mesh, VCoord);
   testPseudoToGeometricOpType<Array1DR8>("Array1DR8", Env, Mesh, VCoord);

   // Test 2D arrays
   testPseudoToGeometricOpType<Array2DR4>("Array2DR4", Env, Mesh, VCoord);
   testPseudoToGeometricOpType<Array2DR8>("Array2DR8", Env, Mesh, VCoord);

   // Test 3D arrays
   testPseudoToGeometricOpType<Array3DR4>("Array3DR4", Env, Mesh, VCoord);
   testPseudoToGeometricOpType<Array3DR8>("Array3DR8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test ExtractRegionOp for a specific array type
template <typename ArrayType>
void testExtractRegionOpType(const std::string &TypeName, const MachEnv *Env,
                             const HorzMesh *Mesh, const VertCoord *VCoord) {

   using Helper       = TestHelper<ArrayType>;
   using ScalarT      = typename Helper::ScalarT;
   constexpr int Rank = Helper::Rank;

   // Get dimensions and create test input field
   auto Dims     = Helper::getDims(Mesh, VCoord);
   auto DimNames = Helper::getDimNames();

   std::string InputFieldName = "TestExtractRegionInput_" + TypeName;
   std::string MaskFieldName  = "TestExtractRegionMask";

   // Create input field with known values
   Helper::createField(InputFieldName, Dims, [](auto... indices) {
      // Sum of indices + 1 to ensure non-zero values
      return static_cast<ScalarT>((indices + ...) + 1);
   });

   // Create mask field (1D horizontal, I4)
   // Mask = 1 for first half of cells, 0 for second half
   auto MaskField = Field::create(MaskFieldName, "Test regional mask", "1", "",
                                  0, 1, 1, {"NCells"});

   Array1DI4 MaskData(MaskFieldName + "_data", Mesh->NCellsSize);
   MaskField->attachData<Array1DI4>(MaskData, false);

   auto MaskDataHost = Kokkos::create_mirror_view(MaskData);
   I4 MidPoint       = Mesh->NCellsOwned / 2;
   for (I4 i = 0; i < Mesh->NCellsSize; ++i) {
      MaskDataHost(i) = (i < MidPoint) ? 1 : 0;
   }
   deepCopy(MaskData, MaskDataHost);

   // Create ExtractRegionOp
   auto Op = AnalysisOpFactory::createOp(
       "ExtractRegion", {InputFieldName},
       makeOpConfig(opParam("MaskName", MaskFieldName)));

   Op->initialize(Env, Mesh, VCoord, makeOpConfig());

   // Compute
   TimeInstant TestTime;
   Op->compute(TestTime);

   // Get result
   auto OutputNames = Op->getOutputFieldNames();
   OMEGA_ASSERT(OutputNames.size() == 1,
                "ExtractRegionOp should have 1 output");

   auto ResultField = Field::get(OutputNames[0]);
   auto ResultData  = ResultField->getDataArray<ArrayType>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Get input for verification
   auto InputField = Field::get(InputFieldName);
   auto InputData  = InputField->getDataArray<ArrayType>();
   auto InputHost  = createHostMirrorCopy(InputData);

   // Verify: where mask=1, result should equal input; where mask=0, result
   // should equal fill value
   bool Passed               = true;
   ScalarT Tol               = Helper::getTolerance();
   int NumVerified           = 0;
   constexpr ScalarT FillVal = FillValue<ScalarT>;

   if constexpr (Rank == 1) {
      // 1D: check all horizontal points
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         I4 MaskVal       = MaskDataHost(i);
         ScalarT Computed = ResultHost(i);
         ScalarT Expected;

         if (MaskVal == 0) {
            // Where mask is 0, expect fill value (not computed)
            Expected = FillVal;
         } else {
            // Where mask is 1, expect input value
            Expected = InputHost(i);
         }

         Real RelErr = std::abs(static_cast<Real>(Computed - Expected)) /
                       std::abs(static_cast<Real>(Expected));
         if (RelErr > Tol) {
            Passed = false;
            LOG_ERROR("  1D mismatch at cell {}: Mask={}, Expected {}, Got {}, "
                      "RelErr={}",
                      i, MaskVal, static_cast<Real>(Expected),
                      static_cast<Real>(Computed), RelErr);
            break;
         }
         NumVerified++;
      }
   } else if constexpr (Rank == 2) {
      // 2D: check all (i,k) points in active layers
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         I4 KMin    = VCoord->MinLayerCellH(i);
         I4 KMax    = VCoord->MaxLayerCellH(i);
         I4 MaskVal = MaskDataHost(i);
         for (I4 k = KMin; k <= KMax; ++k) {
            ScalarT Computed = ResultHost(i, k);
            ScalarT Expected;

            if (MaskVal == 0) {
               // Where mask is 0, expect fill value (not computed)
               Expected = FillVal;
            } else {
               // Where mask is 1, expect input value
               Expected = InputHost(i, k);
            }

            Real RelErr = std::abs(static_cast<Real>(Computed - Expected)) /
                          std::abs(static_cast<Real>(Expected));
            if (RelErr > Tol) {
               Passed = false;
               LOG_ERROR("  2D mismatch at cell {}, layer {}: Mask={}, "
                         "Expected {}, Got {}, RelErr={}",
                         i, k, MaskVal, static_cast<Real>(Expected),
                         static_cast<Real>(Computed), RelErr);
               break;
            }
            NumVerified++;
         }
         if (!Passed)
            break;
      }
   } else if constexpr (Rank == 3) {
      // 3D: check all (i,j,k) points in active layers
      for (I4 i = 0; i < Dims[0]; ++i) {
         for (I4 j = 0; j < Mesh->NCellsOwned; ++j) {
            I4 KMin    = VCoord->MinLayerCellH(j);
            I4 KMax    = VCoord->MaxLayerCellH(j);
            I4 MaskVal = MaskDataHost(j);
            for (I4 k = KMin; k <= KMax; ++k) {
               ScalarT Computed = ResultHost(i, j, k);
               ScalarT Expected;

               if (MaskVal == 0) {
                  // Where mask is 0, expect fill value (not computed)
                  Expected = FillVal;
               } else {
                  // Where mask is 1, expect input value
                  Expected = InputHost(i, j, k);
               }

               Real RelErr = std::abs(static_cast<Real>(Computed - Expected)) /
                             std::abs(static_cast<Real>(Expected));
               if (RelErr > Tol) {
                  Passed = false;
                  LOG_ERROR("  3D mismatch at i={}, cell {}, layer {}: "
                            "Mask={}, Expected {}, Got {}, RelErr={}",
                            i, j, k, MaskVal, static_cast<Real>(Expected),
                            static_cast<Real>(Computed), RelErr);
                  break;
               }
               NumVerified++;
            }
            if (!Passed)
               break;
         }
         if (!Passed)
            break;
      }
   }

   // Verify output field has mask attached
   if (Passed && !ResultField->hasRegionalMask()) {
      Passed = false;
      LOG_ERROR(
          "  ExtractRegionOp output field missing regional mask attachment");
   }

   reportTest("ExtractRegionOp: " + TypeName, Passed);
   checkOutputMetadata("ExtractRegionOp: " + TypeName, OutputNames[0],
                       TestUnits, TestStdName, "");

   if (!Passed) {
      LOG_ERROR("ExtractRegionOp {} test failed after verifying {} points",
                TypeName, NumVerified);
   }

   Field::destroy(MaskFieldName);
}

//------------------------------------------------------------------------------
// Test HorzMeanOp: area-weighted per-layer horizontal mean of a 2D field.
// Sub-test 1 (unmasked): verifies the weighted mean matches a host-computed
// reference using AreaCell weights and CellMask, that cells with an exactly-
// zero input value are excluded from both numerator and denominator (emulating
// an upstream ExtractRegion that zeroes out-of-region cells), and that a
// fully-zeroed level yields a guarded zero result (zero-denominator guard).
// Sub-test 2 (masked): attaches a 1D regional mask to the input Field and
// verifies that per-level means match a host reference computed with the same
// mask — confirming that zero-valued cells inside the region are included and
// non-zero cells outside the region are excluded.
// Sub-test 3 (interface): supplies an interface-dimensioned input
// (NVertLayersP1) and verifies the operator auto-detects the staggering and
// uses the extended active range [MinLayer, MaxLayer+1].
void testHorzMeanOp(const MachEnv *Env, const HorzMesh *Mesh,
                    const VertCoord *VCoord) {

   using Helper = TestHelper<Array2DR8>;

   auto Dims                  = Helper::getDims(Mesh, VCoord);
   std::string InputFieldName = "TestHorzMeanInput";
   const I4 NVert             = VCoord->NVertLayers;

   // Global cell IDs (1-based) so the input values and the excluded region
   // block are identical physical cells regardless of MPI decomposition.
   auto Decomp        = Decomp::getDefault();
   auto CellIDH       = Decomp->CellIDH;
   const I4 NCellsGlb = Mesh->NCellsGlobal;

   // Build a known input keyed on global cell ID. Force one vertical level
   // (ZeroLevel) to be exactly zero at every cell to exercise the
   // zero-denominator guard, and zero out a globally-defined block of cells at
   // another level (RegionLevel) to exercise presence-based exclusion. Both
   // choices are decomposition-independent for MPI correctness.
   const I4 ZeroLevel   = 0;
   const I4 RegionLevel = (NVert > 1) ? 1 : 0;
   const I4 BlockGID    = NCellsGlb / 2; // exclude global IDs >= BlockGID

   Helper::createField(
       InputFieldName, Dims,
       [ZeroLevel, RegionLevel, BlockGID, CellIDH](I4 i, I4 k) {
          I4 GID = CellIDH(i) - 1; // 0-based global id
          if (k == ZeroLevel)
             return 0.0; // whole level zero -> guarded result
          if (k == RegionLevel && GID >= BlockGID)
             return 0.0; // excluded block (emulates out-of-region zeroing)
          return static_cast<double>((GID + 1) * (k + 1));
       });

   // -------------------------------------------------------------------------
   // Sub-test 1: no regional mask attached — legacy Value != 0 behaviour
   // -------------------------------------------------------------------------
   auto Op = AnalysisOpFactory::createOp("HorzMean", {InputFieldName},
                                         makeOpConfig());
   Op->initialize(Env, Mesh, VCoord, makeOpConfig());
   TimeInstant TestTime;
   Op->compute(TestTime);

   auto OutputNames = Op->getOutputFieldNames();
   OMEGA_ASSERT(OutputNames.size() == 1, "HorzMeanOp should have 1 output");
   auto ResultField = Field::get(OutputNames[0]);
   auto ResultData  = ResultField->getDataArray<Array1DReal>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Host references from the same weights/mask the operator uses
   auto InputField = Field::get(InputFieldName);
   auto InputHost = createHostMirrorCopy(InputField->getDataArray<Array2DR8>());
   auto AreaHost  = createHostMirrorCopy(Mesh->AreaCell);
   auto MinLayerH = VCoord->MinLayerCellH;
   auto MaxLayerH = VCoord->MaxLayerCellH;

   bool Passed    = true;
   const Real Tol = 1.0e-8;
   for (I4 k = 0; k < NVert; ++k) {
      // Accumulate this rank's owned contribution, then reduce globally with
      // the same semantics as the operator so the reference is MPI-correct.
      // The input is on mid-layers, so the active range is [MinLayer,
      // MaxLayer].
      Real LocalNum = 0.0, LocalDen = 0.0;
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         const bool Active = (k >= MinLayerH(i) && k <= MaxLayerH(i));
         Real Value        = InputHost(i, k);
         if (Active && Value != 0.0) {
            LocalNum += AreaHost(i) * Value;
            LocalDen += AreaHost(i);
         }
      }
      Real GlobalNum = globalSum(LocalNum, Env->getComm());
      Real GlobalDen = globalSum(LocalDen, Env->getComm());

      Real Expected = (GlobalDen != 0.0) ? (GlobalNum / GlobalDen) : 0.0;
      Real Computed = ResultHost(k);
      Real Err      = std::abs(Computed - Expected);
      if (Err > Tol * (1.0 + std::abs(Expected))) {
         Passed = false;
         LOG_ERROR("  HorzMean (unmasked) mismatch at level {}: "
                   "Expected {}, Got {}",
                   k, Expected, Computed);
         break;
      }
   }

   // Explicitly confirm the fully-zeroed level yielded the guarded zero.
   if (Passed && std::abs(ResultHost(ZeroLevel)) > Tol) {
      Passed = false;
      LOG_ERROR("  HorzMean zero-denominator guard failed at level {}: Got {}",
                ZeroLevel, ResultHost(ZeroLevel));
   }

   reportTest("HorzMeanOp", Passed);
   checkOutputMetadata("HorzMeanOp", OutputNames[0], TestUnits, TestStdName,
                       "area: mean");

   Field::destroy(OutputNames[0]);

   // -------------------------------------------------------------------------
   // Sub-test 2: 1D regional mask attached to input Field.
   //
   // Use a mask that selects global IDs [0, BlockGID) and excludes the rest —
   // the same geographic split used for RegionLevel above, but now expressed
   // as an explicit 1D mask rather than zeroed values. With the mask in place:
   //   - zero-valued cells INSIDE the region are INCLUDED in the mean,
   //   - non-zero cells OUTSIDE the region are EXCLUDED from the mean.
   // This is the opposite of the legacy Value != 0 path, confirming the mask
   // takes full control when attached.
   // -------------------------------------------------------------------------
   Array1DI4 RegMaskData("HorzMeanTestRegMask", Mesh->NCellsSize);
   auto RegMaskHost = Kokkos::create_mirror_view(RegMaskData);
   for (I4 i = 0; i < Mesh->NCellsSize; ++i) {
      I4 GID         = CellIDH(i) - 1; // 0-based global id
      RegMaskHost(i) = (GID < BlockGID) ? 1 : 0;
   }
   Kokkos::deep_copy(RegMaskData, RegMaskHost);
   InputField->setRegionalMask(RegMaskData);

   // Re-create the operator — initialize() must detect the newly attached mask.
   auto Op2 = AnalysisOpFactory::createOp("HorzMean", {InputFieldName},
                                          makeOpConfig());
   Op2->initialize(Env, Mesh, VCoord, makeOpConfig());
   Op2->compute(TestTime);

   auto OutputNames2 = Op2->getOutputFieldNames();
   OMEGA_ASSERT(OutputNames2.size() == 1, "HorzMeanOp should have 1 output");
   auto ResultField2 = Field::get(OutputNames2[0]);
   auto ResultHost2 =
       createHostMirrorCopy(ResultField2->getDataArray<Array1DReal>());

   bool Passed2 = true;
   for (I4 k = 0; k < NVert; ++k) {
      // Host reference: include cell i iff VertCoord mask is active AND the
      // regional mask selects it. No Value != 0 requirement — genuinely zero
      // in-region values participate in the mean.
      Real LocalNum = 0.0, LocalDen = 0.0;
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         const bool Active = (k >= MinLayerH(i) && k <= MaxLayerH(i));
         if (!Active)
            continue;
         if (RegMaskHost(i) == 0)
            continue;
         LocalNum += AreaHost(i) * InputHost(i, k);
         LocalDen += AreaHost(i);
      }
      Real GlobalNum = globalSum(LocalNum, Env->getComm());
      Real GlobalDen = globalSum(LocalDen, Env->getComm());

      Real Expected = (GlobalDen != 0.0) ? (GlobalNum / GlobalDen) : 0.0;
      Real Computed = ResultHost2(k);
      Real Err      = std::abs(Computed - Expected);
      if (Err > Tol * (1.0 + std::abs(Expected))) {
         Passed2 = false;
         LOG_ERROR("  HorzMean (masked) mismatch at level {}: "
                   "Expected {}, Got {}",
                   k, Expected, Computed);
         break;
      }
   }

   reportTest("HorzMeanOp (masked)", Passed2);

   Field::destroy(OutputNames2[0]);

   // -------------------------------------------------------------------------
   // Sub-test 3: interface-dimensioned input (NVertLayersP1).
   //
   // Build an input field on interfaces to confirm the operator auto-detects
   // the staggering and uses the extended active range [MinLayer, MaxLayer+1].
   // No regional mask is attached, so the legacy Value != 0 presence applies.
   // -------------------------------------------------------------------------
   const I4 NVertP1               = VCoord->NVertLayersP1;
   std::string InterfaceFieldName = "TestHorzMeanInterfaceInput";
   auto InterfaceField            = Field::create(
       InterfaceFieldName, "Interface input for HorzMean", "m", "",
       -std::numeric_limits<Real>::max(), std::numeric_limits<Real>::max(), 2,
       {"NCells", "NVertLayersP1"});
   Array2DR8 InterfaceData(InterfaceFieldName + "_data", Mesh->NCellsSize,
                           NVertP1);
   auto InterfaceDataH = Kokkos::create_mirror_view(InterfaceData);
   for (I4 i = 0; i < Mesh->NCellsSize; ++i) {
      I4 GID = CellIDH(i) - 1; // 0-based global id
      for (I4 k = 0; k < NVertP1; ++k)
         InterfaceDataH(i, k) = static_cast<double>((GID + 1) * (k + 1));
   }
   Kokkos::deep_copy(InterfaceData, InterfaceDataH);
   InterfaceField->attachData<Array2DR8>(InterfaceData, false);

   auto Op3 = AnalysisOpFactory::createOp("HorzMean", {InterfaceFieldName},
                                          makeOpConfig());
   Op3->initialize(Env, Mesh, VCoord, makeOpConfig());
   Op3->compute(TestTime);

   auto OutputNames3 = Op3->getOutputFieldNames();
   OMEGA_ASSERT(OutputNames3.size() == 1, "HorzMeanOp should have 1 output");
   auto ResultHost3 = createHostMirrorCopy(
       Field::get(OutputNames3[0])->getDataArray<Array1DReal>());

   bool Passed3 = true;
   for (I4 k = 0; k < NVertP1; ++k) {
      // Interface input: active range is [MinLayer, MaxLayer+1].
      Real LocalNum = 0.0, LocalDen = 0.0;
      for (I4 i = 0; i < Mesh->NCellsOwned; ++i) {
         const bool Active = (k >= MinLayerH(i) && k <= MaxLayerH(i) + 1);
         Real Value        = InterfaceDataH(i, k);
         if (Active && Value != 0.0) {
            LocalNum += AreaHost(i) * Value;
            LocalDen += AreaHost(i);
         }
      }
      Real GlobalNum = globalSum(LocalNum, Env->getComm());
      Real GlobalDen = globalSum(LocalDen, Env->getComm());

      Real Expected = (GlobalDen != 0.0) ? (GlobalNum / GlobalDen) : 0.0;
      Real Computed = ResultHost3(k);
      Real Err      = std::abs(Computed - Expected);
      if (Err > Tol * (1.0 + std::abs(Expected))) {
         Passed3 = false;
         LOG_ERROR("  HorzMean (interface) mismatch at level {}: "
                   "Expected {}, Got {}",
                   k, Expected, Computed);
         break;
      }
   }

   reportTest("HorzMeanOp (interface)", Passed3);

   Field::destroy(InterfaceFieldName);
   Field::destroy(OutputNames3[0]);

   Field::destroy(InputFieldName);
}

//------------------------------------------------------------------------------
// Test ExtractRegionOp with all supported array types
void testExtractRegionOp(const MachEnv *Env, const HorzMesh *Mesh,
                         const VertCoord *VCoord) {

   // Test 1D arrays (all scalar types)
   testExtractRegionOpType<Array1DI4>("Array1DI4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array1DI8>("Array1DI8", Env, Mesh, VCoord);
   testExtractRegionOpType<Array1DR4>("Array1DR4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array1DR8>("Array1DR8", Env, Mesh, VCoord);

   // Test 2D arrays (all scalar types)
   testExtractRegionOpType<Array2DI4>("Array2DI4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array2DI8>("Array2DI8", Env, Mesh, VCoord);
   testExtractRegionOpType<Array2DR4>("Array2DR4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array2DR8>("Array2DR8", Env, Mesh, VCoord);

   // Test 3D arrays (all scalar types)
   testExtractRegionOpType<Array3DI4>("Array3DI4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array3DI8>("Array3DI8", Env, Mesh, VCoord);
   testExtractRegionOpType<Array3DR4>("Array3DR4", Env, Mesh, VCoord);
   testExtractRegionOpType<Array3DR8>("Array3DR8", Env, Mesh, VCoord);
}

//------------------------------------------------------------------------------
// Test TransectAccumulatorOp with edge-based transport
// Tests the accumulation of edge transport across a synthetic transect.
// Uses non-contiguous edge pattern (every 5th edge) to thoroughly test
// MPI reduction; real transects would be spatially contiguous.
void testTransectAccumulatorOp(const MachEnv *Env, const HorzMesh *Mesh,
                               const VertCoord *VCoord) {

   I4 NEdges      = Mesh->NEdgesSize;
   I4 NVertLayers = VCoord->NVertLayers;

   // Create transport field (NEdges × NVertLayers)
   // Value at each point = layer index + 1 (so layer 0→1.0, layer 1→2.0, etc.)
   std::string TransportFieldName = "TestTransectTransport";
   auto TransportField =
       Field::create(TransportFieldName, "Test transport for transect", "m3/s",
                     "", -1e30, 1e30, 2, {"NEdges", "NVertLayers"});

   Array2DR8 TransportData(TransportFieldName + "_data", NEdges, NVertLayers);
   TransportField->attachData<Array2DR8>(TransportData, false);

   auto TransportDataHost = Kokkos::create_mirror_view(TransportData);
   for (I4 i = 0; i < NEdges; ++i) {
      for (I4 k = 0; k < NVertLayers; ++k) {
         TransportDataHost(i, k) =
             static_cast<Real>(k + 1); // Layer-dependent value
      }
   }
   deepCopy(TransportData, TransportDataHost);

   // Create transect mask field (which edges are on transect)
   // Every 5th edge is on transect (mask=1.0), others are off (mask=0.0)
   std::string MaskFieldName = "TestTransectMask";
   auto MaskField = Field::create(MaskFieldName, "Test transect mask", "1", "",
                                  0, 1, 1, {"NEdges"});

   Array1DI4 MaskData(MaskFieldName + "_data", NEdges);
   MaskField->attachData<Array1DI4>(MaskData, false);

   auto MaskDataHost = Kokkos::create_mirror_view(MaskData);
   for (I4 i = 0; i < NEdges; ++i) {
      MaskDataHost(i) = ((i % 5) == 0) ? 1 : 0;
   }
   deepCopy(MaskData, MaskDataHost);

   // Create transect sign field (direction of transport)
   // Alternating +1.0 and -1.0 for edges on transect
   std::string SignFieldName = "TestTransectSign";
   auto SignField = Field::create(SignFieldName, "Test transect sign", "1", "",
                                  -1, 1, 1, {"NEdges"});

   Array1DI4 SignData(SignFieldName + "_data", NEdges);
   SignField->attachData<Array1DI4>(SignData, false);

   auto SignDataHost = Kokkos::create_mirror_view(SignData);
   for (I4 i = 0; i < NEdges; ++i) {
      if ((i % 5) == 0) {
         // Edges on transect: alternate sign based on edge index
         SignDataHost(i) = ((i / 5) % 2 == 0) ? 1 : -1;
      } else {
         SignDataHost(i) = 0; // Doesn't matter, mask is 0
      }
   }
   deepCopy(SignData, SignDataHost);

   // Calculate expected result per layer
   // For each layer k: sum over owned edges of transport(i,k) × mask(i) ×
   // sign(i)
   std::vector<Real> LocalSumPerLayer(NVertLayers, 0.0);

   for (I4 i = 0; i < Mesh->NEdgesOwned; ++i) {
      if (MaskDataHost(i) != 0) { // Edge is on transect
         I4 KMin = VCoord->MinLayerEdgeBotH(i);
         I4 KMax = VCoord->MaxLayerEdgeTopH(i);
         for (I4 k = KMin; k <= KMax; ++k) {
            LocalSumPerLayer[k] +=
                TransportDataHost(i, k) * static_cast<Real>(SignDataHost(i));
         }
      }
   }

   // Global sum across MPI ranks
   std::vector<Real> ExpectedPerLayer(NVertLayers);
   for (I4 k = 0; k < NVertLayers; ++k) {
      ExpectedPerLayer[k] = globalSum(LocalSumPerLayer[k], Env->getComm());
   }

   // Create TransectAccumulator operator
   auto TransectOp = AnalysisOpFactory::createOp(
       "TransectAccumulator",
       {TransportFieldName, MaskFieldName, SignFieldName},
       makeOpConfig(opParam("TransectName", "TestTransect")));
   TransectOp->initialize(Env, Mesh, VCoord, makeOpConfig());

   TimeInstant TestTime;
   TransectOp->compute(TestTime);

   // Get result (should be 1D: NVertLayers)
   // Use operator's OutputNames to get the correct field name
   auto OutputNames = TransectOp->getOutputFieldNames();
   if (OutputNames.empty()) {
      LOG_ERROR("TransectAccumulatorOp has no output names!");
      reportTest("TransectAccumulatorOp", false);
      return;
   }
   auto ResultField = Field::get(OutputNames[0]);
   auto ResultData  = ResultField->getDataArray<Array1DR8>();
   auto ResultHost  = createHostMirrorCopy(ResultData);

   // Verify results for all layers
   bool Passed = true;
   for (I4 k = 0; k < NVertLayers; ++k) {
      Real Computed = ResultHost(k);
      Real Expected = ExpectedPerLayer[k];

      if (std::abs(Computed - Expected) > 1.0e-8) {
         Passed = false;
         LOG_ERROR("  Layer {}: Expected {}, Got {}", k, Expected, Computed);
         break;
      }
   }

   reportTest("TransectAccumulatorOp", Passed);
   checkOutputMetadata("TransectAccumulatorOp", OutputNames[0], "m3/s", "", "");
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

   // Initialize VertAdv
   VertAdv::init();

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

   // Read the initial state so the pseudo-thickness the mass weights use is
   // the one from the mesh file rather than fill values
   if (!IOStream::validateAll())
      ABORT_ERROR("AnalysisOperatorTest: stream validation failed");
   Metadata ReqMeta;
   Error ReadErr = IOStream::read("InitialState", ModelClock, ReqMeta);
   CHECK_ERROR_ABORT(ReadErr,
                     "AnalysisOperatorTest: failed to read initial state");

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

      testWeightedStats(DefEnv, Mesh, VCoord);

      testTimeMeanOp(DefEnv, Mesh, VCoord, ModelClock);

      testScalarMultiplyOp(DefEnv, Mesh, VCoord);

      testBinaryMultiplyOp(DefEnv, Mesh, VCoord);

      testPrefixSumOp(DefEnv, Mesh, VCoord);

      testPrefixSumOpWithBC(DefEnv, Mesh, VCoord);

      testCoordinateBinningOp(DefEnv, Mesh, VCoord);

      testBinnedAccumulatorOp(DefEnv, Mesh, VCoord);

      testPseudoToGeometricOp(DefEnv, Mesh, VCoord);

      testExtractRegionOp(DefEnv, Mesh, VCoord);

      testHorzMeanOp(DefEnv, Mesh, VCoord);

      testTransectAccumulatorOp(DefEnv, Mesh, VCoord);

      testInheritedMetadata(DefEnv, Mesh, VCoord, ModelClock);

      testDerivedUnits(DefEnv, Mesh, VCoord);

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
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return Err;
}
