//===-- Test driver for OMEGA Fields -----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Fields
///
/// This driver tests the capabilities for OMEGA to manage the registration
/// and use of fields and metadata. This test driver will not include IO and
/// only tests the registration and retrieval of fields and metadata.
//
//===-----------------------------------------------------------------------===/

#include "Field.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "Dimension.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "TimeMgr.h"
#include "mpi.h"
#include <vector>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Set some constant reference values for simplicity
const I4 RefI4           = 3;
const I8 RefI8           = 400000000;
const R4 RefR4           = 5.1;
const R8 RefR8           = 6.123456789;
const bool RefBool       = true;
const std::string RefStr = "Reference String";

//------------------------------------------------------------------------------
// A simple test evaluation function
template <typename T>
void TstEval(const std::string &TestName, T TestVal, T ExpectVal, int &Error) {

   if (TestVal == ExpectVal) {
      LOG_INFO("{}: PASS", TestName);
   } else {
      LOG_ERROR("{}: FAIL", TestName);
      ++Error;
   }
}
//------------------------------------------------------------------------------
// Initialization routine to create reference Fields
int initFieldTest() {

   int Err    = 0;
   int ErrRef = 0;

   // Initialize various environments
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   initLogging(DefEnv);
   MPI_Comm DefComm = DefEnv->getComm();
   Err              = IO::init(DefComm);
   if (Err != 0) {
      LOG_ERROR("IO initialization failed");
      return Err;
   }

   // Open config file
   OMEGA::Config("Omega");
   Err = OMEGA::Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("FieldTest: Error reading config file");
      return Err;
   }

   // Initialize decomposition
   Decomp::init();
   Decomp *DefDecomp = Decomp::getDefault();

   // Create offsets for dimension definition
   I4 NCellsSize      = DefDecomp->NCellsSize;
   I4 NCellsOwned     = DefDecomp->NCellsOwned;
   I4 NCellsGlobal    = DefDecomp->NCellsGlobal;
   I4 NEdgesSize      = DefDecomp->NEdgesSize;
   I4 NEdgesOwned     = DefDecomp->NEdgesOwned;
   I4 NEdgesGlobal    = DefDecomp->NEdgesGlobal;
   I4 NVerticesSize   = DefDecomp->NVerticesSize;
   I4 NVerticesOwned  = DefDecomp->NVerticesOwned;
   I4 NVerticesGlobal = DefDecomp->NVerticesGlobal;
   HostArray1DI4 CellOffset("NCellsOffset", NCellsSize);
   HostArray1DI4 EdgeOffset("NEdgesOffset", NEdgesSize);
   HostArray1DI4 VrtxOffset("NVerticesOffset", NVerticesSize);

   for (int N = 0; N < NCellsSize; ++N) {
      if (N < NCellsOwned) {
         CellOffset(N) = DefDecomp->CellIDH(N) - 1; // Offset must be zero-based
      } else {
         CellOffset(N) = -1; // Denotes cells that are not to be used
      }
   }
   for (int N = 0; N < NEdgesSize; ++N) {
      if (N < NEdgesOwned) {
         EdgeOffset(N) = DefDecomp->EdgeIDH(N) - 1; // Offset must be zero-based
      } else {
         EdgeOffset(N) = -1; // Denotes edges that are not to be used
      }
   }
   for (int N = 0; N < NVerticesSize; ++N) {
      if (N < NVerticesOwned) {
         VrtxOffset(N) = DefDecomp->VertexIDH(N) - 1; // zero-based
      } else {
         VrtxOffset(N) = -1; // Denotes edges that are not to be used
      }
   }

   // Define dimensions
   std::shared_ptr<Dimension> CellDim =
       Dimension::create("NCells", NCellsGlobal, NCellsSize, CellOffset);
   std::shared_ptr<Dimension> EdgeDim =
       Dimension::create("NEdges", NEdgesGlobal, NEdgesSize, EdgeOffset);
   std::shared_ptr<Dimension> VrtxDim = Dimension::create(
       "NVertices", NVerticesGlobal, NVerticesSize, VrtxOffset);
   I4 NVertLevels = 100;
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLevels", NVertLevels);
   I4 NTracers                        = 2;
   std::shared_ptr<Dimension> TrcrDim = Dimension::create("NTracers", NTracers);
   I4 NTime                           = 2;
   std::shared_ptr<Dimension> TimeDim = Dimension::create("NTime", NTime);
   I4 NStuff                          = 2;
   std::shared_ptr<Dimension> StuffDim =
       Dimension::create("NStuff", NVertLevels);

   // Create a model clock for time info
   TimeInstant SimStartTime(0001, 1, 1, 0, 0, 0.0);
   TimeInterval TimeStep(2, TimeUnits::Hours);
   Clock *ModelClock = new Clock(SimStartTime, TimeStep);

   // Initialize Field class - creates Model and Sim metadata fields
   int Err1 = Field::init(ModelClock);
   TstEval<int>("Field initialization", Err1, ErrRef, Err);

   // Add some global (Model and Simulation) metadata
   // First using single metadata add functions
   std::shared_ptr<Field> CodeField = Field::get(CodeMeta);
   Err1                             = CodeField->addMetadata("CodeI4", RefI4);
   TstEval<int>("Add metadata I4", Err1, ErrRef, Err);
   Err1 = CodeField->addMetadata("CodeI8", RefI8);
   TstEval<int>("Add metadata I8", Err1, ErrRef, Err);
   Err1 = CodeField->addMetadata("CodeR4", RefR4);
   TstEval<int>("Add metadata R4", Err1, ErrRef, Err);
   Err1 = CodeField->addMetadata("CodeR8", RefR8);
   TstEval<int>("Add metadata R8", Err1, ErrRef, Err);
   Err1 = CodeField->addMetadata("CodeStr", RefStr);
   TstEval<int>("Add metadata String", Err1, ErrRef, Err);
   Err1 = CodeField->addMetadata("CodeBool", true);
   TstEval<int>("Add metadata String", Err1, ErrRef, Err);

   // Now with the multiple add function
   std::shared_ptr<Field> SimField = Field::get(SimMeta);
   Err1                            = SimField->addMetadata(
       {std::make_pair("SimI4", RefI4), std::make_pair("SimI8", RefI8),
                                   std::make_pair("SimR4", RefR4), std::make_pair("SimR8", RefR8),
                                   std::make_pair("SimBool", true), std::make_pair("SimStr", RefStr)});
   TstEval<int>("Add multiple global meta entries", Err1, ErrRef, Err);

   // Define fields for various data types, mesh locations and memory
   // locations.

   std::vector<std::string> DimNames(1);

   // 1D Fields on host

   DimNames[0] = "NCells";
   auto Test1DI4H =
       Field::create("Test1DI4H", "Test 1DI4 field on host", "Units1DI4H",
                     "var_name_1DI4", 0, 100000, 999, 1, DimNames);

   DimNames[0] = "NEdges";
   auto Test1DI8H =
       Field::create("Test1DI8H", "Test 1DI8 field on host", "Units1DI8H",
                     "var_name_1DI8", 0, 100000, 999, 1, DimNames);

   DimNames[0] = "NVertices";
   auto Test1DR4H =
       Field::create("Test1DR4H", "Test 1DR4 field on host", "Units1DR4H",
                     "var_name_1DR4", 0.0, 100000.0, 999, 1, DimNames);

   DimNames[0] = "NVertLevels";
   auto Test1DR8H =
       Field::create("Test1DR8H", "Test 1DR8 field on host", "Units1DR8H",
                     "var_name_1DR8", 0.0, 100000.0, 999, 1, DimNames);

   // 1D Fields on device

   DimNames[0] = "NCells";
   auto Test1DI4 =
       Field::create("Test1DI4", "Test 1DI4 field on device", "Units1DI4",
                     "var_name_1DI4", 0, 100000, 999, 1, DimNames);

   DimNames[0] = "NEdges";
   auto Test1DI8 =
       Field::create("Test1DI8", "Test 1DI8 field on device", "Units1DI8",
                     "var_name_1DI8", 0, 100000, 999, 1, DimNames);

   DimNames[0] = "NVertices";
   auto Test1DR4 =
       Field::create("Test1DR4", "Test 1DR4 field on device", "Units1DR4",
                     "var_name_1DR4", 0.0, 100000.0, 999, 1, DimNames);

   DimNames[0] = "NVertLevels";
   auto Test1DR8 =
       Field::create("Test1DR8", "Test 1DR8 field on device", "Units1DR8",
                     "var_name_1DR8", 0.0, 100000.0, 999, 1, DimNames);

   // 2D Fields on host
   DimNames.resize(2);

   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";
   auto Test2DI4H =
       Field::create("Test2DI4H", "Test 2DI4 field on host", "Units2DI4H",
                     "var_name_2DI4", 0, 100000, 999, 2, DimNames);

   DimNames[0] = "NEdges";
   DimNames[1] = "NVertLevels";
   auto Test2DI8H =
       Field::create("Test2DI8H", "Test 2DI8 field on host", "Units2DI8H",
                     "var_name_2DI8", 0, 100000, 999, 2, DimNames);

   DimNames[0] = "NVertices";
   DimNames[1] = "NVertLevels";
   auto Test2DR4H =
       Field::create("Test2DR4H", "Test 2DR4 field on host", "Units2DR4H",
                     "var_name_2DR4", 0.0, 100000.0, 999, 2, DimNames);

   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";
   auto Test2DR8H =
       Field::create("Test2DR8H", "Test 2DR8 field on host", "Units2DR8H",
                     "var_name_2DR8", 0.0, 100000.0, 999, 2, DimNames);

   // 2D Fields on device

   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";
   auto Test2DI4 =
       Field::create("Test2DI4", "Test 2DI4 field on device", "Units2DI4",
                     "var_name_2DI4", 0, 100000, 999, 2, DimNames);

   DimNames[0] = "NEdges";
   DimNames[1] = "NVertLevels";
   auto Test2DI8 =
       Field::create("Test2DI8", "Test 2DI8 field on device", "Units2DI8",
                     "var_name_2DI8", 0, 100000, 999, 2, DimNames);

   DimNames[0] = "NVertices";
   DimNames[1] = "NVertLevels";
   auto Test2DR4 =
       Field::create("Test2DR4", "Test 2DR4 field on device", "Units2DR4",
                     "var_name_2DR4", 0.0, 100000.0, 999, 2, DimNames);

   DimNames[0] = "NEdges";
   DimNames[1] = "NVertLevels";
   auto Test2DR8 =
       Field::create("Test2DR8", "Test 2DR8 field on device", "Units2DR8",
                     "var_name_2DR8", 0.0, 100000.0, 999, 2, DimNames);

   // Higher dimension fields on device

   DimNames.resize(3);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";
   DimNames[2] = "NTracers";
   auto Test3DI4 =
       Field::create("Test3DI4", "Test 3DI4 field on device", "Units3DI4",
                     "var_name_3DI4", 0, 100000, 999, 3, DimNames);

   DimNames.resize(4);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";
   DimNames[2] = "NTracers";
   DimNames[3] = "NTime";
   auto Test4DI8 =
       Field::create("Test4DI8", "Test 4DI8 field on device", "Units4DI8",
                     "var_name_4DI8", 0, 100000, 999, 4, DimNames);

   DimNames.resize(5);
   DimNames[0] = "NCells";
   DimNames[1] = "NVertLevels";
   DimNames[2] = "NTracers";
   DimNames[3] = "NTime";
   DimNames[4] = "NStuff";
   auto Test5DR4 =
       Field::create("Test5DR4", "Test 5DR4 field on device", "Units5DR4",
                     "var_name_5DR4", 0, 100000, 999, 5, DimNames);

   /// Create two field groups for 1D and 2D fields
   auto FieldGroup1D = FieldGroup::create("FieldGroup1D");
   auto FieldGroup2D = FieldGroup::create("FieldGroup2D");

   // Add fields to each group - for 1d use the member function, for
   // 2d, use the add-by-name interface
   Err1 = FieldGroup1D->addField("Test1DI4");
   TstEval<int>("Add field to group 1DI4", Err1, ErrRef, Err);
   Err1 = FieldGroup1D->addField("Test1DI8");
   TstEval<int>("Add field to group 1DI8", Err1, ErrRef, Err);
   Err1 = FieldGroup1D->addField("Test1DR4");
   TstEval<int>("Add field to group 1DR4", Err1, ErrRef, Err);
   Err1 = FieldGroup1D->addField("Test1DR8");
   TstEval<int>("Add field to group 1DR8", Err1, ErrRef, Err);

   Err1 = FieldGroup::addFieldToGroup("Test2DI4", "FieldGroup2D");
   TstEval<int>("Add field to group 2DI4", Err1, ErrRef, Err);
   Err1 = FieldGroup::addFieldToGroup("Test2DI8", "FieldGroup2D");
   TstEval<int>("Add field to group 2DI8", Err1, ErrRef, Err);
   Err1 = FieldGroup::addFieldToGroup("Test2DR4", "FieldGroup2D");
   TstEval<int>("Add field to group 2DR4", Err1, ErrRef, Err);
   Err1 = FieldGroup::addFieldToGroup("Test2DR8", "FieldGroup2D");
   TstEval<int>("Add field to group 2DR8", Err1, ErrRef, Err);

   // Create data arrays

   HostArray1DI4 Data1DI4H("Test1DI4H", NCellsSize);
   HostArray1DI8 Data1DI8H("Test1DI8H", NEdgesSize);
   HostArray1DR4 Data1DR4H("Test1DR4H", NVerticesSize);
   HostArray1DR8 Data1DR8H("Test1DR8H", NVertLevels);

   HostArray2DI4 Data2DI4H("Test2DI4H", NCellsSize, NVertLevels);
   HostArray2DI8 Data2DI8H("Test2DI8H", NEdgesSize, NVertLevels);
   HostArray2DR4 Data2DR4H("Test2DR4H", NVerticesSize, NVertLevels);
   HostArray2DR8 Data2DR8H("Test2DR8H", NCellsSize, NVertLevels);

   Array1DI4 Data1DI4("Test1DI4", NCellsSize);
   Array1DI8 Data1DI8("Test1DI8", NEdgesSize);
   Array1DR4 Data1DR4("Test1DR4", NVerticesSize);
   Array1DR8 Data1DR8("Test1DR8", NVertLevels);

   Array2DI4 Data2DI4("Test2DI4", NCellsSize, NVertLevels);
   Array2DI8 Data2DI8("Test2DI8", NEdgesSize, NVertLevels);
   Array2DR4 Data2DR4("Test2DR4", NVerticesSize, NVertLevels);
   Array2DR8 Data2DR8("Test2DR8", NEdgesSize, NVertLevels);

   Array3DI4 Data3DI4("Test3DI4", NCellsSize, NVertLevels, NTracers);
   Array4DI8 Data4DI8("Test4DI8", NCellsSize, NVertLevels, NTracers, NTime);
   Array5DR4 Data5DR4("Test5DR4", NCellsSize, NVertLevels, NTracers, NTime,
                      NStuff);

   // Host arrays vertical vector
   for (int K = 0; K < NVertLevels; ++K) {
      Data1DR8H(K) = RefR8 + K;
   }
   // Host arrays on cells
   for (int Cell = 0; Cell < NCellsSize; ++Cell) {
      Data1DI4H(Cell) = RefI4 + Cell;
      for (int K = 0; K < NVertLevels; ++K) {
         Data2DI4H(Cell, K) = RefI4 + Cell + K;
         Data2DR8H(Cell, K) = RefR8 + Cell + K;
      }
   }
   // Host arrays on edges
   for (int Edge = 0; Edge < NEdgesSize; ++Edge) {
      Data1DI8H(Edge) = RefI8 + Edge;
      for (int K = 0; K < NVertLevels; ++K) {
         Data2DI8H(Edge, K) = RefI8 + Edge + K;
      }
   }
   // Host arrays on vertices
   for (int Vrtx = 0; Vrtx < NVerticesSize; ++Vrtx) {
      Data1DR4H(Vrtx) = RefR4 + Vrtx;
      for (int K = 0; K < NVertLevels; ++K) {
         Data2DR4H(Vrtx, K) = RefR4 + Vrtx + K;
      }
   }
   // Device array vertical vector
   parallelFor(
       {NVertLevels}, KOKKOS_LAMBDA(int K) { Data1DR8(K) = RefR8 + K; });
   // Device arrays on cells
   parallelFor(
       {NCellsSize}, KOKKOS_LAMBDA(int Cell) {
          Data1DI4(Cell) = RefI4 + Cell;
          for (int K = 0; K < NVertLevels; ++K) {
             Data2DI4(Cell, K) = RefI4 + Cell + K;
             for (int Trcr = 0; Trcr < NTracers; ++Trcr) {
                Data3DI4(Cell, K, Trcr) = RefI4 + Cell + K + Trcr;
                for (int TimeLvl = 0; TimeLvl < NTime; ++TimeLvl) {
                   Data4DI8(Cell, K, Trcr, TimeLvl) =
                       RefI8 + Cell + K + Trcr + TimeLvl;
                   for (int Stf = 0; Stf < NStuff; ++Stf) {
                      Data5DR4(Cell, K, Trcr, TimeLvl, Stf) =
                          RefR4 + Cell + K + Trcr + TimeLvl + Stf;
                   }
                }
             }
          }
       });
   // Device arrays on edges
   parallelFor(
       {NEdgesSize}, KOKKOS_LAMBDA(int Edge) {
          Data1DI8(Edge) = RefI8 + Edge;
          for (int K = 0; K < NVertLevels; ++K) {
             Data2DI8(Edge, K) = RefI8 + Edge + K;
             Data2DR8(Edge, K) = RefR8 + Edge + K;
          }
       });
   // Device arrays on vertices
   parallelFor(
       {NVerticesSize}, KOKKOS_LAMBDA(int Vrtx) {
          Data1DR4(Vrtx) = RefR4 + Vrtx;
          for (int K = 0; K < NVertLevels; ++K) {
             Data2DR4(Vrtx, K) = RefR4 + Vrtx + K;
          }
       });

   // Attach data arrays
   // Use member function for some and name interface for others

   Err1 = Test1DI4H->attachData<HostArray1DI4>(Data1DI4H);
   TstEval<int>("Attach data to field 1DI4H", Err1, ErrRef, Err);
   Err1 = Test1DI8H->attachData<HostArray1DI8>(Data1DI8H);
   TstEval<int>("Attach data to field 1DI8H", Err1, ErrRef, Err);
   Err1 = Test1DR4H->attachData<HostArray1DR4>(Data1DR4H);
   TstEval<int>("Attach data to field 1DR4H", Err1, ErrRef, Err);
   Err1 = Test1DR8H->attachData<HostArray1DR8>(Data1DR8H);
   TstEval<int>("Attach data to field 1DR8H", Err1, ErrRef, Err);

   Err1 = Field::attachFieldData<HostArray2DI4>("Test2DI4H", Data2DI4H);
   TstEval<int>("Attach data to field 2DI4H", Err1, ErrRef, Err);
   Err1 = Field::attachFieldData<HostArray2DI8>("Test2DI8H", Data2DI8H);
   TstEval<int>("Attach data to field 2DI8H", Err1, ErrRef, Err);
   Err1 = Field::attachFieldData<HostArray2DR4>("Test2DR4H", Data2DR4H);
   TstEval<int>("Attach data to field 2DR4H", Err1, ErrRef, Err);
   Err1 = Field::attachFieldData<HostArray2DR8>("Test2DR8H", Data2DR8H);
   TstEval<int>("Attach data to field 2DR8H", Err1, ErrRef, Err);

   Err1 = Field::attachFieldData("Test1DI4", Data1DI4);
   TstEval<int>("Attach data to field 1DI4", Err1, ErrRef, Err);
   Err1 = Field::attachFieldData("Test1DI8", Data1DI8);
   TstEval<int>("Attach data to field 1DI8", Err1, ErrRef, Err);
   Err1 = Field::attachFieldData("Test1DR4", Data1DR4);
   TstEval<int>("Attach data to field 1DR4", Err1, ErrRef, Err);
   Err1 = Field::attachFieldData("Test1DR8", Data1DR8);
   TstEval<int>("Attach data to field 1DR8", Err1, ErrRef, Err);

   Err1 = Test2DI4->attachData(Data2DI4);
   TstEval<int>("Attach data to field 2DI4", Err1, ErrRef, Err);
   Err1 = Test2DI8->attachData(Data2DI8);
   TstEval<int>("Attach data to field 2DI8", Err1, ErrRef, Err);
   Err1 = Test2DR4->attachData(Data2DR4);
   TstEval<int>("Attach data to field 2DR4", Err1, ErrRef, Err);
   Err1 = Test2DR8->attachData(Data2DR8);
   TstEval<int>("Attach data to field 2DR8", Err1, ErrRef, Err);

   Err1 = Test3DI4->attachData(Data3DI4);
   TstEval<int>("Attach data to field 3DI4", Err1, ErrRef, Err);
   Err1 = Test4DI8->attachData(Data4DI8);
   TstEval<int>("Attach data to field 4DI8", Err1, ErrRef, Err);
   Err1 = Test5DR4->attachData(Data5DR4);
   TstEval<int>("Attach data to field 5DR4", Err1, ErrRef, Err);

   // End of init
   return Err;

} // End initialization Fields

//------------------------------------------------------------------------------
// We will test the IO Field interfaces by defining an IO Field during init
// and then retrieving the field and comparing results. Because IO Field makes
// use of the Metadata class and mostly points to other arrays/classes, we
// do not test all possible combinations - just two data types and host/device
// arrays to make sure generic pointers work fine.

int main(int argc, char **argv) {

   int Err    = 0;
   int Err1   = 0;
   int ErrRef = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {
      // Call initialization to create reference IO field
      Err1 = initFieldTest();
      TstEval<int>("Initialize Fields", Err1, ErrRef, Err);

      // Retrieve some mesh sizes to be used later
      Decomp *DefDecomp = Decomp::getDefault();
      I4 NCellsSize     = DefDecomp->NCellsSize;
      I4 NEdgesSize     = DefDecomp->NEdgesSize;
      I4 NVerticesSize  = DefDecomp->NVerticesSize;
      I4 NVertLevels    = 100;
      I4 NTracers       = 2;
      I4 NTime          = 2;
      I4 NStuff         = 2;

      // Test metadata retrieval functions using the global Code, Sim metadata
      // Use Code to retrieve single metadata entries
      std::shared_ptr<Field> CodeField = Field::get(CodeMeta);

      // Test the existence function for both existence and non-existence
      bool MetaTest = CodeField->hasMetadata("CodeI4") and
                      CodeField->hasMetadata("CodeI8") and
                      CodeField->hasMetadata("CodeR4") and
                      CodeField->hasMetadata("CodeR8") and
                      CodeField->hasMetadata("CodeBool") and
                      CodeField->hasMetadata("CodeStr") and
                      !CodeField->hasMetadata("Junk");
      TstEval<bool>("Metadata existence", MetaTest, true, Err);

      // Retrieve individual metadata entries

      I4 MetaI4 = 0;
      Err1      = CodeField->getMetadata("CodeI4", MetaI4);
      TstEval<I4>("Get I4 Metadata", MetaI4, RefI4, Err);
      I8 MetaI8 = 0;
      Err1      = CodeField->getMetadata("CodeI8", MetaI8);
      TstEval<I8>("Get I8 Metadata", MetaI8, RefI8, Err);
      R4 MetaR4 = 0;
      Err1      = CodeField->getMetadata("CodeR4", MetaR4);
      TstEval<R4>("Get R4 Metadata", MetaR4, RefR4, Err);
      R8 MetaR8 = 0;
      Err1      = CodeField->getMetadata("CodeR8", MetaR8);
      TstEval<R8>("Get R8 Metadata", MetaR8, RefR8, Err);
      bool MetaBool;
      Err1 = CodeField->getMetadata("CodeBool", MetaBool);
      TstEval<bool>("Get Bool Metadata", MetaBool, RefBool, Err);
      std::string MetaStr = " ";
      Err1                = CodeField->getMetadata("CodeStr", MetaStr);
      TstEval<std::string>("Get string Metadata", MetaStr, RefStr, Err);

      // Retrieve all Metadata for the Code Field.
      std::any MetaVal;
      std::shared_ptr<Metadata> CodeMetaAll = CodeField->getAllMetadata();
      MetaVal                               = (*CodeMetaAll)["CodeI4"];
      TstEval<I4>("Get All Metadata I4", std::any_cast<I4>(MetaVal), RefI4,
                  Err);
      MetaVal = (*CodeMetaAll)["CodeI8"];
      TstEval<I8>("Get All Metadata I8", std::any_cast<I8>(MetaVal), RefI8,
                  Err);
      MetaVal = (*CodeMetaAll)["CodeR4"];
      TstEval<R4>("Get All Metadata R4", std::any_cast<R4>(MetaVal), RefR4,
                  Err);
      MetaVal = (*CodeMetaAll)["CodeR8"];
      TstEval<R8>("Get All Metadata R8", std::any_cast<R8>(MetaVal), RefR8,
                  Err);
      MetaVal = (*CodeMetaAll)["CodeBool"];
      TstEval<bool>("Get All Metadata bool", std::any_cast<bool>(MetaVal),
                    RefBool, Err);
      MetaVal = (*CodeMetaAll)["CodeStr"];
      TstEval<std::string>("Get All Metadata string",
                           std::any_cast<std::string>(MetaVal), RefStr, Err);

      // Retrieve all Metadata from SimMeta using the named interface.
      std::shared_ptr<Metadata> SimMetaAll = Field::getFieldMetadata(SimMeta);
      MetaVal                              = (*SimMetaAll)["SimI4"];
      TstEval<I4>("Get All Metadata I4", std::any_cast<I4>(MetaVal), RefI4,
                  Err);
      MetaVal = (*SimMetaAll)["SimI8"];
      TstEval<I8>("Get All Metadata I8", std::any_cast<I8>(MetaVal), RefI8,
                  Err);
      MetaVal = (*SimMetaAll)["SimR4"];
      TstEval<R4>("Get All Metadata R4", std::any_cast<R4>(MetaVal), RefR4,
                  Err);
      MetaVal = (*SimMetaAll)["SimR8"];
      TstEval<R8>("Get All Metadata R8", std::any_cast<R8>(MetaVal), RefR8,
                  Err);
      MetaVal = (*SimMetaAll)["SimBool"];
      TstEval<bool>("Get All Metadata bool", std::any_cast<bool>(MetaVal),
                    RefBool, Err);
      MetaVal = (*SimMetaAll)["SimStr"];
      TstEval<std::string>("Get All Metadata string",
                           std::any_cast<std::string>(MetaVal), RefStr, Err);

      // Test updating a metadata value
      Err1 = CodeField->updateMetadata("CodeI8", MetaI8 + 2);
      Err1 = CodeField->getMetadata("CodeI8", MetaI8);
      TstEval<I8>("Update I8 Metadata", MetaI8, RefI8 + 2, Err);

      // Test removal of a metadata entry with the given name
      Err1 = CodeField->removeMetadata("CodeBool");
      TstEval<int>("Remove metadata call", Err1, ErrRef, Err);
      MetaTest = CodeField->hasMetadata("CodeBool");
      TstEval<bool>("Remove metadata verify", MetaTest, false, Err);

      // Test removal of all metadata entries from a field
      CodeField->removeAllMetadata();
      MetaTest = CodeField->hasMetadata("CodeI4") and
                 CodeField->hasMetadata("CodeI8") and
                 CodeField->hasMetadata("CodeR4") and
                 CodeField->hasMetadata("CodeR8") and
                 CodeField->hasMetadata("CodeBool") and
                 CodeField->hasMetadata("CodeStr");
      TstEval<bool>("Remove all metadata", MetaTest, false, Err);

      // Finished with metadata, now retrieve other field properties
      std::shared_ptr<Field> Test1DI4H = Field::get("Test1DI4H");

      // Retrieve field name
      std::string MyName = Test1DI4H->getName();
      TstEval<std::string>("Retrieve field name", MyName, "Test1DI4H", Err);

      // Retrieve data type of field
      ArrayDataType Type1DI4 = Test1DI4H->getType();
      TstEval<ArrayDataType>("Retrieve field type - member", Type1DI4,
                             ArrayDataType::I4, Err);

      // Retrieve type of a given field by name
      ArrayDataType Type1DR8 = Field::getFieldType("Test1DR8H");
      TstEval<ArrayDataType>("Retrieve type by field name", Type1DR8,
                             ArrayDataType::R8, Err);

      // Retrieve location of field data
      ArrayMemLoc MemLoc1DI4H = Test1DI4H->getMemoryLocation();
      if (MemLoc1DI4H == ArrayMemLoc::Both or
          MemLoc1DI4H == ArrayMemLoc::Host) {
         LOG_INFO("Retrieve mem location - member: PASS");
      } else {
         LOG_ERROR("Retrieve mem location - member: FAIL");
         ++Err;
      }
      bool OnHost = Test1DI4H->isOnHost();
      TstEval<bool>("Retrieve onHost flag", OnHost, true, Err);

      ArrayMemLoc MemLoc1DI4 = Field::getFieldMemoryLocation("Test1DI4");
      if (MemLoc1DI4 == ArrayMemLoc::Both or
          MemLoc1DI4 == ArrayMemLoc::Device) {
         LOG_INFO("Retrieve mem location by name: PASS");
      } else {
         LOG_ERROR("Retrieve mem location by name: FAIL");
         ++Err;
      }
      OnHost = Field::isFieldOnHost("Test1DI4");
      if (MemLoc1DI4 == ArrayMemLoc::Both) {
         TstEval<bool>("Retrieve onHost flag for device", OnHost, true, Err);
      } else {
         TstEval<bool>("Retrieve onHost flag for device", OnHost, false, Err);
      }

      // Check retrieval of dimension information
      int NDims = Test1DI4H->getNumDims();
      TstEval<int>("Retrieve number of dims1", NDims, 1, Err);

      std::vector<std::string> DimNames(NDims);
      Err1 = Test1DI4H->getDimNames(DimNames);
      TstEval<int>("Retrieve dim names 1 function call", Err1, ErrRef, Err);
      TstEval<std::string>("Retrieve dim names 1 result 1", DimNames[0],
                           "NCells", Err);

      std::shared_ptr<Field> Test5DR4 = Field::get("Test5DR4");
      NDims                           = Test5DR4->getNumDims();
      TstEval<int>("Retrieve number of dims5", NDims, 5, Err);

      DimNames.resize(NDims);
      Err1 = Test5DR4->getDimNames(DimNames);
      TstEval<int>("Retrieve dim names 5 function call", Err1, ErrRef, Err);
      TstEval<std::string>("Retrieve dim names 5 result 1", DimNames[0],
                           "NCells", Err);
      TstEval<std::string>("Retrieve dim names 5 result 2", DimNames[1],
                           "NVertLevels", Err);
      TstEval<std::string>("Retrieve dim names 5 result 3", DimNames[2],
                           "NTracers", Err);
      TstEval<std::string>("Retrieve dim names 5 result 4", DimNames[3],
                           "NTime", Err);
      TstEval<std::string>("Retrieve dim names 5 result 5", DimNames[4],
                           "NStuff", Err);

      // Test some field group functions and use to retrieve some remaining
      // fields for data testing

      // Test existence function
      bool GroupTest = FieldGroup::exists("FieldGroup1D") and
                       FieldGroup::exists("FieldGroup2D") and
                       !FieldGroup::exists("Junk");
      TstEval<bool>("Field group existence", GroupTest, true, Err);

      // Retrieve the 1D group and use the 2D group to test the
      // retrieval-by-name interfaces
      std::shared_ptr<FieldGroup> Group1D = FieldGroup::get("FieldGroup1D");

      // Test membership in groups
      GroupTest =
          Group1D->hasField("Test1DI4") and Group1D->hasField("Test1DI8") and
          Group1D->hasField("Test1DR4") and Group1D->hasField("Test1DR8") and
          not Group1D->hasField("Junk");
      TstEval<bool>("Field group members 1D", GroupTest, true, Err);

      GroupTest = FieldGroup::isFieldInGroup("Test2DI4", "FieldGroup2D") and
                  FieldGroup::isFieldInGroup("Test2DI8", "FieldGroup2D") and
                  FieldGroup::isFieldInGroup("Test2DR4", "FieldGroup2D") and
                  FieldGroup::isFieldInGroup("Test2DR8", "FieldGroup2D") and
                  not FieldGroup::isFieldInGroup("Junk", "FieldGroup2D");
      TstEval<bool>("Field group members 2D", GroupTest, true, Err);

      // Similarly test retrieval of field list
      std::set<std::string> List1D = Group1D->getFieldList();
      std::set<std::string> List2D =
          FieldGroup::getFieldListFromGroup("FieldGroup2D");
      for (auto Iter = List1D.begin(); Iter != List1D.end(); ++Iter) {
         std::string FieldName = *Iter;
         GroupTest             = Group1D->hasField(FieldName);
         TstEval<bool>("Get Field list from group 1D", GroupTest, true, Err);
      }
      for (auto Iter = List2D.begin(); Iter != List2D.end(); ++Iter) {
         std::string FieldName = *Iter;
         GroupTest = FieldGroup::isFieldInGroup(FieldName, "FieldGroup2D");
         TstEval<bool>("Get Field list from group 2D", GroupTest, true, Err);
      }

      // Retrieve all group fields to using Group retrievals (skip any fields
      // already retrieved above)

      std::shared_ptr<Field> Test1DI4 = Group1D->getField("Test1DI4");
      std::shared_ptr<Field> Test1DI8 = Group1D->getField("Test1DI8");
      std::shared_ptr<Field> Test1DR4 = Group1D->getField("Test1DR4");
      std::shared_ptr<Field> Test1DR8 = Group1D->getField("Test1DR8");

      // Retrieves a full field from a group with a given name
      std::shared_ptr<Field> Test2DI4 =
          FieldGroup::getFieldFromGroup("Test2DI4", "FieldGroup2D");
      std::shared_ptr<Field> Test2DI8 =
          FieldGroup::getFieldFromGroup("Test2DI8", "FieldGroup2D");
      std::shared_ptr<Field> Test2DR4 =
          FieldGroup::getFieldFromGroup("Test2DR4", "FieldGroup2D");
      std::shared_ptr<Field> Test2DR8 =
          FieldGroup::getFieldFromGroup("Test2DR8", "FieldGroup2D");

      // Retrieve all remaining groups using standard get function
      // 1DI4 already retrieved previously
      std::shared_ptr<Field> Test1DI8H = Field::get("Test1DI8H");
      std::shared_ptr<Field> Test1DR4H = Field::get("Test1DR4H");
      std::shared_ptr<Field> Test1DR8H = Field::get("Test1DR8H");

      std::shared_ptr<Field> Test2DI4H = Field::get("Test2DI4H");
      std::shared_ptr<Field> Test2DI8H = Field::get("Test2DI8H");
      std::shared_ptr<Field> Test2DR4H = Field::get("Test2DR4H");
      std::shared_ptr<Field> Test2DR8H = Field::get("Test2DR8H");

      std::shared_ptr<Field> Test3DI4 = Field::get("Test3DI4");
      std::shared_ptr<Field> Test4DI8 = Field::get("Test4DI8");
      // Test5DR4 already retrieved above

      // Retrieve data arrays and data pointers. As before, use both
      // member functions and get-by-name interfaces.
      HostArray1DI4 Data1DI4H = Test1DI4H->getDataArray<HostArray1DI4>();
      HostArray1DI8 Data1DI8H = Test1DI8H->getDataArray<HostArray1DI8>();
      HostArray1DR4 Data1DR4H = Test1DR4H->getDataArray<HostArray1DR4>();
      HostArray1DR8 Data1DR8H = Test1DR8H->getDataArray<HostArray1DR8>();

      HostArray2DI4 Data2DI4H =
          Field::getFieldDataArray<HostArray2DI4>("Test2DI4H");
      HostArray2DI8 Data2DI8H =
          Field::getFieldDataArray<HostArray2DI8>("Test2DI8H");
      HostArray2DR4 Data2DR4H =
          Field::getFieldDataArray<HostArray2DR4>("Test2DR4H");
      HostArray2DR8 Data2DR8H =
          Field::getFieldDataArray<HostArray2DR8>("Test2DR8H");

      Array1DI4 Data1DI4 = Test1DI4->getDataArray<Array1DI4>();
      Array1DI8 Data1DI8 = Test1DI8->getDataArray<Array1DI8>();
      Array1DR4 Data1DR4 = Test1DR4->getDataArray<Array1DR4>();
      Array1DR8 Data1DR8 = Test1DR8->getDataArray<Array1DR8>();

      Array2DI4 Data2DI4 = Field::getFieldDataArray<Array2DI4>("Test2DI4");
      Array2DI8 Data2DI8 = Field::getFieldDataArray<Array2DI8>("Test2DI8");
      Array2DR4 Data2DR4 = Field::getFieldDataArray<Array2DR4>("Test2DR4");
      Array2DR8 Data2DR8 = Field::getFieldDataArray<Array2DR8>("Test2DR8");

      Array3DI4 Data3DI4 = Test3DI4->getDataArray<Array3DI4>();
      Array4DI8 Data4DI8 = Test4DI8->getDataArray<Array4DI8>();
      Array5DR4 Data5DR4 = Test5DR4->getDataArray<Array5DR4>();

      // Test values for correctness
      // Host arrays vertical vector
      int DataCount1 = 0;
      for (int K = 0; K < NVertLevels; ++K) {
         if (Data1DR8H(K) != RefR8 + K)
            ++DataCount1;
      }
      TstEval<int>("Get data 1DR8H", DataCount1, 0, Err);

      // Host arrays on cells
      DataCount1     = 0;
      int DataCount2 = 0;
      int DataCount3 = 0;
      for (int Cell = 0; Cell < NCellsSize; ++Cell) {
         if (Data1DI4H(Cell) != RefI4 + Cell)
            ++DataCount1;
         for (int K = 0; K < NVertLevels; ++K) {
            I4 I4Ref = RefI4 + Cell + K;
            R8 R8Ref = RefR8 + Cell + K;
            if (Data2DI4H(Cell, K) != RefI4 + Cell + K)
               ++DataCount2;
            if (Data2DR8H(Cell, K) != RefR8 + Cell + K)
               ++DataCount3;
         }
      }
      TstEval<int>("Get data 1DI4H", DataCount1, 0, Err);
      TstEval<int>("Get data 2DI4H", DataCount2, 0, Err);
      TstEval<int>("Get data 2DR8H", DataCount3, 0, Err);

      // Host arrays on edges
      DataCount1 = 0;
      DataCount2 = 0;
      for (int Edge = 0; Edge < NEdgesSize; ++Edge) {
         if (Data1DI8H(Edge) != RefI8 + Edge)
            ++DataCount1;
         for (int K = 0; K < NVertLevels; ++K) {
            if (Data2DI8H(Edge, K) != RefI8 + Edge + K)
               ++DataCount2;
         }
      }
      TstEval<int>("Get data 1DI8H", DataCount1, 0, Err);
      TstEval<int>("Get data 2DI8H", DataCount2, 0, Err);

      // Host arrays on vertices
      DataCount1 = 0;
      DataCount2 = 0;
      for (int Vrtx = 0; Vrtx < NVerticesSize; ++Vrtx) {
         if (Data1DR4H(Vrtx) != RefR4 + Vrtx)
            ++DataCount1;
         for (int K = 0; K < NVertLevels; ++K) {
            if (Data2DR4H(Vrtx, K) != RefR4 + Vrtx + K)
               ++DataCount2;
         }
      }
      TstEval<int>("Get data 1DR4H", DataCount1, 0, Err);
      TstEval<int>("Get data 2DR4H", DataCount2, 0, Err);

      // For device arrays, we need to use reduction loops

      auto DataReducer = Kokkos::Sum<I4>(DataCount1);

      // Device array vertical vector
      parallelReduce(
          {NVertLevels},
          KOKKOS_LAMBDA(int K, I4 &LCount) {
             if (Data1DR8(K) != RefR8 + K)
                ++LCount;
          },
          DataReducer);
      TstEval<int>("Get data 1DR8", DataCount1, 0, Err);

      // To save some effort, consolidate some loops
      // Device arrays on cells
      parallelReduce(
          {NCellsSize},
          KOKKOS_LAMBDA(int Cell, I4 &LCount) {
             if (Data1DI4(Cell) != RefI4 + Cell)
                ++LCount;
             for (int K = 0; K < NVertLevels; ++K) {
                int Add2 = Cell * NVertLevels + K;
                if (Data2DI4(Cell, K) != RefI4 + Cell + K)
                   ++LCount;
                for (int Trcr = 0; Trcr < NTracers; ++Trcr) {
                   int Add3 = Trcr * NCellsSize * NVertLevels + Add2;
                   if (Data3DI4(Cell, K, Trcr) != RefI4 + Cell + K + Trcr)
                      ++LCount;
                   for (int TimeLvl = 0; TimeLvl < NTime; ++TimeLvl) {
                      int Add4 =
                          TimeLvl * NTracers * NCellsSize * NVertLevels + Add3;
                      if (Data4DI8(Cell, K, Trcr, TimeLvl) !=
                          RefI8 + Cell + K + Trcr + TimeLvl)
                         ++LCount;
                      for (int Stf = 0; Stf < NStuff; ++Stf) {
                         int Add5 =
                             Stf * NTime * NTracers * NCellsSize * NVertLevels +
                             Add4;
                         if (Data5DR4(Cell, K, Trcr, TimeLvl, Stf) !=
                             RefR4 + Cell + K + Trcr + TimeLvl + Stf)
                            ++LCount;
                      }
                   }
                }
             }
          },
          DataReducer);
      TstEval<int>("Get data all cell device arrays, ptrs", DataCount1, 0, Err);

      // Device arrays on edges
      parallelReduce(
          {NEdgesSize},
          KOKKOS_LAMBDA(int Edge, I4 &LCount) {
             if (Data1DI8(Edge) != RefI8 + Edge)
                ++LCount;
             for (int K = 0; K < NVertLevels; ++K) {
                int Add = Edge * NVertLevels + K;
                if (Data2DI8(Edge, K) != RefI8 + Edge + K)
                   ++LCount;
                if (Data2DR8(Edge, K) != RefR8 + Edge + K)
                   ++LCount;
             }
          },
          DataReducer);
      TstEval<int>("Get data all edge device arrays, ptrs", DataCount1, 0, Err);

      // Device arrays on vertices
      parallelReduce(
          {NVerticesSize},
          KOKKOS_LAMBDA(int Vrtx, I4 &LCount) {
             if (Data1DR4(Vrtx) != RefR4 + Vrtx)
                ++LCount;
             for (int K = 0; K < NVertLevels; ++K) {
                int Add = Vrtx * NVertLevels + K;
                if (Data2DR4(Vrtx, K) != RefR4 + Vrtx + K)
                   ++LCount;
             }
          },
          DataReducer);
      TstEval<int>("Get data all vertex device arrays, ptrs", DataCount1, 0,
                   Err);

      // Test removing a field from a group
      Err1 = Group1D->removeField("Test1DI4");
      TstEval<int>("Remove field call return 1D", Err1, ErrRef, Err);
      GroupTest = Group1D->hasField("Test1DI4");
      TstEval<bool>("Remove field member result 1D", GroupTest, false, Err);

      Err1 = FieldGroup::removeFieldFromGroup("Test2DI4", "FieldGroup2D");
      TstEval<int>("Remove field call return 2D", Err1, ErrRef, Err);
      GroupTest = FieldGroup::isFieldInGroup("Test2DI4", "FieldGroup2D");
      TstEval<bool>("Remove field member result 2D", GroupTest, false, Err);

      // Remove field groups
      Err1 = FieldGroup::destroy("FieldGroup1D");
      TstEval<int>("Remove field group call return 1D", Err1, ErrRef, Err);
      GroupTest = FieldGroup::exists("FieldGroup1D");
      TstEval<bool>("Remove field group verification", GroupTest, false, Err);

      FieldGroup::clear();
      GroupTest = FieldGroup::exists("FieldGroup2D");
      TstEval<bool>("Remove all field groups", GroupTest, false, Err);

      // Destroy a field and check for removal
      bool ShouldExist = false;
      Field::destroy("Test1DI4H");
      Field::destroy("Test1DI4");
      bool FieldExists = Field::exists("Test1DI4H");
      TstEval<bool>("Destroy field 1DI4H", FieldExists, ShouldExist, Err);
      FieldExists = Field::exists("Test1DI4");
      TstEval<bool>("Destroy field 1DI4", FieldExists, ShouldExist, Err);

      // Clear all fields and check a couple for successful removal
      Field::clear();
      FieldExists = Field::exists("Test2DR8H");
      TstEval<bool>("Clear all fields 2DR8H", FieldExists, ShouldExist, Err);
      FieldExists = Field::exists("Test2DR8");
      TstEval<bool>("Clear all fields 2DR8", FieldExists, ShouldExist, Err);
   }

   // Clean up environments
   Dimension::clear();
   Decomp::clear();
   Kokkos::finalize();
   MPI_Finalize();

   if (Err >= 256)
      Err = 255;

   // End of testing
   return Err;
}
//===--- End test driver for IO Field --------------------------------------===/
