//===-- Test driver for OMEGA Halo -------------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Halo class
///
/// This driver tests the OMEGA model Halo class, which collects and stores
/// everything needed to perform halo exchanges on any supported Kokkos array
/// defined on a mesh in OMEGA with a given parallel decomposition. This
/// unit test driver tests functionality by creating Kokkos arrays of every
/// type, rank, and memory space supported in OMEGA, initializing each array
/// based on global IDs of the mesh elememts, performing halo exchanges, and
/// confirming the exchanged arrays are identical to the initial arrays.
///
//
//===-----------------------------------------------------------------------===/

#include "Halo.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "mpi.h"

using namespace OMEGA;

//------------------------------------------------------------------------------
// This function template performs a single test on a Kokkos array type in a
// given index space. Two Kokkos arrays of the same type and size are input,
// InitArray contains the global IDs of the mesh elements for all the owned and
// halo elements of the array, while TestArray contains the global IDs only in
// the owned elements. The Halo class object, a label describing the test for
// output, and optionally the index space of the input arrays (default is
// OnCell) are also input. An error code is returned as output. A halo exchange
// is performed on TestArray, and then TestArray is compared to InitArray. If
// any elements differ the test is a failure and an error is returned.

template <typename T>
int haloExchangeTest(
    Halo *MyHalo,
    T InitArray,  /// Array initialized based on global IDs of mesh elements
    T &TestArray, /// Array only initialized in owned elements
    const char *Label,            /// Unique label for test
    MeshElement ThisElem = OnCell /// index space, cell by default
) {

   I4 IErr   = 0; // error code
   I4 RetErr = 0; // return error code

   // Set total array size and ensure arrays are of same size
   I4 NTot = InitArray.size();
   if (NTot != TestArray.size()) {
      LOG_ERROR("HaloTest: {} arrays must be of same size", Label);
      LOG_INFO("HaloTest: {} exchange test FAIL", Label);
      RetErr += 1;
      return RetErr;
   }

   // Perform halo exchange
   IErr = MyHalo->exchangeFullArrayHalo(TestArray, ThisElem);
   if (IErr != 0) {
      LOG_ERROR("HaloTest: Error during {} halo exchange", Label);
      LOG_INFO("HaloTest: {} exchange test FAIL", Label);
      RetErr += 1;
      return RetErr;
   }

   auto TestArrayH = createHostMirrorCopy(TestArray);
   auto InitArrayH = createHostMirrorCopy(InitArray);

   // Collapse arrays to 1D for easy iteration
   Kokkos::View<typename T::value_type *, typename T::array_layout,
                HostMemSpace>
       CollapsedInit(InitArrayH.data(), InitArrayH.size());
   Kokkos::View<typename T::value_type *, typename T::array_layout,
                HostMemSpace>
       CollapsedTest(TestArrayH.data(), TestArrayH.size());

   // Confirm all elements are identical, if not set error code
   // and break out of loop
   for (int N = 0; N < NTot; ++N) {
      if (CollapsedInit(N) != CollapsedTest(N)) {
         IErr = -1;
         break;
      }
   }

   if (IErr == 0) {
      LOG_INFO("HaloTest: {} exchange test PASS", Label);
   } else {
      LOG_INFO("HaloTest: {} exchange test FAIL", Label);
      RetErr += 1;
   }

   return RetErr;

} // end haloExchangeTest

//------------------------------------------------------------------------------
// Initialization routine for Halo tests. Calls all the init routines needed
// to create the default Halo.

int initHaloTest() {

   I4 IErr{0};

   // Initialize the machine environment and fetch the default environment
   // pointer and the MPI communicator
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the logging system
   initLogging(DefEnv);

   // Open config file
   Config("Omega");
   IErr = Config::readAll("omega.yml");
   if (IErr != 0) {
      LOG_CRITICAL("HaloTest: Error reading config file");
      return IErr;
   }

   // Initialize the IO system
   IErr = IO::init(DefComm);
   if (IErr != 0)
      LOG_ERROR("HaloTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   IErr = Decomp::init();
   if (IErr != 0)
      LOG_ERROR("HaloTest: error initializing default decomposition");

   // Initialize the default halo
   IErr = Halo::init();
   if (IErr != 0)
      LOG_ERROR("HaloTest: error initializing default halo");

   return IErr;

} // end initHaloTest

//------------------------------------------------------------------------------
// The test driver. Performs halo exchange tests of all index spaces and all
// supported Kokkos array types. For each test, an initial array is set based on
// Global IDs of the mesh elements in the given index space for all owned and
// halo elements, and is copied into the test array. The test array halo
// elements are then set to junk values and a halo exchange is performed, which
// if successful will fetch the proper values from neighboring test arrays such
// that the test array and initial array are equivalent.

int main(int argc, char *argv[]) {

   // Error tracking variables
   I4 TotErr = 0;
   I4 IErr   = 0;

   // Initialize global MPI environment and Kokkos
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {

      // Call Halo test initialization routine
      IErr = initHaloTest();
      if (IErr != 0)
         LOG_ERROR("HaloTest: initHaloTest error");

      // Retrieve pointer to default halo
      Halo *DefHalo = Halo::getDefault();

      // Retrieve pointer to default decomposition
      Decomp *DefDecomp = Decomp::getDefault();

      I4 NumOwned;
      I4 NumAll;

      // Perform 1DI4 array tests for each index space (cell, edge, and vertex)

      HostArray1DI4 Init1DI4CellH("Init1DI4CellH", DefDecomp->NCellsSize);
      HostArray1DI4 Test1DI4CellH("Test1DI4CellH", DefDecomp->NCellsSize);
      Array1DI4 Init1DI4Cell("Init1DI4Cell", DefDecomp->NCellsSize);
      Array1DI4 Test1DI4Cell("Test1DI4Cell", DefDecomp->NCellsSize);

      NumOwned      = DefDecomp->NCellsOwned;
      NumAll        = DefDecomp->NCellsAll;
      Init1DI4CellH = DefDecomp->CellIDH;
      deepCopy(Test1DI4CellH, Init1DI4CellH);

      for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
         Test1DI4CellH(ICell) = -1;
      }

      deepCopy(Init1DI4Cell, Init1DI4CellH);
      deepCopy(Test1DI4Cell, Test1DI4CellH);

      TotErr +=
          haloExchangeTest(DefHalo, Init1DI4CellH, Test1DI4CellH, "1DI4CellH");
      TotErr +=
          haloExchangeTest(DefHalo, Init1DI4Cell, Test1DI4Cell, "1DI4Cell");

      HostArray1DI4 Init1DI4EdgeH("Init1DI4EdgeH", DefDecomp->NEdgesSize);
      HostArray1DI4 Test1DI4EdgeH("Test1DI4EdgeH", DefDecomp->NEdgesSize);
      Array1DI4 Init1DI4Edge("Init1DI4Edge", DefDecomp->NEdgesSize);
      Array1DI4 Test1DI4Edge("Test1DI4Edge", DefDecomp->NEdgesSize);

      NumOwned      = DefDecomp->NEdgesOwned;
      NumAll        = DefDecomp->NEdgesAll;
      Init1DI4EdgeH = DefDecomp->EdgeIDH;
      deepCopy(Test1DI4EdgeH, Init1DI4EdgeH);

      for (int IEdge = NumOwned; IEdge < NumAll; ++IEdge) {
         Test1DI4EdgeH(IEdge) = -1;
      }

      deepCopy(Init1DI4Edge, Init1DI4EdgeH);
      deepCopy(Test1DI4Edge, Test1DI4EdgeH);
      TotErr += haloExchangeTest(DefHalo, Init1DI4EdgeH, Test1DI4EdgeH,
                                 "1DI4EdgeH", OnEdge);
      TotErr += haloExchangeTest(DefHalo, Init1DI4Edge, Test1DI4Edge,
                                 "1DI4Edge", OnEdge);

      HostArray1DI4 Init1DI4VertexH("Init1DI4VertexH",
                                    DefDecomp->NVerticesSize);
      HostArray1DI4 Test1DI4VertexH("Test1DI4VertexH",
                                    DefDecomp->NVerticesSize);
      Array1DI4 Init1DI4Vertex("Init1DI4Vertex", DefDecomp->NVerticesSize);
      Array1DI4 Test1DI4Vertex("Test1DI4Vertex", DefDecomp->NVerticesSize);

      NumOwned        = DefDecomp->NVerticesOwned;
      NumAll          = DefDecomp->NVerticesAll;
      Init1DI4VertexH = DefDecomp->VertexIDH;
      deepCopy(Test1DI4VertexH, Init1DI4VertexH);

      for (int IVertex = NumOwned; IVertex < NumAll; ++IVertex) {
         Test1DI4VertexH(IVertex) = -1;
      }

      deepCopy(Init1DI4Vertex, Init1DI4VertexH);
      deepCopy(Test1DI4Vertex, Test1DI4VertexH);
      TotErr += haloExchangeTest(DefHalo, Init1DI4VertexH, Test1DI4VertexH,
                                 "1DI4VertexH", OnVertex);
      TotErr += haloExchangeTest(DefHalo, Init1DI4Vertex, Test1DI4Vertex,
                                 "1DI4Vertex", OnVertex);

      // Declaration of variables for remaining tests

      // Random dimension sizes
      I4 N2{20};
      I4 N3{10};
      I4 N4{3};
      I4 N5{2};

      I4 NumGlobe = DefDecomp->NCellsGlobal;
      NumOwned    = DefDecomp->NCellsOwned;
      NumAll      = DefDecomp->NCellsAll;

      // Declare init and test arrays for all the remaining array types
      HostArray1DI8 Init1DI8H("Init1DI8H", NumAll);
      HostArray1DR4 Init1DR4H("Init1DR4H", NumAll);
      HostArray1DR8 Init1DR8H("Init1DR8H", NumAll);
      HostArray2DI4 Init2DI4H("Init2DI4H", NumAll, N2);
      HostArray2DI8 Init2DI8H("Init2DI8H", NumAll, N2);
      HostArray2DR4 Init2DR4H("Init2DR4H", NumAll, N2);
      HostArray2DR8 Init2DR8H("Init2DR8H", NumAll, N2);
      HostArray3DI4 Init3DI4H("Init3DI4H", N3, NumAll, N2);
      HostArray3DI8 Init3DI8H("Init3DI8H", N3, NumAll, N2);
      HostArray3DR4 Init3DR4H("Init3DR4H", N3, NumAll, N2);
      HostArray3DR8 Init3DR8H("Init3DR8H", N3, NumAll, N2);
      HostArray4DI4 Init4DI4H("Init4DI4H", N4, N3, NumAll, N2);
      HostArray4DI8 Init4DI8H("Init4DI8H", N4, N3, NumAll, N2);
      HostArray4DR4 Init4DR4H("Init4DR4H", N4, N3, NumAll, N2);
      HostArray4DR8 Init4DR8H("Init4DR8H", N4, N3, NumAll, N2);
      HostArray5DI4 Init5DI4H("Init5DI4H", N5, N4, N3, NumAll, N2);
      HostArray5DI8 Init5DI8H("Init5DI8H", N5, N4, N3, NumAll, N2);
      HostArray5DR4 Init5DR4H("Init5DR4H", N5, N4, N3, NumAll, N2);
      HostArray5DR8 Init5DR8H("Init5DR8H", N5, N4, N3, NumAll, N2);
      HostArray1DI8 Test1DI8H("Test1DI8H", NumAll);
      HostArray1DR4 Test1DR4H("Test1DR4H", NumAll);
      HostArray1DR8 Test1DR8H("Test1DR8H", NumAll);
      HostArray2DI4 Test2DI4H("Test2DI4H", NumAll, N2);
      HostArray2DI8 Test2DI8H("Test2DI8H", NumAll, N2);
      HostArray2DR4 Test2DR4H("Test2DR4H", NumAll, N2);
      HostArray2DR8 Test2DR8H("Test2DR8H", NumAll, N2);
      HostArray3DI4 Test3DI4H("Test3DI4H", N3, NumAll, N2);
      HostArray3DI8 Test3DI8H("Test3DI8H", N3, NumAll, N2);
      HostArray3DR4 Test3DR4H("Test3DR4H", N3, NumAll, N2);
      HostArray3DR8 Test3DR8H("Test3DR8H", N3, NumAll, N2);
      HostArray4DI4 Test4DI4H("Test4DI4H", N4, N3, NumAll, N2);
      HostArray4DI8 Test4DI8H("Test4DI8H", N4, N3, NumAll, N2);
      HostArray4DR4 Test4DR4H("Test4DR4H", N4, N3, NumAll, N2);
      HostArray4DR8 Test4DR8H("Test4DR8H", N4, N3, NumAll, N2);
      HostArray5DI4 Test5DI4H("Test5DI4H", N5, N4, N3, NumAll, N2);
      HostArray5DI8 Test5DI8H("Test5DI8H", N5, N4, N3, NumAll, N2);
      HostArray5DR4 Test5DR4H("Test5DR4H", N5, N4, N3, NumAll, N2);
      HostArray5DR8 Test5DR8H("Test5DR8H", N5, N4, N3, NumAll, N2);

      Array1DI8 Init1DI8("Init1DI8", NumAll);
      Array1DR4 Init1DR4("Init1DR4", NumAll);
      Array1DR8 Init1DR8("Init1DR8", NumAll);
      Array2DI4 Init2DI4("Init2DI4", NumAll, N2);
      Array2DI8 Init2DI8("Init2DI8", NumAll, N2);
      Array2DR4 Init2DR4("Init2DR4", NumAll, N2);
      Array2DR8 Init2DR8("Init2DR8", NumAll, N2);
      Array3DI4 Init3DI4("Init3DI4", N3, NumAll, N2);
      Array3DI8 Init3DI8("Init3DI8", N3, NumAll, N2);
      Array3DR4 Init3DR4("Init3DR4", N3, NumAll, N2);
      Array3DR8 Init3DR8("Init3DR8", N3, NumAll, N2);
      Array4DI4 Init4DI4("Init4DI4", N4, N3, NumAll, N2);
      Array4DI8 Init4DI8("Init4DI8", N4, N3, NumAll, N2);
      Array4DR4 Init4DR4("Init4DR4", N4, N3, NumAll, N2);
      Array4DR8 Init4DR8("Init4DR8", N4, N3, NumAll, N2);
      Array5DI4 Init5DI4("Init5DI4", N5, N4, N3, NumAll, N2);
      Array5DI8 Init5DI8("Init5DI8", N5, N4, N3, NumAll, N2);
      Array5DR4 Init5DR4("Init5DR4", N5, N4, N3, NumAll, N2);
      Array5DR8 Init5DR8("Init5DR8", N5, N4, N3, NumAll, N2);
      Array1DI8 Test1DI8("Test1DI8", NumAll);
      Array1DR4 Test1DR4("Test1DR4", NumAll);
      Array1DR8 Test1DR8("Test1DR8", NumAll);
      Array2DI4 Test2DI4("Test2DI4", NumAll, N2);
      Array2DI8 Test2DI8("Test2DI8", NumAll, N2);
      Array2DR4 Test2DR4("Test2DR4", NumAll, N2);
      Array2DR8 Test2DR8("Test2DR8", NumAll, N2);
      Array3DI4 Test3DI4("Test3DI4", N3, NumAll, N2);
      Array3DI8 Test3DI8("Test3DI8", N3, NumAll, N2);
      Array3DR4 Test3DR4("Test3DR4", N3, NumAll, N2);
      Array3DR8 Test3DR8("Test3DR8", N3, NumAll, N2);
      Array4DI4 Test4DI4("Test4DI4", N4, N3, NumAll, N2);
      Array4DI8 Test4DI8("Test4DI8", N4, N3, NumAll, N2);
      Array4DR4 Test4DR4("Test4DR4", N4, N3, NumAll, N2);
      Array4DR8 Test4DR8("Test4DR8", N4, N3, NumAll, N2);
      Array5DI4 Test5DI4("Test5DI4", N5, N4, N3, NumAll, N2);
      Array5DI8 Test5DI8("Test5DI8", N5, N4, N3, NumAll, N2);
      Array5DR4 Test5DR4("Test5DR4", N5, N4, N3, NumAll, N2);
      Array5DR8 Test5DR8("Test5DR8", N5, N4, N3, NumAll, N2);

      I8 I8Val = 1000000000;
      R4 R4Val = 0.1234;
      R8 R8Val = 0.123456789;

      // Initialize and run remaining 1D tests
      for (int ICell = 0; ICell < NumAll; ++ICell) {
         I4 NewVal        = DefDecomp->CellIDH(ICell);
         Init1DI8H(ICell) = static_cast<I8>(NewVal) * I8Val;
         Init1DR4H(ICell) = static_cast<R4>(NewVal) + R4Val;
         Init1DR8H(ICell) = static_cast<R8>(NewVal) + R8Val;
      }

      deepCopy(Test1DI8H, Init1DI8H);
      deepCopy(Test1DR4H, Init1DR4H);
      deepCopy(Test1DR8H, Init1DR8H);

      for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
         Test1DI8H(ICell) = -1;
         Test1DR4H(ICell) = -1;
         Test1DR8H(ICell) = -1;
      }

      deepCopy(Init1DI8, Init1DI8H);
      deepCopy(Init1DR4, Init1DR4H);
      deepCopy(Init1DR8, Init1DR8H);
      deepCopy(Test1DI8, Test1DI8H);
      deepCopy(Test1DR4, Test1DR4H);
      deepCopy(Test1DR8, Test1DR8H);

      TotErr += haloExchangeTest(DefHalo, Init1DI8H, Test1DI8H, "1DI8H");
      TotErr += haloExchangeTest(DefHalo, Init1DR4H, Test1DR4H, "1DR4H");
      TotErr += haloExchangeTest(DefHalo, Init1DR8H, Test1DR8H, "1DR8H");
      TotErr += haloExchangeTest(DefHalo, Init1DI8, Test1DI8, "1DI8");
      TotErr += haloExchangeTest(DefHalo, Init1DR4, Test1DR4, "1DR4");
      TotErr += haloExchangeTest(DefHalo, Init1DR8, Test1DR8, "1DR8");

      // Initialize and run 2D tests
      for (int ICell = 0; ICell < NumAll; ++ICell) {
         for (int J = 0; J < N2; ++J) {
            I4 NewVal           = DefDecomp->CellIDH(ICell) * N2 + J;
            Init2DI4H(ICell, J) = NewVal;
            Init2DI8H(ICell, J) = static_cast<I8>(NewVal) * I8Val;
            Init2DR4H(ICell, J) = static_cast<R4>(NewVal) + R4Val;
            Init2DR8H(ICell, J) = static_cast<R8>(NewVal) + R8Val;
         }
      }

      deepCopy(Test2DI4H, Init2DI4H);
      deepCopy(Test2DI8H, Init2DI8H);
      deepCopy(Test2DR4H, Init2DR4H);
      deepCopy(Test2DR8H, Init2DR8H);

      for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
         for (int J = 0; J < N2; ++J) {
            Test2DI4H(ICell, J) = -1;
            Test2DI8H(ICell, J) = -1;
            Test2DR4H(ICell, J) = -1;
            Test2DR8H(ICell, J) = -1;
         }
      }

      deepCopy(Init2DI4, Init2DI4H);
      deepCopy(Init2DI8, Init2DI8H);
      deepCopy(Init2DR4, Init2DR4H);
      deepCopy(Init2DR8, Init2DR8H);
      deepCopy(Test2DI4, Test2DI4H);
      deepCopy(Test2DI8, Test2DI8H);
      deepCopy(Test2DR4, Test2DR4H);
      deepCopy(Test2DR8, Test2DR8H);

      TotErr += haloExchangeTest(DefHalo, Init2DI4H, Test2DI4H, "2DI4H");
      TotErr += haloExchangeTest(DefHalo, Init2DI8H, Test2DI8H, "2DI8H");
      TotErr += haloExchangeTest(DefHalo, Init2DR4H, Test2DR4H, "2DR4H");
      TotErr += haloExchangeTest(DefHalo, Init2DR8H, Test2DR8H, "2DR8H");
      TotErr += haloExchangeTest(DefHalo, Init2DI4, Test2DI4, "2DI4");
      TotErr += haloExchangeTest(DefHalo, Init2DI8, Test2DI8, "2DI8");
      TotErr += haloExchangeTest(DefHalo, Init2DR4, Test2DR4, "2DR4");
      TotErr += haloExchangeTest(DefHalo, Init2DR8, Test2DR8, "2DR8");

      // Initialize and run 3D tests
      for (int K = 0; K < N3; ++K) {
         for (int ICell = 0; ICell < NumAll; ++ICell) {
            for (int J = 0; J < N2; ++J) {
               I4 NewVal = (K * NumGlobe + DefDecomp->CellIDH(ICell)) * N2 + J;
               Init3DI4H(K, ICell, J) = NewVal;
               Init3DI8H(K, ICell, J) = static_cast<I8>(NewVal) * I8Val;
               Init3DR4H(K, ICell, J) = static_cast<R4>(NewVal) + R4Val;
               Init3DR8H(K, ICell, J) = static_cast<R8>(NewVal) + R8Val;
            }
         }
      }

      deepCopy(Test3DI4H, Init3DI4H);
      deepCopy(Test3DI8H, Init3DI8H);
      deepCopy(Test3DR4H, Init3DR4H);
      deepCopy(Test3DR8H, Init3DR8H);

      for (int K = 0; K < N3; ++K) {
         for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
            for (int J = 0; J < N2; ++J) {
               Test3DI4H(K, ICell, J) = -1;
               Test3DI8H(K, ICell, J) = -1;
               Test3DR4H(K, ICell, J) = -1;
               Test3DR8H(K, ICell, J) = -1;
            }
         }
      }

      deepCopy(Init3DI4, Init3DI4H);
      deepCopy(Init3DI8, Init3DI8H);
      deepCopy(Init3DR4, Init3DR4H);
      deepCopy(Init3DR8, Init3DR8H);
      deepCopy(Test3DI4, Test3DI4H);
      deepCopy(Test3DI8, Test3DI8H);
      deepCopy(Test3DR4, Test3DR4H);
      deepCopy(Test3DR8, Test3DR8H);

      TotErr += haloExchangeTest(DefHalo, Init3DI4H, Test3DI4H, "3DI4H");
      TotErr += haloExchangeTest(DefHalo, Init3DI8H, Test3DI8H, "3DI8H");
      TotErr += haloExchangeTest(DefHalo, Init3DR4H, Test3DR4H, "3DR4H");
      TotErr += haloExchangeTest(DefHalo, Init3DR8H, Test3DR8H, "3DR8H");
      TotErr += haloExchangeTest(DefHalo, Init3DI4, Test3DI4, "3DI4");
      TotErr += haloExchangeTest(DefHalo, Init3DI8, Test3DI8, "3DI8");
      TotErr += haloExchangeTest(DefHalo, Init3DR4, Test3DR4, "3DR4");
      TotErr += haloExchangeTest(DefHalo, Init3DR8, Test3DR8, "3DR8");

      // Initialize and run 4D tests
      for (int L = 0; L < N4; ++L) {
         for (int K = 0; K < N3; ++K) {
            for (int ICell = 0; ICell < NumAll; ++ICell) {
               for (int J = 0; J < N2; ++J) {
                  I4 NewVal =
                      ((L * N3 + K) * NumGlobe + DefDecomp->CellIDH(ICell)) *
                          N2 +
                      J;
                  Init4DI4H(L, K, ICell, J) = NewVal;
                  Init4DI8H(L, K, ICell, J) = static_cast<I8>(NewVal) * I8Val;
                  Init4DR4H(L, K, ICell, J) = static_cast<R4>(NewVal) + R4Val;
                  Init4DR8H(L, K, ICell, J) = static_cast<R8>(NewVal) + R8Val;
               }
            }
         }
      }

      deepCopy(Test4DI4H, Init4DI4H);
      deepCopy(Test4DI8H, Init4DI8H);
      deepCopy(Test4DR4H, Init4DR4H);
      deepCopy(Test4DR8H, Init4DR8H);

      for (int L = 0; L < N4; ++L) {
         for (int K = 0; K < N3; ++K) {
            for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
               for (int J = 0; J < N2; ++J) {
                  Test4DI4H(L, K, ICell, J) = -1;
                  Test4DI8H(L, K, ICell, J) = -1;
                  Test4DR4H(L, K, ICell, J) = -1;
                  Test4DR8H(L, K, ICell, J) = -1;
               }
            }
         }
      }

      deepCopy(Init4DI4, Init4DI4H);
      deepCopy(Init4DI8, Init4DI8H);
      deepCopy(Init4DR4, Init4DR4H);
      deepCopy(Init4DR8, Init4DR8H);
      deepCopy(Test4DI4, Test4DI4H);
      deepCopy(Test4DI8, Test4DI8H);
      deepCopy(Test4DR4, Test4DR4H);
      deepCopy(Test4DR8, Test4DR8H);

      TotErr += haloExchangeTest(DefHalo, Init4DI4H, Test4DI4H, "4DI4H");
      TotErr += haloExchangeTest(DefHalo, Init4DI8H, Test4DI8H, "4DI8H");
      TotErr += haloExchangeTest(DefHalo, Init4DR4H, Test4DR4H, "4DR4H");
      TotErr += haloExchangeTest(DefHalo, Init4DR8H, Test4DR8H, "4DR8H");
      TotErr += haloExchangeTest(DefHalo, Init4DI4, Test4DI4, "4DI4");
      TotErr += haloExchangeTest(DefHalo, Init4DI8, Test4DI8, "4DI8");
      TotErr += haloExchangeTest(DefHalo, Init4DR4, Test4DR4, "4DR4");
      TotErr += haloExchangeTest(DefHalo, Init4DR8, Test4DR8, "4DR8");

      // Initialize and run 5D tests
      for (int M = 0; M < N5; ++M) {
         for (int L = 0; L < N4; ++L) {
            for (int K = 0; K < N3; ++K) {
               for (int ICell = 0; ICell < NumAll; ++ICell) {
                  for (int J = 0; J < N2; ++J) {
                     I4 NewVal = (((M * N4 + L) * N3 + K) * NumGlobe +
                                  DefDecomp->CellIDH(ICell)) *
                                     N2 +
                                 J;
                     Init5DI4H(M, L, K, ICell, J) = NewVal;
                     Init5DI8H(M, L, K, ICell, J) =
                         static_cast<I8>(NewVal) * I8Val;
                     Init5DR4H(M, L, K, ICell, J) =
                         static_cast<R4>(NewVal) + R4Val;
                     Init5DR8H(M, L, K, ICell, J) =
                         static_cast<R8>(NewVal) + R8Val;
                  }
               }
            }
         }
      }

      deepCopy(Test5DI4H, Init5DI4H);
      deepCopy(Test5DI8H, Init5DI8H);
      deepCopy(Test5DR4H, Init5DR4H);
      deepCopy(Test5DR8H, Init5DR8H);

      for (int M = 0; M < N5; ++M) {
         for (int L = 0; L < N4; ++L) {
            for (int K = 0; K < N3; ++K) {
               for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
                  for (int J = 0; J < N2; ++J) {
                     Test5DI4H(M, L, K, ICell, J) = -1;
                     Test5DI8H(M, L, K, ICell, J) = -1;
                     Test5DR4H(M, L, K, ICell, J) = -1;
                     Test5DR8H(M, L, K, ICell, J) = -1;
                  }
               }
            }
         }
      }

      deepCopy(Init5DI4, Init5DI4H);
      deepCopy(Init5DI8, Init5DI8H);
      deepCopy(Init5DR4, Init5DR4H);
      deepCopy(Init5DR8, Init5DR8H);
      deepCopy(Test5DI4, Test5DI4H);
      deepCopy(Test5DI8, Test5DI8H);
      deepCopy(Test5DR4, Test5DR4H);
      deepCopy(Test5DR8, Test5DR8H);

      TotErr += haloExchangeTest(DefHalo, Init5DI4H, Test5DI4H, "5DI4H");
      TotErr += haloExchangeTest(DefHalo, Init5DI8H, Test5DI8H, "5DI8H");
      TotErr += haloExchangeTest(DefHalo, Init5DR4H, Test5DR4H, "5DR4H");
      TotErr += haloExchangeTest(DefHalo, Init5DR8H, Test5DR8H, "5DR8H");
      TotErr += haloExchangeTest(DefHalo, Init5DI4, Test5DI4, "5DI4");
      TotErr += haloExchangeTest(DefHalo, Init5DI8, Test5DI8, "5DI8");
      TotErr += haloExchangeTest(DefHalo, Init5DR4, Test5DR4, "5DR4");
      TotErr += haloExchangeTest(DefHalo, Init5DR8, Test5DR8, "5DR8");

      // Memory clean up
      Halo::clear();
      Decomp::clear();
      MachEnv::removeAll();

      if (TotErr == 0) {
         LOG_INFO("HaloTest: Successful completion");
      } else {
         LOG_INFO("HaloTest: Failed");
      }
   }
   Kokkos::finalize();
   MPI_Finalize();

   if (TotErr >= 256)
      TotErr = 255;

   return TotErr;

} // end of main
//===-----------------------------------------------------------------------===/
