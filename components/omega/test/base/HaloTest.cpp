//===-- Test driver for OMEGA Halo -------------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA Halo class
///
/// This driver tests the OMEGA model Halo class, which collects and stores
/// everything needed to perform halo exchanges on any supported Kokkos array
/// defined on a mesh in OMEGA with a given parallel decomposition. This
/// unit test driver tests functionality by creating Kokkos arrays of every
/// type and dimensionality supported in OMEGA, initializing each array based
/// on global IDs of the mesh elememts, performing halo exchanges, and
/// confirming the exchanged arrays are identical to the initial arrays.
///
//
//===-----------------------------------------------------------------------===/

#include "Halo.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "mpi.h"

//------------------------------------------------------------------------------
// This function template performs a single test on a Kokkos array type in a
// given index space. Two Kokkos arrays of the same type and size are input,
// InitArray contains the global IDs of the mesh elements for all the owned and
// halo elements of the array, while TestArray contains the global IDs only in
// the owned elements. The Halo class object, a label describing the test for
// output, an integer to accumulate errors, and optionally the index space of
// the input arrays (default is OnCell) are also input. A halo exchange is
// performed on TestArray, and then TestArray is compared to InitArray. If any
// elements differ the test is a failure and an error is returned.

template <typename T>
void haloExchangeTest(
    OMEGA::Halo *MyHalo,
    T InitArray,  /// Array initialized based on global IDs of mesh elements
    T &TestArray, /// Array only initialized in owned elements
    const char *Label,                          /// Unique label for test
    OMEGA::I4 &TotErr,                          /// Integer to track errors
    OMEGA::MeshElement ThisElem = OMEGA::OnCell /// index space, cell by default
) {

   OMEGA::I4 IErr{0}; // error code

   // Set total array size and ensure arrays are of same size
   OMEGA::I4 NTot = InitArray.size();
   if (NTot != TestArray.size()) {
      LOG_ERROR("HaloTest: {} arrays must be of same size", Label);
      LOG_INFO("HaloTest: {} exchange test FAIL", Label);
      TotErr += -1;
      return;
   }

   // Perform halo exchange
   IErr = MyHalo->exchangeFullArrayHalo(TestArray, ThisElem);
   if (IErr != 0) {
      LOG_ERROR("HaloTest: Error during {} halo exchange", Label);
      LOG_INFO("HaloTest: {} exchange test FAIL", Label);
      TotErr += -1;
      return;
   }

   // Collapse arrays to 1D for easy iteration
   Kokkos::View<typename T::value_type *, typename T::array_layout,
                typename T::memory_space>
       CollapsedInit(InitArray.data(), InitArray.size());
   Kokkos::View<typename T::value_type *, typename T::array_layout,
                typename T::memory_space>
       CollapsedTest(TestArray.data(), TestArray.size());

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
      TotErr += -1;
   }

   return;

} // end haloExchangeTest

//------------------------------------------------------------------------------
// Initialization routine for Halo tests. Calls all the init routines needed
// to create the default Halo.

int initHaloTest() {

   OMEGA::I4 IErr{0};

   // Initialize the machine environment and fetch the default environment
   // pointer and the MPI communicator
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   MPI_Comm DefComm       = DefEnv->getComm();

   // Initialize the IO system
   IErr = OMEGA::IO::init(DefComm);
   if (IErr != 0)
      LOG_ERROR("HaloTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   IErr = OMEGA::Decomp::init();
   if (IErr != 0)
      LOG_ERROR("HaloTest: error initializing default decomposition");

   // Initialize the default halo
   IErr = OMEGA::Halo::init();
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
   OMEGA::I4 TotErr{0};
   OMEGA::I4 IErr{0};

   // Initialize global MPI environment and Kokkos
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   {

      // Call Halo test initialization routine
      IErr = initHaloTest();
      if (IErr != 0)
         LOG_ERROR("HaloTest: initHaloTest error");

      // Retrieve pointer to default halo
      OMEGA::Halo *DefHalo = OMEGA::Halo::getDefault();

      // Retrieve pointer to default decomposition
      OMEGA::Decomp *DefDecomp = OMEGA::Decomp::getDefault();

      OMEGA::I4 NumOwned;
      OMEGA::I4 NumAll;

      // Perform 1DI4 array tests for each index space (cell, edge, and vertex)

      OMEGA::HostArray1DI4 Init1DI4Cell("Init1DI4Cell", DefDecomp->NCellsSize);
      OMEGA::HostArray1DI4 Test1DI4Cell("Test1DI4Cell", DefDecomp->NCellsSize);

      NumOwned     = DefDecomp->NCellsOwned;
      NumAll       = DefDecomp->NCellsAll;
      Init1DI4Cell = DefDecomp->CellIDH;
      OMEGA::deepCopy(Test1DI4Cell, Init1DI4Cell);

      for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
         Test1DI4Cell(ICell) = -1;
      }

      haloExchangeTest(DefHalo, Init1DI4Cell, Test1DI4Cell, "1DI4 Cell",
                       TotErr);

      OMEGA::HostArray1DI4 Init1DI4Edge("Init1DI4Edge", DefDecomp->NEdgesSize);
      OMEGA::HostArray1DI4 Test1DI4Edge("Test1DI4Edge", DefDecomp->NEdgesSize);

      NumOwned     = DefDecomp->NEdgesOwned;
      NumAll       = DefDecomp->NEdgesAll;
      Init1DI4Edge = DefDecomp->EdgeIDH;
      OMEGA::deepCopy(Test1DI4Edge, Init1DI4Edge);

      for (int IEdge = NumOwned; IEdge < NumAll; ++IEdge) {
         Test1DI4Edge(IEdge) = -1;
      }

      haloExchangeTest(DefHalo, Init1DI4Edge, Test1DI4Edge, "1DI4 Edge", TotErr,
                       OMEGA::OnEdge);

      OMEGA::HostArray1DI4 Init1DI4Vertex("Init1DI4Vertex",
                                          DefDecomp->NVerticesSize);
      OMEGA::HostArray1DI4 Test1DI4Vertex("Test1DI4Vertex",
                                          DefDecomp->NVerticesSize);

      NumOwned       = DefDecomp->NVerticesOwned;
      NumAll         = DefDecomp->NVerticesAll;
      Init1DI4Vertex = DefDecomp->VertexIDH;
      OMEGA::deepCopy(Test1DI4Vertex, Init1DI4Vertex);

      for (int IVertex = NumOwned; IVertex < NumAll; ++IVertex) {
         Test1DI4Vertex(IVertex) = -1;
      }
      haloExchangeTest(DefHalo, Init1DI4Vertex, Test1DI4Vertex, "1DI4 Vertex",
                       TotErr, OMEGA::OnVertex);

      // Declaration of variables for remaining tests

      // Random dimension sizes
      OMEGA::I4 N2{20};
      OMEGA::I4 N3{10};
      OMEGA::I4 N4{3};
      OMEGA::I4 N5{2};

      NumOwned = DefDecomp->NCellsOwned;
      NumAll   = DefDecomp->NCellsAll;

      // Declare init and test arrays for all the remaining array types
      OMEGA::HostArray1DI8 Init1DI8("Init1DI8", NumAll);
      OMEGA::HostArray1DR4 Init1DR4("Init1DR4", NumAll);
      OMEGA::HostArray1DR8 Init1DR8("Init1DR8", NumAll);
      OMEGA::HostArray2DI4 Init2DI4("Init2DI4", NumAll, N2);
      OMEGA::HostArray2DI8 Init2DI8("Init2DI8", NumAll, N2);
      OMEGA::HostArray2DR4 Init2DR4("Init2DR4", NumAll, N2);
      OMEGA::HostArray2DR8 Init2DR8("Init2DR8", NumAll, N2);
      OMEGA::HostArray3DI4 Init3DI4("Init3DI4", N3, NumAll, N2);
      OMEGA::HostArray3DI8 Init3DI8("Init3DI8", N3, NumAll, N2);
      OMEGA::HostArray3DR4 Init3DR4("Init3DR4", N3, NumAll, N2);
      OMEGA::HostArray3DR8 Init3DR8("Init3DR8", N3, NumAll, N2);
      OMEGA::HostArray4DI4 Init4DI4("Init4DI4", N4, N3, NumAll, N2);
      OMEGA::HostArray4DI8 Init4DI8("Init4DI8", N4, N3, NumAll, N2);
      OMEGA::HostArray4DR4 Init4DR4("Init4DR4", N4, N3, NumAll, N2);
      OMEGA::HostArray4DR8 Init4DR8("Init4DR8", N4, N3, NumAll, N2);
      OMEGA::HostArray5DI4 Init5DI4("Init5DI4", N5, N4, N3, NumAll, N2);
      OMEGA::HostArray5DI8 Init5DI8("Init5DI8", N5, N4, N3, NumAll, N2);
      OMEGA::HostArray5DR4 Init5DR4("Init5DR4", N5, N4, N3, NumAll, N2);
      OMEGA::HostArray5DR8 Init5DR8("Init5DR8", N5, N4, N3, NumAll, N2);
      OMEGA::HostArray1DI8 Test1DI8("Test1DI8", NumAll);
      OMEGA::HostArray1DR4 Test1DR4("Test1DR4", NumAll);
      OMEGA::HostArray1DR8 Test1DR8("Test1DR8", NumAll);
      OMEGA::HostArray2DI4 Test2DI4("Test2DI4", NumAll, N2);
      OMEGA::HostArray2DI8 Test2DI8("Test2DI8", NumAll, N2);
      OMEGA::HostArray2DR4 Test2DR4("Test2DR4", NumAll, N2);
      OMEGA::HostArray2DR8 Test2DR8("Test2DR8", NumAll, N2);
      OMEGA::HostArray3DI4 Test3DI4("Test3DI4", N3, NumAll, N2);
      OMEGA::HostArray3DI8 Test3DI8("Test3DI8", N3, NumAll, N2);
      OMEGA::HostArray3DR4 Test3DR4("Test3DR4", N3, NumAll, N2);
      OMEGA::HostArray3DR8 Test3DR8("Test3DR8", N3, NumAll, N2);
      OMEGA::HostArray4DI4 Test4DI4("Test4DI4", N4, N3, NumAll, N2);
      OMEGA::HostArray4DI8 Test4DI8("Test4DI8", N4, N3, NumAll, N2);
      OMEGA::HostArray4DR4 Test4DR4("Test4DR4", N4, N3, NumAll, N2);
      OMEGA::HostArray4DR8 Test4DR8("Test4DR8", N4, N3, NumAll, N2);
      OMEGA::HostArray5DI4 Test5DI4("Test5DI4", N5, N4, N3, NumAll, N2);
      OMEGA::HostArray5DI8 Test5DI8("Test5DI8", N5, N4, N3, NumAll, N2);
      OMEGA::HostArray5DR4 Test5DR4("Test5DR4", N5, N4, N3, NumAll, N2);
      OMEGA::HostArray5DR8 Test5DR8("Test5DR8", N5, N4, N3, NumAll, N2);

      // Initialize and run remaining 1D tests
      for (int ICell = 0; ICell < NumAll; ++ICell) {
         OMEGA::I4 NewVal = DefDecomp->CellIDH(ICell);
         Init1DI8(ICell)  = static_cast<OMEGA::I8>(NewVal);
         Init1DR4(ICell)  = static_cast<OMEGA::R4>(NewVal);
         Init1DR8(ICell)  = static_cast<OMEGA::R8>(NewVal);
      }

      OMEGA::deepCopy(Test1DI8, Init1DI8);
      OMEGA::deepCopy(Test1DR4, Init1DR4);
      OMEGA::deepCopy(Test1DR8, Init1DR8);

      for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
         Test1DI8(ICell) = -1;
         Test1DR4(ICell) = -1;
         Test1DR8(ICell) = -1;
      }

      haloExchangeTest(DefHalo, Init1DI8, Test1DI8, "1DI8", TotErr);
      haloExchangeTest(DefHalo, Init1DR4, Test1DR4, "1DR4", TotErr);
      haloExchangeTest(DefHalo, Init1DR8, Test1DR8, "1DR8", TotErr);

      // Initialize and run 2D tests
      for (int ICell = 0; ICell < NumAll; ++ICell) {
         for (int J = 0; J < N2; ++J) {
            OMEGA::I4 NewVal   = (J + 1) * DefDecomp->CellIDH(ICell);
            Init2DI4(ICell, J) = NewVal;
            Init2DI8(ICell, J) = static_cast<OMEGA::I8>(NewVal);
            Init2DR4(ICell, J) = static_cast<OMEGA::R4>(NewVal);
            Init2DR8(ICell, J) = static_cast<OMEGA::R8>(NewVal);
         }
      }

      OMEGA::deepCopy(Test2DI4, Init2DI4);
      OMEGA::deepCopy(Test2DI8, Init2DI8);
      OMEGA::deepCopy(Test2DR4, Init2DR4);
      OMEGA::deepCopy(Test2DR8, Init2DR8);

      for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
         for (int J = 0; J < N2; ++J) {
            Test2DI4(ICell, J) = -1;
            Test2DI8(ICell, J) = -1;
            Test2DR4(ICell, J) = -1;
            Test2DR8(ICell, J) = -1;
         }
      }

      haloExchangeTest(DefHalo, Init2DI4, Test2DI4, "2DI4", TotErr);
      haloExchangeTest(DefHalo, Init2DI8, Test2DI8, "2DI8", TotErr);
      haloExchangeTest(DefHalo, Init2DR4, Test2DR4, "2DR4", TotErr);
      haloExchangeTest(DefHalo, Init2DR8, Test2DR8, "2DR8", TotErr);

      // Initialize and run 3D tests
      for (int K = 0; K < N3; ++K) {
         for (int ICell = 0; ICell < NumAll; ++ICell) {
            for (int J = 0; J < N2; ++J) {
               OMEGA::I4 NewVal = (K + 1) * (J + 1) * DefDecomp->CellIDH(ICell);
               Init3DI4(K, ICell, J) = NewVal;
               Init3DI8(K, ICell, J) = static_cast<OMEGA::I8>(NewVal);
               Init3DR4(K, ICell, J) = static_cast<OMEGA::R4>(NewVal);
               Init3DR8(K, ICell, J) = static_cast<OMEGA::R8>(NewVal);
            }
         }
      }

      OMEGA::deepCopy(Test3DI4, Init3DI4);
      OMEGA::deepCopy(Test3DI8, Init3DI8);
      OMEGA::deepCopy(Test3DR4, Init3DR4);
      OMEGA::deepCopy(Test3DR8, Init3DR8);

      for (int K = 0; K < N3; ++K) {
         for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
            for (int J = 0; J < N2; ++J) {
               Test3DI4(K, ICell, J) = -1;
               Test3DI8(K, ICell, J) = -1;
               Test3DR4(K, ICell, J) = -1;
               Test3DR8(K, ICell, J) = -1;
            }
         }
      }

      haloExchangeTest(DefHalo, Init3DI4, Test3DI4, "3DI4", TotErr);
      haloExchangeTest(DefHalo, Init3DI8, Test3DI8, "3DI8", TotErr);
      haloExchangeTest(DefHalo, Init3DR4, Test3DR4, "3DR4", TotErr);
      haloExchangeTest(DefHalo, Init3DR8, Test3DR8, "3DR8", TotErr);

      // Initialize and run 4D tests
      for (int L = 0; L < N4; ++L) {
         for (int K = 0; K < N3; ++K) {
            for (int ICell = 0; ICell < NumAll; ++ICell) {
               for (int J = 0; J < N2; ++J) {
                  OMEGA::I4 NewVal =
                      (L + 1) * (K + 1) * (J + 1) * DefDecomp->CellIDH(ICell);
                  Init4DI4(L, K, ICell, J) = NewVal;
                  Init4DI8(L, K, ICell, J) = static_cast<OMEGA::I8>(NewVal);
                  Init4DR4(L, K, ICell, J) = static_cast<OMEGA::R4>(NewVal);
                  Init4DR8(L, K, ICell, J) = static_cast<OMEGA::R8>(NewVal);
               }
            }
         }
      }

      OMEGA::deepCopy(Test4DI4, Init4DI4);
      OMEGA::deepCopy(Test4DI8, Init4DI8);
      OMEGA::deepCopy(Test4DR4, Init4DR4);
      OMEGA::deepCopy(Test4DR8, Init4DR8);

      for (int L = 0; L < N4; ++L) {
         for (int K = 0; K < N3; ++K) {
            for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
               for (int J = 0; J < N2; ++J) {
                  Test4DI4(L, K, ICell, J) = -1;
                  Test4DI8(L, K, ICell, J) = -1;
                  Test4DR4(L, K, ICell, J) = -1;
                  Test4DR8(L, K, ICell, J) = -1;
               }
            }
         }
      }

      haloExchangeTest(DefHalo, Init4DI4, Test4DI4, "4DI4", TotErr);
      haloExchangeTest(DefHalo, Init4DI8, Test4DI8, "4DI8", TotErr);
      haloExchangeTest(DefHalo, Init4DR4, Test4DR4, "4DR4", TotErr);
      haloExchangeTest(DefHalo, Init4DR8, Test4DR8, "4DR8", TotErr);

      // Initialize and run 5D tests
      for (int M = 0; M < N5; ++M) {
         for (int L = 0; L < N4; ++L) {
            for (int K = 0; K < N3; ++K) {
               for (int ICell = 0; ICell < NumAll; ++ICell) {
                  for (int J = 0; J < N2; ++J) {
                     OMEGA::I4 NewVal = (M + 1) * (L + 1) * (K + 1) * (J + 1) *
                                        DefDecomp->CellIDH(ICell);
                     Init5DI4(M, L, K, ICell, J) = NewVal;
                     Init5DI8(M, L, K, ICell, J) =
                         static_cast<OMEGA::I8>(NewVal);
                     Init5DR4(M, L, K, ICell, J) =
                         static_cast<OMEGA::R4>(NewVal);
                     Init5DR8(M, L, K, ICell, J) =
                         static_cast<OMEGA::R8>(NewVal);
                  }
               }
            }
         }
      }

      OMEGA::deepCopy(Test5DI4, Init5DI4);
      OMEGA::deepCopy(Test5DI8, Init5DI8);
      OMEGA::deepCopy(Test5DR4, Init5DR4);
      OMEGA::deepCopy(Test5DR8, Init5DR8);

      for (int M = 0; M < N5; ++M) {
         for (int L = 0; L < N4; ++L) {
            for (int K = 0; K < N3; ++K) {
               for (int ICell = NumOwned; ICell < NumAll; ++ICell) {
                  for (int J = 0; J < N2; ++J) {
                     Test5DI4(M, L, K, ICell, J) = -1;
                     Test5DI8(M, L, K, ICell, J) = -1;
                     Test5DR4(M, L, K, ICell, J) = -1;
                     Test5DR8(M, L, K, ICell, J) = -1;
                  }
               }
            }
         }
      }

      haloExchangeTest(DefHalo, Init5DI4, Test5DI4, "5DI4", TotErr);
      haloExchangeTest(DefHalo, Init5DI8, Test5DI8, "5DI8", TotErr);
      haloExchangeTest(DefHalo, Init5DR4, Test5DR4, "5DR4", TotErr);
      haloExchangeTest(DefHalo, Init5DR8, Test5DR8, "5DR8", TotErr);

      // Memory clean up
      OMEGA::Decomp::clear();
      OMEGA::MachEnv::removeAll();

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
