//===-- Test driver for OMEGA Dimension class ---------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA Dimension class
///
/// This driver tests the capabilities for OMEGA to manage dimension information
/// for use by Fields and IO.
//
//===-----------------------------------------------------------------------===/

#include "Dimension.h"
#include "Config.h"
#include "DataTypes.h"
#include "Decomp.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "mpi.h"

#include <memory>
#include <vector>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Initialization routine to create dimensions
int initDimensionTest() {

   int Err = 0;

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
   OMEGA::Config("omega");
   Err = OMEGA::Config::readAll("omega.yml");
   if (Err != 0) {
      LOG_CRITICAL("DimensionTest: Error reading config file");
      return Err;
   }

   // Initialize decomposition
   Decomp::init();
   Decomp *DefDecomp = Decomp::getDefault();

   // Create offsets for dimension definition
   I4 NCellsSize   = DefDecomp->NCellsSize;
   I4 NCellsOwned  = DefDecomp->NCellsOwned;
   I4 NCellsGlobal = DefDecomp->NCellsGlobal;
   I4 NEdgesSize   = DefDecomp->NEdgesSize;
   I4 NEdgesOwned  = DefDecomp->NEdgesOwned;
   I4 NEdgesGlobal = DefDecomp->NEdgesGlobal;
   HostArray1DI4 CellOffset("NCellsOffset", NCellsSize);
   HostArray1DI4 EdgeOffset("NEdgesOffset", NEdgesSize);
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

   // Define dimensions
   std::shared_ptr<Dimension> CellDim =
       Dimension::create("NCells", NCellsGlobal, NCellsSize, CellOffset);
   std::shared_ptr<Dimension> EdgeDim =
       Dimension::create("NEdges", NEdgesGlobal, NEdgesSize, EdgeOffset);
   I4 NVertLevels = 100;
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLevels", NVertLevels);

   return Err;

} // End initialization of Dimensions

//------------------------------------------------------------------------------
// Main driver for testing Dimension class interfaces.

int main(int argc, char **argv) {

   int Err = 0;

   // Initialize the global MPI environment
   // We do not actually use message passing but need to test the
   // array types and behavior within the distributed environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();

   {
      int Err1 = 0;
      // Call initialization to create sample dimensions
      Err1 = initDimensionTest();
      if (Err1 != 0) {
         LOG_ERROR("Error initializating dimensions");
         ++Err;
      }

      // Retrieve the default domain decomposition and set reference values
      Decomp *DefDecomp = Decomp::getDefault();

      I4 NCellsLocRef = DefDecomp->NCellsSize;
      I4 NCellsGlbRef = DefDecomp->NCellsGlobal;
      I4 NCellsOwned  = DefDecomp->NCellsOwned;
      I4 NEdgesLocRef = DefDecomp->NEdgesSize;
      I4 NEdgesGlbRef = DefDecomp->NEdgesGlobal;
      I4 NEdgesOwned  = DefDecomp->NEdgesOwned;
      I4 NVertLvlsRef = 100;

      // Retrieve the number of defined dimensions
      I4 NDimsRef = 3;
      I4 NDims    = Dimension::getNumDefinedDims();
      if (NDims == NDimsRef) {
         LOG_INFO("Retrieve number of dimensions: PASS");
      } else {
         LOG_ERROR("Retrieve number of dimensions: FAIL");
         ++Err;
      }

      // Check to see if expected dimensions exist (and a non-existent one
      // doesn't)
      if (Dimension::exists("NCells") and Dimension::exists("NEdges") and
          Dimension::exists("NVertLevels") and !Dimension::exists("Garbage")) {
         LOG_INFO("Test dimension existence function: PASS");
      } else {
         LOG_ERROR("Test dimension existence function: FAIL");
      }

      // Test length retrieval by name

      I4 NCellsGlb = Dimension::getDimLengthGlobal("NCells");
      I4 NCellsLoc = Dimension::getDimLengthLocal("NCells");
      if (NCellsGlb == NCellsGlbRef and NCellsLoc == NCellsLocRef) {
         LOG_INFO("Length retrievals for NCells: PASS");
      } else {
         LOG_ERROR("Length retrievals for NCells: FAIL");
         ++Err;
      }

      I4 NEdgesGlb = Dimension::getDimLengthGlobal("NEdges");
      I4 NEdgesLoc = Dimension::getDimLengthLocal("NEdges");
      if (NEdgesGlb == NEdgesGlbRef and NEdgesLoc == NEdgesLocRef) {
         LOG_INFO("Length retrievals for NEdges: PASS");
      } else {
         LOG_ERROR("Length retrievals for NEdges: FAIL");
         ++Err;
      }

      I4 NVertLvlsGlb = Dimension::getDimLengthGlobal("NVertLevels");
      I4 NVertLvlsLoc = Dimension::getDimLengthLocal("NVertLevels");
      if (NVertLvlsGlb == NVertLvlsRef and NVertLvlsLoc == NVertLvlsRef) {
         LOG_INFO("Length retrievals for NVertLevels: PASS");
      } else {
         LOG_ERROR("Length retrievals for NVertLevels: FAIL");
         ++Err;
      }

      // Test distributed property by name
      if (Dimension::isDistributedDim("NCells") and
          Dimension::isDistributedDim("NEdges") and
          !Dimension::isDistributedDim("NVertLevels")) {
         LOG_INFO("Test distributed property by name: PASS");
      } else {
         LOG_ERROR("Test distributed property by name: FAIL");
         ++Err;
      }

      // Test get offset array by dim name
      HostArray1DI4 OffsetCell = Dimension::getDimOffset("NCells");
      HostArray1DI4 OffsetEdge = Dimension::getDimOffset("NEdges");
      HostArray1DI4 OffsetVert = Dimension::getDimOffset("NVertLevels");
      I4 Count                 = 0;
      for (int N = 0; N < NCellsLocRef; ++N) {
         if (N < NCellsOwned) {
            if (OffsetCell(N) != DefDecomp->CellIDH(N) - 1)
               ++Count;
         } else {
            if (OffsetCell(N) != -1)
               ++Count;
         }
      }
      for (int N = 0; N < NEdgesLocRef; ++N) {
         if (N < NEdgesOwned) {
            if (OffsetEdge(N) != DefDecomp->EdgeIDH(N) - 1)
               ++Count;
         } else {
            if (OffsetEdge(N) != -1)
               ++Count;
         }
      }
      for (int N = 0; N < NVertLvlsRef; ++N) {
         if (OffsetVert(N) != N)
            ++Count;
      }
      if (Count == 0) {
         LOG_INFO("Offset retrieval by name: PASS");
      } else {
         LOG_ERROR("Offset retrieval by name: FAIL");
         ++Err;
      }

      // Test iterators and also retrieval by instance
      for (auto Iter = Dimension::begin(); Iter != Dimension::end(); ++Iter) {
         std::string ThisName               = Iter->first;
         std::shared_ptr<Dimension> ThisDim = Iter->second;
         std::string MyName                 = ThisDim->getName();
         I4 LengthLoc                       = ThisDim->getLengthLocal();
         I4 LengthGlb                       = ThisDim->getLengthGlobal();
         bool Distrib                       = ThisDim->isDistributed();
         HostArray1DI4 OffsetTest           = ThisDim->getOffset();

         bool ScalarPass = false;
         bool OffsetPass = false;
         if (MyName == "NCells") {
            if (LengthLoc == NCellsLocRef and LengthGlb == NCellsGlbRef and
                Distrib)
               ScalarPass = true;
            Count = 0;
            for (int N = 0; N < NCellsLocRef; ++N) {
               if (N < NCellsOwned) {
                  if (OffsetTest(N) != DefDecomp->CellIDH(N) - 1)
                     ++Count;
               } else {
                  if (OffsetTest(N) != -1)
                     ++Count;
               }
            }
            if (Count == 0)
               OffsetPass = true;
         } else if (MyName == "NEdges") {
            if (LengthLoc == NEdgesLocRef and LengthGlb == NEdgesGlbRef and
                Distrib)
               ScalarPass = true;
            Count = 0;
            for (int N = 0; N < NEdgesLocRef; ++N) {
               if (N < NEdgesOwned) {
                  if (OffsetTest(N) != DefDecomp->EdgeIDH(N) - 1)
                     ++Count;
               } else {
                  if (OffsetTest(N) != -1)
                     ++Count;
               }
            }
            if (Count == 0)
               OffsetPass = true;
         } else if (MyName == "NVertLevels") {
            if (LengthLoc == NVertLvlsRef and LengthGlb == NVertLvlsRef and
                !Distrib)
               ScalarPass = true;
            Count = 0;
            for (int N = 0; N < NVertLvlsRef; ++N) {
               if (OffsetTest(N) != N)
                  ++Count;
            }
            if (Count == 0)
               OffsetPass = true;
         } else {
            LOG_ERROR("Unknown dimension name in iteration loop: FAIL");
            ++Err;
         }
         if (ScalarPass and OffsetPass) {
            LOG_INFO("Retrieval by instance in loop: PASS");
         } else {
            LOG_ERROR("Retrieval by instance in loop: FAIL");
            ++Err;
         }
      } // end iteration over dimensions

      // Check retrieval of full dimension by name - just check scalars since
      // other tests should have picked up offset errors
      std::shared_ptr<Dimension> TestDim = Dimension::get("NCells");
      NCellsGlb                          = TestDim->getLengthGlobal();
      NCellsLoc                          = TestDim->getLengthLocal();
      if (NCellsGlb == NCellsGlbRef and NCellsLoc == NCellsLocRef and
          TestDim->isDistributed()) {
         LOG_INFO("Retrieval of full dim: PASS");
      } else {
         LOG_ERROR("Retrieval of full dim: FAIL");
         ++Err;
      }

      // Destroy a dimension
      Dimension::destroy("NCells");
      if (Dimension::exists("NCells")) {
         LOG_ERROR("Dimension destroy test: FAIL");
         ++Err;
      } else {
         LOG_INFO("Dimension destroy test: PASS");
      }

      // Test removal of all dims
      Dimension::clear();
      NDims = Dimension::getNumDefinedDims();
      if (NDims == 0) {
         LOG_INFO("Dimension clear test: PASS");
      } else {
         LOG_ERROR("Dimension clear test: FAIL");
         ++Err;
      }
   }

   // Clean up environments
   Decomp::clear();
   Kokkos::finalize();
   MPI_Finalize();

   if (Err >= 256)
      Err = 255;

   // End of testing
   return Err;
}
//===--- End test driver for Dimension
//--------------------------------------===/
