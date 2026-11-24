//===-- Test driver for OMEGA Dimension class --------------------*- C++ -*-===/
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
#include "Error.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "OmegaKokkos.h"
#include "Pacer.h"
#include "mpi.h"

#include <memory>
#include <vector>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Initialization routine to create dimensions
void initDimensionTest() {

   // Initialize various environments
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = MachEnv::getDefault();
   initLogging(DefEnv);
   MPI_Comm DefComm = DefEnv->getComm();
   Pacer::initialize(DefComm);
   Pacer::setPrefix("Omega:");
   LOG_INFO("------ Dimension unit tests ------");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Initialize parallel IO
   IO::init(DefComm);

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
   I4 NVertLayers = 100;
   std::shared_ptr<Dimension> VertDim =
       Dimension::create("NVertLayers", NVertLayers);

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
      // Call initialization to create sample dimensions
      initDimensionTest();

      // Retrieve the default domain decomposition and set reference values
      Decomp *DefDecomp = Decomp::getDefault();

      I4 NCellsLocRef = DefDecomp->NCellsSize;
      I4 NCellsGlbRef = DefDecomp->NCellsGlobal;
      I4 NCellsOwned  = DefDecomp->NCellsOwned;
      I4 NEdgesLocRef = DefDecomp->NEdgesSize;
      I4 NEdgesGlbRef = DefDecomp->NEdgesGlobal;
      I4 NEdgesOwned  = DefDecomp->NEdgesOwned;
      I4 NVertLyrsRef = 100;

      // Retrieve the number of defined dimensions
      I4 NDimsRef = 3;
      I4 NDims    = Dimension::getNumDefinedDims();
      if (NDims != NDimsRef)
         ABORT_ERROR("DimensionTest: Retrieve number of dimensions - FAIL");

      // Check to see if expected dimensions exist (and a non-existent one
      // doesn't)
      if (!Dimension::exists("NCells") or !Dimension::exists("NEdges") or
          !Dimension::exists("NVertLayers") or Dimension::exists("Garbage"))
         ABORT_ERROR("DimensionTest: dimension existence function - FAIL");

      // Test length retrieval by name

      I4 NCellsGlb = Dimension::getDimLengthGlobal("NCells");
      I4 NCellsLoc = Dimension::getDimLengthLocal("NCells");
      if (NCellsGlb != NCellsGlbRef or NCellsLoc != NCellsLocRef)
         ABORT_ERROR("DimensionTest: Length retrieval NCells FAIL {} {} {} {}",
                     NCellsGlb, NCellsGlbRef, NCellsLoc, NCellsLocRef);

      I4 NEdgesGlb = Dimension::getDimLengthGlobal("NEdges");
      I4 NEdgesLoc = Dimension::getDimLengthLocal("NEdges");
      if (NEdgesGlb != NEdgesGlbRef or NEdgesLoc != NEdgesLocRef)
         ABORT_ERROR("DimensionTest: Length retrieval NEdges FAIL {} {} {} {}",
                     NEdgesGlb, NEdgesGlbRef, NEdgesLoc, NEdgesLocRef);

      I4 NVertLyrsGlb = Dimension::getDimLengthGlobal("NVertLayers");
      I4 NVertLyrsLoc = Dimension::getDimLengthLocal("NVertLayers");
      if (NVertLyrsGlb != NVertLyrsRef or NVertLyrsLoc != NVertLyrsRef)
         ABORT_ERROR(
             "DimensionTest: Length retrieval NVertLayers FAIL {} {} {} {}",
             NVertLyrsGlb, NVertLyrsRef, NVertLyrsLoc, NVertLyrsRef)

      // Test distributed property by name
      if (!Dimension::isDistributedDim("NCells") or
          !Dimension::isDistributedDim("NEdges") or
          Dimension::isDistributedDim("NVertLayers"))
         ABORT_ERROR("DimensionTest: distributed property by name - FAIL");

      // Test get offset array by dim name
      HostArray1DI4 OffsetCell = Dimension::getDimOffset("NCells");
      HostArray1DI4 OffsetEdge = Dimension::getDimOffset("NEdges");
      HostArray1DI4 OffsetVert = Dimension::getDimOffset("NVertLayers");
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
      for (int N = 0; N < NVertLyrsRef; ++N) {
         if (OffsetVert(N) != N)
            ++Count;
      }
      if (Count > 0)
         ABORT_ERROR("DimensionTest: Offset retrieval by name - FAIL");

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
         } else if (MyName == "NVertLayers") {
            if (LengthLoc == NVertLyrsRef and LengthGlb == NVertLyrsRef and
                !Distrib)
               ScalarPass = true;
            Count = 0;
            for (int N = 0; N < NVertLyrsRef; ++N) {
               if (OffsetTest(N) != N)
                  ++Count;
            }
            if (Count == 0)
               OffsetPass = true;
         } else {
            ABORT_ERROR(
                "DimensionTest: Unknown dimension {} in iteration loop: FAIL",
                MyName);
         }
         if (!ScalarPass or !OffsetPass)
            ABORT_ERROR("DimensionTest: Retrieval by instance in loop: FAIL");
      } // end iteration over dimensions

      // Check retrieval of full dimension by name - just check scalars since
      // other tests should have picked up offset errors
      std::shared_ptr<Dimension> TestDim = Dimension::get("NCells");
      NCellsGlb                          = TestDim->getLengthGlobal();
      NCellsLoc                          = TestDim->getLengthLocal();
      if (NCellsGlb != NCellsGlbRef or NCellsLoc != NCellsLocRef or
          !(TestDim->isDistributed()))
         ABORT_ERROR("DimensionTest: Retrieval of full dim - FAIL");

      // Destroy a dimension
      Dimension::destroy("NCells");
      if (Dimension::exists("NCells"))
         ABORT_ERROR("DimensionTest: destroy test - FAIL");

      // Test removal of all dims
      Dimension::clear();
      NDims = Dimension::getNumDefinedDims();
      if (NDims > 0)
         ABORT_ERROR("DimensionTest: remove all test - FAIL");
   }

   LOG_INFO("------ Dimension unit tests successful ------");
   // Clean up environments
   Decomp::clear();
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   // End of testing
   return Err;
}
//===--- End test driver for Dimension -------------------------------------===/
