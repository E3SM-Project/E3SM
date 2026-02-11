//===-- Test driver for OMEGA Decomp -----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA domain decomposition (Decomp)
///
/// This driver tests the OMEGA domain decomposition, decomposing the
/// horizontal domain and creating a number of index-space arrays for
/// locating and describing mesh locations within a parallel distributed
/// memory.
///
//
//===-----------------------------------------------------------------------===/

#include "Decomp.h"
#include "Config.h"
#include "DataTypes.h"
#include "Error.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "Pacer.h"
#include "Reductions.h"
#include "mpi.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// The initialization routine for Decomp testing. It calls various
// init routines, including the creation of the default decomposition.

void initDecompTest() {

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   initLogging(DefEnv);
   LOG_INFO("---- Starting Decomp unit tests ----");

   // Open config file
   Config("Omega");
   Config::readAll("omega.yml");

   // Initialize the IO system
   IO::init(DefComm);

   // Create the default decomposition (initializes the decomposition)
   Decomp::init();

   return;
}

//------------------------------------------------------------------------------
// The test driver for Decomp. This tests the decomposition of a sample
// horizontal domain and verifies the mesh is decomposed correctly.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   Kokkos::initialize();
   Pacer::initialize(MPI_COMM_WORLD);
   Pacer::setPrefix("Omega:");
   {
      // Call initialization routine to create the default decomposition
      initDecompTest();

      // Get MPI vars if needed
      MachEnv *DefEnv = MachEnv::getDefault();
      MPI_Comm Comm   = DefEnv->getComm();
      I4 MyTask       = DefEnv->getMyTask();
      I4 NumTasks     = DefEnv->getNumTasks();
      bool IsMaster   = DefEnv->isMasterTask();

      // Test retrieval of the default decomposition
      Decomp *DefDecomp = Decomp::getDefault();

      // Extract some index info from Decomp
      I4 NCellsGlobal    = DefDecomp->NCellsGlobal;
      I4 NCellsOwned     = DefDecomp->NCellsOwned;
      I4 NEdgesGlobal    = DefDecomp->NEdgesGlobal;
      I4 NEdgesOwned     = DefDecomp->NEdgesOwned;
      I4 NVerticesGlobal = DefDecomp->NVerticesGlobal;
      I4 NVerticesOwned  = DefDecomp->NVerticesOwned;

      // Test that all Cells, Edges, Vertices are accounted for by
      // summing the list of owned values by all tasks. The result should
      // be the sum of the integers from 1 to NCellsGlobal (or edges, vertices).
      I4 RefSumCells    = 0;
      I4 RefSumEdges    = 0;
      I4 RefSumVertices = 0;
      for (int n = 0; n < NCellsGlobal; ++n)
         RefSumCells += n + 1;
      for (int n = 0; n < NEdgesGlobal; ++n)
         RefSumEdges += n + 1;
      for (int n = 0; n < NVerticesGlobal; ++n)
         RefSumVertices += n + 1;
      I4 LocSumCells          = 0;
      I4 LocSumEdges          = 0;
      I4 LocSumVertices       = 0;
      HostArray1DI4 CellIDH   = DefDecomp->CellIDH;
      HostArray1DI4 EdgeIDH   = DefDecomp->EdgeIDH;
      HostArray1DI4 VertexIDH = DefDecomp->VertexIDH;
      for (int n = 0; n < NCellsOwned; ++n)
         LocSumCells += CellIDH(n);
      for (int n = 0; n < NEdgesOwned; ++n)
         LocSumEdges += EdgeIDH(n);
      for (int n = 0; n < NVerticesOwned; ++n)
         LocSumVertices += VertexIDH(n);
      I4 SumCells    = globalSum(LocSumCells, Comm);
      I4 SumEdges    = globalSum(LocSumEdges, Comm);
      I4 SumVertices = globalSum(LocSumVertices, Comm);

      if (SumCells != RefSumCells)
         ABORT_ERROR("DecompTest: Sum cell ID test FAIL {} {}", SumCells,
                     RefSumCells);
      if (SumEdges != RefSumEdges)
         ABORT_ERROR("DecompTest: Sum edge ID test FAIL {} {}", SumEdges,
                     RefSumEdges);
      if (SumVertices != RefSumVertices)
         ABORT_ERROR("DecompTest: Sum vertex ID test FAIL {} {}", SumVertices,
                     RefSumVertices);

      // Clean up
      Decomp::clear();
      MachEnv::removeAll();

      LOG_INFO("---- DecompTest: Successful completion ----");
   }
   Pacer::finalize();
   Kokkos::finalize();
   MPI_Finalize();

   // if we made it to the end, return success
   return 0;

} // end of main
//===-----------------------------------------------------------------------===/
