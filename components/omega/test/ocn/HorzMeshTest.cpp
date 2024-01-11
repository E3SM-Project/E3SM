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
#include "HorzMesh.h"
#include "DataTypes.h"
#include "IO.h"
#include "Logging.h"
#include "MachEnv.h"
#include "mpi.h"


#include <iostream>

//------------------------------------------------------------------------------
// The initialization routine for Decomp testing. It calls various
// init routines, including the creation of the default decomposition.

int initHorzMeshTest() {

   int Err = 0;

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv. Then retrieve the default environment and
   // some needed data members.
   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   MPI_Comm DefComm       = DefEnv->getComm();

   // Initialize the IO system
   Err = OMEGA::IOInit(DefComm);
   if (Err != 0)
      LOG_ERROR("HorzMeshTest: error initializing parallel IO");

   // Create the default decomposition (initializes the decomposition)
   Err = OMEGA::Decomp::init();
   if (Err != 0)
      LOG_ERROR("HorzMeshTest: error initializing default decomposition");

   return Err;
}

OMEGA::R8 distance(OMEGA::R8 x, OMEGA::R8 y, OMEGA::R8 z) {
   
   OMEGA::R8 dist;

   dist = sqrt(x*x + y*y + z*z);

   return dist;

}

//------------------------------------------------------------------------------
// The test driver for Decomp. This tests the decomposition of a sample
// horizontal domain and verifies the mesh is decomposed correctly.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   yakl::init();

   // Call initialization routine to create the default decomposition
   int Err = initHorzMeshTest();
   if (Err != 0)
      LOG_CRITICAL("HorzMeshTest: Error initializing");

   // Get MPI vars if needed
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();
   MPI_Comm Comm          = DefEnv->getComm();
   OMEGA::I4 MyTask       = DefEnv->getMyTask();
   OMEGA::I4 NumTasks     = DefEnv->getNumTasks();
   bool IsMaster          = DefEnv->isMasterTask();

   // Test retrieval of the default decomposition
   OMEGA::Decomp *DefDecomp = OMEGA::Decomp::getDefault();
   if (DefDecomp) { // true if non-null ptr
      LOG_INFO("HorzMeshTest: Default decomp retrieval PASS");
   } else {
      LOG_INFO("HorzMeshTest: Default decomp retrieval FAIL");
      return -1;
   }

   OMEGA::HorzMesh Mesh(DefDecomp);
   OMEGA::I4 SumCells;
   OMEGA::I4 LocCells;
   LocCells = Mesh.NCellsOwned;
   Err = MPI_Allreduce(&LocCells, &SumCells, 1, MPI_INT32_T, MPI_SUM, Comm);


   if (SumCells == DefDecomp->NCellsGlobal){
     LOG_INFO("HorzMeshTest: Sum cell ID test PASS");
   } else {
      LOG_INFO("HorzMeshTest: Sum cell ID test FAIL {} {}", SumCells,
              DefDecomp->NCellsGlobal);
   }

   OMEGA::R8 sphere_radius = distance(Mesh.XCellH(0), Mesh.YCellH(0), Mesh.ZCellH(0));
   OMEGA::R8 tol = 1e-6;
   OMEGA::R8 dist;
   OMEGA::I4 count = 0;
   for (int Cell = 0; Cell < LocCells; Cell++) {
       dist = distance(Mesh.XCellH(Cell), Mesh.YCellH(Cell), Mesh.ZCellH(Cell));
       if ( abs(sphere_radius - dist) > tol)
          count++;
   }

   if (count > 0) {
     LOG_INFO("HorzMeshTest: Cell sphere radius test FAIL");
   } else {
     LOG_INFO("HorzMeshTest: Cell sphere radius test PASS");
   }


   // Test that device arrays are identical

   // MPI_Status status;
   if (Err == 0)
      LOG_INFO("HorzMeshTest: Successful completion");
   yakl::finalize();
   MPI_Finalize();

} // end of main
//===-----------------------------------------------------------------------===/
