//===-- Test driver for OMEGA MachEnv ----------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA MachEnv
///
/// This driver tests the OMEGA model MachEnv module that sets up various
/// machine parameters, including message passing (MPI) variables, threading
/// and other potential variables related to the underlying machine
/// architecture. This unit test driver primarily tests that the quantities
/// and retrieval functions are working correctly.
///
//
//===-----------------------------------------------------------------------===/

#include "MachEnv.h"
#include "Error.h"
#include "Logging.h"
#include "mpi.h"

using namespace OMEGA;

//------------------------------------------------------------------------------
//
/// \brief Initialization for OMEGA MachEnv tests
///
/// This initialization routine initializes several environments in a setting
/// other than the test driver to make sure the environments are persistent
/// across subroutine calls.

void initMachEnvs() {

   // Initialize several environments in reverse order that they
   // are tested.  Use the default environment as the parent
   MachEnv *DefEnv = MachEnv::getDefault();

   // Initialize the Logging system
   initLogging(DefEnv);
   LOG_INFO("------ MachEnv Unit Tests ------");

   // Initialize general subset environment
   int InclSize     = 4;
   int InclTasks[4] = {1, 2, 5, 7};
   MachEnv::create("Subset", DefEnv, InclSize, InclTasks);

   // Initialize strided environment
   MachEnv::create("Stride", DefEnv, 4, 1, 2);

   // Initialize contiguous subset environment
   MachEnv::create("Contig", DefEnv, 4);

   // Initialize contiguous subset environment but different master task
   MachEnv::create("Contig2", DefEnv, 4, 2);

} // end of InitMachEnvs

//------------------------------------------------------------------------------
// The test driver for MachEnv. This tests the values stored in the Default
// Environment and three other based on the three subsetting options.  All
// current values and get routines are tested.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);

   // Create reference values based on MPI_COMM_WORLD
   int WorldTask;
   int WorldSize;
   MPI_Comm_rank(MPI_COMM_WORLD, &WorldTask);
   MPI_Comm_size(MPI_COMM_WORLD, &WorldSize);
   int WorldMaster = 0;

   // The subset environments create 4-task sub-environments so
   // make sure the unit test is run with at least 8 to adequately
   // test all subsets.
   if (WorldSize < 8)
      ABORT_ERROR("MachEnvTest: FAIL must run unit test with at least 8 tasks");

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv
   MachEnv::init(MPI_COMM_WORLD);

   // Initialize the test environments. We do this in a separate routine
   // to make sure environments created during the OMEGA init phase
   // persist to the later run stage.
   initMachEnvs();

   // Verify retrieved values of the Default environment match the
   // expected reference values
   MachEnv *DefEnv = MachEnv::getDefault();

   int MyTask = DefEnv->getMyTask();
   if (MyTask != WorldTask)
      ABORT_ERROR("MachEnvTest: DefaultEnv task test: FAIL "
                  "MyTask = {}, WorldTask = {}",
                  MyTask, WorldTask);

   int MySize = DefEnv->getNumTasks();
   if (MySize != WorldSize)
      ABORT_ERROR("MachEnvTest: DefaultEnv NumTasks test: FAIL "
                  "MySize = {}, WorldSize = {}",
                  MySize, WorldSize);

   int MyMaster = DefEnv->getMasterTask();
   if (MyMaster != WorldMaster)
      ABORT_ERROR("MachEnvTest: DefaultEnv master task test: FAIL "
                  "MyMaster = {}, WorldMaster = {}",
                  MyMaster, WorldMaster);

   bool IsMyMaster = DefEnv->isMasterTask();
   if (MyTask == MyMaster) {
      if (!IsMyMaster)
         ABORT_ERROR("MachEnvTest: DefaultEnv is master task test: FAIL");
   } else {
      if (IsMyMaster)
         ABORT_ERROR("MachEnvTest: DefaultEnv is master task test: PASS");
   }

   // Test setting a new master task

   DefEnv->setMasterTask(2);

   MyMaster = DefEnv->getMasterTask();
   if (MyMaster != 2)
      ABORT_ERROR("MachEnvTest: DefaultEnv set master task test: FAIL "
                  "MyMaster = {}",
                  MyMaster);

   IsMyMaster = DefEnv->isMasterTask();
   if (MyTask == 2) {
      if (!IsMyMaster)
         ABORT_ERROR("MachEnvTest: DefaultEnv isMaster after setMaster: FAIL");
   } else {
      if (IsMyMaster)
         ABORT_ERROR("MachEnvTest: DefaultEnv isMaster after setMaster: FAIL");
   }

   //---------------------------------------------------------------------------
   // Test contiguous subset environment (first four tasks of default)

   // Test retrieval of the contiguous environment
   MachEnv *ContigEnv = MachEnv::get("Contig");

   // Check membership in the new communicator
   if (MyTask < 4) {
      if (!(ContigEnv->isMember()))
         ABORT_ERROR("MachEnvTest: contiguous member test: FAIL MyTask {}",
                     MyTask);
   } else {
      if (ContigEnv->isMember())
         ABORT_ERROR("MachEnvTest: contiguous non-member test: FAIL MyTask {}",
                     MyTask);
   }

   // Perform standard checks on new communicator
   if (ContigEnv->isMember()) {
      int ContigTask = ContigEnv->getMyTask();
      if (ContigTask != WorldTask)
         ABORT_ERROR("MachEnvTest: contiguous task test: FAIL "
                     "ContigTask = {}, WorldTask = {}",
                     ContigTask, WorldTask);

      int ContigSize = ContigEnv->getNumTasks();
      if (ContigSize != 4)
         ABORT_ERROR("MachEnvTest: contiguous NumTasks test: FAIL "
                     "ContigSize = {} (should be 4)",
                     ContigSize);

      int ContigMaster = ContigEnv->getMasterTask();
      if (ContigMaster != WorldMaster)
         ABORT_ERROR("MachEnvTest: contiguous master task test: FAIL "
                     "MyMaster = {}, WorldMaster = {}",
                     MyMaster, WorldMaster);

      bool IsContigMaster = ContigEnv->isMasterTask();
      if (ContigTask == ContigMaster) {
         if (!IsContigMaster)
            ABORT_ERROR("MachEnvTest: contiguous is master task test: FAIL");
      } else {
         if (IsContigMaster)
            ABORT_ERROR("MachEnvTest: contiguous is master task test: FAIL");
      }

   } // end if member of Contiguous subset

   //---------------------------------------------------------------------------
   // Test a similar contiguous subset environment (first four tasks of default)
   // in which the master task has been initialized to task 2

   // Test retrieval of the contiguous environment
   MachEnv *Contig2Env = MachEnv::get("Contig2");

   // Check membership in the new communicator
   if (MyTask < 4) {
      if (!(Contig2Env->isMember()))
         ABORT_ERROR("MachEnvTest: contiguous2 member test: FAIL MyTask {}",
                     MyTask);
   } else {
      if (Contig2Env->isMember())
         ABORT_ERROR("MachEnvTest: contiguous2 non-member test: FAIL MyTask {}",
                     MyTask);
   }

   // Perform standard checks on new communicator
   if (Contig2Env->isMember()) {
      int Contig2Task = Contig2Env->getMyTask();
      if (Contig2Task != WorldTask)
         ABORT_ERROR("MachEnvTest: contiguous2 task test: FAIL "
                     "Contig2Task = {}, WorldTask = {}",
                     Contig2Task, WorldTask);

      int Contig2Size = Contig2Env->getNumTasks();
      if (Contig2Size != 4)
         ABORT_ERROR("MachEnvTest: contiguous2 NumTasks test: FAIL "
                     "Contig2Size = {} (should be 4)",
                     Contig2Size);

      int Contig2Master = Contig2Env->getMasterTask();
      if (Contig2Master != 2)
         ABORT_ERROR("MachEnvTest: contiguous2 master task test: FAIL "
                     "MyMaster = {}, WorldMaster = {}",
                     MyMaster, WorldMaster);

      bool IsContig2Master = Contig2Env->isMasterTask();
      if (Contig2Task == Contig2Master) {
         if (!IsContig2Master)
            ABORT_ERROR("MachEnvTest: contiguous2 is master task test: FAIL");
      } else {
         if (IsContig2Master)
            ABORT_ERROR("MachEnvTest: contiguous2 is master task test: FAIL");
      }

   } // end if member of Contig2 subset

   //---------------------------------------------------------------------------
   // Test the strided constructor with only odd-numbered tasks

   // Test retrieval
   MachEnv *StrideEnv = MachEnv::get("Stride");

   // Check membership in the new communicator
   if (MyTask % 2 == 1) {
      if (!(StrideEnv->isMember()))
         ABORT_ERROR("MachEnvTest: strided member test: FAIL MyTask {}",
                     MyTask);
   } else {
      if (StrideEnv->isMember())
         ABORT_ERROR("MachEnvTest: strided non-member test: FAIL MyTask {}",
                     MyTask);
   }

   // Perform standard checks on new communicator
   if (StrideEnv->isMember()) {
      int StrideTask = StrideEnv->getMyTask();
      if (StrideTask != WorldTask / 2)
         ABORT_ERROR("MachEnvTest: strided task test: FAIL "
                     "StrideTask = {}, WorldTask = {}",
                     StrideTask, WorldTask);

      int StrideSize = StrideEnv->getNumTasks();
      if (StrideSize != 4)
         ABORT_ERROR("MachEnvTest: strided NumTasks test: FAIL "
                     "StrideSize = {} (should be 4)",
                     StrideSize);

      int StrideMaster = StrideEnv->getMasterTask();
      if (StrideMaster != 0)
         ABORT_ERROR("MachEnvTest: strided master task test: FAIL master = {}",
                     StrideMaster);

      bool IsStrideMaster = StrideEnv->isMasterTask();
      if (StrideTask == StrideMaster) {
         if (!IsStrideMaster)
            ABORT_ERROR("MachEnvTest: strided is master task test: FAIL");
      } else {
         if (IsStrideMaster)
            ABORT_ERROR("MachEnvTest: strided is master task test: FAIL");
      }

   } // end if member of strided subset

   //---------------------------------------------------------------------------
   // Test general subset constructor using tasks 1,2,5,7

   int InclSize     = 4;
   int InclTasks[4] = {1, 2, 5, 7};

   // Test retrieval
   MachEnv *SubsetEnv = MachEnv::get("Subset");

   // Check membership in the new communicator
   MyTask          = SubsetEnv->getMyTask();
   bool TaskInList = false;
   int NewTask     = -1;
   for (int i = 0; i < InclSize; ++i) {
      ++NewTask;
      if (WorldTask == InclTasks[i]) {
         TaskInList = true;
         break;
      }
   }

   if (TaskInList) {
      if (!(SubsetEnv->isMember()))
         ABORT_ERROR("MachEnvTest: subset member test: FAIL Task {}", MyTask);

   } else {
      if (SubsetEnv->isMember())
         ABORT_ERROR("MachEnvTest: subset non-member test: FAIL Task {}",
                     MyTask);
   }

   // Perform standard checks on new communicator
   if (SubsetEnv->isMember()) {
      int SubsetTask = SubsetEnv->getMyTask();
      if (SubsetTask != NewTask)
         ABORT_ERROR("MachEnvTest: subset task test: FAIL "
                     "SubsetTask = {}, NewTask = {}",
                     SubsetTask, NewTask);

      int SubsetSize = SubsetEnv->getNumTasks();
      if (SubsetSize != InclSize)
         ABORT_ERROR("MachEnvTest: subset size test: FAIL "
                     "SubsetSize = {}, InclSize  = {}",
                     SubsetSize, InclSize);

      int SubsetMaster = SubsetEnv->getMasterTask();
      if (SubsetMaster != 0)
         ABORT_ERROR("MachEnvTest: subset master task test: FAIL master = {}",
                     SubsetMaster);

      bool IsSubsetMaster = SubsetEnv->isMasterTask();
      if (SubsetTask == SubsetMaster) {
         if (!IsSubsetMaster)
            ABORT_ERROR("MachEnvTest: subset is master task test: FAIL");
      } else {
         if (IsSubsetMaster)
            ABORT_ERROR("MachEnvTest: subset is master task test: FAIL");
      }

   } // end if member of general subset env

   //---------------------------------------------------------------------------
   // Test setting of compile-time vector length

#ifdef OMEGA_VECTOR_LENGTH
   if (VecLength != OMEGA_VECTOR_LENGTH)
      ABORT_ERROR("MachEnvTest: MPI vector length test: FAIL");
#endif

   // finalize and clean up environments (test both removal functions)
   MachEnv::removeEnv("Contig");
   MachEnv::removeAll();

   LOG_INFO("------ MachEnv Unit Tests Successful ------");

   // MPI_Status status;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0; // if we made it here, we were successful

} // end of main
//===-----------------------------------------------------------------------===/
