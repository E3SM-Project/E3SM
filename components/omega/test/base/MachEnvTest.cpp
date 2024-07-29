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
#include "Logging.h"
#include "mpi.h"

#include <iostream>

//------------------------------------------------------------------------------
//
/// \brief Initialization for OMEGA MachEnv tests
///
/// This initialization routine initializes several environments in a setting
/// other than the test driver to make sure the environments are persistent
/// across subroutine calls.

void InitMachEnvs() {

   // Initialize several environments in reverse order that they
   // are tested.  Use the default environment as the parent
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();

   // Initialize the Logging system
   OMEGA::initLogging(DefEnv);

   // Initialize general subset environment
   int InclSize     = 4;
   int InclTasks[4] = {1, 2, 5, 7};
   OMEGA::MachEnv::create("Subset", DefEnv, InclSize, InclTasks);

   // Initialize strided environment
   OMEGA::MachEnv::create("Stride", DefEnv, 4, 1, 2);

   // Initialize contiguous subset environment
   OMEGA::MachEnv::create("Contig", DefEnv, 4);

   // Initialize contiguous subset environment but different master task
   OMEGA::MachEnv::create("Contig2", DefEnv, 4, 2);

} // end of InitMachEnvs

//------------------------------------------------------------------------------
// The test driver for MachEnv. This tests the values stored in the Default
// Environment and three other based on the three subsetting options.  All
// current values and get routines are tested.
//
int main(int argc, char *argv[]) {

   int RetVal = 0;

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);

   // Create reference values based on MPI_COMM_WORLD
   int WorldTask;
   int WorldSize;
   MPI_Comm_rank(MPI_COMM_WORLD, &WorldTask);
   MPI_Comm_size(MPI_COMM_WORLD, &WorldSize);
   int WorldMaster = 0;
   bool IsWorldMaster;
   if (WorldTask == WorldMaster) {
      IsWorldMaster = true;
   } else {
      IsWorldMaster = false;
   }

   // The subset environments create 4-task sub-environments so
   // make sure the unit test is run with at least 8 to adequately
   // test all subsets.
   if (WorldSize < 8) {
      std::cerr << "Please run unit test with at least 8 tasks" << std::endl;
      std::cout << "MachEnv unit test: FAIL" << std::endl;
      return -1;
   }

   // Initialize the Machine Environment class - this also creates
   // the default MachEnv
   OMEGA::MachEnv::init(MPI_COMM_WORLD);

   // Initialize the test environments. We do this in a separate routine
   // to make sure environments created during the OMEGA init phase
   // persist to the later run stage.
   InitMachEnvs();

   // Verify retrieved values of the Default environment match the
   // expected reference values
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();

   int MyTask = DefEnv->getMyTask();
   if (MyTask == WorldTask)
      std::cout << "DefaultEnv task test: PASS" << std::endl;
   else {
      RetVal += 1;
      std::cout << "DefaultEnv task test: FAIL "
                << "MyTask, WorldTask = " << MyTask << WorldTask << std::endl;
   }

   int MySize = DefEnv->getNumTasks();
   if (MySize == WorldSize)
      std::cout << "DefaultEnv NumTasks test: PASS" << std::endl;
   else {
      RetVal += 1;
      std::cout << "DefaultEnv NumTasks test: FAIL "
                << "MySize, WorldSize = " << MySize << " " << WorldSize
                << std::endl;
   }

   int MyMaster = DefEnv->getMasterTask();
   if (MyMaster == WorldMaster)
      std::cout << "DefaultEnv master task test: PASS" << std::endl;
   else {
      RetVal += 1;
      std::cout << "DefaultEnv master task test: FAIL "
                << "MyMaster, WorldMaster = " << MyMaster << " " << WorldMaster
                << std::endl;
   }

   bool IsMyMaster = DefEnv->isMasterTask();
   if (MyTask == MyMaster) {
      if (IsMyMaster)
         std::cout << "DefaultEnv is master task test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "DefaultEnv is master task test: FAIL" << std::endl;
      }
   } else {
      if (IsMyMaster) {
         RetVal += 1;
         std::cout << "DefaultEnv is master task test: FAIL" << std::endl;
      } else
         std::cout << "DefaultEnv is master task test: PASS" << std::endl;
   }

   // Test setting a new master task

   DefEnv->setMasterTask(2);

   MyMaster = DefEnv->getMasterTask();
   if (MyMaster == 2)
      std::cout << "DefaultEnv set master task test: PASS" << std::endl;
   else {
      RetVal += 1;
      std::cout << "DefaultEnv set master task test: FAIL "
                << "MyMaster = " << MyMaster << std::endl;
   }

   IsMyMaster = DefEnv->isMasterTask();
   if (MyTask == 2) {
      if (IsMyMaster)
         std::cout << "DefaultEnv isMaster after setMaster: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "DefaultEnv isMaster after setMaster: FAIL" << std::endl;
      }
   } else {
      if (IsMyMaster) {
         RetVal += 1;
         std::cout << "DefaultEnv isMaster after setMaster: FAIL" << std::endl;
      } else
         std::cout << "DefaultEnv isMaster after setMaster: PASS" << std::endl;
   }

   //---------------------------------------------------------------------------
   // Test contiguous subset environment (first four tasks of default)

   // Test retrieval of the contiguous environment
   OMEGA::MachEnv *ContigEnv = OMEGA::MachEnv::getEnv("Contig");

   // Check membership in the new communicator
   if (MyTask < 4) {
      if (ContigEnv->isMember())
         std::cout << "contiguous member test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "contiguous member test: FAIL " << "MyTask " << MyTask
                   << std::endl;
      }
   } else {
      if (ContigEnv->isMember()) {
         RetVal += 1;
         std::cout << "contiguous member test: FAIL " << "MyTask " << MyTask
                   << std::endl;
      } else
         std::cout << "contiguous member test: PASS" << std::endl;
   }

   // Perform standard checks on new communicator
   if (ContigEnv->isMember()) {
      int ContigTask = ContigEnv->getMyTask();
      if (ContigTask == WorldTask)
         std::cout << "contiguous task test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "contiguous task test: FAIL "
                   << "ContigTask, WorldTask = " << ContigTask << " "
                   << WorldTask << std::endl;
      }

      int ContigSize = ContigEnv->getNumTasks();
      if (ContigSize == 4)
         std::cout << "contiguous NumTasks test: PASS " << std::endl;
      else {
         RetVal += 1;
         std::cout << "contiguous NumTasks test: FAIL "
                   << "ContigSize (should be 4)  = " << ContigSize << std::endl;
      }

      int ContigMaster = ContigEnv->getMasterTask();
      if (ContigMaster == WorldMaster)
         std::cout << "contiguous master task test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "contiguous master task test: FAIL "
                   << "MyMaster, WorldMaster = " << MyMaster << " "
                   << WorldMaster << std::endl;
      }

      bool IsContigMaster = ContigEnv->isMasterTask();
      if (ContigTask == ContigMaster) {
         if (IsContigMaster)
            std::cout << "contiguous is master task test: PASS" << std::endl;
         else {
            RetVal += 1;
            std::cout << "contiguous is master task test: FAIL" << std::endl;
         }
      } else {
         if (IsContigMaster) {
            RetVal += 1;
            std::cout << "contiguous is master task test: FAIL" << std::endl;
         } else
            std::cout << "contiguous is master task test: PASS" << std::endl;
      }

   } // end if member of Contiguous subset

   //---------------------------------------------------------------------------
   // Test a similar contiguous subset environment (first four tasks of default)
   // in which the master task has been initialized to task 2

   // Test retrieval of the contiguous environment
   OMEGA::MachEnv *Contig2Env = OMEGA::MachEnv::getEnv("Contig2");

   // Check membership in the new communicator
   if (MyTask < 4) {
      if (Contig2Env->isMember())
         std::cout << "contiguous2 member test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "contiguous2 member test: FAIL " << "MyTask " << MyTask
                   << std::endl;
      }
   } else {
      if (Contig2Env->isMember()) {
         RetVal += 1;
         std::cout << "contiguous2 member test: FAIL " << "MyTask " << MyTask
                   << std::endl;
      } else
         std::cout << "contiguous2 member test: PASS" << std::endl;
   }

   // Perform standard checks on new communicator
   if (Contig2Env->isMember()) {
      int Contig2Task = Contig2Env->getMyTask();
      if (Contig2Task == WorldTask)
         std::cout << "contiguous2 task test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "contiguous2 task test: FAIL "
                   << "Contig2Task, WorldTask = " << Contig2Task << " "
                   << WorldTask << std::endl;
      }

      int Contig2Size = Contig2Env->getNumTasks();
      if (Contig2Size == 4)
         std::cout << "contiguous2 NumTasks test: PASS " << std::endl;
      else {
         RetVal += 1;
         std::cout << "contiguous2 NumTasks test: FAIL "
                   << "Contig2Size (should be 4)  = " << Contig2Size
                   << std::endl;
      }

      int Contig2Master = Contig2Env->getMasterTask();
      if (Contig2Master == 2)
         std::cout << "contiguous2 master task test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "contiguous2 master task test: FAIL "
                   << "MyMaster, WorldMaster = " << MyMaster << " "
                   << WorldMaster << std::endl;
      }

      bool IsContig2Master = Contig2Env->isMasterTask();
      if (Contig2Task == Contig2Master) {
         if (IsContig2Master)
            std::cout << "contiguous2 is master task test: PASS" << std::endl;
         else {
            RetVal += 1;
            std::cout << "contiguous2 is master task test: FAIL" << std::endl;
         }
      } else {
         if (IsContig2Master) {
            RetVal += 1;
            std::cout << "contiguous2 is master task test: FAIL" << std::endl;
         } else
            std::cout << "contiguous2 is master task test: PASS" << std::endl;
      }

   } // end if member of Contig2 subset

   //---------------------------------------------------------------------------
   // Test the strided constructor with only odd-numbered tasks

   // Test retrieval
   OMEGA::MachEnv *StrideEnv = OMEGA::MachEnv::getEnv("Stride");

   // Check membership in the new communicator
   if (MyTask % 2 == 1) {
      if (StrideEnv->isMember())
         std::cout << "strided member test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "strided member test: FAIL " << "MyTask " << MyTask
                   << std::endl;
      }
   } else {
      if (StrideEnv->isMember()) {
         RetVal += 1;
         std::cout << "strided member test: FAIL " << "MyTask " << MyTask
                   << std::endl;
      } else
         std::cout << "strided member test: PASS" << std::endl;
   }

   // Perform standard checks on new communicator
   if (StrideEnv->isMember()) {
      int StrideTask = StrideEnv->getMyTask();
      if (StrideTask == WorldTask / 2)
         std::cout << "strided task test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "strided task test: FAIL "
                   << "StrideTask, WorldTask = " << StrideTask << " "
                   << WorldTask << std::endl;
      }

      int StrideSize = StrideEnv->getNumTasks();
      if (StrideSize == 4)
         std::cout << "strided NumTasks test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "strided NumTasks test: FAIL "
                   << "StrideSize (should be 4)  = " << StrideSize << std::endl;
      }

      int StrideMaster = StrideEnv->getMasterTask();
      if (StrideMaster == 0)
         std::cout << "strided master task test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "strided master task test: FAIL "
                   << "master = " << StrideMaster << std::endl;
      }

      bool IsStrideMaster = StrideEnv->isMasterTask();
      if (StrideTask == StrideMaster) {
         if (IsStrideMaster)
            std::cout << "strided is master task test: PASS" << std::endl;
         else {
            RetVal += 1;
            std::cout << "strided is master task test: FAIL" << std::endl;
         }
      } else {
         if (IsStrideMaster) {
            RetVal += 1;
            std::cout << "strided is master task test: FAIL" << std::endl;
         } else
            std::cout << "strided is master task test: PASS" << std::endl;
      }

   } // end if member of strided subset

   //---------------------------------------------------------------------------
   // Test general subset constructor using tasks 1,2,5,7

   int InclSize     = 4;
   int InclTasks[4] = {1, 2, 5, 7};

   // Test retrieval
   OMEGA::MachEnv *SubsetEnv = OMEGA::MachEnv::getEnv("Subset");

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
      if (SubsetEnv->isMember())
         std::cout << "subset member test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "subset member test: FAIL" << std::endl;
      }

   } else {
      if (SubsetEnv->isMember()) {
         RetVal += 1;
         std::cout << "subset non-member test: FAIL" << " MyTask " << MyTask
                   << " InclTasks ";
         for (int i = 0; i < InclSize; ++i) {
            std::cout << InclTasks[i];
         }
         std::cout << std::endl;
      } else {
         std::cout << "subset non-member test: PASS " << std::endl;
      }
   }

   // Perform standard checks on new communicator
   if (SubsetEnv->isMember()) {
      int SubsetTask = SubsetEnv->getMyTask();
      if (SubsetTask == NewTask)
         std::cout << "subset task test: PASS " << std::endl;
      else {
         RetVal += 1;
         std::cout << "subset task test: FAIL "
                   << "SubsetTask, NewTask = " << SubsetTask << " " << NewTask
                   << std::endl;
      }

      int SubsetSize = SubsetEnv->getNumTasks();
      if (SubsetSize == InclSize)
         std::cout << "subset size test: PASS " << std::endl;
      else {
         RetVal += 1;
         std::cout << "subset size test: FAIL "
                   << "SubsetSize, InclSize  = " << SubsetSize << " "
                   << InclSize << std::endl;
      }

      int SubsetMaster = SubsetEnv->getMasterTask();
      if (SubsetMaster == 0)
         std::cout << "subset master task test: PASS" << std::endl;
      else {
         RetVal += 1;
         std::cout << "subset master task test: FAIL" << std::endl;
         std::cout << "master = " << SubsetMaster << std::endl;
      }

      bool IsSubsetMaster = SubsetEnv->isMasterTask();
      if (SubsetTask == SubsetMaster) {
         if (IsSubsetMaster)
            std::cout << "subset is master task test: PASS" << std::endl;
         else {
            RetVal += 1;
            std::cout << "subset is master task test: FAIL" << std::endl;
         }
      } else {
         if (IsSubsetMaster) {
            RetVal += 1;
            std::cout << "subset is master task test: FAIL" << std::endl;
         } else
            std::cout << "subset is master task test: PASS" << std::endl;
      }

   } // end if member of general subset env

   //---------------------------------------------------------------------------
   // Test setting of compile-time vector length

#ifdef OMEGA_VECTOR_LENGTH
   if (OMEGA::VecLength == 16)
      std::cout << "MPI vector length test: PASS" << std::endl;
   else {
      RetVal += 1;
      std::cout << "MPI vector length test: FAIL" << std::endl;
      std::cout << "Was test driver built with -D OMEGA_VECTOR_LENGTH=16 ?"
                << std::endl;
   }
#endif

   // finalize and clean up environments (test both removal functions)
   OMEGA::MachEnv::removeEnv("Contig");
   OMEGA::MachEnv::removeAll();

   // MPI_Status status;
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;
} // end of main
//===-----------------------------------------------------------------------===/
