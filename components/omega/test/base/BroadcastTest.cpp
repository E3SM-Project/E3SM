//===-- Test driver for OMEGA Broadcast ----------------------------*- C++
//-*-===/
//
/// \file
/// \brief Test driver for OMEGA Broadcast
///
/// This driver tests the OMEGA model Broadcast module that wraps MPI broadcast
/// functions for Omega data types.
///
//
//===-----------------------------------------------------------------------===/

#include "Broadcast.h"
#include "mpi.h"

#include <iostream>

template <class MyType>
void TestBroadcast(OMEGA::MachEnv *Env, std::string TypeName, int *RetVal) {

   const int MyTask     = Env->getMyTask();
   const int IsMyMaster = Env->isMasterTask();
   const int RootTask   = 2;

   MyType MyVal, FromVal, ToVal;

   if constexpr (std::is_same_v<MyType, bool>) {
      FromVal = true;
      ToVal   = false;

   } else if constexpr (std::is_same_v<MyType, std::string>) {
      FromVal = "a";
      ToVal   = "b";

   } else {
      FromVal = 1;
      ToVal   = -1;
   }

   if (IsMyMaster)
      MyVal = FromVal;
   else
      MyVal = ToVal;

   // Test broadcasting value from master task
   OMEGA::Broadcast(MyVal);

   if (!IsMyMaster) {
      if (MyVal == FromVal)
         std::cout << TypeName << " scalar broadcast from a default task: PASS"
                   << std::endl;
      else {
         std::cout << TypeName << " scalar broadcast from a default task: FAIL"
                   << std::endl;
         *RetVal += 1;
      }
   }

   if (MyTask == RootTask)
      MyVal = FromVal;
   else
      MyVal = ToVal;

   // Test broadcasting value from a non-master task
   OMEGA::Broadcast(MyVal, Env, RootTask);

   if (IsMyMaster) {
      if (MyVal == FromVal)
         std::cout << TypeName
                   << " scalar broadcast from a non-master task: PASS"
                   << std::endl;
      else {
         std::cout << TypeName
                   << " scalar broadcast from a non-master task: FAIL"
                   << std::endl;
         *RetVal += 1;
      }
   }

   if (MyTask == RootTask)
      MyVal = FromVal;
   else
      MyVal = ToVal;

   // Test broadcasting value from a non-master task
   OMEGA::Broadcast(MyVal, RootTask);

   if (IsMyMaster) {
      if (MyVal == FromVal)
         std::cout << TypeName
                   << " scalar broadcast from the default env.: PASS"
                   << std::endl;
      else {
         std::cout << TypeName
                   << " scalar broadcast from the default env.: FAIL"
                   << std::endl;
         *RetVal += 1;
      }
   }

   // length of vector<string> is not fixed
   // elements of vector<bool> seems to be non-addressable
   if constexpr (!std::is_same_v<MyType, std::string> &&
                 !std::is_same_v<MyType, bool>) {

      std::vector<MyType> MyVector;

      for (int i = 1; i <= 5; i++) {
         if (MyTask == RootTask)
            MyVector.push_back(FromVal);
         else
            MyVector.push_back(ToVal);
      }

      // Test broadcasting value from a non-master task
      OMEGA::Broadcast(MyVector, RootTask);

      if (IsMyMaster) {
         if (std::all_of(MyVector.cbegin(), MyVector.cend(),
                         [&](MyType v) { return v == FromVal; }))
            std::cout << TypeName
                      << " vector broadcast from the default env.: PASS"
                      << std::endl;
         else {
            std::cout << TypeName
                      << " vector broadcast from the default env.: FAIL"
                      << std::endl;
            *RetVal += 1;
         }
      }
   }
}

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
   int WorldSize;
   MPI_Comm_size(MPI_COMM_WORLD, &WorldSize);

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

   // Get the Default environment
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefaultEnv();

   // I4 Broadcast tests
   TestBroadcast<OMEGA::I4>(DefEnv, "I4", &RetVal);

   // I8 Broadcast tests
   TestBroadcast<OMEGA::I8>(DefEnv, "I8", &RetVal);

   // R4 Broadcast tests
   TestBroadcast<OMEGA::R4>(DefEnv, "R4", &RetVal);

   // R8 Broadcast tests
   TestBroadcast<OMEGA::R8>(DefEnv, "R8", &RetVal);

   // Real Broadcast tests
   TestBroadcast<OMEGA::Real>(DefEnv, "Real", &RetVal);

   // boolean Broadcast tests
   TestBroadcast<bool>(DefEnv, "bool", &RetVal);

   // string Broadcast tests
   TestBroadcast<std::string>(DefEnv, "string", &RetVal);

   // Initialize general subset environment
   int InclSize     = 4;
   int InclTasks[4] = {1, 2, 5, 7};
   OMEGA::MachEnv *SubsetEnv =
       OMEGA::MachEnv::create("Subset", DefEnv, InclSize, InclTasks);
   OMEGA::I4 MyVal;

   const int MyTask = DefEnv->getMyTask();

   if (MyTask == InclTasks[0])
      MyVal = 1;
   else
      MyVal = -1;

   // Test broadcasting value within a subset-env

   OMEGA::I4 *pos =
       std::find(std::begin(InclTasks), std::end(InclTasks), MyTask);

   if (pos != std::end(InclTasks))
      OMEGA::Broadcast(MyVal, SubsetEnv);

   if (pos != std::end(InclTasks)) {
      if (MyVal == 1)
         std::cout << "I4" << " sub-group broadcast at rank " << MyTask
                   << " : PASS" << std::endl;
      else {
         std::cout << "I4" << " sub-group broadcast at rank " << MyTask
                   << " : FAIL" << std::endl;
         RetVal += 1;
      }
   } else {
      if (MyVal == -1)
         std::cout << "I4" << " sub-group broadcast at rank " << MyTask
                   << " : PASS" << std::endl;
      else {
         std::cout << "I4" << " sub-group broadcast at rank " << MyTask
                   << " : FAIL" << std::endl;
         RetVal += 1;
      }
   }

   OMEGA::MachEnv::removeEnv("Subset");

   // MPI_Status status;
   MPI_Finalize();

   if (RetVal >= 256)
      RetVal = 255;

   return RetVal;

} // end of main
//===-----------------------------------------------------------------------===/
