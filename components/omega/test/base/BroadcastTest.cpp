//===-- Test driver for OMEGA Broadcast --------------------------*- C++ -*-===/
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
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "mpi.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Test function templated by type. Tests both default and an alternate MachEnv
template <class MyType>
void TestBroadcast(MachEnv *DefEnv, MachEnv *AltEnv, std::string TypeName) {

   // Set test values based on type
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

   // Test broadcasting value from master task in default environment
   int MyTask      = DefEnv->getMyTask();
   bool IsMyMaster = DefEnv->isMasterTask();
   if (IsMyMaster)
      MyVal = FromVal;
   else
      MyVal = ToVal;

   Broadcast(MyVal);

   if (MyVal != FromVal)
      ABORT_ERROR("BroadcastTest: scalar {} broadcast from master task: FAIL",
                  TypeName);

   // Test broadcasting from non-root task with default env
   int RootTask = 2;
   if (MyTask == RootTask)
      MyVal = FromVal;
   else
      MyVal = ToVal;

   Broadcast(MyVal, RootTask);

   if (MyVal != FromVal)
      ABORT_ERROR(
          "BroadcastTest: scalar {} broadcast from non-master task: FAIL",
          TypeName);

   // Test broadcasting value from a non-master task in default env using full
   // interface
   if (MyTask == RootTask)
      MyVal = FromVal;
   else
      MyVal = ToVal;

   Broadcast(MyVal, DefEnv, RootTask);

   if (MyVal != FromVal)
      ABORT_ERROR(
          "BroadcastTest: scalar {} broadcast from non-master task in full "
          "interface: FAIL",
          TypeName);

   // Test scalar broadcast from master task in the alternative MachEnv
   MyTask     = AltEnv->getMyTask();
   IsMyMaster = AltEnv->isMasterTask();
   if (IsMyMaster)
      MyVal = FromVal;
   else
      MyVal = ToVal;

   Broadcast(MyVal, AltEnv);

   if (AltEnv->isMember()) {
      if (MyVal != FromVal)
         ABORT_ERROR(
             "BroadcastTest: scalar {} broadcast in alt environment: FAIL",
             TypeName);
   } else { // if not a member of environment - should still have ToVal
      if (MyVal == FromVal)
         ABORT_ERROR(
             "BroadcastTest: scalar {} broadcast in alt environment: FAIL"
             " for non-member task",
             TypeName);
   } // endif member of AltEnv

   // Test scalar broadcast from non-master task in the alternative MachEnv
   if (MyTask == RootTask)
      MyVal = FromVal;
   else
      MyVal = ToVal;

   Broadcast(MyVal, AltEnv, RootTask);

   if (AltEnv->isMember()) {
      if (MyVal != FromVal)
         ABORT_ERROR(
             "BroadcastTest: scalar {} broadcast in alt environment: FAIL",
             TypeName);
   } else { // if not a member of environment - should still have ToVal
      if (MyVal == FromVal)
         ABORT_ERROR(
             "BroadcastTest: scalar {} broadcast in alt environment: FAIL"
             " for non-member task",
             TypeName);
   } // endif member of AltEnv

   // Test vector broadcasts for all but strings - use default env, non-master
   //   - cannot broadcast strings since length of vector<string> is not fixed
   MyTask     = DefEnv->getMyTask();
   IsMyMaster = DefEnv->isMasterTask();

   if constexpr (!std::is_same_v<MyType, std::string>) {

      // create and fill a vector with appropriate type
      std::vector<MyType> MyVector;

      for (int i = 1; i <= 5; i++) {
         if (MyTask == RootTask)
            MyVector.push_back(FromVal);
         else
            MyVector.push_back(ToVal);
      }

      // Test broadcasting vector from a non-master task
      Broadcast(MyVector, RootTask);

      int Icount = 0;
      for (int i = 0; i < 5; i++) {
         if (MyVector[i] != FromVal)
            ++Icount;
      }
      if (Icount > 0)
         ABORT_ERROR("BroadcastTest: Vector {} broadcast FAIL", TypeName);

   } // endif not string

} // End test routine

//------------------------------------------------------------------------------
// The test driver for MachEnv. This tests the values stored in the Default
// Environment and three other based on the three subsetting options.  All
// current values and get routines are tested.
//
int main(int argc, char *argv[]) {

   // Initialize the global MPI environment
   MPI_Init(&argc, &argv);
   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv  = MachEnv::getDefault();
   MPI_Comm DefComm = DefEnv->getComm();

   // Initialize the Logging system
   initLogging(DefEnv);
   LOG_INFO("------ Broadcast Unit Tests ------");

   // The subset environments create 4-task sub-environments so
   // make sure the unit test is run with at least 8 to adequately
   // test all subsets.
   int WorldSize;
   MPI_Comm_size(DefComm, &WorldSize);
   if (WorldSize < 8)
      ABORT_ERROR("BroadcastTest: FAIL Must run with at least 8 tasks");

   // Initialize general subset environment
   int InclSize       = 4;
   int InclTasks[4]   = {1, 2, 5, 7};
   MachEnv *SubsetEnv = MachEnv::create("Subset", DefEnv, InclSize, InclTasks);

   // I4 Broadcast tests
   TestBroadcast<I4>(DefEnv, SubsetEnv, "I4");

   // I8 Broadcast tests
   TestBroadcast<I8>(DefEnv, SubsetEnv, "I8");

   // R4 Broadcast tests
   TestBroadcast<R4>(DefEnv, SubsetEnv, "R4");

   // R8 Broadcast tests
   TestBroadcast<R8>(DefEnv, SubsetEnv, "R8");

   // Real Broadcast tests
   TestBroadcast<Real>(DefEnv, SubsetEnv, "Real");

   // boolean Broadcast tests
   TestBroadcast<bool>(DefEnv, SubsetEnv, "bool");

   // string Broadcast tests
   TestBroadcast<std::string>(DefEnv, SubsetEnv, "string");

   // Cleanup
   MachEnv::removeEnv("Subset");
   LOG_INFO("------ Broadcast Unit Tests Successful ------");

   // MPI_Status status;
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();

   return 0; // if we made it here, return a success code

} // end of main
//===-----------------------------------------------------------------------===/
