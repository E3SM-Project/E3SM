//===-- Test driver for OMEGA error handler ----------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA error handler
///
/// This driver tests the error handling capabilities for the OMEGA
/// model. In particular, it tests error checking for various error types
/// and the writing of stack traces for critical errors.
///
//
//===-----------------------------------------------------------------------===/

#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"

#include <iostream>

using namespace OMEGA;

//------------------------------------------------------------------------------
// Create a call tree for less critical errors to test passing return codes
// up the chain

// At lowest level, create and return a warning message
Error returnErrorRoutine3(void) {
   Error RetVal; // init to success
   int Level = 3;
   RETURN_ERROR(RetVal, ErrorCode::Warn, "Test warning at level {}", Level);
}

// At second level, add another warning message
Error returnErrorRoutine2(void) {
   Error RetVal; // init to success
   int Level = 2;
   RetVal += returnErrorRoutine3(); // tests accumulator
   RETURN_ERROR(RetVal, ErrorCode::Warn, "Test warning at level {}", Level);
}

Error returnErrorRoutine1(void) {
   Error ErrorTmp1; // init to success
   int Level       = 1;
   Error ErrorTmp2 = returnErrorRoutine2();
   Error RetVal    = ErrorTmp1 + ErrorTmp2;
   // upgrade from warning to fail
   RETURN_ERROR(RetVal, ErrorCode::Fail, "Test error at level {}", Level);
}

//------------------------------------------------------------------------------
// Create critical error call chain to test an abort and a correct stack trace

void passError2(Error &Err) {
   R8 MsgArgReal         = 1.23456789012;
   I4 MsgArgInt          = 5;
   std::string MsgArgStr = "Random";
   if (!Err.isSuccess())
      ABORT_ERROR("Expected critical error with args: {} {} {}", MsgArgStr,
                  MsgArgInt, MsgArgReal);
   // Alternative critical failure tests
   // CHECK_ERROR_ABORT(Err, "Expected critical error with args: {} {} {}",
   //                   MsgArgStr, MsgArgInt, MsgArgReal);
   // OMEGA_REQUIRE(Err.isSuccess(),
   //               "Expected Error with args: {} {} {}", MsgArgStr, MsgArgInt,
   //                MsgArgReal);
   // This one requires a debug build
   // OMEGA_ASSERT(Err.isSuccess(),
   //              "Expected Error with args: {} {} {}", MsgArgStr, MsgArgInt,
   //              MsgArgReal);
}

void passError1(Error &Err) { passError2(Err); }

//------------------------------------------------------------------------------
// Main test driver
int main(int argc, char **argv) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *DefEnv = OMEGA::MachEnv::getDefault();
   initLogging(DefEnv);

   LOG_INFO("------ Error handler test ------");
   // Test default constructor and success check
   Error ErrDefault;
   if (ErrDefault.isSuccess()) {
      LOG_INFO("Default error constructor: PASS");
   } else {
      ++RetVal;
      LOG_ERROR("Default error constructor: FAIL");
   }

   // Test constructors with no message and boolean tests
   Error ErrorTestSuccess2(ErrorCode::Success);
   if (ErrorTestSuccess2.isSuccess() and !ErrorTestSuccess2.isFail()) {
      LOG_INFO("No message error constructor success: PASS");
   } else {
      ++RetVal;
      LOG_ERROR("No message error constructor success: FAIL");
   }
   Error ErrorTestWarn2(ErrorCode::Warn);
   if (!ErrorTestWarn2.isSuccess() and !ErrorTestWarn2.isFail() and
       ErrorTestWarn2.Code == ErrorCode::Warn) {
      LOG_INFO("No message error constructor warn: PASS");
   } else {
      ++RetVal;
      LOG_ERROR("No message error constructor warn: FAIL");
   }
   Error ErrorTestFail2(ErrorCode::Fail);
   if (ErrorTestFail2.isFail() and !ErrorTestFail2.isSuccess()) {
      LOG_INFO("No message error constructor fail: PASS");
   } else {
      ++RetVal;
      LOG_ERROR("No message error constructor fail: FAIL");
   }

   // Test constructor with message
   // Create a total error to accumulate codes and messages
   std::string LineTxt   = std::to_string(__LINE__ + 3);
   std::string ExpectMsg = "   [error] [ErrorTest.cpp:" + LineTxt +
                           "] Initialized to fail String 1\n";
   Error TotalError(ErrorCode::Fail, __LINE__, __FILE__,
                    "Initialized to fail {} {}", "String", 1);
   if (TotalError.isFail() and TotalError.Msg == ExpectMsg) {
      LOG_INFO("Error constructor with message: PASS");
   } else {
      ++RetVal;
      LOG_ERROR("Error constructor with message: FAIL");
      LOG_ERROR(" Error constructor with message expected: {}", ExpectMsg);
      LOG_ERROR(" Error constructor with message actual:   {}", TotalError.Msg);
   }

   // Test constructor with an argument that itself contains braces, which
   // must be inserted verbatim and not be taken for a placeholder
   {
      std::string BraceLine   = std::to_string(__LINE__ + 3);
      std::string BraceExpect = "   [error] [ErrorTest.cpp:" + BraceLine +
                                "] Bad units m^{-2} in field\n";
      Error BraceError(ErrorCode::Fail, __LINE__, __FILE__,
                       "Bad units {} in field", "m^{-2}");
      if (BraceError.isFail() and BraceError.Msg == BraceExpect) {
         LOG_INFO("Error constructor with braces in argument: PASS");
      } else {
         ++RetVal;
         LOG_ERROR("Error constructor with braces in argument: FAIL");
         LOG_ERROR(" expected: {}", BraceExpect);
         LOG_ERROR(" actual:   {}", BraceError.Msg);
      }
   }

   // Test reset of error code
   TotalError.reset();
   if (TotalError.isSuccess()) {
      LOG_INFO("Reset of error code: PASS");
   } else {
      ++RetVal;
      CHECK_ERROR(TotalError, "Reset of error code: FAIL");
   }

   // Call routines to check error return functions
   TotalError += returnErrorRoutine1();

   // Test ERROR_CHECK for adding and printing error message
   CHECK_ERROR(TotalError, "Encountered error in returnErrorRoutine1");

   // Test ERROR_CHECK for adding warning and error message
   CHECK_ERROR_WARN(TotalError, "Encountered error in returnErrorRoutine1");

   // The above calls have reset the error so create a new failure
   TotalError += Error(ErrorCode::Fail, __LINE__, __FILE__, "Error in driver");

   // If we've encountered any test failures to this point, return a success
   // code (since unit test is expecting a failure)
   if (RetVal > 0) {
      return 0;

      // Otherwise, test critical error which should abort
   } else {
      passError1(TotalError);
   }

   // If code reaches here, the critical error test has failed so return a
   // success code
   MPI_Barrier(MPI_COMM_WORLD);
   MPI_Finalize();
   return 0;
}
//------------------------------------------------------------------------------
