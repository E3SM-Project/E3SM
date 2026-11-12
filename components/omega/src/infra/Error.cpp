//===-- Error.cpp - Omega error handling implementation ---------*- C++ -*-===//
//
/// \file
/// \brief Utilities for trapping and managing exceptions
///
/// This file contains the implementation of the error handling functions for
/// Omega. The Macros and some templated functions are implemented in the
/// header file Error.h
//
//===----------------------------------------------------------------------===//

#include "Error.h"
#include "DataTypes.h"
#include "Logging.h"
#include "mpi.h"
#include <cpptrace/cpptrace.hpp>
#include <cstdio>
#include <cstdlib>
#include <memory>
#include <spdlog/spdlog.h>
#include <string>

namespace OMEGA {
//------------------------------------------------------------------------------
// Default constructor
Error::Error() : Code{ErrorCode::Success}, Msg{""} {}

//------------------------------------------------------------------------------
// Constructor with no message
Error::Error(ErrorCode ErrCode) // [in] error code to assign
    : Code{ErrCode}, Msg{""} {}

//------------------------------------------------------------------------------
// Constructor with msg is defined in header

//------------------------------------------------------------------------------
// Resets an error to a success with no message
void Error::reset() {
   Code = ErrorCode::Success;
   Msg  = "";
}

//------------------------------------------------------------------------------
// Aborts the simulation by calling the MPI abort command
void Error::abort() {
   LOG_CRITICAL("Omega aborting");
   spdlog::apply_all([](std::shared_ptr<spdlog::logger> L) { L->flush(); });
   int ErrCode = static_cast<int>(ErrorCode::Critical);

   // With a single task, Cray MPICH's MPI_Abort calls exit(), which runs the
   // destructors of static objects (e.g. the VertCoord registry). Those free
   // Kokkos Views after Kokkos is no longer usable, which crashes or spins
   // forever and leaves the job step hung. Exit without running them instead.
   int MPIInitialized = 0;
   int MPIFinalized   = 0;
   MPI_Initialized(&MPIInitialized);
   MPI_Finalized(&MPIFinalized);
   int NTasks = 1;
   if (MPIInitialized && !MPIFinalized)
      MPI_Comm_size(MPI_COMM_WORLD, &NTasks);
   if (NTasks > 1)
      MPI_Abort(MPI_COMM_WORLD, ErrCode);

   std::fflush(nullptr);
   std::_Exit(ErrCode);
}

//------------------------------------------------------------------------------
// Query error for success - true if the error code is success
bool Error::isSuccess() { return (Code == ErrorCode::Success); }

//------------------------------------------------------------------------------
// Query error for failure. True only if error code is fail so this
// is not equivalent to not success.
bool Error::isFail() { return (Code == ErrorCode::Fail); }

//------------------------------------------------------------------------------
// Templated function createMsg is defined in header

//------------------------------------------------------------------------------
// Operator for adding two Errors. The result is the maximum of the
// two error codes and the concantenation of error messages.
Error Error::operator+(const Error &Err) const {

   Error Result(std::max(Code, Err.Code));

   Result.Msg = Msg + Err.Msg;

   return Result;
}

//------------------------------------------------------------------------------
// Increment operator for adding two Errors. The result is the maximum of the
// two error codes and the concantenation of error messages. This form
// guarantes the new message is prepended so that the latest message appears
// first.

Error &Error::operator+=(const Error &Err) {

   // reuse (+) operator defined above
   *this = Err + *this;
   return *this;

} // end Error increment operator

//------------------------------------------------------------------------------
// Utility function to create an error message by inserting arguments into
// placeholder locations and adding a prefix. A newline character is also
// added to the end.
std::string Error::createErrMsg(
    const std::string &Prefix,       // [in] prefix with info (can be empty)
    const std::string &InMsg,        // [in] error message with placeholders
    std::vector<std::string> MsgArgs // [in] additional arguments for msg
) {

   // Initialize the result string with a prefix (can be an empty string)
   // and the input message string with placeholders
   std::string Result = Prefix + InMsg;

   // Look for placeholder {} pairs to replace with each argument and count to
   // make sure the number of expected and actual arguments match
   int EndPos   = 0;
   int ArgIndex = 0;
   int StartPos = 0;
   while ((StartPos = Result.find('{', EndPos)) != std::string::npos) {
      EndPos = Result.find('}', StartPos);
      if (EndPos != std::string::npos) {
         // Replace {} with argument
         int Length = EndPos - StartPos + 1;
         if (ArgIndex < MsgArgs.size()) {
            Result.replace(StartPos, Length, MsgArgs[ArgIndex]);
            // Resume scanning after the inserted argument so that braces
            // inside it are not taken for placeholders
            EndPos = StartPos + MsgArgs[ArgIndex].size();
            ArgIndex++;
         } else {
            LOG_ERROR("Not enough arguments for placeholders in error msg");
            break;
         }

      } else {
         LOG_ERROR("Unpaired curly bracket in error message format");
         break;
      }
   }

   if (ArgIndex < MsgArgs.size())
      LOG_ERROR("Too many arguments for placeholders in error msg {} {}",
                ArgIndex, MsgArgs.size());

   // Add newline to end and return final message
   Result += "\n";
   return Result;

} // end createErrMsg

} // namespace OMEGA
//===----------------------------------------------------------------------===//
