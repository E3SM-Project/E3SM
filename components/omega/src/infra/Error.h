#ifndef OMEGA_ERROR_H
#define OMEGA_ERROR_H
//===-- Error.h - Omega error handling facility -----------------*- C++ -*-===//
//
/// \file
/// \brief Utilities for trapping and managing exceptions
///
/// This header defines variables, routines and macros for trapping exceptions
/// in Omega and logging those errors with stack traces if needed. It is built
/// on both the Omega logging facility with spdlog and the cpptrace library for
/// writing stack traces.
//
//===----------------------------------------------------------------------===//

#include "DataTypes.h"
#include "Logging.h"
#include <cpptrace/cpptrace.hpp>
#include <spdlog/spdlog.h>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace OMEGA {

//------------------------------------------------------------------------------
/// Standard error codes. The order must be increasing severity.
enum class ErrorCode {
   Success = 0, /// Standard success code
   Warn,        /// Standard warning code
   Fail,        /// Standard failure code
   Critical     /// Standard critical error code
};

//------------------------------------------------------------------------------
/// Error class for managing error and return codes and messages
class Error {

 public:
   ErrorCode Code;  /// error code to return
   std::string Msg; /// message to return with error code (can be blank)

   /// Default constructor creates an instance with success code and empty msg
   Error();

   /// Constructor with no message
   Error(ErrorCode ErrCode ///< [in] error code to assign
   );

   /// Constructor with message and prefix
   template <typename... Args>
   Error(ErrorCode ErrCode,        ///< [in] error code to assign
         int Line,                 ///< [in] line where err occurs (cpp LINE)
         const char *File,         ///< [in] file where err occurs (cpp FILE)
         const std::string &InMsg, ///< [in] error message
         Args &&...MsgArgs         ///< [in] additional arguments for error msg
   );

   /// Constructor with message and no prefix
   template <typename... Args>
   Error(ErrorCode ErrCode,        ///< [in] error code to assign
         const std::string &InMsg, ///< [in] error message
         Args &&...MsgArgs         ///< [in] additional arguments for error msg
   );

   /// Query error for success. True only if the error code is success and is
   /// false for warn and fail.
   bool isSuccess();

   /// Query error for failure. True only if error code is fail and is false
   /// for both warn and success.
   bool isFail();

   /// Resets an error code to success with an empty message
   void reset();

   /// Aborts the simulation by calling the MPI abort command
   static void abort();

   /// Operator for adding two Errors. The result is the maximum of the
   /// two error codes and the concantenation of error messages (with a
   /// new line between them).
   Error operator+(const Error &) const;

   /// Increment operator for adding two Errors. The result is the maximum of
   /// the two error codes and the concantenation of error messages (with a
   /// new line between them). The new message is prepending to the existing
   /// message.
   Error &operator+=(const Error &);

 private:
   /// Utility function to create error message by inserting arguments into
   /// placeholder locations and adding a prefix.
   static std::string createErrMsg(
       const std::string &Prefix, ///< [in] prefix with info (can be empty)
       const std::string &InMsg,  ///< [in] error message with placeholders
       std::vector<std::string>
           MsgArgs ///< [in] additional arguments for error msg
   );

}; // end Error class

//------------------------------------------------------------------------------
// Implementation for templated functions

//------------------------------------------------------------------------------
// Utility function to convert arguments to strings and add them to a
// vector of string arguments
template <typename T>
void buildErrorArgList(std::vector<std::string> &ArgVector, T Arg) {
   std::ostringstream TempStr;
   TempStr << Arg;
   ArgVector.push_back(TempStr.str());
}

//------------------------------------------------------------------------------
// Constructor with error message and prefix. The prefix is based on the error
// level and file/line location. The error message can contain {} placeholders
// for the additional arguments.
template <typename... Args>
Error::Error(ErrorCode ErrCode, // [in] error code to assign
             int Line,          // [in] line where err occurs (cpp __LINE__)
             const char *File,  // [in] file where err occurs (cpp __FILE__)
             const std::string &InMsg, // [in] error message
             Args &&...MsgArgs // [in] additional arguments for error msg
) {

   // Convert input variadic arguments to strings and store in a vector of
   // strings. This uses a fold expression to loop over the arg list.
   std::vector<std::string> ArgStrings;
   (buildErrorArgList(ArgStrings, std::forward<Args>(MsgArgs)), ...);

   // Convert file argument into a string for easier conversion
   std::string FileWPath(File);

   // Strip pathname from file argument
   std::string Filename;
   auto Pos = FileWPath.find_last_of("\\/");
   if (Pos != std::string::npos) {
      Filename = FileWPath.substr(Pos + 1);
   } else {
      Filename = FileWPath; // no path detected so use full string
   }

   // Begin constructing the message prefix based on the error code level.
   // Also add an indent since this message will appear under any final error
   // message.

   std::string MsgPrefix = "   [";
   switch (ErrCode) {
   case ErrorCode::Critical:
      MsgPrefix += "critical] ";
      break;
   case ErrorCode::Fail:
      MsgPrefix += "error] ";
      break;
   case ErrorCode::Warn:
      MsgPrefix += "warn] ";
      break;
   case ErrorCode::Success:
      MsgPrefix = "";
      break;
   default:
      MsgPrefix = "";
   }

   // Add the filename and line number in square brackets
   MsgPrefix += "[" + Filename + ":" + std::to_string(Line) + "] ";

   // Create the final message based on input arguments.
   Msg = Error::createErrMsg(MsgPrefix, InMsg, ArgStrings);

   // Store the code
   Code = ErrCode;
}

//------------------------------------------------------------------------------
// Constructor with error message and no prefix. The error message can contain
// {} placeholders for the additional arguments.
template <typename... Args>
Error::Error(ErrorCode ErrCode,        // [in] error code to assign
             const std::string &InMsg, // [in] error message
             Args &&...MsgArgs // [in] additional arguments for error msg
) {

   // Convert input variadic arguments to strings and store in a vector of
   // strings. This uses a fold expression to loop over the arg list.
   std::vector<std::string> ArgStrings;
   (buildErrorArgList(ArgStrings, std::forward<Args>(MsgArgs)), ...);

   // Use blank prefix
   std::string MsgPrefix = "";

   // Create the final message string based on input arguments.
   Msg = Error::createErrMsg(MsgPrefix, InMsg, ArgStrings);

   // Store the code
   Code = ErrCode;
}

} // namespace OMEGA
//------------------------------------------------------------------------------
// Define Macros for Error handling interfaces. Macros are used to obtain
// and preserve the line number at which an error occurs.

///  The critical error macro prints an error message and aborts the simulation
#define ABORT_ERROR(_Msg, ...)            \
   {                                      \
      LOG_CRITICAL(_Msg, ##__VA_ARGS__);  \
      cpptrace::generate_trace().print(); \
      OMEGA::Error::abort();              \
   }

/// The assert macro checks for a required condition and exits if not met
/// It is only active for debug builds.
#ifdef OMEGA_DEBUG
#define OMEGA_ASSERT(_Condition, _ErrMsg, ...) \
   if (!(_Condition)) {                        \
      LOG_CRITICAL(_ErrMsg, ##__VA_ARGS__);    \
      cpptrace::generate_trace().print();      \
      OMEGA::Error::abort();                   \
   }
#else
#define OMEGA_ASSERT(_Condition, _ErrMsg, ...)
#endif

/// This macro checks for a required condition and exits if not met
/// It is always evaluated (unlike ASSERT for debug builds)
#define OMEGA_REQUIRE(_Condition, _ErrMsg, ...) \
   if (!(_Condition)) {                         \
      LOG_CRITICAL(_ErrMsg, ##__VA_ARGS__);     \
      cpptrace::generate_trace().print();       \
      OMEGA::Error::abort();                    \
   }

/// This macro adds an error to the return code and returns this code to the
/// calling routine.
#define RETURN_ERROR(_ReturnErr, _ErrCode, _ErrMsg, ...)          \
   {                                                              \
      OMEGA::Error _NewErr(_ErrCode, __LINE__, __FILE__, _ErrMsg, \
                           ##__VA_ARGS__);                        \
      return (_ReturnErr += _NewErr);                             \
   }

/// This macro checks an existing error code and if it is not success, it
/// prints the error message, resets the code and continues the simulation.
#define CHECK_ERROR(_ErrChk, _NewMsg, ...)                              \
   if (!_ErrChk.isSuccess()) {                                          \
      _ErrChk += Error(OMEGA::ErrorCode::Fail, _NewMsg, ##__VA_ARGS__); \
      LOG_ERROR(_ErrChk.Msg);                                           \
      _ErrChk.reset();                                                  \
   }

/// This macro checks an existing error code and if it is not success, it
/// prints a warning message, resets the code and continues the simulation.
#define CHECK_ERROR_WARN(_ErrChk, _NewMsg, ...)                         \
   if (!_ErrChk.isSuccess()) {                                          \
      _ErrChk += Error(OMEGA::ErrorCode::Warn, _NewMsg, ##__VA_ARGS__); \
      LOG_WARN(_ErrChk.Msg);                                            \
      _ErrChk.reset();                                                  \
   }

/// This macro checks an existing error code and if it is a fail, it
/// aborts with a critical error
#define CHECK_ERROR_ABORT(_ErrChk, _NewMsg, ...)                        \
   if (_ErrChk.isFail()) {                                              \
      _ErrChk += Error(OMEGA::ErrorCode::Fail, _NewMsg, ##__VA_ARGS__); \
      LOG_CRITICAL(_ErrChk.Msg);                                        \
      cpptrace::generate_trace().print();                               \
      OMEGA::Error::abort();                                            \
   }

//===----------------------------------------------------------------------===//
#endif // OMEGA_ERROR_H
