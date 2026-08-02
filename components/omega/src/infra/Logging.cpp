//===-- infra/Logging.cpp - Implements Omega Logging functions --*- C++ -*-===//
//
/// \file
/// \brief Implements Omega logging functions
///
/// This implements Omega logging initialization and includes some utility
/// functions for formatting message prefixes.
//
//===----------------------------------------------------------------------===//
#include "Logging.h"
#include "MachEnv.h"
#include <algorithm>
#include <cctype>
#include <cstdio>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iostream>
#include <limits>
#include <set>
#include <spdlog/sinks/basic_file_sink.h>
#include <sstream>
#include <vector>

namespace OMEGA {

//------------------------------------------------------------------------------
// Utility function to create a log message with prefix
std::string
_PackLogMsg(const char *file, // [in] file from where log called (cpp __FILE__)
            int line,         // [in] src code line where called (cpp __LINE__)
            const std::string &msg // [in] message text
) {
   // find position after path where filename starts
   std::string path(file);
   size_t pos = path.find_last_of("\\/");
   if (pos != std::string::npos) {
      // add prefix of form [file:line] to message string
      return "[" + path.substr(pos + 1) + ":" + std::to_string(line) + "] " +
             msg;
   }
   // no path separator found, use the full file string
   return "[" + path + ":" + std::to_string(line) + "] " + msg;
}

//------------------------------------------------------------------------------
// Trim surrounding whitespace from a string
static std::string trimStr(const std::string &Str) {
   const std::string WhiteSpace = " \t\n\r\f\v";
   std::size_t Begin            = Str.find_first_not_of(WhiteSpace);
   if (Begin == std::string::npos)
      return "";
   std::size_t End = Str.find_last_not_of(WhiteSpace);
   return Str.substr(Begin, End - Begin + 1);
}

//------------------------------------------------------------------------------
// Parse a string as a non-negative int. Returns false if the string is empty,
// contains any non-digit character, or overflows an int.
static bool toNonNegInt(const std::string &Str, int &Out) {
   if (Str.empty())
      return false;
   for (char C : Str) {
      if (!std::isdigit(static_cast<unsigned char>(C)))
         return false;
   }
   try {
      std::size_t Pos = 0;
      long Val        = std::stol(Str, &Pos);
      if (Pos != Str.size() || Val < 0 ||
          Val > static_cast<long>(std::numeric_limits<int>::max()))
         return false;
      Out = static_cast<int>(Val);
      return true;
   } catch (...) {
      return false;
   }
}

//------------------------------------------------------------------------------
// Resolve a logging task selector string into the list of ranks that should log
std::vector<int> _selectLogTasks(const std::string &Selector,
                                 OMEGA::I4 NumTasks, OMEGA::I4 MasterTask,
                                 bool &Valid) {

   Valid = true;

   std::string Sel = trimStr(Selector);

   // Lowercase copy for case-insensitive keyword matching
   std::string Lower = Sel;
   std::transform(Lower.begin(), Lower.end(), Lower.begin(),
                  [](unsigned char C) { return std::tolower(C); });

   // Standalone keyword: all ranks in the sub-communicator
   if (Lower == "*") {
      std::vector<int> All;
      for (int I = 0; I < NumTasks; ++I)
         All.push_back(I);
      return All;
   }

   // Standalone keyword: master rank only
   if (Lower == "m" || Lower == "master")
      return std::vector<int>{MasterTask};

   // Empty selector is malformed
   if (Sel.empty()) {
      Valid = false;
      return std::vector<int>{MasterTask};
   }

   // Numeric expression: comma-separated single ranks and/or "a-b" ranges.
   // A std::set keeps the result sorted and deduplicated.
   std::set<int> Selected;
   std::stringstream Ss(Sel);
   std::string Token;
   while (std::getline(Ss, Token, ',')) {
      Token = trimStr(Token);

      std::size_t DashPos = Token.find('-');
      if (DashPos == std::string::npos) {
         int Rank;
         if (!toNonNegInt(Token, Rank)) {
            Valid = false;
            return std::vector<int>{MasterTask};
         }
         Selected.insert(Rank);
      } else {
         int Lo, Hi;
         if (!toNonNegInt(trimStr(Token.substr(0, DashPos)), Lo) ||
             !toNonNegInt(trimStr(Token.substr(DashPos + 1)), Hi) || Lo > Hi) {
            Valid = false;
            return std::vector<int>{MasterTask};
         }
         for (int I = Lo; I <= Hi; ++I)
            Selected.insert(I);
      }
   }

   return std::vector<int>(Selected.begin(), Selected.end());
}

//------------------------------------------------------------------------------
// Source the logging task selector from the OMEGA_LOG_TASKS environment
// variable, falling back to the master rank when unset/empty.
static std::string getLogTaskSelector() {
   const char *Env = std::getenv("OMEGA_LOG_TASKS");
   if (Env != nullptr && Env[0] != '\0')
      return std::string(Env);
   return std::string("master");
}

//------------------------------------------------------------------------------
// Determine the log file name to use when the caller does not supply one.
// The OMEGA_LOG_FILE environment variable takes precedence so that a harness
// (eg CTest, which gives each unit test its own log) can name the file at run
// time; otherwise the built-in default name is used.
static std::string getDefaultLogFilePath() {
   const char *Env = std::getenv("OMEGA_LOG_FILE");
   if (Env != nullptr && Env[0] != '\0')
      return std::string(Env);
   return OmegaDefaultLogfile;
}

//------------------------------------------------------------------------------
// Insert the MPI task id into a log file name so that each logging rank gets
// its own file, eg omega.log -> omega_3.log. When the name has no extension,
// the task id is simply appended.
static std::string addTaskIdToPath(const std::string &LogFilePath,
                                   OMEGA::I4 TaskId) {

   const std::string Suffix = "_" + std::to_string(TaskId);

   std::size_t DotPos = LogFilePath.find_last_of('.');
   std::size_t SepPos = LogFilePath.find_last_of("\\/");

   // only treat a dot as an extension when it is part of the file name and
   // not part of a directory in the path
   bool HasExtension = (DotPos != std::string::npos) &&
                       (SepPos == std::string::npos || DotPos > SepPos);

   if (HasExtension)
      return LogFilePath.substr(0, DotPos) + Suffix +
             LogFilePath.substr(DotPos);

   return LogFilePath + Suffix;
}

//------------------------------------------------------------------------------
// Resolve which sub-communicator ranks should log from the runtime selector,
// emitting master-rank diagnostics for invalid or empty selections. Sets
// NumLogging to the count of selected ranks that exist in this communicator.
static std::vector<int> resolveLogTasks(const OMEGA::MachEnv *DefEnv,
                                        int &NumLogging) {

   OMEGA::I4 NumTasks   = DefEnv->getNumTasks();
   OMEGA::I4 MasterTask = DefEnv->getMasterTask();

   std::string Selector = getLogTaskSelector();
   bool Valid           = false;
   std::vector<int> Tasks =
       _selectLogTasks(Selector, NumTasks, MasterTask, Valid);

   if (!Valid && DefEnv->isMasterTask()) {
      std::cerr << "[Omega Logging] Invalid OMEGA_LOG_TASKS selector \""
                << Selector << "\"; falling back to master rank only."
                << std::endl;
   }

   // Count selected ranks that actually exist in this communicator
   NumLogging = 0;
   for (int Rank : Tasks) {
      if (Rank >= 0 && Rank < NumTasks)
         ++NumLogging;
   }
   if (Valid && NumLogging == 0 && DefEnv->isMasterTask()) {
      std::cerr << "[Omega Logging] OMEGA_LOG_TASKS selector \"" << Selector
                << "\" matches no ranks in this communicator; logging disabled."
                << std::endl;
   }

   return Tasks;
}

//------------------------------------------------------------------------------
// Initialize logging for case of a pre-defined custom logger
// Return code: 1->enabled, 0->disabled, negative values->errors
int initLogging(
    const OMEGA::MachEnv *DefEnv,          // [in] MachEnv for MPI task info
    std::shared_ptr<spdlog::logger> Logger // [in] custom logger to use
) {

   int RetVal = 0;

   OMEGA::I4 TaskId = DefEnv->getMyTask();

   int NumLogging;
   std::vector<int> Tasks = resolveLogTasks(DefEnv, NumLogging);

   bool ThisTaskLogs =
       (std::find(Tasks.begin(), Tasks.end(), TaskId) != Tasks.end());

   if (ThisTaskLogs) {

      spdlog::set_default_logger(Logger);

      // set prefix format - here n is logger name, l is level and v is msg txt
      spdlog::set_pattern("[%n %l] %v");
      // set default log level
      spdlog::set_level(
          static_cast<spdlog::level::level_enum>(SPDLOG_ACTIVE_LEVEL));
      // flush output buffers for levels above warn
      spdlog::flush_on(spdlog::level::warn);

      RetVal = 1; // log enabled

   } else {
      spdlog::set_level(spdlog::level::off);
      RetVal = 0; // log disabled
   }

   return RetVal;
}

//------------------------------------------------------------------------------
// Initialize logging for a default logger and provide log file name
// return code: 1->enabled, 0->disabled, negative values->errors
int initLogging(
    const OMEGA::MachEnv *DefEnv,  //[in] default environment for MPI info
    std::string const &LogFilePath //[in] log file name with path
) {

   int RetVal = 0;

   OMEGA::I4 TaskId = DefEnv->getMyTask();

   std::string NewLogFilePath;

   int NumLogging;
   std::vector<int> Tasks = resolveLogTasks(DefEnv, NumLogging);

   bool ThisTaskLogs =
       (std::find(Tasks.begin(), Tasks.end(), TaskId) != Tasks.end());

   if (ThisTaskLogs) {
      try {
         // create log file name/path and set default (*) logger. When more than
         // one rank logs, append the rank id to keep per-rank files distinct.
         // An empty path means the caller wants the default name, which may be
         // supplied by the environment.
         NewLogFilePath =
             LogFilePath.empty() ? getDefaultLogFilePath() : LogFilePath;

         if (NumLogging > 1)
            NewLogFilePath = addTaskIdToPath(NewLogFilePath, TaskId);

         // Start the run with an empty log file so that it only ever holds
         // messages from this run. Both the spdlog sink created below and the
         // stdout/stderr stream opened at the end of this function append to
         // the file, so that their writes interleave correctly rather than
         // overwriting each other; the file is therefore truncated here, ahead
         // of both, rather than by either of them.
         std::ofstream Truncate(NewLogFilePath, std::ios::trunc);
         Truncate.close();

         // Create default logger
         spdlog::set_default_logger(
             spdlog::basic_logger_mt("*", NewLogFilePath));

         // Set the prefix for all messages l is log level and v is msg txt
         spdlog::set_pattern("[%l] %v");
         // Set the default log level based on cpp input
         spdlog::set_level(
             static_cast<spdlog::level::level_enum>(SPDLOG_ACTIVE_LEVEL));
         // Flush the message buffer for all messages above warning level
         spdlog::flush_on(spdlog::level::warn);

         RetVal = 1; // log enabled

      } catch (spdlog::spdlog_ex const &Ex) {
         std::cout << "Log init failed: " << Ex.what() << std::endl;
         RetVal = -1; // error occured
      }

   } else {
      spdlog::set_level(spdlog::level::off);
      RetVal = 0; // log disabled
   }

   // If logging is successful, also redirect stdout and stderr to log file
   if (RetVal == 1) {
      // Open an output filestream to the new log file defined above
      LogFileStream.open(NewLogFilePath, std::ios::app);

      // Set the stdout and stderr buffers to the file streambuffer
      std::cout.rdbuf(LogFileStream.rdbuf());
      std::cerr.rdbuf(LogFileStream.rdbuf());

      // Open the log with a banner giving the file and the time it was
      // opened, so the log can always be tied to the run that wrote it
      std::time_t Now = std::time(nullptr);
      char TimeStr[64];
      if (std::strftime(TimeStr, sizeof(TimeStr), "%Y-%m-%d %H:%M:%S",
                        std::localtime(&Now)) == 0)
         TimeStr[0] = '\0';
      spdlog::info(fmt::runtime("----- Omega log " + NewLogFilePath +
                                " opened " + std::string(TimeStr) +
                                " on task " + std::to_string(TaskId) +
                                " -----"));
   }

   return RetVal;
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//
