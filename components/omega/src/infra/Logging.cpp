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
#include <cstdio>
#include <fstream>
#include <iostream>
#include <spdlog/sinks/basic_file_sink.h>

#define _OMEGA_STRINGIFY(x) #x
#define _OMEGA_TOSTRING(x)  _OMEGA_STRINGIFY(x)

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
// Utility function to determine which tasks perform logging based on
// input variable string containing lists of tasks or other options
std::vector<int>
splitTasks(const std::string &str,  //[in] option for which tasks should write
           const OMEGA::I4 NumTasks //[in] total number of tasks
) {

   std::vector<int> Tasks;
   std::stringstream ss(str);
   std::string Task;
   int Start, Stop;
   char Dash;

   // If all tasks should write, create the explicit list of all tasks
   if (str == "ALL") {
      for (int i = 0; i < NumTasks; ++i) {
         Tasks.push_back(i);
      }
      return Tasks;
   }

   // Parse the task string. If there is a single task, extract the task id.
   // If there is a range (beg:end), create a list with that range
   while (std::getline(ss, Task, ':')) {
      std::istringstream iss(Task);
      iss >> Start;

      if (iss.eof()) { // no range found, so extract single task
         Tasks.push_back(Start);

      } else {
         iss >> Dash;
         iss >> Stop;

         for (int i = Start; i <= Stop; ++i) {
            Tasks.push_back(i);
         }
      }
   }

   // Default to all tasks if a -1 task found in list
   if (std::find(Tasks.begin(), Tasks.end(), -1) != Tasks.end()) {
      Tasks.clear();
      for (int i = 0; i < NumTasks; ++i) {
         Tasks.push_back(i);
      }
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

   OMEGA::I4 TaskId   = DefEnv->getMyTask();
   OMEGA::I4 NumTasks = DefEnv->getNumTasks();

   // determine which tasks will write log files
   std::vector<int> Tasks =
       splitTasks(_OMEGA_TOSTRING(OMEGA_LOG_TASKS), NumTasks);

   if (Tasks.size() > 0 &&
       (std::find(Tasks.begin(), Tasks.end(), TaskId) != Tasks.end())) {

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

   OMEGA::I4 TaskId   = DefEnv->getMyTask();
   OMEGA::I4 NumTasks = DefEnv->getNumTasks();
   std::string NewLogFilePath;

   // Determine which tasks will write log files
   std::vector<int> Tasks =
       splitTasks(_OMEGA_TOSTRING(OMEGA_LOG_TASKS), NumTasks);

   if (Tasks.size() > 0 &&
       (std::find(Tasks.begin(), Tasks.end(), TaskId) != Tasks.end())) {

      try {
         std::size_t dotPos = LogFilePath.find_last_of('.');

         // create log file name/path and set default (*) logger
         if (Tasks.size() > 1 && dotPos != std::string::npos) {
            NewLogFilePath = LogFilePath.substr(0, dotPos) + "_" +
                             std::to_string(TaskId) +
                             LogFilePath.substr(dotPos);
         } else {
            NewLogFilePath = LogFilePath;
         }

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
   }

   return RetVal;
}

} // namespace OMEGA
//===----------------------------------------------------------------------===//
