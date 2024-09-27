//===-- infra/Logging.cpp - Implements Omega Logging functions --*- C++ -*-===//
//
/// \file
/// \brief implements Omega logging functions
///
/// This implements Omega logging initialization.
//
//===----------------------------------------------------------------------===//
#include "Logging.h"
#include "MachEnv.h"
#include <iostream>
#include <spdlog/sinks/basic_file_sink.h>

#define _OMEGA_STRINGIFY(x) #x
#define _OMEGA_TOSTRING(x)  _OMEGA_STRINGIFY(x)

namespace OMEGA {

// Function to pack log message
std::string _PackLogMsg(const char *file, int line, const std::string &msg) {

   // check if MachEnv is initialized
   std::string path(file);
   size_t pos = path.find_last_of("\\/");
   if (pos != std::string::npos) {
      return "[" + path.substr(pos + 1) + ":" + std::to_string(line) + "] " +
             msg;
   }
   return "[" + path + ":" + std::to_string(line) + "] " + msg;
}

std::vector<int> splitTasks(const std::string &str, const OMEGA::I4 NumTasks) {
   std::vector<int> Tasks;
   std::stringstream ss(str);
   std::string Task;
   int Start, Stop;
   char Dash;

   if (str == "ALL") {
      for (int i = 0; i < NumTasks; ++i) {
         Tasks.push_back(i);
      }
      return Tasks;
   }

   while (std::getline(ss, Task, ':')) {
      std::istringstream iss(Task);
      iss >> Start;

      if (iss.eof()) {
         Tasks.push_back(Start);

      } else {
         iss >> Dash;
         iss >> Stop;

         for (int i = Start; i <= Stop; ++i) {
            Tasks.push_back(i);
         }
      }
   }

   if (std::find(Tasks.begin(), Tasks.end(), -1) != Tasks.end()) {
      Tasks.clear();
      for (int i = 0; i < NumTasks; ++i) {
         Tasks.push_back(i);
      }
   }

   return Tasks;
}

// return code: 1->enabled, 0->disabled, negative values->errors
int initLogging(const OMEGA::MachEnv *DefEnv,
                std::shared_ptr<spdlog::logger> Logger) {

   int RetVal = 0;

   OMEGA::I4 TaskId   = DefEnv->getMyTask();
   OMEGA::I4 NumTasks = DefEnv->getNumTasks();

   std::vector<int> Tasks =
       splitTasks(_OMEGA_TOSTRING(OMEGA_LOG_TASKS), NumTasks);

   if (Tasks.size() > 0 &&
       (std::find(Tasks.begin(), Tasks.end(), TaskId) != Tasks.end())) {

      spdlog::set_default_logger(Logger);

      spdlog::set_pattern("[%n %l] %v");
      spdlog::set_level(
          static_cast<spdlog::level::level_enum>(SPDLOG_ACTIVE_LEVEL));
      spdlog::flush_on(spdlog::level::warn);

      RetVal = 1; // log enabled

   } else {
      spdlog::set_level(spdlog::level::off);
      RetVal = 0; // log disabled
   }

   return RetVal;
}

// return code: 1->enabled, 0->disabled, negative values->errors
int initLogging(const OMEGA::MachEnv *DefEnv, std::string const &LogFilePath) {

   int RetVal = 0;

   OMEGA::I4 TaskId   = DefEnv->getMyTask();
   OMEGA::I4 NumTasks = DefEnv->getNumTasks();

   std::vector<int> Tasks =
       splitTasks(_OMEGA_TOSTRING(OMEGA_LOG_TASKS), NumTasks);

   if (Tasks.size() > 0 &&
       (std::find(Tasks.begin(), Tasks.end(), TaskId) != Tasks.end())) {

      try {
         std::size_t dotPos = LogFilePath.find_last_of('.');

         if (Tasks.size() > 1 && dotPos != std::string::npos) {
            auto NewLogFilePath = LogFilePath.substr(0, dotPos) + "_" +
                                  std::to_string(TaskId) +
                                  LogFilePath.substr(dotPos);
            spdlog::set_default_logger(
                spdlog::basic_logger_mt("*", NewLogFilePath));
         } else {
            spdlog::set_default_logger(
                spdlog::basic_logger_mt("*", LogFilePath));
         }

         spdlog::set_pattern("[%n %l] %v");
         spdlog::set_level(
             static_cast<spdlog::level::level_enum>(SPDLOG_ACTIVE_LEVEL));
         spdlog::flush_on(spdlog::level::warn);

         RetVal = 1; // log enalbed

      } catch (spdlog::spdlog_ex const &Ex) {
         std::cout << "Log init failed: " << Ex.what() << std::endl;
         RetVal = -1; // error occured
      }

   } else {
      spdlog::set_level(spdlog::level::off);
      RetVal = 0; // log disabled
   }

   return RetVal;
}

} // namespace OMEGA
