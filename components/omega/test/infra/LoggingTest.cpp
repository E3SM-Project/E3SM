//===-- Test driver for OMEGA logging -------------------------*- C++ -*-===/
//
/// \file
/// \brief Test driver for OMEGA logging
///
/// This driver tests the logging capabilities for the OMEGA
/// model. In particular, it tests creating a log file according to
/// log levels and supporting Kokkos data types.
///
//
//===-----------------------------------------------------------------------===/

#include <iostream>

#define _OMEGA_STRINGIFY(x) #x
#define _OMEGA_TOSTRING(x)  _OMEGA_STRINGIFY(x)

#include "Logging.h"

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/ringbuffer_sink.h"

using namespace OMEGA;

enum CheckType { EndsWith, StartsWith, Contains };

auto TestSink = std::make_shared<spdlog::sinks::ringbuffer_sink_mt>(1);

bool hasEnding(std::string const &fullString, std::string const &ending) {

   if (fullString.length() < ending.length()) {
      return false;
   } else {
      return (0 == fullString.compare(fullString.length() - ending.length(),
                                      ending.length(), ending));
   }
}

bool hasSubstring(std::string const &fullString, std::string const &subString) {
   return (fullString.find(subString) != std::string::npos);
}

int outputTestResult(std::string const &TestName, std::string const &Expected,
                     CheckType Type) {

   int RetVal = 0;

   std::vector<std::string> Msgs = TestSink->last_formatted();

   if (Expected.length() == 0) {
      if (Msgs.size() == 0) {
         std::cout << TestName << ": PASS" << std::endl;
      } else {
         std::cout << TestName << ": FAIL" << std::endl;
         RetVal = 1;
      }
   } else {
      if (Msgs.size() == 0) {
         std::cout << TestName << ": FAIL" << std::endl;
         RetVal = 1;
      } else {
         std::string Actual = Msgs[0];
         bool pf            = false;

         if (Type == EndsWith) {
            std::string NewExpected(Expected + "\n");
            pf = hasEnding(Actual, NewExpected);

         } else if (Type == Contains) {
            pf = hasSubstring(Actual, Expected);
         }

         if (pf) {
            std::cout << TestName << ": PASS" << std::endl;
         } else {
            std::cout << TestName << ": FAIL" << std::endl;
            RetVal = 1;
         }
      }
   }

   return RetVal;
}

int testDefaultLogLevel(bool LogEnabled) {

   int RetVal = 0;

   // NOTE: default log level is INFO
   LOG_DEBUG("This shouldn't be logged.");

   // check if no debug log
   if (LogEnabled)
      RetVal +=
          outputTestResult("Default log level 1", std::string(""), EndsWith);

   const std::string TMSG2 = "This should be logged.";
   LOG_INFO(TMSG2);

   // check if this is info log
   if (LogEnabled)
      RetVal += outputTestResult("Default log level 2", TMSG2, EndsWith);

   return RetVal;
}

int testKokkosDataTypes(bool LogEnabled) {

   int RetVal       = 0;
   int constexpr d1 = 2;
   int constexpr d2 = 3;

   Kokkos::initialize();
   {
      HostArray1DReal test1d("test1d", d1);
      HostArray2DReal test2d("test2d", d1, d2);

      LOG_INFO("1d var {}", test1d);

      // check if HostArray1DReal is detected
      if (LogEnabled)
         RetVal += outputTestResult("Kokkos data type 1", "HostArray1DReal",
                                    Contains);

      LOG_INFO("2d var {}", test2d);

      // check if HostArray2DReal is detected
      if (LogEnabled)
         RetVal += outputTestResult("Kokkos data type 2", "HostArray2DReal",
                                    Contains);
   }
   Kokkos::finalize();

   return RetVal;
}

int main(int argc, char **argv) {

   int RetVal = 0;

   MPI_Init(&argc, &argv);

   OMEGA::MachEnv::init(MPI_COMM_WORLD);
   OMEGA::MachEnv *DefEnv = OMEGA::MachEnv::getDefault();
   OMEGA::I4 TaskId       = DefEnv->getMyTask();

   std::string TasksStr = _OMEGA_TOSTRING(OMEGA_LOG_TASKS);

   if (TasksStr.find(std::to_string(TaskId)) != std::string::npos) {
      try {

         std::string LogFilePath = "tmplog_" + std::to_string(TaskId) + ".log";
         std::remove(LogFilePath.c_str());

         // "sinks" hold pointers to "spdlog" sinks.
         std::vector<spdlog::sink_ptr> sinks;
         // adds the first sink of basic file sink mt
         sinks.push_back(
             std::make_shared<spdlog::sinks::basic_file_sink_mt>(LogFilePath));
         // adds the second sink for this unit testing
         sinks.push_back(TestSink);

         // creates a logger that sends log messages to multiple sinks
         auto Logger = std::make_shared<spdlog::logger>(
             "unit", std::begin(sinks), std::end(sinks));

         // initialize Omega logging with the logger
         auto LogEnabled = initLogging(Logger) == 1;

         spdlog::set_pattern("[%n %l] %v");
         spdlog::set_level(
             static_cast<spdlog::level::level_enum>(SPDLOG_ACTIVE_LEVEL));
         spdlog::flush_on(spdlog::level::warn);

         RetVal += testDefaultLogLevel(LogEnabled);
         RetVal += testKokkosDataTypes(LogEnabled);

      } catch (const std::exception &Ex) {
         std::cout << Ex.what() << ": FAIL" << std::endl;
         RetVal += 1;
      } catch (...) {
         std::cout << "Unknown: FAIL" << std::endl;
         RetVal += 1;
      }

   } else {
      spdlog::set_level(spdlog::level::off);
      RetVal = 0;
   }

   // Finalize environments
   MPI_Finalize();

   return RetVal;
}
