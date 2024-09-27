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

#include "Logging.h"
#include "MachEnv.h"

#include "spdlog/sinks/basic_file_sink.h"
#include "spdlog/sinks/ringbuffer_sink.h"

#define _OMEGA_STRINGIFY(x) #x
#define _OMEGA_TOSTRING(x)  _OMEGA_STRINGIFY(x)

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

   int RetVal              = 0;
   bool TestRun            = true;
   const std::string TMSG1 = "This shouldn't be logged.";
   const std::string TMSG2 = "This should be logged.";

   if (OMEGA_LOG_LEVEL == 1) {
      LOG_TRACE(TMSG1);
   } else if (OMEGA_LOG_LEVEL == 2) {
      LOG_DEBUG(TMSG1);
   } else if (OMEGA_LOG_LEVEL == 3) {
      LOG_INFO(TMSG1);
   } else if (OMEGA_LOG_LEVEL == 4) {
      LOG_WARN(TMSG1);
   } else if (OMEGA_LOG_LEVEL == 5) {
      LOG_ERROR(TMSG1);
   } else {
      TestRun = false; // trace or off
   }

   if (LogEnabled && TestRun)
      RetVal +=
          outputTestResult("Default log level 1", std::string(""), EndsWith);

   TestRun = true;

   if (OMEGA_LOG_LEVEL == 0) {
      LOG_TRACE(TMSG2);
   } else if (OMEGA_LOG_LEVEL == 1) {
      LOG_DEBUG(TMSG2);
   } else if (OMEGA_LOG_LEVEL == 2) {
      LOG_INFO(TMSG2);
   } else if (OMEGA_LOG_LEVEL == 3) {
      LOG_WARN(TMSG2);
   } else if (OMEGA_LOG_LEVEL == 4) {
      LOG_ERROR(TMSG2);
   } else if (OMEGA_LOG_LEVEL == 5) {
      LOG_CRITICAL(TMSG2);
   } else {
      TestRun = false; // off
   }

   if (LogEnabled && TestRun)
      RetVal += outputTestResult("Default log level 2", TMSG2, EndsWith);

   return RetVal;
}

int testKokkosDataTypes(bool LogEnabled) {

   int RetVal       = 0;
   int constexpr d1 = 2;
   int constexpr d2 = 3;
   bool TestRun;

   Kokkos::initialize();
   {
      HostArray1DReal test1d("test1d", d1);
      HostArray2DReal test2d("test2d", d1, d2);

      TestRun = true;

      if (OMEGA_LOG_LEVEL == 0) {
         LOG_INFO("1d var {}", test1d);
      } else if (OMEGA_LOG_LEVEL == 1) {
         LOG_DEBUG("1d var {}", test1d);
      } else if (OMEGA_LOG_LEVEL == 2) {
         LOG_INFO("1d var {}", test1d);
      } else if (OMEGA_LOG_LEVEL == 3) {
         LOG_WARN("1d var {}", test1d);
      } else if (OMEGA_LOG_LEVEL == 4) {
         LOG_ERROR("1d var {}", test1d);
      } else if (OMEGA_LOG_LEVEL == 5) {
         LOG_CRITICAL("1d var {}", test1d);
      } else {
         TestRun = false; // off
      }

      // check if HostArray1DReal is detected
      if (LogEnabled && TestRun)
         RetVal += outputTestResult("Kokkos data type 1", "test1d", Contains);

      TestRun = true;

      if (OMEGA_LOG_LEVEL == 0) {
         LOG_INFO("2d var {}", test2d);
      } else if (OMEGA_LOG_LEVEL == 1) {
         LOG_DEBUG("2d var {}", test2d);
      } else if (OMEGA_LOG_LEVEL == 2) {
         LOG_INFO("2d var {}", test2d);
      } else if (OMEGA_LOG_LEVEL == 3) {
         LOG_WARN("2d var {}", test2d);
      } else if (OMEGA_LOG_LEVEL == 4) {
         LOG_ERROR("2d var {}", test2d);
      } else if (OMEGA_LOG_LEVEL == 5) {
         LOG_CRITICAL("2d var {}", test2d);
      } else {
         TestRun = false; // off
      }

      // check if HostArray2DReal is detected
      if (LogEnabled && TestRun)
         RetVal += outputTestResult("Kokkos data type 2", "test2d", Contains);
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
      auto Logger = std::make_shared<spdlog::logger>("unit", std::begin(sinks),
                                                     std::end(sinks));

      // initialize Omega logging with the logger
      bool LogEnabled = (initLogging(DefEnv, Logger) == 1);

      RetVal += testDefaultLogLevel(LogEnabled);
      RetVal += testKokkosDataTypes(LogEnabled);

   } catch (const std::exception &Ex) {
      std::cout << Ex.what() << ": FAIL" << std::endl;
      RetVal += 1;
   } catch (...) {
      std::cout << "Unknown: FAIL" << std::endl;
      RetVal += 1;
   }

   // Finalize environments
   MPI_Finalize();

   return RetVal;
}
