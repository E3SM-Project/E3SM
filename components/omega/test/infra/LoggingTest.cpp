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

// #include "DataTypes.h"
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
         RetVal = -1;
      }
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
         RetVal = -1;
      }
   }

   return RetVal;
}

int testDefaultLogLevel() {

   int RetVal = 0;

   // NOTE: default log level is INFO
   LOG_DEBUG("This shouldn't be logged.");

   // check if no debug log
   RetVal -= outputTestResult("Default log level 1", std::string(""), EndsWith);

   const std::string TMSG2 = "This should be logged.";
   LOG_INFO(TMSG2);

   // check if this is info log
   RetVal -= outputTestResult("Default log level 2", TMSG2, EndsWith);

   return RetVal;
}

int testKokkosDataTypes() {

   int RetVal       = 0;
   int constexpr d1 = 2;
   int constexpr d2 = 3;

   Kokkos::initialize();
   {
      HostArray1DReal test1d("test1d", d1);
      HostArray2DReal test2d("test2d", d1, d2);

      LOG_INFO("1d var {}", test1d);

      // check if HostArray1DReal is detected
      RetVal -=
          outputTestResult("Kokkos data type 1", "HostArray1DReal", Contains);

      LOG_INFO("2d var {}", test2d);

      // check if HostArray2DReal is detected
      RetVal -=
          outputTestResult("Kokkos data type 2", "HostArray2DReal", Contains);
   }
   Kokkos::finalize();

   return RetVal;
}

int main(int argc, char **argv) {

   int RetVal                    = 0;
   const std::string LogFilePath = "tmplog.log";

   try {

      std::remove(LogFilePath.c_str());

      std::vector<spdlog::sink_ptr> sinks;
      sinks.push_back(
          std::make_shared<spdlog::sinks::basic_file_sink_mt>(LogFilePath));
      sinks.push_back(TestSink);

      auto logger = std::make_shared<spdlog::logger>("unit", std::begin(sinks),
                                                     std::end(sinks));

      initLogging(logger);

      RetVal -= testDefaultLogLevel();
      RetVal -= testKokkosDataTypes();

      // std::remove(LogFilePath.c_str());

   } catch (const std::exception &Ex) {
      std::cout << Ex.what() << ": FAIL" << std::endl;
      RetVal -= -1;
   } catch (...) {
      std::cout << "Unknown: FAIL" << std::endl;
      RetVal -= -1;
   }

   return RetVal;
}
