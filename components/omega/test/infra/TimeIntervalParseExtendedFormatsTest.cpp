//===-- Test driver for OMEGA TimeInterval extended string formats -*- C++
//-*-===/
//
// This test encodes the *desired* behavior for parsing time interval strings
// from config (e.g., TimeIntegration::TimeStep and RunDuration).
//
// It intentionally checks "short" forms like HH:MM:SS(.sss), MM:SS(.sss), and
// SS(.sss). Today, the implementation assumes the string always begins with
// days (DDDD_HH:MM:SS(.sss)), so these cases reproduce the bug.
//
// This test is expected to fail until parsing is improved.
//
//===----------------------------------------------------------------------===/

#include "DataTypes.h"
#include "Error.h"
#include "Logging.h"
#include "MachEnv.h"
#include "TimeMgr.h"
#include "mpi.h"

#include <cmath>
#include <string>

using namespace OMEGA;

namespace {

struct ExpectedParts {
   I8 days{0};
   I8 hours{0};
   I8 minutes{0};
   I8 secondsWhole{0};
   R8 secondsFrac{0.0};
};

bool nearlyEqual(R8 a, R8 b, R8 tol) { return std::fabs(a - b) <= tol; }

int checkSeconds(const std::string &label, const std::string &input,
                 const ExpectedParts &exp, R8 tol = 1e-12) {

   const R8 expectedSeconds =
       static_cast<R8>(exp.days) * static_cast<R8>(SECONDS_PER_DAY) +
       static_cast<R8>(exp.hours) * static_cast<R8>(SECONDS_PER_HOUR) +
       static_cast<R8>(exp.minutes) * static_cast<R8>(SECONDS_PER_MINUTE) +
       static_cast<R8>(exp.secondsWhole) + exp.secondsFrac;

   std::string s = input; // ctor takes non-const ref
   TimeInterval interval(s);

   R8 seconds{0.0};
   interval.get(seconds, TimeUnits::Seconds);

   if (!nearlyEqual(seconds, expectedSeconds, tol)) {
      LOG_ERROR("{}: '{}' parsed seconds mismatch: got {}, expected {}", label,
                input, seconds, expectedSeconds);
      return 1;
   }

   return 0;
}

} // namespace

int main(int argc, char **argv) {

   MPI_Init(&argc, &argv);

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *defEnv = MachEnv::getDefault();
   initLogging(defEnv);

   int errCount = 0;

   // Desired: HH:MM:SS(.sss)
   errCount += checkSeconds("Extended TimeInterval parse", "01:23:45",
                            ExpectedParts{0, 1, 23, 45, 0.0});
   errCount += checkSeconds("Extended TimeInterval parse", "01:23:45.678",
                            ExpectedParts{0, 1, 23, 45, 0.678});

   // Desired: MM:SS(.sss)
   errCount += checkSeconds("Extended TimeInterval parse", "23:45.678",
                            ExpectedParts{0, 0, 23, 45, 0.678});

   // Desired: SS(.sss)
   errCount += checkSeconds("Extended TimeInterval parse", "45.678",
                            ExpectedParts{0, 0, 0, 45, 0.678});
   errCount += checkSeconds("Extended TimeInterval parse", "45.6",
                            ExpectedParts{0, 0, 0, 45, 0.6});
   errCount += checkSeconds("Extended TimeInterval parse", "45.000001",
                            ExpectedParts{0, 0, 0, 45, 0.000001});

   if (errCount != 0) {
      LOG_ERROR("TimeIntervalParseExtendedFormatsTest: FAIL ({} errors)",
                errCount);
      MPI_Finalize();
      return 1;
   }

   LOG_INFO("TimeIntervalParseExtendedFormatsTest: PASS");

   MPI_Finalize();
   return 0;
}
