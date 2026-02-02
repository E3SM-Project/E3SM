//===-- Test driver for OMEGA TimeInterval string parsing ---------*- C++
//-*-===/
//
// This test exercises the TimeInterval(std::string&) constructor used by
// TimeIntegration options like TimeStep and RunDuration.
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

ExpectedParts splitInterval(const TimeInterval &interval) {
   I8 whole{0}, numer{0}, denom{1};
   interval.get(whole, numer, denom);

   // Decompose whole seconds into D/H/M/S (non-negative assumed for this test)
   ExpectedParts parts;
   parts.days         = whole / SECONDS_PER_DAY;
   I8 rem             = whole % SECONDS_PER_DAY;
   parts.hours        = rem / SECONDS_PER_HOUR;
   rem                = rem % SECONDS_PER_HOUR;
   parts.minutes      = rem / SECONDS_PER_MINUTE;
   parts.secondsWhole = rem % SECONDS_PER_MINUTE;

   parts.secondsFrac = static_cast<R8>(numer) / static_cast<R8>(denom);
   return parts;
}

bool nearlyEqual(R8 a, R8 b, R8 tol) { return std::fabs(a - b) <= tol; }

int checkCase(const std::string &label, const std::string &input,
              const ExpectedParts &exp, R8 tol = 1e-12) {

   std::string s = input; // ctor takes non-const ref
   TimeInterval interval(s);

   R8 seconds{0.0};
   interval.get(seconds, TimeUnits::Seconds);

   const R8 expectedSeconds =
       static_cast<R8>(exp.days) * static_cast<R8>(SECONDS_PER_DAY) +
       static_cast<R8>(exp.hours) * static_cast<R8>(SECONDS_PER_HOUR) +
       static_cast<R8>(exp.minutes) * static_cast<R8>(SECONDS_PER_MINUTE) +
       static_cast<R8>(exp.secondsWhole) + exp.secondsFrac;

   auto parts = splitInterval(interval);

   int err = 0;

   if (!nearlyEqual(seconds, expectedSeconds, tol)) {
      ++err;
      LOG_ERROR("{}: '{}' seconds mismatch: got {}, expected {}", label, input,
                seconds, expectedSeconds);
   }

   if (parts.days != exp.days || parts.hours != exp.hours ||
       parts.minutes != exp.minutes || parts.secondsWhole != exp.secondsWhole) {
      ++err;
      LOG_ERROR(
          "{}: '{}' parts mismatch: got {}d {:02d}h {:02d}m {:02d}s, expected "
          "{}d {:02d}h {:02d}m {:02d}s",
          label, input, parts.days, static_cast<int>(parts.hours),
          static_cast<int>(parts.minutes), static_cast<int>(parts.secondsWhole),
          exp.days, static_cast<int>(exp.hours), static_cast<int>(exp.minutes),
          static_cast<int>(exp.secondsWhole));
   }

   if (!nearlyEqual(parts.secondsFrac, exp.secondsFrac, tol)) {
      ++err;
      LOG_ERROR("{}: '{}' fractional seconds mismatch: got {}, expected {}",
                label, input, parts.secondsFrac, exp.secondsFrac);
   }

   return err;
}

} // namespace

int main(int argc, char **argv) {

   MPI_Init(&argc, &argv);

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *defEnv = MachEnv::getDefault();
   initLogging(defEnv);

   int errCount = 0;

   // Currently documented/implemented format is:
   //   DDDD_HH:MM:SS(.sss...) with arbitrary day/second widths and any
   //   single-character separators.

   errCount += checkCase("TimeInterval parse", "0000_01:23:45",
                         ExpectedParts{0, 1, 23, 45, 0.0});

   errCount += checkCase("TimeInterval parse", "0000_01:23:45.678",
                         ExpectedParts{0, 1, 23, 45, 0.678});

   errCount += checkCase("TimeInterval parse", "0000_01:23:45.6789",
                         ExpectedParts{0, 1, 23, 45, 0.6789});

   // Mixed separators are allowed by the implementation.
   errCount += checkCase("TimeInterval parse", "0012-01.23/45.5",
                         ExpectedParts{12, 1, 23, 45, 0.5});

   if (errCount != 0) {
      LOG_ERROR("TimeIntervalParseTest: FAIL ({} errors)", errCount);
      MPI_Finalize();
      return 1;
   }

   LOG_INFO("TimeIntervalParseTest: PASS");

   MPI_Finalize();
   return 0;
}
