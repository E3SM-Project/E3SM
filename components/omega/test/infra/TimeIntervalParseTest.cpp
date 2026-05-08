//===-- Test driver for OMEGA TimeInterval string parsing   ------*- C++ -*-===/
//
// This test exercises the TimeInterval(std::string&) constructor used by
// TimeIntegration options like TimeStep and RunDuration.
//
//===-----------------------------------------------------------------------===/

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
   I8 Days{0};
   I8 Hours{0};
   I8 Minutes{0};
   I8 SecondsWhole{0};
   R8 SecondsFrac{0.0};
};

ExpectedParts splitInterval(const TimeInterval &Interval) {
   I8 Whole{0}, Numer{0}, Denom{1};
   Interval.get(Whole, Numer, Denom);

   // Decompose whole seconds into D/H/M/S (non-negative assumed for this test)
   ExpectedParts Parts;
   Parts.Days         = Whole / SECONDS_PER_DAY;
   I8 Rem             = Whole % SECONDS_PER_DAY;
   Parts.Hours        = Rem / SECONDS_PER_HOUR;
   Rem                = Rem % SECONDS_PER_HOUR;
   Parts.Minutes      = Rem / SECONDS_PER_MINUTE;
   Parts.SecondsWhole = Rem % SECONDS_PER_MINUTE;

   Parts.SecondsFrac = static_cast<R8>(Numer) / static_cast<R8>(Denom);
   return Parts;
}

bool nearlyEqual(R8 A, R8 B, R8 Tol) { return std::fabs(A - B) <= Tol; }

void checkCase(const std::string &Label, const std::string &Input,
               const ExpectedParts &Exp, R8 Tol = 1e-12) {

   std::string S = Input; // ctor takes non-const ref
   TimeInterval Interval(S);

   R8 Seconds{0.0};
   Interval.get(Seconds, TimeUnits::Seconds);

   const R8 ExpectedSeconds =
       static_cast<R8>(Exp.Days) * static_cast<R8>(SECONDS_PER_DAY) +
       static_cast<R8>(Exp.Hours) * static_cast<R8>(SECONDS_PER_HOUR) +
       static_cast<R8>(Exp.Minutes) * static_cast<R8>(SECONDS_PER_MINUTE) +
       static_cast<R8>(Exp.SecondsWhole) + Exp.SecondsFrac;

   auto Parts = splitInterval(Interval);

   if (!nearlyEqual(Seconds, ExpectedSeconds, Tol)) {
      ABORT_ERROR("{}: '{}' seconds mismatch: got {}, expected {}", Label,
                  Input, Seconds, ExpectedSeconds);
   }

   if (Parts.Days != Exp.Days || Parts.Hours != Exp.Hours ||
       Parts.Minutes != Exp.Minutes || Parts.SecondsWhole != Exp.SecondsWhole) {
      ABORT_ERROR(
          "{}: '{}' parts mismatch: got {}d {:02d}h {:02d}m {:02d}s, expected "
          "{}d {:02d}h {:02d}m {:02d}s",
          Label, Input, Parts.Days, static_cast<int>(Parts.Hours),
          static_cast<int>(Parts.Minutes), static_cast<int>(Parts.SecondsWhole),
          Exp.Days, static_cast<int>(Exp.Hours), static_cast<int>(Exp.Minutes),
          static_cast<int>(Exp.SecondsWhole));
   }

   if (!nearlyEqual(Parts.SecondsFrac, Exp.SecondsFrac, Tol)) {
      ABORT_ERROR("{}: '{}' fractional seconds mismatch: got {}, expected {}",
                  Label, Input, Parts.SecondsFrac, Exp.SecondsFrac);
   }

   return;
}

} // namespace

int main(int argc, char **argv) {

   MPI_Init(&argc, &argv);

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *defEnv = MachEnv::getDefault();
   initLogging(defEnv);

   LOG_INFO("----- TimeIntervalParseTest -----");

   // Currently documented/implemented format is:
   //   DDDD_HH:MM:SS(.sss...) with arbitrary day/second widths and any
   //   single-character separators.

   checkCase("TimeInterval parse", "0000_01:23:45",
             ExpectedParts{0, 1, 23, 45, 0.0});

   checkCase("TimeInterval parse", "0000_01:23:45.678",
             ExpectedParts{0, 1, 23, 45, 0.678});

   checkCase("TimeInterval parse", "0000_01:23:45.6789",
             ExpectedParts{0, 1, 23, 45, 0.6789});

   // Mixed separators are allowed by the implementation.
   checkCase("TimeInterval parse", "0012-01.23/45.5",
             ExpectedParts{12, 1, 23, 45, 0.5});

   // If we made it here, we are successful
   LOG_INFO("----- TimeIntervalParseTest Successful -----");

   MPI_Finalize();
   return 0;
}
