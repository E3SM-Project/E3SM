//===-- Test driver for OMEGA TimeInterval extended string formats -*- C++ -*-=/
//
// This test encodes the *desired* behavior for parsing time interval strings
// from config (e.g., TimeIntegration::TimeStep and RunDuration).
//
// It intentionally checks "short" forms like HH:MM:SS(.sss), MM:SS(.sss), and
// SS(.sss).
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

bool nearlyEqual(R8 A, R8 B, R8 Tol) { return std::fabs(A - B) <= Tol; }

void checkSeconds(const std::string &Label, const std::string &Input,
                  const ExpectedParts &Exp, R8 Tol = 1e-12) {

   const R8 ExpectedSeconds =
       static_cast<R8>(Exp.Days) * static_cast<R8>(SECONDS_PER_DAY) +
       static_cast<R8>(Exp.Hours) * static_cast<R8>(SECONDS_PER_HOUR) +
       static_cast<R8>(Exp.Minutes) * static_cast<R8>(SECONDS_PER_MINUTE) +
       static_cast<R8>(Exp.SecondsWhole) + Exp.SecondsFrac;

   std::string S = Input; // ctor takes non-const ref
   TimeInterval Interval(S);

   R8 Seconds{0.0};
   Interval.get(Seconds, TimeUnits::Seconds);

   if (!nearlyEqual(Seconds, ExpectedSeconds, Tol))
      ABORT_ERROR("{}: '{}' parsed seconds mismatch: got {}, expected {}",
                  Label, Input, Seconds, ExpectedSeconds);

   return;
}

} // namespace

int main(int argc, char **argv) {

   MPI_Init(&argc, &argv);

   MachEnv::init(MPI_COMM_WORLD);
   MachEnv *defEnv = MachEnv::getDefault();
   initLogging(defEnv);

   LOG_INFO("----- TimeIntervalParseExtendedFormatsTest -----");

   // Desired: HH:MM:SS(.sss)
   checkSeconds("Extended TimeInterval parse", "01:23:45",
                ExpectedParts{0, 1, 23, 45, 0.0});
   checkSeconds("Extended TimeInterval parse", "01:23:45.678",
                ExpectedParts{0, 1, 23, 45, 0.678});

   // Desired: MM:SS(.sss)
   checkSeconds("Extended TimeInterval parse", "23:45.678",
                ExpectedParts{0, 0, 23, 45, 0.678});

   // Desired: SS(.sss)
   checkSeconds("Extended TimeInterval parse", "45.678",
                ExpectedParts{0, 0, 0, 45, 0.678});
   checkSeconds("Extended TimeInterval parse", "45.6",
                ExpectedParts{0, 0, 0, 45, 0.6});
   checkSeconds("Extended TimeInterval parse", "45.000001",
                ExpectedParts{0, 0, 0, 45, 0.000001});

   // If we made it here, the tests were successful
   LOG_INFO("----- TimeIntervalParseExtendedFormatsTest Successful -----");

   MPI_Finalize();
   return 0;
}
