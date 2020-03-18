#include <catch2/catch.hpp>

#include "share/scream_assert.hpp"
#include <csignal>
#include <unistd.h>
#include  <csetjmp>
#include <iostream>

namespace {

jmp_buf JumpBuffer;

static volatile std::sig_atomic_t gSignalStatus = 0;

void signal_handler (int /* signum */) {
  gSignalStatus = 1;
  std::longjmp(JumpBuffer,gSignalStatus);
}

void run_fpe_tests () {
  int mask = scream::get_enabled_fpes();

  std::cout << "tests mask: " << mask << "\n";
  // Get the current fenv.
  // Note: for some reason, each time SIGFPE is thrown or getenv
  //       is called, the fenv is reset to 0. So remember to
  //       set the fenv before each test.
  std::fenv_t fenv;
  feholdexcept(&fenv);
  fesetenv(&fenv);

  const auto has_fe_divbyzero = [](int mask)->bool{
    return mask & FE_DIVBYZERO;
  };
  const auto has_fe_overflow = [](int mask)->bool{
    return mask & FE_OVERFLOW;
  };
  const auto has_fe_invalid = [](int mask)->bool{
    return mask & FE_INVALID;
  };

  std::cout << " has FE_DIVBYZERO: " << (has_fe_divbyzero(mask) ? "yes" : "no") << "\n";
  std::cout << " has FE_INVALID: " << (has_fe_invalid(mask) ? "yes" : "no") << "\n";
  std::cout << " has FE_OVERFLOW: " << (has_fe_overflow(mask) ? "yes" : "no") << "\n";

  double one = 1.0;
  double zero = 0.0;
  double inf, nan, ovfl;

  // Run the tests.
  // Note: sometimes a FPE is not thrown when the bad number
  //       is generated, but rather the next time is used.
  //       Therefore, sometimes we multiply the result
  //       by 1.0 before testing that the FPE was thrown.
  
  // Test 1/0
  setjmp(JumpBuffer);
  inf = one/zero;
  REQUIRE (gSignalStatus==(has_fe_divbyzero(mask) ? 1 : 0));
  gSignalStatus = 0;
  fesetenv(&fenv);

  // Test 0/0
  setjmp(JumpBuffer);
  nan = zero/zero;
  REQUIRE (gSignalStatus==(has_fe_invalid(mask) ? 1 : 0));
  gSignalStatus = 0;
  fesetenv(&fenv);

  // Test invalid arg
  setjmp(JumpBuffer);
  nan = std::sqrt(-1.0);
  nan *= 1.0;
  REQUIRE (gSignalStatus==(has_fe_invalid(mask) ? 1 : 0));
  gSignalStatus = 0;
  fesetenv(&fenv);

  // Test overflow
  setjmp(JumpBuffer);
  ovfl = exp(710.0);
  ovfl *= 1.0;
  REQUIRE (gSignalStatus==(has_fe_overflow(mask) ? 1 : 0));
  gSignalStatus = 0;

  (void) inf;
  (void) nan;
  (void) ovfl;
}

TEST_CASE ("fpes","") {
  using namespace scream;

  // Set a handler, which simply sets the global gSignalStatus,
  // so we can check with catch whether the FPE was raised.
  struct sigaction sa;
  sa.sa_handler = signal_handler;
  sa.sa_flags = SA_NODEFER;
  sigemptyset(&sa.sa_mask);
  sigaction(SIGFPE, &sa, NULL);

  SECTION ("default_fpes") {
    run_fpe_tests ();
  }

  SECTION ("user-requested fpe") {
    disable_all_fpes();
    enable_fpes(FE_DIVBYZERO);
    run_fpe_tests ();
  }

  SECTION ("user-requested fpe") {
    enable_fpes(FE_ALL_EXCEPT);
    disable_fpes(FE_DIVBYZERO);
    run_fpe_tests ();
  }
}

TEST_CASE ("asserts") {
  auto test_req_msg = [](const bool test, const std::string& msg) {
    scream_require_msg(test,msg);
  };
  auto test_err_msg = [](const std::string& msg) {
    scream_error_msg(msg);
  };
  REQUIRE_THROWS (test_req_msg(1>3,"Uh? I wonder what Sharkowsky would have to say about this...\n"));

  REQUIRE_NOTHROW (test_req_msg(3>1,"Uh? I wonder what Sharkowsky would have to say about this...\n"));

  REQUIRE_THROWS (test_err_msg("Hello world!\n"));
}

} // anonymous namespace
