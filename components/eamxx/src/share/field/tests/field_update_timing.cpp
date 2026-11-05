#include <catch2/catch.hpp>

#include "share/field/field_identifier.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/core/eamxx_setup_random_test.hpp"

#include <chrono>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <limits>

// Clock the Field update/scale/max/min methods (plain, fill-aware, masked).
// This is NOT a correctness test: it is meant to compare different implementations
// of the update kernels (e.g., compile-time vs runtime CombineMode).
// The problem size can be tuned via env vars:
//   EAMXX_FIELD_TIMING_NCOL  (default 16384)
//   EAMXX_FIELD_TIMING_NLEV  (default 128)
//   EAMXX_FIELD_TIMING_NREPS (default 50)
namespace {

int get_env_int (const char* name, const int def) {
  const char* v = std::getenv(name);
  return v ? std::atoi(v) : def;
}

// Returns best-of-nreps time in microseconds (Field update methods are blocking)
double clock_it (const std::function<void()>& f, const int nreps) {
  using clock = std::chrono::steady_clock;
  f(); // warmup
  double best = std::numeric_limits<double>::max();
  for (int r=0; r<nreps; ++r) {
    auto t0 = clock::now();
    f();
    auto t1 = clock::now();
    best = std::min(best,std::chrono::duration<double,std::micro>(t1-t0).count());
  }
  return best;
}

TEST_CASE ("field_update_timing") {
  using namespace scream;
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  int seed = get_random_test_seed();

  const int ncol  = get_env_int("EAMXX_FIELD_TIMING_NCOL", 16384);
  const int nlev  = get_env_int("EAMXX_FIELD_TIMING_NLEV", 128);
  const int nreps = get_env_int("EAMXX_FIELD_TIMING_NREPS", 50);

  std::vector<FieldTag> tags = {COL,LEV};
  std::vector<int>      dims = {ncol,nlev};

  FieldIdentifier fid_d ("fd", {tags,dims}, kg, "some_grid", DataType::DoubleType);
  FieldIdentifier fid_f ("ff", {tags,dims}, kg, "some_grid", DataType::FloatType);
  FieldIdentifier fid_i ("fi", {tags,dims}, kg, "some_grid", DataType::IntType);

  Field y_d (fid_d), x_d (fid_d), x_f (fid_f), y_i (fid_i), x_i (fid_i);
  for (auto* f : {&y_d,&x_d,&x_f,&y_i,&x_i}) f->allocate_view();

  // Keep values away from 0 (we do divisions) and from the fill value
  randomize_uniform (y_d,seed++,1.0,2.0);
  randomize_uniform (x_d,seed++,1.0,2.0);
  randomize_uniform (x_f,seed++,1.0,2.0);
  randomize_uniform (y_i,seed++,1,100);
  randomize_uniform (x_i,seed++,1,100);

  // Fill-aware RHS (same data, just flagged as possibly containing fill values)
  Field x_d_fill = x_d.clone("x_d_fill");
  x_d_fill.deep_copy(x_d);
  x_d_fill.get_header().set_may_be_filled(true);

  Field mask (FieldIdentifier("mask",{tags,dims},kg,"some_grid",DataType::IntType));
  mask.allocate_view();
  mask.deep_copy(1);
  mask.subfield(COL,0).deep_copy(0);

  struct Case { std::string name; std::function<void()> f; };
  std::vector<Case> cases = {
    // Plain (same data type)
    {"dd update      ", [&]{ y_d.update(x_d,2.0,0.5); }},
    {"dd update(+g)  ", [&]{ y_d.update(x_d,2.0,0.5,1.0); }},
    {"dd scale       ", [&]{ y_d.scale(x_d); }},
    {"dd scale_inv   ", [&]{ y_d.scale_inv(x_d); }},
    {"dd max         ", [&]{ y_d.max(x_d); }},
    {"dd min         ", [&]{ y_d.min(x_d); }},
    {"dd scale(scal) ", [&]{ y_d.scale(1.0); }},
    {"ii update      ", [&]{ y_i.update(x_i,2,1); }},
    // Mixed data types
    {"df update      ", [&]{ y_d.update(x_f,2.0,0.5); }},
    // Fill-aware
    {"dd fill update ", [&]{ y_d.update(x_d_fill,2.0,0.5); }},
    {"dd fill scale  ", [&]{ y_d.scale(x_d_fill); }},
    {"dd fill max    ", [&]{ y_d.max(x_d_fill); }},
    // Masked
    {"dd mask update ", [&]{ y_d.update(x_d,2.0,0.5,mask); }},
    {"dd mask scale  ", [&]{ y_d.scale(x_d,mask); }},
    {"dd mask max    ", [&]{ y_d.max(x_d,mask); }},
  };

  if (comm.am_i_root()) {
    std::cout << "\n=== Field update timing (ncol=" << ncol << ", nlev=" << nlev
              << ", best of " << nreps << ") ===\n";
  }
  for (const auto& c : cases) {
    // Reset y between cases to keep values bounded (not timed)
    y_d.deep_copy(1.5);
    y_i.deep_copy(10);
    const double t = clock_it(c.f,nreps);
    if (comm.am_i_root()) {
      std::cout << "  " << c.name << " : " << std::fixed << std::setprecision(2) << std::setw(10) << t << " us\n";
    }
  }
  REQUIRE (true);
}

} // anonymous namespace
