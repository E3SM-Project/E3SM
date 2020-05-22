#include <catch2/catch.hpp>

#include "ekat/util/scream_kokkos_utils.hpp"
#include "ekat/util/scream_arch.hpp"

#include <thread>

namespace {

TEST_CASE("data_type", "[kokkos_utils]") {
  using namespace scream::util;

  // Check meta-util to get underlying scalar type of a raw MD array
  REQUIRE(std::is_same<ValueType<double**&>::type,double>::value);
  REQUIRE(std::is_same<ValueType<double*[3]>::type,double>::value);
  REQUIRE(std::is_same<ValueType<double[2][3]>::type,double>::value);

  // Check meta-util to get rank and dynamic rank of a raw MD array
  REQUIRE(GetRanks<double[2][3]>::rank==2);
  REQUIRE(GetRanks<double[2][3]>::rank_dynamic==0);
  REQUIRE(GetRanks<double*[2][3]>::rank==3);
  REQUIRE(GetRanks<double*[2][3]>::rank_dynamic==1);
  REQUIRE(GetRanks<double**[2][3]>::rank==4);
  REQUIRE(GetRanks<double**[2][3]>::rank_dynamic==2);

  // Check meta-util that allows to reshape a view
  Kokkos::View<double*> v1d("",100);
  auto v2d = reshape<double[2][50]>(v1d);
  REQUIRE(v2d.size()==100);

  auto v3d = reshape<double*[5][5]>(v2d,4);
  REQUIRE (v3d.size()==100);
}

TEST_CASE("team_policy", "[kokkos_utils]") {
  using namespace scream::util;
  using namespace scream;

  using Device = DefaultDevice;
  using ExeSpace = typename KokkosTypes<Device>::ExeSpace;

  for (int nk: {128, 122, 255, 42}) {
    const int ni = 1000;
    const auto p = ExeSpaceUtils<ExeSpace>::get_default_team_policy(ni, nk);
    REQUIRE(p.league_size() == ni);
    if (OnGpu<ExeSpace>::value) {
      if (nk == 42) {
        REQUIRE(p.team_size() == 64);
      } else {
        REQUIRE(p.team_size() == 128);
      }
    }
    else {
#if defined SCREAM_MIMIC_GPU && defined KOKKOS_ENABLE_OPENMP
      REQUIRE((Kokkos::OpenMP::concurrency() == 1 || p.team_size() > 1));
#endif
    }
  }
}

TEST_CASE("team_utils_omp", "[kokkos_utils]")
{
#ifdef KOKKOS_ENABLE_OPENMP
  using namespace scream::util;
  using namespace scream;

  using Device = DefaultDevice;
  using ExeSpace = typename KokkosTypes<Device>::ExeSpace;
  using MemberType = typename KokkosTypes<Device>::MemberType;
  const int n = omp_get_max_threads();
  // test will not work with more than 16 threads
  if (n > 16) {
    WARN("Skipped because this test doesn't support more than 16 threads");
    return;
  }

  const int ni = n*5;
  for (int s = 1; s <= n; ++s) {
    const auto p = ExeSpaceUtils<ExeSpace>::get_team_policy_force_team_size(ni, s);
    TeamUtils<ExeSpace> tu(p);
    const int c = tu.get_num_concurrent_teams();
    typename KokkosTypes<Device>::template view_2d<int> ws_idxs("ws_idxs", ni, s);
#if 0
    const int real_ts = omp_get_max_threads() / c;
    std::cout << "thrds " << n << " teamsizeV " << s << " teamsizeR " << real_ts << " ni " << ni << " conc " << c <<  std::endl;
#endif
    int kernel_errors = 0;
    Kokkos::parallel_reduce("unittest_team_utils", p, KOKKOS_LAMBDA(MemberType team_member, int& total_errs) {
      int nerrs_local = 0;
      const int i  = team_member.league_rank();
      const int wi = tu.get_workspace_idx(team_member);

#if 0
      const int thread_num = omp_get_thread_num();
      for (int j = 0; j < n; ++j) {
        if (j == thread_num) {
          if (j == 0) {
            std::cout << "===================================" << std::endl;
          }
          std::cout << " For total_threads: " << n << " league size " << team_member.league_size() << " and team_size: " << s << ", team: " << i << ", team_rank=" << team_member.team_rank() << ", thread: " << thread_num << " , conc: " << c << ", idx: " << wi << std::endl;
        }
        std::this_thread::sleep_for(std::chrono::milliseconds(100));
        team_member.team_barrier();
      }
#endif
      ws_idxs(i, team_member.team_rank()) = wi+1; // never zero
      if (wi >= c) ++nerrs_local;

      total_errs += nerrs_local;
    }, kernel_errors);
#if 0
    std::cout << "===================== DONE ==========================" << std::endl;
#endif
    REQUIRE(kernel_errors == 0);

    // post processing
    const int teams_per_idx = (ni + c - 1) / c;
    for (int i = 0; i < ni; ++i) {
      int exp_wi = 0;
      // all threads in a team should share wsidx
      for (int t = 0; t < s; ++t) {
        int curr_wi = ws_idxs(i, t);
#if 0
        std::cout << "idxs(" << i << ", " << t << ") = " << curr_wi << std::endl;
#endif
        if (t == 0) exp_wi = curr_wi;
        REQUIRE(curr_wi != 0);
        REQUIRE(curr_wi == exp_wi);
      }
    }

    // Check that each wsidx used correct number of times
    for (int ci = 1; ci <= c; ++ci) {
      for (int t = 0; t < s; ++t) {
        int cnt = 0;
        for (int i = 0; i < ni; ++i) {
          if (ws_idxs(i,t) == ci) ++cnt;
        }
        REQUIRE(cnt <= teams_per_idx);
      }
    }
  }
#endif
}

void test_utils_large_ni(const double saturation_multiplier)
{
  using namespace scream::util;
  using namespace scream;

  using Device = DefaultDevice;
  using ExeSpace = typename KokkosTypes<Device>::ExeSpace;
  using MemberType = typename KokkosTypes<Device>::MemberType;

  const int nk = 128;
  const Real overprov_factor = 1.5;
  const auto temp_policy = ExeSpaceUtils<ExeSpace>::get_default_team_policy(1, nk);
  TeamUtils<ExeSpace> tu_temp(temp_policy);
  const int num_conc = tu_temp.get_max_concurrent_threads() / temp_policy.team_size();

  int ni = num_conc*saturation_multiplier;
  if (ni == 0) ni = 1;
  const auto p = ExeSpaceUtils<ExeSpace>::get_default_team_policy(ni, nk);
  TeamUtils<ExeSpace> tu(p, overprov_factor);

  REQUIRE(p.league_size() == ni);
  if (saturation_multiplier <= 1.0) {
    REQUIRE(tu.get_num_ws_slots() == ni);
  }
  else if (!OnGpu<ExeSpace>::value) {
    REQUIRE(tu.get_num_ws_slots() == num_conc);
  }
  else {
    REQUIRE(tu.get_num_ws_slots() == num_conc*overprov_factor);
  }

  int max_workspace_idx = 0;
  typename KokkosTypes<Device>::template view_1d<int> test_data("test_data", tu.get_num_ws_slots());
  Kokkos::parallel_reduce("unique_token_check", p, KOKKOS_LAMBDA(MemberType team_member, int& max_ws_idx) {
    const int wi = tu.get_workspace_idx(team_member);

    if (wi > max_ws_idx) { max_ws_idx = wi; }

    Kokkos::single(Kokkos::PerTeam(team_member), [&] () {
      int volatile* const data = &test_data(wi);
      *data += 1;
    });

    tu.release_workspace_idx(team_member, wi);
  }, Kokkos::Max<int>(max_workspace_idx));

  const auto test_data_h = Kokkos::create_mirror_view(test_data);
  Kokkos::deep_copy(test_data_h, test_data);

  int sum = 0;
  for(int i = 0; i < tu.get_num_ws_slots(); ++i) {
    sum += test_data_h(i);
  }

  REQUIRE(sum == ni);
}

TEST_CASE("team_utils_large_ni", "[kokkos_utils]")
{
  test_utils_large_ni(10);
  test_utils_large_ni(1);
  test_utils_large_ni(.5);
}

} // anonymous namespace
