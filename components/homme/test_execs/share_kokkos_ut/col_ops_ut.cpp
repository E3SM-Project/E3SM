#include <catch2/catch.hpp>

#include "ColumnOps.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

#include "Types.hpp"

#include <random>
#include <iostream>

using namespace Homme;


// ============= COLUMN OPS ================ //
constexpr int NUM_PTS = 4;
constexpr int num_elems = 10;

TEST_CASE("col_ops_interpolation", "interpolation") {

  constexpr auto LAST_LEV_P = ColInfo<NUM_INTERFACE_LEV>::LastPack;
  constexpr auto LAST_LEV   = ColInfo<NUM_PHYSICAL_LEV>::LastPack;
  constexpr auto LAST_LEV_PACK_IDX  = ColInfo<NUM_PHYSICAL_LEV>::LastPackEnd;

  ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV]>   d_midpoints_field_in  ("",num_elems);
  ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV_P]> d_interface_field_in  ("",num_elems);
  ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV]>   d_midpoints_field_out ("",num_elems);
  ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV_P]> d_interface_field_out ("",num_elems);

  auto h_midpoints_field_in  = Kokkos::create_mirror_view(d_midpoints_field_in);
  auto h_interface_field_in  = Kokkos::create_mirror_view(d_interface_field_in);
  auto h_midpoints_field_out = Kokkos::create_mirror_view(d_midpoints_field_out);
  auto h_interface_field_out = Kokkos::create_mirror_view(d_interface_field_out);

  // Fill input fields columns with 0,1,2,3,...
  for (int ie=0; ie<num_elems; ++ie) {
    for (int igp=0; igp<NUM_PTS; ++igp) {
      for (int jgp=0; jgp<NUM_PTS; ++jgp) {
        for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
          const int ilev = k / VECTOR_SIZE;
          const int ivec = k % VECTOR_SIZE;

          h_midpoints_field_in(ie,igp,jgp,ilev)[ivec] = k;
          h_interface_field_in(ie,igp,jgp,ilev)[ivec] = k;
        }
        // Fill last interface
        h_interface_field_in(ie,igp,jgp,LAST_LEV_P)[(NUM_INTERFACE_LEV-1) % VECTOR_SIZE] = NUM_INTERFACE_LEV-1;
      }
    }
  }
  Kokkos::deep_copy(d_midpoints_field_in,h_midpoints_field_in);
  Kokkos::deep_copy(d_interface_field_in,h_interface_field_in);

  SECTION ("deltas") {
    // Compute deltas
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        ColumnOps::compute_midpoint_delta(kv,Homme::subview(d_interface_field_in,kv.ie,igp,jgp),
                                           Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        ColumnOps::compute_interface_delta(kv,Homme::subview(d_midpoints_field_in,kv.ie,igp,jgp),
                                            Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    // Check answers
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));

          // First midpoint should be fine, but first interface should be 0
          REQUIRE (mid_out(0) == 1.0);
          REQUIRE (int_out(0) == 0.0);

          for (int k=1; k<NUM_PHYSICAL_LEV; ++k) {
            REQUIRE (mid_out(k) == 1.0);
            REQUIRE (int_out(k) == 1.0);
          }

          // Last interface should also be 0
          REQUIRE(int_out(NUM_INTERFACE_LEV-1) == 0.0);
        }
      }
    }
  }

  SECTION("interpolation") {

    // Interpolate
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        ColumnOps::compute_midpoint_values(kv,Homme::subview(d_interface_field_in,kv.ie,igp,jgp),
                                            Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        ColumnOps::compute_interface_values(kv,Homme::subview(d_midpoints_field_in,kv.ie,igp,jgp),
                                             Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });
    Kokkos::fence();

    // Check answers
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          // First midpoint should be fine, but first interface should be equal to the first midpoint
          REQUIRE (mid_out(0) == 0.5);
          REQUIRE (int_out(0) == 0.0);

          for (int k=1; k<NUM_PHYSICAL_LEV; ++k) {

            // Note: the difference of sign in k +/- 0.5 comes from the fact that midpoint k
            //       is between interfaces k and k+1, while interface k is between
            //       midpoints k-1 and k.
            REQUIRE (mid_out(k) == k+0.5);
            REQUIRE (int_out(k) == k-0.5);
          }

          // Last interface should be equal to the last midpoint
          REQUIRE(int_out(NUM_INTERFACE_LEV-1) == h_midpoints_field_in(ie,igp,jgp,LAST_LEV)[LAST_LEV_PACK_IDX]);
        }
      }
    }
  }

  SECTION("update_with_product") {

    Kokkos::deep_copy (d_midpoints_field_out,1.0);

    // Interpolate
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;
        auto int_in = Homme::subview(d_interface_field_in,kv.ie,igp,jgp);

        auto prod_provider = [&](const int ilev)->Scalar {
          return int_in(ilev)*int_in(ilev);
        };
        ColumnOps::compute_midpoint_values<CombineMode::Add>(
                  kv, prod_provider,
                  Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp)
        );
      });
    });

    // Check answers: x_m += average(x_i*x_i), so we should have
    //   x_m(k) = 1.0 + (k*k + (k+1)*(k+1))/2
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {

            REQUIRE (mid_out(k) == (1.0 + (k*k + (k+1)*(k+1))/2.0));
          }
        }
      }
    }
  }
}

TEST_CASE("col_ops_reduction", "packed_reduction") {

  if (!OnGpu<ExecSpace>::value) {
    ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV]>   d_midpoints_field_in  ("",num_elems);
    ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV_P]> d_interface_field_in  ("",num_elems);
    ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV]>   d_midpoints_field_out ("",num_elems);
    ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV_P]> d_interface_field_out ("",num_elems);

    auto h_midpoints_field_in  = Kokkos::create_mirror_view(d_midpoints_field_in);
    auto h_interface_field_in  = Kokkos::create_mirror_view(d_interface_field_in);
    auto h_midpoints_field_out = Kokkos::create_mirror_view(d_midpoints_field_out);
    auto h_interface_field_out = Kokkos::create_mirror_view(d_interface_field_out);

    std::random_device rd;
    using rngAlg = std::mt19937_64;
    const unsigned int catchRngSeed = Catch::rngSeed();
    const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
    std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
    rngAlg engine(seed);
    std::uniform_real_distribution<Real> pdf(0.01, 1.0);
    genRandArray(h_midpoints_field_in, engine, pdf);
    genRandArray(h_interface_field_in, engine, pdf);
    Kokkos::deep_copy(d_midpoints_field_in,h_midpoints_field_in);
    Kokkos::deep_copy(d_interface_field_in,h_interface_field_in);
    
    // To use with std library for checking results
    std::vector<double> mid_data(NUM_PHYSICAL_LEV);
    std::vector<double> int_data(NUM_INTERFACE_LEV);

    // Forward and inclusive
    SECTION("packed") {
      using CO = ColumnOps;
      constexpr int NPL = NUM_PHYSICAL_LEV;
      constexpr int NIL = NUM_INTERFACE_LEV;
      Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                           KOKKOS_LAMBDA(const TeamMember& team) {
        KernelVariables kv(team);
        Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                             [&](const int idx) {
          const int igp = idx / NUM_PTS;
          const int jgp = idx % NUM_PTS;

          // Providers
          auto provide_field_mid = [&] (const int ilev)->Scalar {
            return d_midpoints_field_in(kv.ie,igp,jgp,ilev);
          };

          auto provide_field_int = [&] (const int ilev)->Scalar {
            return d_interface_field_in(kv.ie,igp,jgp,ilev);
          };

          using pfm_type = decltype(provide_field_mid);
          using pfi_type = decltype(provide_field_int);
          Real& mid_sum = d_midpoints_field_out(kv.ie,igp,jgp,0)[0];
          Real& int_sum = d_interface_field_out(kv.ie,igp,jgp,0)[0];
          CO::column_reduction<NPL,pfm_type,true>(kv, provide_field_mid, mid_sum);
          CO::column_reduction<NIL,pfi_type,true>(kv, provide_field_int, int_sum);
        });
      });

      Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
      Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

      for (int ie=0; ie<num_elems; ++ie) {
        for (int igp=0; igp<NUM_PTS; ++igp) {
          for (int jgp=0; jgp<NUM_PTS; ++jgp) {
            // Compute using std functions

            // Copy from host view to std vector
            auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
            auto int_in = viewAsReal(Homme::subview(h_interface_field_in,ie,igp,jgp));
            std::copy_n(mid_in.data(),NPL,mid_data.begin());
            std::copy_n(int_in.data(),NIL,int_data.begin());

            // Manuallly do the reduction, in the right order
            Real sum_mid = 0.0;
            Real sum_int = 0.0;

            Real tmp_m[VECTOR_SIZE] = {0};
            Real tmp_i[VECTOR_SIZE] = {0};
            for (int ivec=0; ivec<VECTOR_SIZE; ++ivec) {
              for (int ilev=0; ilev<NUM_LEV-1; ++ilev) {
                tmp_m[ivec] += mid_data[ilev*VECTOR_SIZE+ivec];
              }
              for (int ilev=0; ilev<NUM_LEV_P-1; ++ilev) {
                tmp_i[ivec] += int_data[ilev*VECTOR_SIZE+ivec];
              }
            }
            for (int ivec=0; ivec<ColInfo<NUM_PHYSICAL_LEV>::LastPackLen; ++ivec) {
              tmp_m[ivec] += mid_data[(NUM_LEV-1)*VECTOR_SIZE+ivec];
            }
            for (int ivec=0; ivec<ColInfo<NUM_INTERFACE_LEV>::LastPackLen; ++ivec) {
              tmp_i[ivec] += int_data[(NUM_LEV_P-1)*VECTOR_SIZE+ivec];
            }

            for (int ivec=0; ivec<VECTOR_SIZE; ++ivec) {
              sum_mid += tmp_m[ivec];
              sum_int += tmp_i[ivec];
            }

            // Check answer
            auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
            auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
            REQUIRE (mid_out(0) == sum_mid);
            REQUIRE (int_out(0) == sum_int);
          }
        }
      }
    }
  }
}

TEST_CASE("col_ops_scan_sum", "scan_sum") {

  ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV]>   d_midpoints_field_in  ("",num_elems);
  ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV_P]> d_interface_field_in  ("",num_elems);
  ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV]>   d_midpoints_field_out ("",num_elems);
  ExecViewManaged<Scalar*[NUM_PTS][NUM_PTS][NUM_LEV_P]> d_interface_field_out ("",num_elems);

  auto h_midpoints_field_in  = Kokkos::create_mirror_view(d_midpoints_field_in);
  auto h_interface_field_in  = Kokkos::create_mirror_view(d_interface_field_in);
  auto h_midpoints_field_out = Kokkos::create_mirror_view(d_midpoints_field_out);
  auto h_interface_field_out = Kokkos::create_mirror_view(d_interface_field_out);

  std::random_device rd;
  using rngAlg = std::mt19937_64;
  const unsigned int catchRngSeed = Catch::rngSeed();
  const unsigned int seed = catchRngSeed==0 ? rd() : catchRngSeed;
  std::cout << "seed: " << seed << (catchRngSeed==0 ? " (catch rng seed was 0)\n" : "\n");
  rngAlg engine(seed);
  std::uniform_real_distribution<Real> pdf(0.01, 1.0);
  genRandArray(h_midpoints_field_in, engine, pdf);
  genRandArray(h_interface_field_in, engine, pdf);
  Kokkos::deep_copy(d_midpoints_field_in,h_midpoints_field_in);
  Kokkos::deep_copy(d_interface_field_in,h_interface_field_in);
  
  // To use with std library for checking results
  std::vector<double> mid_data(NUM_PHYSICAL_LEV);
  std::vector<double> int_data(NUM_INTERFACE_LEV);
  std::vector<double> mid_sums(NUM_PHYSICAL_LEV);
  std::vector<double> int_sums(NUM_INTERFACE_LEV);

  // Forward and inclusive
  SECTION("fwd_inclusive") {
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        // Providers
        auto provide_field_mid = [&] (const int ilev)->Scalar {
          return d_midpoints_field_in(kv.ie,igp,jgp,ilev);
        };

        auto provide_field_int = [&] (const int ilev)->Scalar {
          return d_interface_field_in(kv.ie,igp,jgp,ilev);
        };

        ColumnOps::column_scan<true,true,NUM_PHYSICAL_LEV>(kv, provide_field_mid,
                                       Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        ColumnOps::column_scan<true,true,NUM_INTERFACE_LEV>(kv, provide_field_int,
                                       Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          // Compute using std functions

          // Copy from host view to std vector
          auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
          auto int_in = viewAsReal(Homme::subview(h_interface_field_in,ie,igp,jgp));
          std::copy_n(mid_in.data(),NUM_PHYSICAL_LEV,mid_data.begin());
          std::copy_n(int_in.data(),NUM_INTERFACE_LEV,int_data.begin());

          // scan sum the std vectors
          std::fill(mid_sums.begin(),mid_sums.end(),0.0);
          std::fill(int_sums.begin(),int_sums.end(),0.0);
          std::partial_sum(mid_data.begin(),mid_data.end(),mid_sums.begin());
          std::partial_sum(int_data.begin(),int_data.end(),int_sums.begin());

          // Check answer
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            REQUIRE (mid_out(k) == mid_sums[k]);
            REQUIRE (int_out(k) == int_sums[k]);
          }
          REQUIRE (int_out(NUM_INTERFACE_LEV-1) == int_sums[NUM_INTERFACE_LEV-1]);
        }
      }
    }
  }

  // Backward and inclusive
  SECTION("bwd_inclusive") {
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        // Providers
        auto provide_field_mid = [&] (const int ilev)->Scalar {
          return d_midpoints_field_in(kv.ie,igp,jgp,ilev);
        };

        auto provide_field_int = [&] (const int ilev)->Scalar {
          return d_interface_field_in(kv.ie,igp,jgp,ilev);
        };

        ColumnOps::column_scan<false,true,NUM_PHYSICAL_LEV>(kv, provide_field_mid,
                                        Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        ColumnOps::column_scan<false,true,NUM_INTERFACE_LEV>(kv, provide_field_int,
                                        Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          // Compute using std functions

          // Copy from host view to std vector
          auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
          auto int_in = viewAsReal(Homme::subview(h_interface_field_in,ie,igp,jgp));
          std::copy_n(mid_in.data(),NUM_PHYSICAL_LEV,mid_data.begin());
          std::copy_n(int_in.data(),NUM_INTERFACE_LEV,int_data.begin());

          // scan sum the std vectors
          std::fill(mid_sums.begin(),mid_sums.end(),0.0);
          std::fill(int_sums.begin(),int_sums.end(),0.0);
          std::partial_sum(mid_data.rbegin(),mid_data.rend(),mid_sums.rbegin());
          std::partial_sum(int_data.rbegin(),int_data.rend(),int_sums.rbegin());

          // Check answer
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            REQUIRE (mid_out(k) == mid_sums[k]);
            REQUIRE (int_out(k) == int_sums[k]);
          }
          REQUIRE (int_out(NUM_INTERFACE_LEV-1) == int_sums[NUM_INTERFACE_LEV-1]);
        }
      }
    }
  }

  // Forward and exclusive
  SECTION("fwd_exclusive") {
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        // Providers
        auto provide_field_mid = [&] (const int ilev)->Scalar {
          return d_midpoints_field_in(kv.ie,igp,jgp,ilev);
        };

        auto provide_field_int = [&] (const int ilev)->Scalar {
          return d_interface_field_in(kv.ie,igp,jgp,ilev);
        };

        ColumnOps::column_scan<true,false,NUM_PHYSICAL_LEV>(kv, provide_field_mid,
                                        Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        ColumnOps::column_scan<true,false,NUM_INTERFACE_LEV>(kv, provide_field_int,
                                        Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          // Compute using std functions

          // Copy from host view to std vector
          auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
          auto int_in = viewAsReal(Homme::subview(h_interface_field_in,ie,igp,jgp));
          std::copy_n(mid_in.data(),NUM_PHYSICAL_LEV,mid_data.begin());
          std::copy_n(int_in.data(),NUM_INTERFACE_LEV,int_data.begin());

          // scan sum the std vectors
          // std::exclusive_sum is only in c++17, not c++11. To emulate that, do a standard partial sum,
          // but store results one element later. Obiously, you also have to stop one element sooner
          std::fill(mid_sums.begin(),mid_sums.end(),0.0);
          std::fill(int_sums.begin(),int_sums.end(),0.0);
          std::partial_sum(mid_data.begin(),--mid_data.end(),++mid_sums.begin());
          std::partial_sum(int_data.begin(),--int_data.end(),++int_sums.begin());
    
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            REQUIRE (mid_out(k) == mid_sums[k]);
            REQUIRE (int_out(k) == int_sums[k]);
          }
          REQUIRE (int_out(NUM_INTERFACE_LEV-1) == int_sums[NUM_INTERFACE_LEV-1]);
        }
      }
    }
  }

  // Backward and exclusive
  SECTION("bwd_exclusive") {
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        // Providers
        auto provide_field_mid = [&] (const int ilev)->Scalar {
          return d_midpoints_field_in(kv.ie,igp,jgp,ilev);
        };

        auto provide_field_int = [&] (const int ilev)->Scalar {
          return d_interface_field_in(kv.ie,igp,jgp,ilev);
        };

        ColumnOps::column_scan<false,false,NUM_PHYSICAL_LEV>(kv, provide_field_mid,
                                         Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        ColumnOps::column_scan<false,false,NUM_INTERFACE_LEV>(kv, provide_field_int,
                                         Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          // Compute using std functions

          // Copy from host view to std vector
          auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
          auto int_in = viewAsReal(Homme::subview(h_interface_field_in,ie,igp,jgp));
          std::copy_n(mid_in.data(),NUM_PHYSICAL_LEV,mid_data.begin());
          std::copy_n(int_in.data(),NUM_INTERFACE_LEV,int_data.begin());

          // scan sum the std vectors
          // std::exclusive_sum is only in c++17, not c++11. To emulate that, do a standard partial sum,
          // but store results one element later. Obiously, you also have to stop one element sooner
          std::fill(mid_sums.begin(),mid_sums.end(),0.0);
          std::fill(int_sums.begin(),int_sums.end(),0.0);
          std::partial_sum(mid_data.rbegin(),--mid_data.rend(),++mid_sums.rbegin());
          std::partial_sum(int_data.rbegin(),--int_data.rend(),++int_sums.rbegin());
    
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            REQUIRE (mid_out(k) == mid_sums[k]);
            REQUIRE (int_out(k) == int_sums[k]);
          }
          REQUIRE (int_out(NUM_INTERFACE_LEV-1) == int_sums[NUM_INTERFACE_LEV-1]);
        }
      }
    }
  }

  // Forward midpoints to interfaces
  SECTION("fwd_mid_to_int_field") {
    const Real s0 = pdf(engine);
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        auto in  = Homme::subview(d_midpoints_field_in,kv.ie,igp,jgp);
        auto out = Homme::subview(d_interface_field_out,kv.ie,igp,jgp);

        out(0)[0] = s0;
        ColumnOps::column_scan_mid_to_int<true>(kv, in, out);
      });
    });

    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          // Compute manually

          // Copy from host view to std vector
          auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
          std::copy_n(mid_in.data(),NUM_PHYSICAL_LEV,mid_data.begin());

          std::fill(int_sums.begin(),int_sums.end(),0.0);
          int_sums[0] = s0; // Init the first entry
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            int_sums[k+1] = int_sums[k] + mid_in[k];
          }

          // Check answer
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            REQUIRE (int_out(k) == int_sums[k]);
          }
          REQUIRE (int_out(NUM_INTERFACE_LEV-1) == int_sums[NUM_INTERFACE_LEV-1]);
        }
      }
    }
  }

  // Backward midpoints to interfaces
  SECTION("bwd_mid_to_int_field") {
    const Real s0 = pdf(engine);

    using Info = ColInfo<NUM_INTERFACE_LEV>;
    constexpr int LAST_PACK     = Info::LastPack;
    constexpr int LAST_PACK_END = Info::LastPackEnd;
    Kokkos::deep_copy(d_interface_field_out,0.0);
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        auto in  = Homme::subview(d_midpoints_field_in,kv.ie,igp,jgp);
        auto out = Homme::subview(d_interface_field_out,kv.ie,igp,jgp);

        out(LAST_PACK)[LAST_PACK_END] = s0;
        ColumnOps::column_scan_mid_to_int<false>(kv, in, out);
      });
    });

    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          // Compute manually

          // Copy from host view to std vector
          auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
          std::copy_n(mid_in.data(),NUM_PHYSICAL_LEV,mid_data.begin());

          std::fill(int_sums.begin(),int_sums.end(),0.0);
          int_sums[NUM_INTERFACE_LEV-1] = s0; // Init the last entry
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            int kk = NUM_INTERFACE_LEV - 1 - k;
            int_sums[kk-1] = int_sums[kk] + mid_in[kk-1];
          }

          // Check answer
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=NUM_INTERFACE_LEV-1; k>=0; --k) {
            REQUIRE (int_out(k) == int_sums[k]);
          }
        }
      }
    }
  }

  // Forward midpoints to interfaces
  SECTION("fwd_mid_to_int_provider") {
    const Real s0 = pdf(engine);
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        auto in  = Homme::subview(d_midpoints_field_in,kv.ie,igp,jgp);
        auto provider = [&](const int ilev) -> Scalar{
          return in(ilev);
        };
        auto out = Homme::subview(d_interface_field_out,kv.ie,igp,jgp);

        out(0)[0] = s0;
        ColumnOps::column_scan_mid_to_int<true>(kv, provider, out);
      });
    });

    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          // Compute manually

          // Copy from host view to std vector
          auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
          std::copy_n(mid_in.data(),NUM_PHYSICAL_LEV,mid_data.begin());

          std::fill(int_sums.begin(),int_sums.end(),0.0);
          int_sums[0] = s0; // Init the first entry
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            int_sums[k+1] = int_sums[k] + mid_in[k];
          }

          // Check answer
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_INTERFACE_LEV; ++k) {
            REQUIRE (int_out(k) == int_sums[k]);
          }
        }
      }
    }
  }

  // Backward midpoints to interfaces
  SECTION("bwd_mid_to_int_provider") {
    const Real s0 = pdf(engine);
    using Info = ColInfo<NUM_INTERFACE_LEV>;
    constexpr int LAST_PACK     = Info::LastPack;
    constexpr int LAST_PACK_END = Info::LastPackEnd;
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NUM_PTS*NUM_PTS),
                           [&](const int idx) {
        const int igp = idx / NUM_PTS;
        const int jgp = idx % NUM_PTS;

        auto in  = Homme::subview(d_midpoints_field_in,kv.ie,igp,jgp);
        auto provider = [&](const int ilev) -> Scalar{
          return in(ilev);
        };
        auto out = Homme::subview(d_interface_field_out,kv.ie,igp,jgp);

        out(LAST_PACK)[LAST_PACK_END] = s0;
        ColumnOps::column_scan_mid_to_int<false>(kv, provider, out);
      });
    });

    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NUM_PTS; ++igp) {
        for (int jgp=0; jgp<NUM_PTS; ++jgp) {
          // Compute manually

          // Copy from host view to std vector
          auto mid_in = viewAsReal(Homme::subview(h_midpoints_field_in,ie,igp,jgp));
          std::copy_n(mid_in.data(),NUM_PHYSICAL_LEV,mid_data.begin());

          std::fill(int_sums.begin(),int_sums.end(),0.0);
          int_sums[NUM_INTERFACE_LEV-1] = s0; // Init the last entry
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            int kk = NUM_INTERFACE_LEV - 1 - k;
            int_sums[kk-1] = int_sums[kk] + mid_in[kk-1];
          }

          // Check answer
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_INTERFACE_LEV; ++k) {
            REQUIRE (int_out(k) == int_sums[k]);
          }
        }
      }
    }
  }
}
