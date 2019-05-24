#include <catch2/catch.hpp>

#include <iostream>
#include <random>

#include "ElementOps.hpp"
#include "Types.hpp"

#include "utilities/TestUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"
#include "utilities/ViewUtils.hpp"

using namespace Homme;


// ============= ELEMENT OPS ================ //

TEST_CASE("elem_ops", "elem_ops") {

  constexpr int num_elems = 10;

  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   d_midpoints_field_in  ("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> d_interface_field_in  ("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV]>   d_midpoints_field_out ("",num_elems);
  ExecViewManaged<Scalar*[NP][NP][NUM_LEV_P]> d_interface_field_out ("",num_elems);

  auto h_midpoints_field_in  = Kokkos::create_mirror_view(d_midpoints_field_in);
  auto h_interface_field_in  = Kokkos::create_mirror_view(d_interface_field_in);
  auto h_midpoints_field_out = Kokkos::create_mirror_view(d_midpoints_field_out);
  auto h_interface_field_out = Kokkos::create_mirror_view(d_interface_field_out);

  ElementOps elem_ops;

  // Fill input fields columns with 0,1,2,3,...
  for (int ie=0; ie<num_elems; ++ie) {
    for (int igp=0; igp<NP; ++igp) {
      for (int jgp=0; jgp<NP; ++jgp) {
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
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        elem_ops.compute_midpoint_delta(kv,Homme::subview(d_interface_field_in,kv.ie,igp,jgp),
                                           Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        elem_ops.compute_interface_delta(kv,Homme::subview(d_midpoints_field_in,kv.ie,igp,jgp),
                                            Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    // Check answers
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
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
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        elem_ops.compute_midpoint_values(kv,Homme::subview(d_interface_field_in,kv.ie,igp,jgp),
                                            Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        elem_ops.compute_interface_values(kv,Homme::subview(d_midpoints_field_in,kv.ie,igp,jgp),
                                             Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    // Check answers
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
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
          REQUIRE(int_out(NUM_INTERFACE_LEV-1) == h_midpoints_field_in(ie,igp,jgp,LAST_LEV)[LAST_MIDPOINT_VEC_IDX]);
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
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        elem_ops.update_midpoint_values_with_product(
                  kv,1.0,
                  Homme::subview(d_interface_field_in,kv.ie,igp,jgp),
                  Homme::subview(d_interface_field_in,kv.ie,igp,jgp),
                  Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp)
        );
      });
    });

    // Check answers: x_m += average(x_i*x_i), so we should have
    //   x_m(k) = 1.0 + (k*k + (k+1)*(k+1))/2
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {

            REQUIRE (mid_out(k) == (1.0 + (k*k + (k+1)*(k+1))/2.0));
          }
        }
      }
    }
  }

  // Precompute squares.
  auto square = [](const double x) ->double { return x*x; };

  std::vector<double> mid_squares(NUM_PHYSICAL_LEV);
  std::vector<double> int_squares(NUM_INTERFACE_LEV);

  std::iota(mid_squares.begin(),mid_squares.end(), 1.0);
  std::iota(int_squares.begin(),int_squares.end(), 1.0);

  std::transform(mid_squares.begin(),mid_squares.end(),mid_squares.begin(),square);
  std::transform(int_squares.begin(),int_squares.end(),int_squares.begin(),square);

  // Copy to host/device views
  auto h_mid_in_data = reinterpret_cast<Real*>(Homme::subview(h_midpoints_field_in,0,0,0).data());
  auto h_int_in_data = reinterpret_cast<Real*>(Homme::subview(h_interface_field_in,0,0,0).data());

  std::copy(mid_squares.begin(),mid_squares.end(),h_mid_in_data);
  std::copy(int_squares.begin(),int_squares.end(),h_int_in_data);

  Kokkos::deep_copy(d_midpoints_field_in,h_midpoints_field_in);
  Kokkos::deep_copy(d_interface_field_in,h_interface_field_in);

  // Providers
  auto provide_square_m = KOKKOS_LAMBDA (const int ilev)->Scalar {
    return d_midpoints_field_in(0,0,0,ilev);
  };

  auto provide_square_i = KOKKOS_LAMBDA (const int ilev)->Scalar {
    return d_interface_field_in(0,0,0,ilev);
  };

  std::vector<double> mid_sums(NUM_PHYSICAL_LEV);
  std::vector<double> int_sums(NUM_INTERFACE_LEV);

  SECTION("column_scan_fwd_inclusive") {

    std::partial_sum(mid_squares.begin(),mid_squares.end(),mid_sums.begin());
    std::partial_sum(int_squares.begin(),int_squares.end(),int_sums.begin());
    
    // Scan sum forward
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        elem_ops.column_scan<true,true,NUM_PHYSICAL_LEV>(kv, provide_square_m,
                                       Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        elem_ops.column_scan<true,true,NUM_INTERFACE_LEV>(kv, provide_square_i,
                                       Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    // Check answers: out(n) = sum_{k=1}^n k*k
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
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

  SECTION("column_scan_bwd_inclusive") {
    std::fill(mid_sums.begin(),mid_sums.end(),0.0);
    std::fill(int_sums.begin(),int_sums.end(),0.0);

    std::partial_sum(mid_squares.rbegin(),mid_squares.rend(),mid_sums.rbegin());
    std::partial_sum(int_squares.rbegin(),int_squares.rend(),int_sums.rbegin());
    
    // Scan sum forward
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        elem_ops.column_scan<false,true,NUM_PHYSICAL_LEV>(kv, provide_square_m,
                                        Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        elem_ops.column_scan<false,true,NUM_INTERFACE_LEV>(kv, provide_square_i,
                                        Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    // Check answers: out(n) = sum_{k=1}^n k*k
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          auto mid_out = viewAsReal(Homme::subview(h_midpoints_field_out,ie,igp,jgp));
          auto int_out = viewAsReal(Homme::subview(h_interface_field_out,ie,igp,jgp));
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
// printf("sum {%d^2,...,%d^2} = %3.12f\n",NUM_PHYSICAL_LEV+1,k+1,mid_out(k));

            REQUIRE (mid_out(k) == mid_sums[k]);
            REQUIRE (int_out(k) == int_sums[k]);
          }
          REQUIRE (int_out(NUM_INTERFACE_LEV-1) == int_sums[NUM_INTERFACE_LEV-1]);
        }
      }
    }
  }

  SECTION("column_scan_fwd_exclusive") {
    std::fill(mid_sums.begin(),mid_sums.end(),0.0);
    std::fill(int_sums.begin(),int_sums.end(),0.0);

    // std::exclusive_sum is only in c++17, not c++11. To emulate that, do a standard partial sum,
    // but store results one element later. Obiously, you also have to stop one element sooner
    std::partial_sum(mid_squares.begin(),--mid_squares.end(),++mid_sums.begin());
    std::partial_sum(int_squares.begin(),--int_squares.end(),++int_sums.begin());
    
    // Scan sum forward
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        elem_ops.column_scan<true,false,NUM_PHYSICAL_LEV>(kv, provide_square_m,
                                        Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        elem_ops.column_scan<true,false,NUM_INTERFACE_LEV>(kv, provide_square_i,
                                        Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    // Check answers: out(n) = sum_{k=1}^n k*k
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
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

  SECTION("column_scan_bwd_exclusive") {
    std::fill(mid_sums.begin(),mid_sums.end(),0.0);
    std::fill(int_sums.begin(),int_sums.end(),0.0);

    // std::exclusive_sum is only in c++17, not c++11. To emulate that, do a standard partial sum,
    // but store results one element later. Obiously, you also have to stop one element sooner
    std::partial_sum(mid_squares.rbegin(),--mid_squares.rend(),++mid_sums.rbegin());
    std::partial_sum(int_squares.rbegin(),--int_squares.rend(),++int_sums.rbegin());
    
    // Scan sum forward
    Kokkos::parallel_for(Homme::get_default_team_policy<ExecSpace>(num_elems),
                         KOKKOS_LAMBDA(const TeamMember& team) {
      KernelVariables kv(team);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team,NP*NP),
                           [&](const int idx) {
        const int igp = idx / NP;
        const int jgp = idx % NP;

        elem_ops.column_scan<false,false,NUM_PHYSICAL_LEV>(kv, provide_square_m,
                                         Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
        elem_ops.column_scan<false,false,NUM_INTERFACE_LEV>(kv, provide_square_i,
                                         Homme::subview(d_interface_field_out,kv.ie,igp,jgp));
      });
    });

    // Check answers: out(n) = sum_{k=1}^n k*k
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);

    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
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
}
