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

  constexpr int last_midpoint_vec_idx  = (NUM_PHYSICAL_LEV + VECTOR_SIZE - 1) % VECTOR_SIZE;
  constexpr int last_interface_vec_idx = (last_midpoint_vec_idx+1) % VECTOR_SIZE;
  constexpr int LAST_LEV_P = NUM_LEV_P-1;

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
          // First midpoint should be fine, but first interface should be 0
          REQUIRE (h_midpoints_field_out(ie,igp,jgp,0)[0] == 1.0);
          REQUIRE (h_interface_field_out(ie,igp,jgp,0)[0] == 0.0);

          for (int k=1; k<NUM_PHYSICAL_LEV; ++k) {
            const int ilev = k / VECTOR_SIZE;
            const int ivec = k % VECTOR_SIZE;

            REQUIRE (h_midpoints_field_out(ie,igp,jgp,ilev)[ivec] == 1.0);
            REQUIRE (h_interface_field_out(ie,igp,jgp,ilev)[ivec] == 1.0);
          }

          // Last interface should also be 0
          REQUIRE(h_interface_field_out(ie,igp,jgp,LAST_LEV_P)[last_interface_vec_idx] == 0.0);
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
          // First midpoint should be fine, but first interface should be equal to the first midpoint
          REQUIRE (h_midpoints_field_out(ie,igp,jgp,0)[0] == 0.5);
          REQUIRE (h_interface_field_out(ie,igp,jgp,0)[0] == 0.0);

          for (int k=1; k<NUM_PHYSICAL_LEV; ++k) {
            const int ilev = k / VECTOR_SIZE;
            const int ivec = k % VECTOR_SIZE;

            // Note: the difference of sign in k +/- 0.5 comes from the fact that midpoint k
            //       is between interfaces k and k+1, while interface k is between
            //       midpoints k-1 and k.
            REQUIRE (h_midpoints_field_out(ie,igp,jgp,ilev)[ivec] == k+0.5);
            REQUIRE (h_interface_field_out(ie,igp,jgp,ilev)[ivec] == k-0.5);
          }

          // Last interface should be equal to the last midpoint
          REQUIRE(h_interface_field_out(ie,igp,jgp,LAST_LEV_P)[last_interface_vec_idx] == h_midpoints_field_in(ie,igp,jgp,NUM_LEV-1)[last_midpoint_vec_idx]);
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

        elem_ops.update_midpoint_values_with_product(kv,1.0,
                                                        Homme::subview(d_interface_field_in,kv.ie,igp,jgp),
                                                        Homme::subview(d_interface_field_in,kv.ie,igp,jgp),
                                                        Homme::subview(d_midpoints_field_out,kv.ie,igp,jgp));
      });
    });

    // Check answers: x_m += average(x_i*x_i), so we should have
    //   x_m(k) = 1.0 + (k*k + (k+1)*(k+1))/2
    Kokkos::deep_copy(h_midpoints_field_out,d_midpoints_field_out);
    Kokkos::deep_copy(h_interface_field_out,d_interface_field_out);
    for (int ie=0; ie<num_elems; ++ie) {
      for (int igp=0; igp<NP; ++igp) {
        for (int jgp=0; jgp<NP; ++jgp) {
          for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
            const int ilev = k / VECTOR_SIZE;
            const int ivec = k % VECTOR_SIZE;

            REQUIRE (h_midpoints_field_out(ie,igp,jgp,ilev)[ivec] == (1.0 + (k*k + (k+1)*(k+1))/2.0));
          }
        }
      }
    }
  }
}
