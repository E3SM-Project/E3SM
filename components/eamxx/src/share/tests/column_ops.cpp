#include <catch2/catch.hpp>

#include "ekat/kokkos/ekat_kokkos_types.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"
#include "ekat/util/ekat_arch.hpp"
#include "share/util/scream_column_ops.hpp"

namespace {

TEST_CASE ("combine_ops") {
  using namespace scream;
  using pack_type = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  constexpr auto Replace     = CombineMode::Replace;
  constexpr auto Scale       = CombineMode::Scale;
  constexpr auto Update      = CombineMode::Update;
  constexpr auto ScaleUpdate = CombineMode::ScaleUpdate;
  constexpr auto ScaleAdd    = CombineMode::ScaleAdd;
  constexpr auto Add         = CombineMode::Add;
  constexpr auto Multiply    = CombineMode::Multiply;
  constexpr auto Divide      = CombineMode::Divide;

  pack_type two (2.0);
  pack_type four (4.0);
  pack_type six (6.0);
  pack_type x;

  x = two;
  combine<Replace>(two,x);
  REQUIRE ( (x==two).all() );

  x = two;
  combine<Scale>(two,x,3.0);
  REQUIRE ( (x==six).all() );

  x = two;
  combine<Update>(two,x,1.0,2.0);
  REQUIRE ( (x==six).all() );

  x = two;
  combine<ScaleUpdate>(two,x,2.0,1.0);
  REQUIRE ( (x==six).all() );

  x = two;
  combine<ScaleAdd>(two,x,2.0);
  REQUIRE ( (x==six).all() );

  x = two;
  combine<Add>(two,x);
  REQUIRE ( (x==four).all() );

  x = two;
  combine<Multiply>(two,x);
  REQUIRE ( (x==four).all() );

  x = four;
  combine<Divide>(two,x);
  REQUIRE ( (x==two).all() );
}

TEST_CASE("column_ops_ps_1") {
  using namespace scream;
  using device_type = DefaultDevice;
  using KT = ekat::KokkosTypes<device_type>;
  using pack_type = ekat::Pack<Real,1>;
  using exec_space = typename device_type::execution_space;

  using view_2d_type  = KT::view_2d<pack_type>;
  using policy_type   = KT::TeamPolicy;
  using member_type   = KT::MemberType;
  using column_ops    = ColumnOps<device_type,Real>;

  constexpr int num_cols = 1;
  constexpr int num_levs = 16;

  policy_type policy(num_cols,std::min(num_levs,exec_space::concurrency()));

  view_2d_type v_int("",num_cols,num_levs+1);
  view_2d_type v_mid("",num_cols,num_levs);
  view_2d_type dv_mid("",num_cols,num_levs);
  view_2d_type dz_mid("",num_cols,num_levs);
  auto v_int_h = Kokkos::create_mirror_view(v_int);
  auto v_mid_h = Kokkos::create_mirror_view(v_mid);
  auto dv_mid_h = Kokkos::create_mirror_view(dv_mid);
  auto dz_mid_h = Kokkos::create_mirror_view(dz_mid);

  SECTION ("int_to_mid") {

    // Fill v_int with odd numbers, so that v_mid should contain even ones.
    for (int k=0; k<num_levs+1; ++k) {
      v_int_h(0,k)[0] = 2*k+1;
    }
    Kokkos::deep_copy(v_int,v_int_h);

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto v_m = ekat::subview(v_mid,icol);

      column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_mid_h,v_mid);
    for (int k=0; k<num_levs; ++k) {
      REQUIRE (v_mid_h(0,k)[0] == 2*(k+1) );
    }

    // Re-do with lambda as provider
    Kokkos::deep_copy(v_mid,pack_type(0));
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      
      auto v_i = [&](const int k)->pack_type {
        return v_int(icol,k);
      };
      auto v_m = ekat::subview(v_mid,icol);

      column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_mid_h,v_mid);
    for (int k=0; k<num_levs; ++k) {
      REQUIRE (v_mid_h(0,k)[0] == 2*(k+1) );
    }
  }

  SECTION ("mid_to_int_linear") {

    // Fill v_mid with even numbers and set dz=1, so that v_int should contain odd ones.
    for (int k=0; k<num_levs; ++k) {
      v_mid_h(0,k)[0] = 2*(k+1);
      dz_mid_h(0,k)[0] = 1;
    }
    Kokkos::deep_copy(v_mid,v_mid_h);
    Kokkos::deep_copy(dz_mid,dz_mid_h);

    // Use boundary conditions s.t. v_int becomes monotone increasing odd numbers
    const Real bc_top = 1;
    const Real bc_bot = 2*num_levs+1;

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto v_m = ekat::subview(v_mid,icol);
      auto dz  = ekat::subview(dz_mid,icol);

      column_ops::compute_interface_values_linear(team,num_levs,v_m,dz,bc_top,bc_bot,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1 );
    }

    // Re-do with lambda as provider
    Kokkos::deep_copy(v_int,pack_type(0));
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();

      auto v_m = [&](const int k)->pack_type {
        return v_mid(icol,k);
      };
      auto dz = [&](const int k)->pack_type {
        return dz_mid(icol,k);
      };
      auto v_i = ekat::subview(v_int,icol);

      column_ops::compute_interface_values_linear(team,num_levs,v_m,dz,bc_top,bc_bot,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1 );
    }
  }

  SECTION ("int_to_mid_to_int_linear") {

    // Fill v_mid with even numbers and set dz=1, so that v_int should contain odd ones.
    for (int k=0; k<num_levs+1; ++k) {
      v_int_h(0,k)[0] = 2*k+1;
      if (k < num_levs) dz_mid_h(0,k)[0] = 1;
    }
    Kokkos::deep_copy(v_int,v_int_h);
    Kokkos::deep_copy(dz_mid,dz_mid_h);

    // Use boundary conditions s.t. v_int becomes monotone increasing odd numbers
    const Real bc_top = v_int_h(0,0)[0];
    const Real bc_bot = v_int_h(0,num_levs)[0];

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto v_m = ekat::subview(v_mid,icol);
      auto dz  = ekat::subview(dz_mid,icol);

      column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
      team.team_barrier();
      column_ops::compute_interface_values_linear(team,num_levs,v_m,dz,bc_top,bc_bot,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1);
    }
  }

  SECTION ("mid_to_int_compatible_fix_top") {

    // Fill v_mid with even numbers, so that v_int should contain odd ones.
    for (int k=0; k<num_levs; ++k) {
      v_mid_h(0,k)[0] = 2*(k+1);
    }
    Kokkos::deep_copy(v_mid,v_mid_h);

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto v_m = ekat::subview(v_mid,icol);

      column_ops::compute_interface_values_compatible<true>(team,num_levs,v_m,1,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1 );
    }

    // Re-do with lambda as provider
    Kokkos::deep_copy(v_int,pack_type(0));
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      
      auto v_m = [&](const int k)->pack_type {
        return v_mid(icol,k);
      };
      auto v_i = ekat::subview(v_int,icol);

      column_ops::compute_interface_values_compatible<true>(team,num_levs,v_m,1,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1 );
    }
  }

  SECTION ("mid_to_int_compatible_fix_bot") {

    // Fill v_mid with even numbers, so that v_int should contain odd ones.
    for (int k=0; k<num_levs; ++k) {
      v_mid_h(0,k)[0] = 2*(k+1);
    }
    Kokkos::deep_copy(v_mid,v_mid_h);

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto v_m = ekat::subview(v_mid,icol);

      column_ops::compute_interface_values_compatible<false>(team,num_levs,v_m,2*num_levs+1,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1 );
    }

    // Re-do with lambda as provider
    Kokkos::deep_copy(v_int,pack_type(0));
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      
      auto v_m = [&](const int k)->pack_type {
        return v_mid(icol,k);
      };
      auto v_i = ekat::subview(v_int,icol);

      column_ops::compute_interface_values_compatible<false>(team,num_levs,v_m,2*num_levs+1,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1 );
    }
  }

  SECTION ("int_to_mid_to_int_compatible_fix_top") {

    // Fill v_int with odd numbers, so that v_mid should contain even ones.
    for (int k=0; k<num_levs+1; ++k) {
      v_int_h(0,k)[0] = 2*k+1;
    }
    Kokkos::deep_copy(v_int,v_int_h);

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto v_m = ekat::subview(v_mid,icol);

      column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
      team.team_barrier();
      column_ops::compute_interface_values_compatible<true>(team,num_levs,v_m,1,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1 );
    }
  }

  SECTION ("int_to_mid_to_int_compatible_fix_bot") {

    // Fill v_int with odd numbers, so that v_mid should contain even ones.
    for (int k=0; k<num_levs+1; ++k) {
      v_int_h(0,k)[0] = 2*k+1;
    }
    Kokkos::deep_copy(v_int,v_int_h);

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto v_m = ekat::subview(v_mid,icol);

      column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
      team.team_barrier();
      column_ops::compute_interface_values_compatible<false>(team,num_levs,v_m,2*num_levs+1,v_i);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == 2*k+1 );
    }
  }

  SECTION ("delta") {

    // Fill v_int with odd numbers, so that dv_mid=2 on all levels
    for (int k=0; k<num_levs+1; ++k) {
      v_int_h(0,k)[0] = 2*k+1;
    }
    Kokkos::deep_copy(v_int,v_int_h);

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto dv_m = ekat::subview(dv_mid,icol);

      column_ops::compute_midpoint_delta(team,num_levs,v_i,dv_m);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(dv_mid_h,dv_mid);
    for (int k=0; k<num_levs; ++k) {
      REQUIRE (dv_mid_h(0,k)[0] == 2 );
    }

    // Re-do with lambda as provider
    Kokkos::deep_copy(dv_mid,pack_type(0));
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      
      auto v_i = [&](const int k)->pack_type {
        return v_int(icol,k);
      };
      auto dv_m = ekat::subview(dv_mid,icol);

      column_ops::compute_midpoint_delta(team,num_levs,v_i,dv_m);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(dv_mid_h,dv_mid);
    for (int k=0; k<num_levs; ++k) {
      REQUIRE (dv_mid_h(0,k)[0] == 2 );
    }
  }

  SECTION ("scan_from_top") {

    const Real s0 = -2.0; // If it worskfor s0!=0, it works for s0=0.

    // Fill v_mid with integers, then use Gauss formula to check the sum
    for (int k=0; k<num_levs; ++k) {
      dv_mid_h(0,k)[0] = k+1;
    }
    Kokkos::deep_copy(dv_mid,dv_mid_h);

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto dv_m = ekat::subview(dv_mid,icol);

      column_ops::column_scan<true>(team,num_levs,dv_m,v_i,s0);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == s0+k*(k+1)/2.0);
    }

    // Re-do with lambda as provider
    Kokkos::deep_copy(v_int,pack_type(0));
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      
      auto dv_m = [&](const int k)->pack_type {
        return dv_mid(icol,k);
      };
      auto v_i = ekat::subview(v_int,icol);

      column_ops::column_scan<true>(team,num_levs,dv_m,v_i,s0);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      REQUIRE (v_int_h(0,k)[0] == s0+k*(k+1)/2.0);
    }
  }

  SECTION ("scan_from_bot") {

    const Real s0 = -2.0; // If it worskfor s0!=0, it works for s0=0.

    // Fill v_mid with integers, then use Gauss formula to check the sum
    for (int k=0; k<num_levs; ++k) {
      const auto k_bwd = num_levs - k - 1;
      dv_mid_h(0,k_bwd)[0] = k+1;
    }
    Kokkos::deep_copy(dv_mid,dv_mid_h);

    // Run column kernel
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      auto v_i = ekat::subview(v_int,icol);
      auto dv_m = ekat::subview(dv_mid,icol);

      column_ops::column_scan<false>(team,num_levs,dv_m,v_i,s0);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      const auto k_bwd = num_levs - k;
      REQUIRE (v_int_h(0,k_bwd)[0] == s0+k*(k+1)/2.0);
    }

    // Re-do with lambda as provider
    Kokkos::deep_copy(v_int,pack_type(0));
    Kokkos::parallel_for(policy,
                         KOKKOS_LAMBDA(const member_type& team){
      const int icol = team.league_rank();
      
      auto dv_m = [&](const int k)->pack_type {
        return dv_mid(icol,k);
      };
      auto v_i = ekat::subview(v_int,icol);

      column_ops::column_scan<false>(team,num_levs,dv_m,v_i,s0);
    });
    Kokkos::fence();

    // Check answer
    Kokkos::deep_copy(v_int_h,v_int);
    for (int k=0; k<num_levs+1; ++k) {
      const auto k_bwd = num_levs - k;
      REQUIRE (v_int_h(0,k_bwd)[0] == s0+k*(k+1)/2.0);
    }
  }
}

TEST_CASE("column_ops_ps_N") {
  // No point in re-running test for a larger pack size
  if (!ekat::OnGpu<ekat::DefaultDevice::execution_space>::value) {
    using namespace scream;
    using device_type = DefaultDevice;
    using KT = ekat::KokkosTypes<device_type>;
    constexpr int ps = SCREAM_PACK_SIZE;
    using pack_type = ekat::Pack<Real,ps>;
    using exec_space = typename device_type::execution_space;

    using view_2d_type  = KT::view_2d<pack_type>;
    using policy_type   = KT::TeamPolicy;
    using member_type   = KT::MemberType;
    using column_ops    = ColumnOps<device_type,Real>;

    constexpr int num_cols = 1;
    // Test both the case where num_mid_packs==num_int_packs, and
    // the case where num_int_packs=num_mid_packs+1.
    for (int num_levs : {2*ps, 2*ps+1} ) {
      std::cout << "num_levs: " << num_levs << "\n";
      const int num_mid_packs = ekat::PackInfo<ps>::num_packs(num_levs);
      const int num_int_packs = ekat::PackInfo<ps>::num_packs(num_levs+1);

      view_2d_type v_mid("",num_cols,num_mid_packs);
      view_2d_type v_int("",num_cols,num_int_packs);
      view_2d_type dv_mid("",num_cols,num_mid_packs);
      view_2d_type dz_mid("",num_cols,num_mid_packs);
      auto v_int_h = Kokkos::create_mirror_view(v_int);
      auto v_mid_h = Kokkos::create_mirror_view(v_mid);
      auto dv_mid_h = Kokkos::create_mirror_view(dv_mid);
      auto dz_mid_h = Kokkos::create_mirror_view(dz_mid);

      policy_type policy(num_cols,std::min(num_mid_packs,exec_space::concurrency()));

      SECTION ("int_to_mid") {

        // Fill v_int with odd numbers, so that v_mid should contain even ones.
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          v_int_h(0,ipack)[ivec] = 2*k+1;
        }
        Kokkos::deep_copy(v_int,v_int_h);

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto v_m = ekat::subview(v_mid,icol);

          column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_mid_h,v_mid);
        for (int k=0; k<num_levs; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_mid_h(0,ipack)[ivec] == 2*(k+1) );
        }

        // Re-do with lambda as provider
        Kokkos::deep_copy(v_mid,pack_type(0));
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          
          auto v_i = [&](const int k)->pack_type {
            return v_int(icol,k);
          };
          auto v_m = ekat::subview(v_mid,icol);

          column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_mid_h,v_mid);
        for (int k=0; k<num_levs; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_mid_h(0,ipack)[ivec] == 2*(k+1) );
        }
      }

      SECTION ("mid_to_int_linear") {

        // Fill v_mid with even numbers and set dz=1, so that v_int should contain odd ones.
        for (int k=0; k<num_levs; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          v_mid_h(0,ipack)[ivec] = 2*(k+1);
          dz_mid_h(0,ipack)[ivec] = 1;
        }
        Kokkos::deep_copy(v_mid,v_mid_h);
        Kokkos::deep_copy(dz_mid,dz_mid_h);

        // Use boundary conditions s.t. v_int becomes monotone increasing odd numbers
        const Real bc_top = 1;
        const Real bc_bot = 2*num_levs+1;

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto v_m = ekat::subview(v_mid,icol);
          auto dz  = ekat::subview(dz_mid,icol);

          column_ops::compute_interface_values_linear(team,num_levs,v_m,dz,bc_top,bc_bot,v_i);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == 2*k+1 );
        }
      }

      SECTION ("int_to_mid_to_int_linear") {

        // Fill v_mid with even numbers and set dz=1, so that v_int should contain odd ones.
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          v_int_h(0,ipack)[ivec] = 2*k+1;
          if (k < num_levs) dz_mid_h(0,ipack)[ivec] = 1;
        }
        Kokkos::deep_copy(v_int,v_int_h);
        Kokkos::deep_copy(dz_mid,dz_mid_h);

        // Use boundary conditions s.t. v_int becomes monotone increasing odd numbers
        const Real bc_top = v_int_h(0,0)[0];
        const Real bc_bot = v_int_h(0,num_levs/ps)[num_levs%ps];

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto v_m = ekat::subview(v_mid,icol);
          auto dz  = ekat::subview(dz_mid,icol);

          column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
          team.team_barrier();
          column_ops::compute_interface_values_linear(team,num_levs,v_m,dz,bc_top,bc_bot,v_i);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == 2*k+1);
        }
      }

      SECTION ("mid_to_int_compatible_fix_top") {

        // Fill v_mid with even numbers, so that v_int should contain odd ones.
        for (int k=0; k<num_levs; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          v_mid_h(0,ipack)[ivec] = 2*(k+1);
        }
        Kokkos::deep_copy(v_mid,v_mid_h);

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto v_m = ekat::subview(v_mid,icol);

          column_ops::compute_interface_values_compatible<true>(team,num_levs,v_m,1,v_i);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == 2*k+1 );
        }

        // Re-do with lambda as provider
        Kokkos::deep_copy(v_int,pack_type(0));
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          
          auto v_m = [&](const int k)->pack_type {
            return v_mid(icol,k);
          };
          auto v_i = ekat::subview(v_int,icol);

          column_ops::compute_interface_values_compatible<true>(team,num_levs,v_m,1,v_i);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == 2*k+1 );
        }
      }

      SECTION ("mid_to_int_compatible_fix_bot") {

        // Fill v_mid with even numbers, so that v_int should contain odd ones.
        for (int k=0; k<num_levs; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          v_mid_h(0,ipack)[ivec] = 2*(k+1);
        }
        Kokkos::deep_copy(v_mid,v_mid_h);

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto v_m = ekat::subview(v_mid,icol);

          column_ops::compute_interface_values_compatible<false>(team,num_levs,v_m,2*num_levs+1,v_i);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == 2*k+1 );
        }

        // Re-do with lambda as provider
        Kokkos::deep_copy(v_int,pack_type(0));
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          
          auto v_m = [&](const int k)->pack_type {
            return v_mid(icol,k);
          };
          auto v_i = ekat::subview(v_int,icol);

          column_ops::compute_interface_values_compatible<false>(team,num_levs,v_m,2*num_levs+1,v_i);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == 2*k+1 );
        }
      }

      SECTION ("int_to_mid_to_int_compatible_fix_top") {

        // Fill v_int with odd numbers, so that v_mid should contain even ones.
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          v_int_h(0,ipack)[ivec] = 2*k+1;
        }
        Kokkos::deep_copy(v_int,v_int_h);

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto v_m = ekat::subview(v_mid,icol);

          column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
          team.team_barrier();
          column_ops::compute_interface_values_compatible<true>(team,num_levs,v_m,1,v_i);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == 2*k+1 );
        }
      }

      SECTION ("int_to_mid_to_int_compatible_fix_bot") {

        // Fill v_int with odd numbers, so that v_mid should contain even ones.
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          v_int_h(0,ipack)[ivec] = 2*k+1;
        }
        Kokkos::deep_copy(v_int,v_int_h);

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto v_m = ekat::subview(v_mid,icol);

          column_ops::compute_midpoint_values(team,num_levs,v_i,v_m);
          team.team_barrier();
          column_ops::compute_interface_values_compatible<false>(team,num_levs,v_m,2*num_levs+1,v_i);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == 2*k+1 );
        }
      }

      SECTION ("delta") {

        // Fill v_int with odd numbers, so that dv_mid=2 on all levels
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          v_int_h(0,ipack)[ivec] = 2*k+1;
        }
        Kokkos::deep_copy(v_int,v_int_h);

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto dv_m = ekat::subview(dv_mid,icol);

          column_ops::compute_midpoint_delta(team,num_levs,v_i,dv_m);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(dv_mid_h,dv_mid);
        for (int k=0; k<num_levs; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (dv_mid_h(0,ipack)[ivec] == 2 );
        }

        // Re-do with lambda as provider
        Kokkos::deep_copy(dv_mid,pack_type(0));
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          
          auto v_i = [&](const int k)->pack_type {
            return v_int(icol,k);
          };
          auto dv_m = ekat::subview(dv_mid,icol);

          column_ops::compute_midpoint_delta(team,num_levs,v_i,dv_m);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(dv_mid_h,dv_mid);
        for (int k=0; k<num_levs; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (dv_mid_h(0,ipack)[ivec] == 2 );
        }
      }

      SECTION ("scan_from_top") {

        const Real s0 = -2.0; // If it worskfor s0!=0, it works for s0=0.

        // Fill v_mid with integers, then use Gauss formula to check the sum
        for (int k=0; k<num_levs; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          dv_mid_h(0,ipack)[ivec] = k+1;
        }
        Kokkos::deep_copy(dv_mid,dv_mid_h);

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto dv_m = ekat::subview(dv_mid,icol);

          column_ops::column_scan<true>(team,num_levs,dv_m,v_i,s0);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == s0+k*(k+1)/2.0);
        }

        // Re-do with lambda as provider
        Kokkos::deep_copy(v_int,pack_type(0));
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          
          auto dv_m = [&](const int k)->pack_type {
            return dv_mid(icol,k);
          };
          auto v_i = ekat::subview(v_int,icol);

          column_ops::column_scan<true>(team,num_levs,dv_m,v_i,s0);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const int ipack = k / ps;
          const int ivec  = k % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == s0+k*(k+1)/2.0);
        }
      }

      SECTION ("scan_from_bot") {

        const Real s0 = -2.0; // If it worskfor s0!=0, it works for s0=0.

        // Fill v_mid with integers, then use Gauss formula to check the sum
        for (int k=0; k<num_levs; ++k) {
          const auto k_bwd = num_levs - k - 1;
          const int ipack = k_bwd / ps;
          const int ivec  = k_bwd % ps;
          dv_mid_h(0,ipack)[ivec] = k+1;
        }
        Kokkos::deep_copy(dv_mid,dv_mid_h);

        // Run column kernel
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          auto v_i = ekat::subview(v_int,icol);
          auto dv_m = ekat::subview(dv_mid,icol);

          column_ops::column_scan<false>(team,num_levs,dv_m,v_i,s0);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const auto k_bwd = num_levs - k;
          const int ipack = k_bwd / ps;
          const int ivec  = k_bwd % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == s0+k*(k+1)/2.0);
        }

        // Re-do with lambda as provider
        Kokkos::deep_copy(v_int,pack_type(0));
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const member_type& team){
          const int icol = team.league_rank();
          
          auto dv_m = [&](const int k)->pack_type {
            return dv_mid(icol,k);
          };
          auto v_i = ekat::subview(v_int,icol);

          column_ops::column_scan<false>(team,num_levs,dv_m,v_i,s0);
        });
        Kokkos::fence();

        // Check answer
        Kokkos::deep_copy(v_int_h,v_int);
        for (int k=0; k<num_levs+1; ++k) {
          const auto k_bwd = num_levs - k;
          const int ipack = k_bwd / ps;
          const int ivec  = k_bwd % ps;
          REQUIRE (v_int_h(0,ipack)[ivec] == s0+k*(k+1)/2.0);
        }
      }
    }
  }
}

} // anonymous namespace
