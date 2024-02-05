#include "catch2/catch.hpp"

#include "physics/nudging/eamxx_nudging_process_interface.hpp"

#include "nudging_tests_helpers.hpp"

#include "share/field/field_utils.hpp"

using namespace scream;

std::shared_ptr<Nudging>
create_nudging (const ekat::Comm& comm,
                const ekat::ParameterList& params,
                const std::shared_ptr<FieldManager>& fm,
                const std::shared_ptr<GridsManager>& gm,
                const util::TimeStamp& t0)
{
  auto nudging = std::make_shared<Nudging>(comm,params);
  nudging->set_grids(gm);
  for (const auto& req : nudging->get_required_field_requests()) {
    auto f = fm->get_field(req.fid.name());
    nudging->set_required_field(f);
  }
  for (const auto& req : nudging->get_computed_field_requests()) {
    auto f = fm->get_field(req.fid.name());
    nudging->set_computed_field(f);
  }
  nudging->initialize(t0,RunType::Initial);

  return nudging;
}

TEST_CASE("nudging_tests") {
  using strvec_t = std::vector<std::string>;

  ekat::Comm comm(MPI_COMM_WORLD);

  // Init scorpio
  scorpio::eam_init_pio_subsystem(comm);

  const std::string nudging_data        = "nudging_data.INSTANT.nsteps_x1." + get_t0().to_string() + ".nc";
  const std::string nudging_data_filled = "nudging_data_filled.INSTANT.nsteps_x1." + get_t0().to_string() + ".nc";

  // A refined grid, with one extra node in between each of the coarse ones
  const int ngcols_fine = 2*ngcols_data - 1;
  const int nlevs_fine  = 2*nlevs_data -1;

  // For grids managers, depending on whether ncols/nlevs match the (coarse)
  // values used to generate the data or are finer
  auto gm_data   = create_gm (comm,ngcols_data,nlevs_data);
  auto gm_fine_h = create_gm (comm,ngcols_fine,nlevs_data);
  auto gm_fine_v = create_gm (comm,ngcols_data,nlevs_fine);
  auto gm_fine   = create_gm (comm,ngcols_fine,nlevs_fine);

  auto grid_data   = gm_data->get_grid("Point Grid");
  auto grid_fine_h = gm_fine_h->get_grid("Point Grid");
  auto grid_fine_v = gm_fine_v->get_grid("Point Grid");
  auto grid_fine   = gm_fine->get_grid("Point Grid");
  
  const int ncols_data = grid_data->get_num_local_dofs();
  // const int ncols_fine = grid_fine->get_num_local_dofs();

  // First section tests nudging when there is no horiz-vert interp
  SECTION ("no-horiz-no-vert") {
    ekat::ParameterList params;
    params.set<strvec_t>("nudging_filename",{nudging_data});
    params.set<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
    params.set<strvec_t>("nudging_fields",{"U","V"});
    params.get<std::string>("log_level","warn");

    // Create fm. Init p_mid, since it's constant in this file
    auto fm = create_fm(grid_data);
    update_field(fm->get_field("p_mid"),get_t0(),0);

    auto U = fm->get_field("U");
    auto V = fm->get_field("V");
    auto p = fm->get_field("p_mid");

    // Test case where model times coincide with input data times
    SECTION ("same-time") {
      std::cout << " -> Testing same time/horiz/vert grid as data ...........\n";

      // Create and init nudging process
      auto nudging = create_nudging(comm,params,fm,gm_data,get_t0());

      // Same space-time grid, should return the same data
      auto fm_tgt = create_fm(grid_data);
      auto time = get_t0();

      auto U_tgt = fm_tgt->get_field("U");
      auto V_tgt = fm_tgt->get_field("V");
      for (int n=0; n<nsteps_data; ++n) {
        time += dt_data;

        // Run nudging
        nudging->run(dt_data);

        // Recompute original data
        update_fields(fm_tgt,time,0);

        // Since all values are integers, we should have no rounding
        REQUIRE (views_are_equal(U,U_tgt));
        REQUIRE (views_are_equal(V,V_tgt));
      }
      std::cout << " -> Testing same time/horiz/vert grid as data ........... PASS\n";
    }

    // Test case where model times are in the middle of input data time intervals
    SECTION ("half-time") {
      std::cout << " -> Testing same horiz/vert grid, different time grid ...\n";

      // Init time as t0-dt/2, so we're "half way" between data slices
      auto time = get_t0() - dt_data/2;

      // Create and init nudging process
      auto nudging = create_nudging(comm,params,fm,gm_data,time);

      auto tmp1 = U.clone("");
      auto tmp2 = U.clone("");

      auto check_f = [&](const Field& f,
                         const util::TimeStamp& t_prev,
                         const util::TimeStamp& t_next) {
        update_field(tmp1,t_prev,0);
        update_field(tmp2,t_next,0);
        tmp1.update(tmp2,0.5,0.5);

        // Since input data are integers, and the time-interp coeff is 0.5,
        // we should be getting the exact answer
        REQUIRE (views_are_equal(f,tmp1));
      };

      auto t_prev = get_t0();
      for (int n=1; n<nsteps_data; ++n) {
        auto t_next = t_prev+dt_data;
        time += dt_data;

        // Run nudging
        nudging->run(dt_data);

        // Compare the two. Since we're exactly half way, we should get exact fp representation
        check_f(fm->get_field("U"),t_prev,t_next);

        t_prev = t_next;
      }
      std::cout << " -> Testing same horiz/vert grid, different time grid ... PASS\n";
    }
  }

  // Now test the case where we do have vertical interp.
  SECTION ("no-horiz-yes-vert") {
    std::cout << " -> Testing same time/horiz grid, different vert grid ...\n";
    const auto Pa = ekat::units::Pa;

    // Helper lambda, to compute f on the "fine" vert grid from f on the data vert grid
    auto manual_interp = [&](const Field& data, const Field& fine) {
      auto fine_h = fine.get_view<Real**,Host>();
      auto data_h = data.get_view<Real**,Host>();
      for (int icol=0; icol<ncols_data; ++icol) {
        // Even entries match original data
        for (int ilev=0; ilev<nlevs_data; ++ilev) {
          fine_h(icol,2*ilev) = data_h(icol,ilev);
        }
        // Odd entries are avg of the two adjacent even entries
        for (int ilev=0; ilev<nlevs_data-1; ++ilev) {
          fine_h(icol,2*ilev+1) = (fine_h(icol,2*ilev)+fine_h(icol,2*ilev+2))/2;
        }
      }
      fine.sync_to_dev();
    };

    ekat::ParameterList params;
    params.set<strvec_t>("nudging_filename",{nudging_data});
    params.set<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
    params.set<strvec_t>("nudging_fields",{"U","V"});
    params.get<std::string>("log_level","warn");

    // Create fm
    auto fm = create_fm(grid_fine_v);
    auto U = fm->get_field("U");
    auto V = fm->get_field("V");
    auto p_mid = fm->get_field("p_mid");

    // Create and init nudging process
    auto nudging = create_nudging(comm,params,fm,gm_fine_v,get_t0());

    // Compute pmid on data grid
    auto layout_data = grid_data->get_3d_scalar_layout(true);
    Field p_mid_data(FieldIdentifier("p_mid",layout_data,Pa,grid_data->name()));
    p_mid_data.allocate_view();
    update_field(p_mid_data,get_t0(),0);

    manual_interp(p_mid_data,p_mid);

    auto time = get_t0();
    Field tmp_data = p_mid_data.clone("tmp data");
    Field tmp_fine = p_mid.clone("tmp fine");
    for (int n=0; n<nsteps_data; ++n) {
      // Run nudging
      nudging->run(dt_data);

      // Compute data on fine grid, by manually interpolating
      // (recall that nudging runs at t+dt)
      update_field(tmp_data,time+dt_data,0);
      manual_interp(tmp_data,tmp_fine);

      REQUIRE (views_are_equal(tmp_fine,fm->get_field("U")));
      time += dt_data;
    }
    std::cout << " -> Testing same time/horiz grid, different vert grid ... PASS\n";
  }

  // Clean up scorpio
  scorpio::eam_pio_finalize();
}
