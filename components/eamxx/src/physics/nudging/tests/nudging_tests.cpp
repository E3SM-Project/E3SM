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
  auto& catch_capture = Catch::getResultCapture();

  using strvec_t = std::vector<std::string>;

  ekat::Comm comm(MPI_COMM_WORLD);

  auto root_print = [&](const std::string& msg) {
    if (comm.am_i_root()) {
      printf("%s",msg.c_str());
    }
  };

  // Init scorpio
  scorpio::eam_init_pio_subsystem(comm);

  // A refined grid, with one extra node in between each of the coarse ones
  const int ngcols_fine = 2*ngcols_data - 1;
  const int nlevs_fine  = 2*nlevs_data -1;

  // Files names
  auto postfix = ".INSTANT.nsteps_x1." + get_t0().to_string() + ".nc";
  auto nudging_data        = "nudging_data" + postfix;
  auto nudging_data_filled = "nudging_data_filled" + postfix;
  auto map_file = "map_ncol" + std::to_string(ngcols_data)
                + "_to_"     + std::to_string(ngcols_fine) + ".nc";

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

  // First section tests nudging when there is no horiz-vert interp
  SECTION ("no-horiz-no-vert") {
    ekat::ParameterList params;
    params.set<strvec_t>("nudging_filenames_patterns",{nudging_data});
    params.set<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
    params.set<strvec_t>("nudging_fields",{"U"});
    params.get<std::string>("log_level","warn");

    // Create fm. Init p_mid, since it's constant in this file
    auto fm = create_fm(grid_data);
    compute_field(fm->get_field("p_mid"),get_t0(),comm,0);

    auto U = fm->get_field("U");
    SECTION ("same-time") {
      std::string msg = " -> Testing same time/horiz/vert grid as data ...........";
      root_print (msg + "\n");
      bool ok = true;

      // Create and init nudging process
      auto nudging = create_nudging(comm,params,fm,gm_data,get_t0());

      auto time = get_t0();

      auto U_tgt = U.clone("U_tgt");
      for (int n=0; ok and n<nsteps_data; ++n) {
        time += dt_data;

        // Run nudging
        nudging->run(dt_data);

        // Recompute original data
        compute_field(U_tgt,time,comm,0);

        // Since all values are integers, we should have no rounding
        CHECK (views_are_equal(U,U_tgt));
        ok &= catch_capture.lastAssertionPassed();
      }
      root_print (msg + (ok ? " PASS\n" : " FAIL\n"));
    }

    // Test case where model times are in the middle of input data time intervals
    SECTION ("half-time") {
      std::string msg = " -> Testing same horiz/vert grid, different time grid ...";
      root_print (msg + "\n");
      bool ok = true;

      // Init time as t0-dt/2, so we're "half way" between data slices
      auto time = get_t0() - dt_data/2;

      // Create and init nudging process
      auto nudging = create_nudging(comm,params,fm,gm_data,time);

      auto tmp1 = U.clone("");
      auto tmp2 = U.clone("");

      auto check_f = [&](const Field& f,
                         const util::TimeStamp& t_prev,
                         const util::TimeStamp& t_next) {
        compute_field(tmp1,t_prev,comm,0);
        compute_field(tmp2,t_next,comm,0);
        tmp1.update(tmp2,0.5,0.5);

        // Since input data are integers, and the time-interp coeff is 0.5,
        // we should be getting the exact answer
        CHECK (views_are_equal(f,tmp1));
        ok &= catch_capture.lastAssertionPassed();
      };

      auto t_prev = get_t0();
      for (int n=1; ok and n<nsteps_data; ++n) {
        auto t_next = t_prev+dt_data;
        time += dt_data;

        // Run nudging
        nudging->run(dt_data);

        // Compare the two. Since we're exactly half way, we should get exact fp representation
        check_f(U,t_prev,t_next);

        t_prev = t_next;
      }
      root_print (msg + (ok ? " PASS\n" : " FAIL\n"));
    }
  }

  // Now test the case where we do have vertical interp.
  SECTION ("no-horiz-yes-vert") {
    const auto Pa = ekat::units::Pa;

    // Helper lambda, to compute f on the "fine" vert grid from f on the data vert grid
    // If in_bounds=false, top/bot entries are extrapolated:
    //  top: f_out(0) = f_in(1) / 2
    //  bot: f_out(bot) = f_in(bot-1)

    auto manual_interp = [&](const Field& data, const Field& fine, const bool in_bounds) {
      auto fine_h = fine.get_view<Real**,Host>();
      auto data_h = data.get_view<Real**,Host>();
      const bool is_pmid = data.name()=="p_mid";
      for (int icol=0; icol<ncols_data; ++icol) {
        // Even entries match original data
        for (int ilev=0; ilev<nlevs_data; ++ilev) {
          fine_h(icol,2*ilev) = data_h(icol,ilev);
        }
        // Odd entries are avg of the two adjacent even entries
        for (int ilev=0; ilev<nlevs_data-1; ++ilev) {
          fine_h(icol,2*ilev+1) = (fine_h(icol,2*ilev)+fine_h(icol,2*ilev+2))/2;
        }
        if (not in_bounds) {
          const int top = 0;
          const int bot = 2*nlevs_data - 1;
          fine_h(icol,top) *= 0.5;
          fine_h(icol,bot) *= is_pmid ? 2 : 1;
        }
      }
      fine.sync_to_dev();
    };

    ekat::ParameterList params;
    params.set<strvec_t>("nudging_filenames_patterns",{nudging_data});
    params.set<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
    params.set<strvec_t>("nudging_fields",{"U"});
    params.get<std::string>("log_level","warn");

    std::string msg = " -> Testing same time/horiz grid, different vert grid";
    root_print (msg + "\n");
    SECTION ("pmid_in_bounds") {
      std::string msg = "   -> Target pressure within source pressure bounds .....";
      root_print (msg + "\n");
      bool ok = true;

      // Create fm
      auto fm = create_fm(grid_fine_v);
      auto U = fm->get_field("U");
      auto p_mid = fm->get_field("p_mid");

      // Create and init nudging process
      auto nudging = create_nudging(comm,params,fm,gm_fine_v,get_t0());

      // Compute pmid on data grid
      auto layout_data = grid_data->get_3d_scalar_layout(true);
      Field p_mid_data(FieldIdentifier("p_mid",layout_data,Pa,grid_data->name()));
      p_mid_data.allocate_view();
      compute_field(p_mid_data,get_t0(),comm,0);

      manual_interp(p_mid_data,p_mid,true);

      auto time = get_t0();
      Field tmp_data = p_mid_data.clone("tmp data");
      Field tmp_fine = p_mid.clone("tmp fine");
      for (int n=0; ok and n<nsteps_data; ++n) {
        // Run nudging
        nudging->run(dt_data);

        // Compute data on fine grid, by manually interpolating
        // (recall that nudging runs at t+dt)
        compute_field(tmp_data,time+dt_data,comm,0);
        manual_interp(tmp_data,tmp_fine,true);

        CHECK (views_are_equal(tmp_fine,U));
        ok &= catch_capture.lastAssertionPassed();
        time += dt_data;
      }
      root_print (msg + (ok ? " PASS\n" : " FAIL\n"));
    }

    SECTION ("pmid_out_of_bounds") {
      std::string msg = "   -> Target pressure outside source pressure bounds ....";
      root_print (msg + "\n");
      bool ok = true;

      // Create fm
      auto fm = create_fm(grid_fine_v);
      auto U = fm->get_field("U");
      auto p_mid = fm->get_field("p_mid");

      // Create and init nudging process
      auto nudging = create_nudging(comm,params,fm,gm_fine_v,get_t0());

      // Compute pmid on data grid
      auto layout_data = grid_data->get_3d_scalar_layout(true);
      Field p_mid_data(FieldIdentifier("p_mid",layout_data,Pa,grid_data->name()));
      p_mid_data.allocate_view();
      compute_field(p_mid_data,get_t0(),comm,0);

      manual_interp(p_mid_data,p_mid,false);

      auto time = get_t0();
      Field tmp_data = p_mid_data.clone("tmp data");
      Field tmp_fine = p_mid.clone("tmp fine");
      for (int n=0; ok and n<nsteps_data; ++n) {
        // Run nudging
        nudging->run(dt_data);

        // Compute data on fine grid, by manually interpolating
        // (recall that nudging runs at t+dt)
        compute_field(tmp_data,time+dt_data,comm,0);
        manual_interp(tmp_data,tmp_fine,false);

        CHECK (views_are_equal(tmp_fine,U));
        ok &= catch_capture.lastAssertionPassed();
        time += dt_data;
      }
      root_print (msg + (ok ? " PASS\n" : " FAIL\n"));
    }
  }

  SECTION ("horiz-remap") {
    std::string msg = " -> Testing different horiz grids .......................";
    root_print (msg + "\n");
    bool ok = true;

    const auto Pa = ekat::units::Pa;

    // Helper lambda, to compute f on the "fine" horiz grid from f on the data
    // Recall that the fine grid is a 1d grid, with all coarse grid cols and
    // a new col in between each two coarse grids. The mapping weights are
    // such that coarse cols are copied, and the new cols are the avg of the
    // two coarse cols left and right.

    auto create_global_f = [&] (const Field& f) {
      const auto& fid = f.get_header().get_identifier();
      int ncols = fid.get_layout().dim(0);
      comm.all_reduce(&ncols,1,MPI_SUM);

      FieldLayout glb_layout = fid.get_layout().clone().reset_dim(0,ncols);
      FieldIdentifier glb_fid(fid.name(),glb_layout,fid.get_units(),fid.get_grid_name());
      Field glb(glb_fid);
      glb.allocate_view();
      return glb;
    };
    auto gather_global_f = [&] (const Field& f) {
      const auto& fid = f.get_header().get_identifier();
      const int ncols = fid.get_layout().dim(0);
      const int nlevs = fid.get_layout().dim(1);

      auto glb = create_global_f(f);
      f.sync_to_host();

      int offset=0;
      auto view_h = f.get_view<Real**,Host>();
      auto glb_view_h = glb.get_view<Real**,Host>();
      for (int rank=0; rank<comm.size(); ++rank) {
        if (rank==comm.rank()) {
          for (int i=0; i<ncols; ++i) {
            auto col_h = ekat::subview(view_h,i);
            auto glb_col_h = ekat::subview(glb_view_h,i+offset);
            std::copy(col_h.data(),col_h.data()+nlevs,glb_col_h.data());
          }
        }
        auto data = glb_view_h.data()+offset*nlevs;
        int rank_size = view_h.size();
        comm.broadcast(&rank_size,1,rank);
        comm.broadcast(data,rank_size,rank);

        int rank_ncols = ncols;
        comm.broadcast(&rank_ncols,1,rank);
        offset += rank_ncols;

        glb.sync_to_dev();
      }
      glb.sync_to_dev();

      return glb;
    };

    auto manual_interp = [&](const Field& data, const Field& fine) {
      auto glb_fine = create_global_f(fine);
      auto glb_data = gather_global_f(data);
      auto glb_fine_h = glb_fine.get_view<Real**,Host>();
      auto glb_data_h = glb_data.get_view<Real**,Host>();
      for (int icol=0; icol<ngcols_data; ++icol) {
        for (int ilev=0; ilev<nlevs_data; ++ilev) {
          glb_fine_h(icol,ilev) = glb_data_h(icol,ilev);
          if (icol<ngcols_data-1) {
            glb_fine_h(ngcols_data+icol,ilev) = (glb_data_h(icol,ilev) + glb_data_h(icol+1,ilev))/2;
          }
        }
      }

      auto ncols_fine = grid_fine->get_num_local_dofs();
      auto offset = ncols_fine;
      comm.scan(&offset,1,MPI_SUM);
      offset -= ncols_fine;

      auto fine_h = fine.get_view<Real**,Host>();
      for (int i=0; i<ncols_fine; ++i) {
        for (int k=0; k<nlevs_data; ++k) {
          fine_h(i,k) = glb_fine_h(offset+i,k);
        }
      }
    };

    ekat::ParameterList params;
    params.set<strvec_t>("nudging_filenames_patterns",{nudging_data});
    params.set<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
    params.set<std::string>("nudging_refine_remap_mapfile",map_file);
    params.set<strvec_t>("nudging_fields",{"U"});
    params.get<std::string>("log_level","warn");

    // Create fm
    auto fm = create_fm(grid_fine_h);
    auto U = fm->get_field("U");
    auto p_mid = fm->get_field("p_mid");

    // Create and init nudging process
    auto nudging = create_nudging(comm,params,fm,gm_fine_h,get_t0());

    // Compute pmid on data grid
    auto layout_data = grid_data->get_3d_scalar_layout(true);
    Field p_mid_data(FieldIdentifier("p_mid",layout_data,Pa,grid_data->name()));
    p_mid_data.allocate_view();
    compute_field(p_mid_data,get_t0(),comm,0);

    manual_interp(p_mid_data,p_mid);

    auto time = get_t0();
    Field tmp_data = p_mid_data.clone("tmp data");
    Field tmp_fine = p_mid.clone("tmp fine");
    for (int n=0; ok and n<nsteps_data; ++n) {
      // Run nudging
      nudging->run(dt_data);

      // Compute data on fine grid, by manually interpolating
      // (recall that nudging runs at t+dt)
      compute_field(tmp_data,time+dt_data,comm,0);
      manual_interp(tmp_data,tmp_fine);

      CHECK (views_are_equal(tmp_fine,U));
      ok &= catch_capture.lastAssertionPassed();
      time += dt_data;
    }
    root_print (msg + (ok ? " PASS\n" : " FAIL\n"));
  }

  SECTION ("filled-data") {
    std::string msg = " -> Testing data with top/bot levels filled .............";
    root_print (msg + "\n");
    bool ok = true;

    // Helper lambda, to manually cure filled levels
    auto manual_cure = [&](const Field& f) {
      auto f_h = f.get_view<Real**,Host>();
      auto ncols = f.get_header().get_identifier().get_layout().dim(0);
      auto nlevs = f.get_header().get_identifier().get_layout().dim(1);
      auto first_good = nlevs_filled;
      auto last_good  = nlevs - nlevs_filled - 1;
      for (int icol=0; icol<ncols; ++icol) {
        for (int ilev=0; ilev<nlevs_filled; ++ilev) {
          f_h(icol,ilev) = f_h(icol,first_good);
        }
        for (int ilev=last_good+1; ilev<nlevs; ++ilev) {
          f_h(icol,ilev) = f_h(icol,last_good);
        }
      }
      f.sync_to_dev();
    };

    ekat::ParameterList params;
    params.set<strvec_t>("nudging_filenames_patterns",{nudging_data_filled});
    params.set<std::string>("source_pressure_type","TIME_DEPENDENT_3D_PROFILE");
    params.set<strvec_t>("nudging_fields",{"U"});
    params.get<std::string>("log_level","warn");

    // Create fm. Init p_mid, since it's constant in this file
    auto fm = create_fm(grid_data);
    auto U = fm->get_field("U");
    auto p_mid = fm->get_field("p_mid");
    compute_field(p_mid,get_t0(),comm,0);

    // Create and init nudging process
    auto nudging = create_nudging(comm,params,fm,gm_data,get_t0());

    auto time = get_t0();
    Field tmp = p_mid.clone("tmp");
    for (int n=0; ok and n<nsteps_data; ++n) {
      // Run nudging
      nudging->run(dt_data);

      // Compute data on fine grid, by manually interpolating
      // (recall that nudging runs at t+dt)
      compute_field(tmp,time+dt_data,comm,0);
      manual_cure(tmp);

      CHECK (views_are_equal(tmp,U));
      ok &= catch_capture.lastAssertionPassed();
      time += dt_data;
    }
    root_print (msg + (ok ? " PASS\n" : " FAIL\n"));
  }

  // Clean up scorpio
  scorpio::eam_pio_finalize();
}
