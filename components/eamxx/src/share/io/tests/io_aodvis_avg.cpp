#include <catch2/catch.hpp>

#include "share/data_managers/field_manager.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field.hpp"
#include "share/field/field_reader.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/util/eamxx_utils.hpp"

#include <ekat_comm.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_units.hpp>

#include <fstream>
#include <string>
#include <vector>

/*
 * AODVis is only valid in sunlit columns. Check that its average is taken over
 * the sunlit steps only, and that a column that is never sunlit comes out as
 * fill_value.
 */

namespace scream {

namespace {

constexpr int ncols = 4;
constexpr int nlevs = 3;
constexpr int nsteps = 4;

constexpr Real FILL = constants::fill_value<Real>;

// Column g is sunlit on step n iff (n+g) is even, except the last column, which
// never is. tau=n in every band and level, so aod=n*nlevs on a sunlit step.
bool sunlit_on (const int n, const int g) { return g<ncols-1 and (n+g)%2==0; }

void calc_fields (Field& tau, Field& sunlit, const AbstractGrid& grid, const int n)
{
  using gid_type = AbstractGrid::gid_type;
  auto gids_h = grid.get_dofs_gids().get_view<const gid_type*,Host>();
  auto tau_h  = tau.get_view<Real***,Host>();
  auto sun_h  = sunlit.get_view<int*,Host>();
  const int nl = grid.get_num_local_dofs();
  for (int i=0; i<nl; ++i) {
    sun_h(i) = sunlit_on(n,gids_h(i)) ? 1 : 0;
    for (int b=0; b<eamxx_swbands(); ++b) {
      for (int k=0; k<nlevs; ++k) {
        tau_h(i,b,k) = n;
      }
    }
  }
  tau.sync_to_dev();
  sunlit.sync_to_dev();
}

} // anonymous namespace

TEST_CASE ("averaged_aodvis")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  register_diagnostics();

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ncols);
  gm->build_grids();
  auto grid = gm->get_grid("point_grid");
  const auto gname = grid->name();

  util::TimeStamp t0({2023,1,1},{0,0,0});

  auto fm = std::make_shared<FieldManager>(gm);
  Field tau(FieldIdentifier("aero_tau_sw",
                            grid->get_3d_vector_layout(LEV,eamxx_swbands(),"swband"),
                            none,gname));
  Field sunlit(FieldIdentifier("sunlit_mask",grid->get_2d_scalar_layout(),
                               none,gname,DataType::IntType));
  for (auto* f : {&tau,&sunlit}) {
    f->allocate_view();
    f->get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(*f);
  }
  calc_fields(tau,sunlit,*grid,0);

  const std::string prefix = "io_aodvis_avg";
  const std::string aod_name = "AerosolOpticalDepth550nm";
  const int dt = 1;

  ekat::ParameterList params;
  params.set<std::string>("filename_prefix",prefix);
  params.set<std::string>("averaging_type","AVERAGE");
  params.set<std::string>("floating_point_precision","real");
  auto& f_pl = params.sublist("fields").sublist(gname);
  f_pl.set<std::vector<std::string>>("field_names",{aod_name});
  auto& ctrl_pl = params.sublist("output_control");
  ctrl_pl.set<std::string>("frequency_units","nsteps");
  ctrl_pl.set<int>("frequency",nsteps);
  ctrl_pl.set<bool>("save_grid_data",false);

  auto t = t0;
  {
    OutputManager om;
    om.initialize(comm,params,t0,false);
    om.setup(fm,gm->get_grid_names());

    for (int n=1; n<=nsteps; ++n) {
      om.init_timestep(t,dt);
      t += dt;
      auto tau_f = fm->get_field("aero_tau_sw");
      auto sun_f = fm->get_field("sunlit_mask");
      calc_fields(tau_f,sun_f,*grid,n);
      for (auto f : {tau_f,sun_f}) {
        f.get_header().get_tracking().update_time_stamp(t);
      }
      om.run(t);
    }
    om.finalize();
  }

  const auto filename = prefix + ".AVERAGE.nsteps_x" + std::to_string(nsteps) +
                        ".np" + std::to_string(comm.size()) + "." +
                        t0.to_string() + ".nc";
  std::ifstream file_check(filename);
  REQUIRE (file_check.good());
  file_check.close();

  scorpio::register_file(filename,scorpio::Read);
  REQUIRE (scorpio::has_var(filename,aod_name));
  REQUIRE (scorpio::has_var(filename,"avg_count_" + aod_name + "_ncol"));
  CHECK (scorpio::has_attribute(filename,aod_name,"_FillValue"));
  scorpio::release_file(filename);

  {
    Field aod(FieldIdentifier(aod_name,grid->get_2d_scalar_layout(),none,gname));
    aod.allocate_view();
    aod.get_header().get_tracking().update_time_stamp(t0);

    FieldReader reader;
    reader.set_file_specs(filename);
    reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
    reader.set_fields({aod});
    reader.read(0);
    aod.sync_to_host();
    auto aod_h = aod.get_view<const Real*,Host>();

    using gid_type = AbstractGrid::gid_type;
    auto gids_h = grid->get_dofs_gids().get_view<const gid_type*,Host>();
    const int nlocal = grid->get_num_local_dofs();
    bool saw_partial = false, saw_empty = false;
    for (int i=0; i<nlocal; ++i) {
      const int g = gids_h(i);
      Real sum = 0; int cnt = 0;
      for (int n=1; n<=nsteps; ++n) {
        if (sunlit_on(n,g)) { sum += n*nlevs; ++cnt; }
      }
      if (cnt==0) {
        REQUIRE (aod_h(i)==FILL);
        saw_empty = true;
      } else {
        REQUIRE (aod_h(i)==Approx(sum/cnt));
        REQUIRE (cnt<nsteps);
        saw_partial = true;
      }
    }
    // Make sure both branches were exercised
    if (comm.size()==1) {
      REQUIRE (saw_partial);
      REQUIRE (saw_empty);
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
