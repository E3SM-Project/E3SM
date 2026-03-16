#include <catch2/catch.hpp>

#include "share/atm_process/atmosphere_diagnostic.hpp"

#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_input.hpp"

#include "share/data_managers/mesh_free_grids_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field.hpp"
#include "share/data_managers/field_manager.hpp"

#include "share/core/eamxx_setup_random_test.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/core/eamxx_types.hpp"

#include <ekat_units.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_assert.hpp>
#include <ekat_comm.hpp>

#include <iomanip>
#include <memory>

namespace scream {

// A multi-output diagnostic that takes two 3D input fields and produces:
//   "sum_f1_f2" = f1 + f2
//   "diff_f1_f2" = f1 - f2
//   "prod_f1_f2" = f1 * f2
// This exercises the multi-output infrastructure without Fortran dependencies.
class MultiOutDiag : public AtmosphereDiagnostic
{
public:
  MultiOutDiag (const ekat::Comm& comm, const ekat::ParameterList& params)
    : AtmosphereDiagnostic(comm,params)
  {
    // Nothing to do
  }

  std::string name() const override { return "MultiOutDiag"; }

  std::vector<std::string> get_diagnostic_names () const override {
    return {"sum_f1_f2", "diff_f1_f2", "prod_f1_f2"};
  }

  void create_requests () override {
    const auto grid = m_grids_manager->get_grid("point_grid");
    const auto& grid_name = grid->name();
    auto units = ekat::units::Units::nondimensional();
    auto layout = grid->get_3d_scalar_layout(true);

    add_field<Required>("f1", layout, units, grid_name);
    add_field<Required>("f2", layout, units, grid_name);

    // Allocate multi-output fields
    for (const auto& oname : get_diagnostic_names()) {
      Field f(FieldIdentifier(oname, layout, units, grid_name));
      f.allocate_view();
      m_diagnostic_outputs[oname] = f;
    }
  }

  int get_num_evaluations () const { return m_num_evaluations; }

protected:

  void compute_diagnostic_impl () override {
    const auto& f1 = get_field_in("f1");
    const auto& f2 = get_field_in("f2");

    auto f1_v = f1.get_view<const Real**>();
    auto f2_v = f2.get_view<const Real**>();
    auto sum_v  = m_diagnostic_outputs.at("sum_f1_f2").get_view<Real**>();
    auto diff_v = m_diagnostic_outputs.at("diff_f1_f2").get_view<Real**>();
    auto prod_v = m_diagnostic_outputs.at("prod_f1_f2").get_view<Real**>();

    const int ncols = f1.get_header().get_identifier().get_layout().dim(0);
    const int nlevs = f1.get_header().get_identifier().get_layout().dim(1);

    Kokkos::parallel_for(Kokkos::MDRangePolicy<Kokkos::Rank<2>>({0,0},{ncols,nlevs}),
      KOKKOS_LAMBDA(const int icol, const int ilev) {
        sum_v(icol,ilev)  = f1_v(icol,ilev) + f2_v(icol,ilev);
        diff_v(icol,ilev) = f1_v(icol,ilev) - f2_v(icol,ilev);
        prod_v(icol,ilev) = f1_v(icol,ilev) * f2_v(icol,ilev);
    });
    Kokkos::fence();

    ++m_num_evaluations;
  }

  int m_num_evaluations = 0;
};

// ---- Helpers ----

util::TimeStamp get_t0 () {
  return util::TimeStamp({2023,2,17},{0,0,0});
}

std::shared_ptr<const GridsManager>
get_gm (const ekat::Comm& comm)
{
  const int nlcols = 3;
  const int nlevs = 4;
  const int ngcols = nlcols*comm.size();
  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ngcols);
  gm->build_grids();
  return gm;
}

std::shared_ptr<FieldManager>
get_fm (const std::shared_ptr<const AbstractGrid>& grid,
        const util::TimeStamp& t0, int seed)
{
  using FL  = FieldLayout;
  using FID = FieldIdentifier;
  using namespace ShortFieldTagsNames;

  const int nlcols = grid->get_num_local_dofs();
  const int nlevs  = grid->get_num_vertical_levels();

  auto fm = std::make_shared<FieldManager>(grid);
  const auto units = ekat::units::Units::nondimensional();
  FL fl ({COL,LEV}, {nlcols,nlevs});

  // Create two input fields
  for (const auto& name : {"f1","f2"}) {
    FID fid(name, fl, units, grid->name());
    Field f(fid);
    f.allocate_view();
    randomize(f, seed++);
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }

  return fm;
}

// ---- Tests ----

TEST_CASE ("multi_output_diag_basic") {
  // Test the multi-output diagnostic interface directly (no IO)
  ekat::Comm comm(MPI_COMM_WORLD);
  auto t0 = get_t0();

  auto gm = get_gm(comm);
  auto grid = gm->get_grid("point_grid");

  auto seed = get_random_test_seed(&comm);
  auto fm = get_fm(grid, t0, seed);

  // Register and create diagnostic
  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("MultiOutDiag", &create_atmosphere_diagnostic<MultiOutDiag>);

  ekat::ParameterList params("multi_out_test");
  auto diag = diag_factory.create("MultiOutDiag", comm, params);
  diag->set_grids(gm);

  // Verify multi-output interface
  REQUIRE(diag->is_multi_output());
  auto names = diag->get_diagnostic_names();
  REQUIRE(names.size() == 3);
  REQUIRE(names[0] == "sum_f1_f2");
  REQUIRE(names[1] == "diff_f1_f2");
  REQUIRE(names[2] == "prod_f1_f2");

  // get_diagnostic() (no args) should throw for multi-output
  REQUIRE_THROWS(diag->get_diagnostic());

  // Set required fields
  diag->set_required_field(fm->get_field("f1").get_const());
  diag->set_required_field(fm->get_field("f2").get_const());

  // Initialize and compute
  diag->initialize(t0, RunType::Initial);
  diag->compute_diagnostic();

  // Verify all outputs have valid timestamps
  for (const auto& oname : names) {
    auto out = diag->get_diagnostic(oname);
    REQUIRE(out.is_allocated());
    REQUIRE(out.get_header().get_tracking().get_time_stamp().is_valid());
  }

  // Verify computed values
  auto f1_h = fm->get_field("f1").get_view<const Real**, Host>();
  auto f2_h = fm->get_field("f2").get_view<const Real**, Host>();

  diag->get_diagnostic("sum_f1_f2").sync_to_host();
  diag->get_diagnostic("diff_f1_f2").sync_to_host();
  diag->get_diagnostic("prod_f1_f2").sync_to_host();

  auto sum_h  = diag->get_diagnostic("sum_f1_f2").get_view<const Real**, Host>();
  auto diff_h = diag->get_diagnostic("diff_f1_f2").get_view<const Real**, Host>();
  auto prod_h = diag->get_diagnostic("prod_f1_f2").get_view<const Real**, Host>();

  const int nlcols = grid->get_num_local_dofs();
  const int nlevs  = grid->get_num_vertical_levels();

  for (int i = 0; i < nlcols; ++i) {
    for (int j = 0; j < nlevs; ++j) {
      REQUIRE(sum_h(i,j)  == f1_h(i,j) + f2_h(i,j));
      REQUIRE(diff_h(i,j) == f1_h(i,j) - f2_h(i,j));
      REQUIRE(prod_h(i,j) == f1_h(i,j) * f2_h(i,j));
    }
  }

  // Verify recomputation skips when timestamps haven't changed
  auto d = std::dynamic_pointer_cast<MultiOutDiag>(diag);
  REQUIRE(d->get_num_evaluations() == 1);
  diag->compute_diagnostic();
  REQUIRE(d->get_num_evaluations() == 1); // skipped, same timestamp

  // Update input timestamp and recompute
  fm->get_field("f1").get_header().get_tracking().update_time_stamp(t0 + 10.0);
  fm->get_field("f2").get_header().get_tracking().update_time_stamp(t0 + 10.0);
  diag->compute_diagnostic();
  REQUIRE(d->get_num_evaluations() == 2);

  // Invalid output name should throw
  REQUIRE_THROWS(diag->get_diagnostic("nonexistent_field"));

  diag->finalize();
}

TEST_CASE ("multi_output_diag_io") {
  // Test multi-output diagnostic through the full IO pipeline:
  // write output fields via OutputManager, then read back and verify
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto& diag_factory = AtmosphereDiagnosticFactory::instance();
  diag_factory.register_product("MultiOutDiag", &create_atmosphere_diagnostic<MultiOutDiag>);

  auto t0 = get_t0();
  auto seed = get_random_test_seed(&comm);

  auto gm = get_gm(comm);
  auto grid = gm->get_grid("point_grid");
  auto fm = get_fm(grid, t0, seed);

  // Request all multi-output fields plus a regular field
  std::vector<std::string> fnames = {"f1", "sum_f1_f2", "diff_f1_f2", "prod_f1_f2"};

  // Create output params
  ekat::ParameterList om_pl;
  om_pl.set("filename_prefix", std::string("io_multi_diags"));
  om_pl.set("field_names", fnames);
  om_pl.set("averaging_type", std::string("instant"));
  auto& ctrl_pl = om_pl.sublist("output_control");
  ctrl_pl.set("frequency_units", std::string("nsteps"));
  ctrl_pl.set("frequency", 1);
  ctrl_pl.set("save_grid_data", false);

  // Write
  {
    OutputManager om;
    om.initialize(comm, om_pl, t0, false);
    om.setup(fm, gm->get_grid_names());

    // Advance fields by one step
    const double dt = 10.0;
    auto t1 = t0 + dt;
    fm->get_field("f1").get_header().get_tracking().update_time_stamp(t1);
    fm->get_field("f2").get_header().get_tracking().update_time_stamp(t1);

    om.init_timestep(t0, dt);
    om.run(t1);
    om.finalize();
  }

  // Read back and verify
  {
    auto filename = "io_multi_diags.INSTANT.nsteps_x1.np"
                  + std::to_string(comm.size())
                  + "." + t0.to_string() + ".nc";

    // Create fields to read into
    const int nlcols = grid->get_num_local_dofs();
    const int nlevs = grid->get_num_vertical_levels();
    auto units = ekat::units::Units::nondimensional();
    FieldLayout fl({ShortFieldTagsNames::COL, ShortFieldTagsNames::LEV}, {nlcols, nlevs});

    auto read_fm = std::make_shared<FieldManager>(grid);
    for (const auto& name : fnames) {
      FieldIdentifier fid(name, fl, units, grid->name());
      Field f(fid);
      f.allocate_view();
      read_fm->add_field(f);
    }

    ekat::ParameterList reader_pl;
    reader_pl.set("filename", filename);
    reader_pl.set("field_names", fnames);
    AtmosphereInput reader(reader_pl, read_fm);
    reader.read_variables();

    // Verify: sum = f1 + f2, etc.
    auto f1_h = read_fm->get_field("f1").get_view<const Real**, Host>();
    auto sum_h = read_fm->get_field("sum_f1_f2").get_view<const Real**, Host>();
    auto diff_h = read_fm->get_field("diff_f1_f2").get_view<const Real**, Host>();
    auto prod_h = read_fm->get_field("prod_f1_f2").get_view<const Real**, Host>();

    for (int i = 0; i < nlcols; ++i) {
      for (int j = 0; j < nlevs; ++j) {
        // f2 was also in the FM when the diagnostic ran, so we can
        // reconstruct f2 = sum - f1
        auto f2_val = sum_h(i,j) - f1_h(i,j);
        REQUIRE(sum_h(i,j)  == f1_h(i,j) + f2_val);
        REQUIRE(diff_h(i,j) == f1_h(i,j) - f2_val);
        REQUIRE(prod_h(i,j) == f1_h(i,j) * f2_val);
      }
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
