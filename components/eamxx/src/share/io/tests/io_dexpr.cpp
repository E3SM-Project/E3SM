#include <catch2/catch.hpp>

#include "share/data_managers/field_manager.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/diagnostics/register_diagnostics.hpp"
#include "share/field/field.hpp"
#include "share/field/field_reader.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_comm.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_units.hpp>

#include <fstream>
#include <string>
#include <vector>

/*
 * End-to-end checks for names resolved by the expression front end.
 *
 * create_diag.cpp covers the translation. Only here can we check that the field
 * survives an output stream, which looks fields up by the REQUESTED name while
 * a diag derives its own from its params -- not the same thing for an
 * expression, nor for the _atm_backtend alias.
 */

namespace scream {

namespace {

constexpr int ncols = 6;
constexpr int nlevs = 4;

// qv[i][k] = n + i*(k+1), T_mid[i][k] = n - i*(k+1): simple enough that
// expected values can be written down by hand.
void calc_fields (Field& qv, Field& T_mid, const int n)
{
  auto qv_h = qv.get_view<Real**,Host>();
  auto T_h  = T_mid.get_view<Real**,Host>();
  const int nl = qv.get_header().get_identifier().get_layout().dim(0);
  for (int i=0; i<nl; ++i) {
    for (int k=0; k<nlevs; ++k) {
      qv_h(i,k) = n + i*(k+1);
      T_h (i,k) = n - i*(k+1);
    }
  }
  qv.sync_to_dev();
  T_mid.sync_to_dev();
}

// NOTE: built on the grids manager, not the grid: the 'fields->$grid' output
//       syntax needs a non-const grid.
std::shared_ptr<FieldManager>
create_test_fm (const std::shared_ptr<const GridsManager>& gm,
                const std::shared_ptr<const AbstractGrid>& grid,
                const util::TimeStamp& t0)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto layout3d = grid->get_3d_scalar_layout(LEV);
  auto fm = std::make_shared<FieldManager>(gm);

  Field qv    (FieldIdentifier("qv",   layout3d, kg/kg, grid->name()));
  Field T_mid (FieldIdentifier("T_mid",layout3d, K,     grid->name()));
  qv.allocate_view();
  T_mid.allocate_view();
  qv.get_header().get_tracking().update_time_stamp(t0);
  T_mid.get_header().get_tracking().update_time_stamp(t0);

  fm->add_field(qv);
  fm->add_field(T_mid);

  calc_fields(qv,T_mid,0);

  return fm;
}

// NOTE: the 'aliases' section only exists in the full fields->$grid syntax
ekat::ParameterList output_params (const std::string& prefix,
                                   const std::string& grid_name,
                                   const std::vector<std::string>& field_names,
                                   const std::vector<std::string>& aliases = {})
{
  ekat::ParameterList params;
  params.set<std::string>("filename_prefix",prefix);
  params.set<std::string>("averaging_type","INSTANT");
  params.set<std::string>("floating_point_precision","real");

  auto& f_pl = params.sublist("fields").sublist(grid_name);
  f_pl.set("field_names",field_names);
  if (not aliases.empty()) {
    f_pl.set("aliases",aliases);
  }

  auto& ctrl_pl = params.sublist("output_control");
  ctrl_pl.set<std::string>("frequency_units","nsteps");
  ctrl_pl.set<int>("frequency",1);
  ctrl_pl.set<bool>("save_grid_data",false);
  // No snapshot at t0: a tendency has no meaning before the first step
  ctrl_pl.set<bool>("skip_t0_output",true);

  return params;
}

} // anonymous namespace

TEST_CASE ("io_with_expressions")
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
  auto fm = create_test_fm(gm,grid,t0);

  const std::string prefix = "io_dexpr";
  const int dt = 1;
  const int nsteps = 3;

  auto params = output_params(prefix,gname,
      {
        // An expression must be given a name, since field names go into the
        // nc file verbatim
        "prod := qv*T_mid",
        // Grouping is explicit, unlike in the underscore spelling, where it
        // falls out of regex greediness
        "shifted := (qv+qv)*T_mid",
        // Built on the intermediate below, by name
        "half := qvT/2",
        // A method call
        "col_qv := qv.mean('lev')",
        // Regression: expands to T_mid_minus_T_mid_prev_over_dt, so the diag
        // names its field after the expansion. Without an explicit output
        // name the stream throws during setup.
        "T_mid_atm_backtend",
      },
      {
        // An intermediate: resolved and usable by name, but not written
        "qvT := qv*T_mid",
      });

  // NOTE: scoped, so the stream (and the diags it owns) are gone before Kokkos
  //       is finalized
  auto t = t0;
  {
  OutputManager om;
  om.initialize(comm,params,t0,false);
  om.setup(fm,gm->get_grid_names());

  for (int n=1; n<=nsteps; ++n) {
    om.init_timestep(t,dt);
    t += dt;

    auto qv    = fm->get_field("qv");
    auto T_mid = fm->get_field("T_mid");
    calc_fields(qv,T_mid,n);
    qv.get_header().get_tracking().update_time_stamp(t);
    T_mid.get_header().get_tracking().update_time_stamp(t);

    om.run(t);
  }
  om.finalize();
  }

  // NOTE: the file is named after the FIRST snapshot in it, which is t0+dt
  //       here, since we skip the t0 output
  const auto filename = prefix + ".INSTANT.nsteps_x1.np" +
                        std::to_string(comm.size()) + "." + (t0+dt).to_string() + ".nc";
  std::ifstream file_check(filename);
  REQUIRE (file_check.good());
  file_check.close();

  scorpio::register_file(filename,scorpio::Read);
  REQUIRE (scorpio::has_var(filename,"prod"));
  REQUIRE (scorpio::has_var(filename,"shifted"));
  REQUIRE (scorpio::has_var(filename,"half"));
  REQUIRE (scorpio::has_var(filename,"col_qv"));
  // The one that used to throw during setup
  REQUIRE (scorpio::has_var(filename,"T_mid_atm_backtend"));
  // An 'aliases' entry is an intermediate: resolved, used, but not written
  REQUIRE_FALSE (scorpio::has_var(filename,"qvT"));
  // The file records where the numbers came from
  REQUIRE (scorpio::get_attribute<std::string>(filename,"prod","alias_of")=="qv*T_mid");
  REQUIRE (scorpio::get_attribute<std::string>(filename,"half","alias_of")=="qvT/2");
  scorpio::release_file(filename);

  // Values
  // NOTE: scoped, so the fields are gone before Kokkos is finalized
  {
  const auto layout3d = grid->get_3d_scalar_layout(LEV);
  const auto layout2d = grid->get_2d_scalar_layout();
  Field prod    (FieldIdentifier("prod",   layout3d, kg/kg*K, gname));
  Field shifted (FieldIdentifier("shifted",layout3d, kg/kg*K, gname));
  Field half    (FieldIdentifier("half",   layout3d, kg/kg*K, gname));
  Field col_qv  (FieldIdentifier("col_qv", layout2d, kg/kg,   gname));
  Field backtend(FieldIdentifier("T_mid_atm_backtend",layout3d, K/s, gname));
  for (auto* f : {&prod,&shifted,&half,&col_qv,&backtend}) {
    f->allocate_view();
    f->get_header().get_tracking().update_time_stamp(t0);
  }

  FieldReader reader;
  reader.set_file_specs(filename);
  reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
  reader.set_fields({prod,shifted,half,col_qv,backtend});

  const int nlocal = grid->get_num_local_dofs();
  for (int n=1; n<=nsteps; ++n) {
    reader.read(n-1);
    prod.sync_to_host();
    shifted.sync_to_host();
    half.sync_to_host();
    col_qv.sync_to_host();
    backtend.sync_to_host();
    auto prod_h     = prod.get_view<const Real**,Host>();
    auto shifted_h  = shifted.get_view<const Real**,Host>();
    auto half_h     = half.get_view<const Real**,Host>();
    auto col_qv_h   = col_qv.get_view<const Real*,Host>();
    auto backtend_h = backtend.get_view<const Real**,Host>();

    for (int i=0; i<nlocal; ++i) {
      Real qv_col_sum = 0;
      for (int k=0; k<nlevs; ++k) {
        const Real qv_v = n + i*(k+1);
        const Real T_v  = n - i*(k+1);
        qv_col_sum += qv_v;

        REQUIRE (prod_h(i,k)==qv_v*T_v);
        // (qv+qv)*T_mid, not qv+(qv*T_mid): '*' binds tighter than '+'
        REQUIRE (shifted_h(i,k)==(qv_v+qv_v)*T_v);
        REQUIRE (half_h(i,k)==qv_v*T_v/2);
        // T_mid gains exactly -i*(k+1)+... i.e. 1 per step, over dt=1
        REQUIRE (backtend_h(i,k)==1);
      }
      REQUIRE (col_qv_h(i)==qv_col_sum/nlevs);
    }
  }
  reader.clean_up();
  }

  // An expression in field_names would become an nc variable of that name.
  // It needs a name -- decided by how it RESOLVED, not by its characters.
  {
    auto bad = output_params("io_dexpr_unnamed",gname,{"qv*T_mid"});
    OutputManager om_bad;
    om_bad.initialize(comm,bad,t0,false);
    REQUIRE_THROWS (om_bad.setup(fm,gm->get_grid_names()));
  }

  // ...so a legacy name with a '.' is still fine unaliased, which a
  // character-based rule would have broken.
  {
    auto legacy = output_params("io_dexpr_legacy",gname,{"T_mid_where_qv_gt_0.01"});
    OutputManager om_legacy;
    om_legacy.initialize(comm,legacy,t0,false);
    REQUIRE_NOTHROW (om_legacy.setup(fm,gm->get_grid_names()));
    om_legacy.finalize();
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
