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

#include <ekat_comm.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_units.hpp>

#include <set>
#include <string>
#include <vector>

/*
 * How far back a shift looks, through an output stream.
 *
 * FieldPrev copies its operand in init_timestep, and the stream calls
 * init_timestep in evaluation order. So the shift of a shift copies a value
 * that was already updated in this step, and lags one step rather than two.
 * That is rejected at setup; with an operation between the two shifts, as in
 * (X.shift(time=1)*1.0).shift(time=1), the lag is two steps, which we check
 * against hand-computed values.
 */

namespace scream {

namespace {

constexpr int ncols = 6;
constexpr int nlevs = 3;
constexpr int nsteps = 5;

constexpr Real FILL = constants::fill_value<Real>;

// T_mid on step n. Quadratic in n, so a lag of one step and of two steps give
// different values, and so do the first and second differences.
Real T_val (const int n, const int i, const int k) {
  return n*n + i*(k+1);
}

void set_T (Field& T_mid, const int n)
{
  auto T_h = T_mid.get_view<Real**,Host>();
  const int nl = T_mid.get_header().get_identifier().get_layout().dim(0);
  for (int i=0; i<nl; ++i) {
    for (int k=0; k<nlevs; ++k) {
      T_h(i,k) = T_val(n,i,k);
    }
  }
  T_mid.sync_to_dev();
}

// By default, an instant stream that writes every step, the only kind of
// instant stream where the shift of a diagnostic is allowed
ekat::ParameterList output_params (const std::string& prefix,
                                   const std::string& grid_name,
                                   const std::vector<std::string>& field_names,
                                   const std::vector<std::string>& aliases = {},
                                   const std::string& avg = "INSTANT")
{
  ekat::ParameterList params;
  params.set<std::string>("filename_prefix",prefix);
  params.set<std::string>("averaging_type",avg);
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
  ctrl_pl.set<bool>("skip_t0_output",true);

  return params;
}

// The names of the diags in the repo that all output streams share
struct DiagRepo : public AtmosphereOutput {
  static std::set<std::string> names () {
    std::set<std::string> n;
    for (const auto& it : m_diag_repo) {
      n.insert(it.first);
    }
    return n;
  }
};

// What a written field must hold on step n. The shift of a diagnostic has no
// value on the first step: that diag has not been computed when the step starts.
struct Expected {
  std::string name;
  int first_valid;
  Real (*value) (int n, int i, int k);
};

} // anonymous namespace

TEST_CASE ("io_dexpr_lag")
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
  const auto layout3d = grid->get_3d_scalar_layout(LEV);

  util::TimeStamp t0({2023,1,1},{0,0,0});
  const int dt = 1;

  auto fm = std::make_shared<FieldManager>(gm);
  {
    Field T_mid (FieldIdentifier("T_mid",layout3d,K,gname));
    T_mid.allocate_view();
    T_mid.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(T_mid);
    set_T(T_mid,0);
  }

  // A stream whose setup throws must say why, and must not leave its diags
  // in the repo shared by all streams
  auto require_setup_throws = [&](const ekat::ParameterList& params, const std::string& msg) {
    const auto repo = DiagRepo::names();
    OutputManager om_fail;
    om_fail.initialize(comm,params,t0,false);
    REQUIRE_THROWS_WITH (om_fail.setup(fm,gm->get_grid_names()), Catch::Contains(msg));
    REQUIRE (DiagRepo::names()==repo);
  };

  // The shift of a shift is rejected: directly, through an alias (also an
  // alias of an alias, and an 'aliases' intermediate), inside tend(), and in
  // the underscore spelling. Whatever the averaging type.
  const std::string nested_msg = "shift(time=1) of a shift(time=1) is not available";
  for (const std::string avg : {"INSTANT","AVERAGE"}) {
    for (const auto& names : std::vector<std::vector<std::string>>{
            {"s2 := T_mid.shift(time=1).shift(time=1)"},
            {"l1 := T_mid.shift(time=1)", "l2 := l1.shift(time=1)"},
            {"l1 := T_mid.shift(time=1)", "l1b := l1", "l2 := l1b.shift(time=1)"},
            {"st := T_mid.shift(time=1).tend()"},
            {"T_mid_prev_prev"}}) {
      require_setup_throws(output_params("io_dexpr_lag_nested",gname,names,{},avg),nested_msg);
    }
    require_setup_throws(output_params("io_dexpr_lag_nested",gname,
                                       {"l2 := l1.shift(time=1)"},{"l1 := T_mid.shift(time=1)"},avg),
                         nested_msg);
  }

  // The shift of anything else is fine, e.g., of an alias of a diag, of a
  // tendency, or of a model field
  for (const auto& names : std::vector<std::vector<std::string>>{
          {"a := T_mid*1.0", "ap := a.shift(time=1)"},
          {"tp := T_mid.tend().shift(time=1)"},
          {"T_mid_prev", "T_mid_atm_backtend"}}) {
    auto ok = output_params("io_dexpr_lag_ok",gname,names,{},"AVERAGE");
    OutputManager om_ok;
    om_ok.initialize(comm,ok,t0,false);
    REQUIRE_NOTHROW (om_ok.setup(fm,gm->get_grid_names()));
    om_ok.finalize();
  }

  // Values, including the shift of an operation on a shift, which is two steps back
  const std::vector<Expected> expected = {
    {"s1",   1, [](int n, int i, int k) { return T_val(n-1,i,k); }},
    {"d1",   2, [](int n, int i, int k) { return T_val(n-1,i,k); }},
    {"a2",   2, [](int n, int i, int k) { return T_val(n-2,i,k); }},
    {"bt1",  1, [](int n, int i, int k) { return T_val(n,i,k)-T_val(n-1,i,k); }},
    {"bt1p", 2, [](int n, int i, int k) { return T_val(n-1,i,k)-T_val(n-2,i,k); }},
    {"bt2",  2, [](int n, int i, int k) { return T_val(n,i,k)-2*T_val(n-1,i,k)+T_val(n-2,i,k); }},
  };
  const std::vector<std::string> requests = {
    "s1 := T_mid.shift(time=1)",
    // the shift of a derived quantity
    "d1 := (T_mid*1.0).shift(time=1)",
    // two steps back: an operation between the shifts
    "a2 := (T_mid.shift(time=1)*1.0).shift(time=1)",
    // the backtend chain of the surface-oscillation diagnostics
    "bt1 := T_mid-T_mid.shift(time=1)",
    "bt1p := bt1.shift(time=1)",
    "bt2 := bt1-bt1.shift(time=1)",
  };

  const std::string prefix = "io_dexpr_lag";
  auto t = t0;
  {
    OutputManager om;
    om.initialize(comm,output_params(prefix,gname,requests),t0,false);
    om.setup(fm,gm->get_grid_names());
    for (int n=1; n<=nsteps; ++n) {
      om.init_timestep(t,dt);
      t += dt;
      auto T_mid = fm->get_field("T_mid");
      set_T(T_mid,n);
      T_mid.get_header().get_tracking().update_time_stamp(t);
      om.run(t);
    }
    om.finalize();
  }

  // NOTE: scoped, so the fields are gone before Kokkos is finalized
  {
  const auto filename = prefix + ".INSTANT.nsteps_x1.np" +
                        std::to_string(comm.size()) + "." + (t0+dt).to_string() + ".nc";
  std::vector<Field> fields;
  for (const auto& e : expected) {
    Field f (FieldIdentifier(e.name,layout3d,K,gname));
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    fields.push_back(f);
  }
  FieldReader reader;
  reader.set_file_specs(filename);
  reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
  reader.set_fields(fields);

  const int nlocal = grid->get_num_local_dofs();
  for (int n=1; n<=nsteps; ++n) {
    reader.read(n-1);
    for (size_t j=0; j<expected.size(); ++j) {
      const auto& e = expected[j];
      fields[j].sync_to_host();
      auto v = fields[j].get_view<const Real**,Host>();
      for (int i=0; i<nlocal; ++i) {
        for (int k=0; k<nlevs; ++k) {
          INFO ("field " << e.name << ", step " << n << ", col " << i << ", lev " << k);
          if (n<e.first_valid) {
            REQUIRE (v(i,k)==FILL);
          } else {
            REQUIRE (v(i,k)==e.value(n,i,k));
          }
        }
      }
    }
  }
  reader.clean_up();
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
