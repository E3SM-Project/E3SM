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
#include "share/physics/physics_constants.hpp"

#include <ekat_comm.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_units.hpp>

#include <cmath>
#include <cstdio>
#include <fstream>
#include <map>
#include <set>
#include <string>
#include <vector>

/*
 * End-to-end checks for names resolved by the expression front end.
 *
 * create_diag.cpp covers the translation. Only here can we check that the field
 * survives an output stream, which looks fields up by the REQUESTED name while
 * a diag derives its own from its params -- not the same thing for an
 * expression.
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
                                   const std::vector<std::string>& aliases = {},
                                   const int freq = 1)
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
  ctrl_pl.set<int>("frequency",freq);
  ctrl_pl.set<bool>("save_grid_data",false);
  // No snapshot at t0: a tendency has no meaning before the first step
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
        // Grouping is explicit
        "shifted := (qv+qv)*T_mid",
        // Built on the intermediate below, by name
        "half := qvT/2",
        // A method call
        "col_qv := qv.mean('lev')",
        // A method call on a method call, bringing the column mean back to a 3d layout
        "col_qv_3d := qv.mean('lev').broadcast_like(qv)",
        // A shorthand that expands to (T_mid-T_mid.shift(time=1)).over_dt()
        "T_mid_atm_backtend := T_mid.tend()",
        // A negative literal is a number, on either side of the operator
        "neg_l := -2*T_mid",
        "neg_r := T_mid*(-2)",
        // ...and a legal threshold
        "T_pos := T_mid.where(T_mid>-1)",
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
  REQUIRE (scorpio::has_var(filename,"col_qv_3d"));
  REQUIRE (scorpio::has_var(filename,"T_mid_atm_backtend"));
  REQUIRE (scorpio::has_var(filename,"neg_l"));
  REQUIRE (scorpio::has_var(filename,"neg_r"));
  REQUIRE (scorpio::has_var(filename,"T_pos"));
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
  Field col_qv_3d (FieldIdentifier("col_qv_3d", layout3d, kg/kg, gname));
  Field backtend(FieldIdentifier("T_mid_atm_backtend",layout3d, K/s, gname));
  Field neg_l   (FieldIdentifier("neg_l",  layout3d, K, gname));
  Field neg_r   (FieldIdentifier("neg_r",  layout3d, K, gname));
  Field T_pos   (FieldIdentifier("T_pos",  layout3d, K, gname));
  for (auto* f : {&prod,&shifted,&half,&col_qv,&col_qv_3d,&backtend,&neg_l,&neg_r,&T_pos}) {
    f->allocate_view();
    f->get_header().get_tracking().update_time_stamp(t0);
  }

  FieldReader reader;
  reader.set_file_specs(filename);
  reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
  reader.set_fields({prod,shifted,half,col_qv,col_qv_3d,backtend,neg_l,neg_r,T_pos});

  const int nlocal = grid->get_num_local_dofs();
  for (int n=1; n<=nsteps; ++n) {
    reader.read(n-1);
    prod.sync_to_host();
    shifted.sync_to_host();
    half.sync_to_host();
    col_qv.sync_to_host();
    col_qv_3d.sync_to_host();
    backtend.sync_to_host();
    neg_l.sync_to_host();
    neg_r.sync_to_host();
    T_pos.sync_to_host();
    auto neg_l_h    = neg_l.get_view<const Real**,Host>();
    auto neg_r_h    = neg_r.get_view<const Real**,Host>();
    auto T_pos_h    = T_pos.get_view<const Real**,Host>();
    auto prod_h     = prod.get_view<const Real**,Host>();
    auto shifted_h  = shifted.get_view<const Real**,Host>();
    auto half_h     = half.get_view<const Real**,Host>();
    auto col_qv_h   = col_qv.get_view<const Real*,Host>();
    auto col_qv_3d_h = col_qv_3d.get_view<const Real**,Host>();
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
        REQUIRE (neg_l_h(i,k)==-2*T_v);
        REQUIRE (neg_r_h(i,k)==-2*T_v);
        // Kept where T_mid>-1, missing elsewhere; T_v takes both signs
        if (T_v>-1) {
          REQUIRE (T_pos_h(i,k)==T_v);
        } else {
          REQUIRE (T_pos_h(i,k)==constants::fill_value<Real>);
        }
      }
      REQUIRE (col_qv_h(i)==qv_col_sum/nlevs);
      for (int k=0; k<nlevs; ++k) {
        REQUIRE (col_qv_3d_h(i,k)==col_qv_h(i));
      }
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

  // A stream whose setup throws must say why, and must not leave its diags
  // in the repo shared by all streams
  auto require_setup_throws = [&](const ekat::ParameterList& params, const std::string& msg) {
    const auto repo = DiagRepo::names();
    OutputManager om_fail;
    om_fail.initialize(comm,params,t0,false);
    REQUIRE_THROWS_WITH (om_fail.setup(fm,gm->get_grid_names()), Catch::Contains(msg));
    REQUIRE (DiagRepo::names()==repo);
  };

  // shift(time=1) of a diag is not allowed in an instant stream, unless it writes every step
  {
    const std::string lag_msg = "shift(time=1) of a diagnostic is not available";
    const std::vector<std::string> lagged = {"d := (T_mid*1.0)-(T_mid*1.0).shift(time=1)"};
    require_setup_throws(output_params("io_dexpr_lag_sparse",gname,lagged,{},3),lag_msg);
    // ...also for the shift of an alias of a diag, for tend(), and for a named diag
    for (const auto& names : std::vector<std::vector<std::string>>{
            {"a := T_mid*1.0", "d := a.shift(time=1)"},
            {"bt := (T_mid*1.0).tend()"},
            {"vp := VirtualTemperature.shift(time=1)"},
            {"vbt := VirtualTemperature.tend()"}}) {
      require_setup_throws(output_params("io_dexpr_lag_sparse",gname,names,{},3),lag_msg);
    }

    // Writing every step is fine
    auto dense = output_params("io_dexpr_lag_dense",gname,lagged,{},1);
    OutputManager om_dense;
    om_dense.initialize(comm,dense,t0,false);
    REQUIRE_NOTHROW (om_dense.setup(fm,gm->get_grid_names()));
    om_dense.finalize();

    // ...however the frequency units are spelled
    auto dense2 = output_params("io_dexpr_lag_dense2",gname,lagged,{},1);
    dense2.sublist("output_control").set<std::string>("frequency_units","nstep");
    OutputManager om_dense2;
    om_dense2.initialize(comm,dense2,t0,false);
    REQUIRE_NOTHROW (om_dense2.setup(fm,gm->get_grid_names()));
    om_dense2.finalize();

    // ...and so is an averaged stream, whatever its frequency
    auto avg = output_params("io_dexpr_lag_avg",gname,lagged,{},3);
    avg.set<std::string>("averaging_type","AVERAGE");
    OutputManager om_avg;
    om_avg.initialize(comm,avg,t0,false);
    REQUIRE_NOTHROW (om_avg.setup(fm,gm->get_grid_names()));
    om_avg.finalize();

    // The shift of a model field, or of an alias of one, is always fine
    for (const auto& names : std::vector<std::vector<std::string>>{
            {"d := T_mid-T_mid.shift(time=1)"},
            {"a := T_mid", "d := a.shift(time=1)"}}) {
      auto plain = output_params("io_dexpr_lag_field",gname,names,{},3);
      OutputManager om_plain;
      om_plain.initialize(comm,plain,t0,false);
      REQUIRE_NOTHROW (om_plain.setup(fm,gm->get_grid_names()));
      om_plain.finalize();
    }
  }

  // An int field (e.g., a mask) cannot be written, in either stream type...
  {
    const std::string int_msg = "Output fields must be real-valued";
    require_setup_throws(output_params("io_dexpr_intmask",gname,{"m := mask.where(qv>1)"}),int_msg);

    auto imask_avg = output_params("io_dexpr_intmask_avg",gname,{"m := mask.where(qv>1)"});
    imask_avg.set<std::string>("averaging_type","AVERAGE");
    require_setup_throws(imask_avg,int_msg);

    // ...unless converted to real
    auto rmask = output_params("io_dexpr_realmask",gname,{"m := mask.where(qv>1)*1.0"});
    OutputManager om_r;
    om_r.initialize(comm,rmask,t0,false);
    REQUIRE_NOTHROW (om_r.setup(fm,gm->get_grid_names()));
    om_r.finalize();
  }

  // At t0, a diag whose input is not valid yet is set to fill_value.
  // That must work also for an int diag (the mask)
  {
    Field late (FieldIdentifier("late",grid->get_3d_scalar_layout(LEV),K,gname));
    late.allocate_view();
    fm->add_field(late);

    auto mask_t0 = output_params("io_dexpr_mask_t0",gname,{"m := mask.where(late>0)*1.0"});
    mask_t0.sublist("output_control").set<bool>("skip_t0_output",false);
    OutputManager om_t0;
    om_t0.initialize(comm,mask_t0,t0,false);
    // NOTE: the t0 output is written by setup
    REQUIRE_NOTHROW (om_t0.setup(fm,gm->get_grid_names()));
    om_t0.finalize();
  }

  // An operand that is both a constant and a field is an error, whether the
  // field is there when the operation is created (alias of a field) or not (of a diag)
  const std::string both_msg = "An operand names both a physical constant and a field";
  for (const std::string target : {"T_mid", "T_mid*1.0"}) {
    require_setup_throws(output_params("io_dexpr_collision",gname,{"gravit := " + target,"x := T_mid*gravit"}),both_msg);
  }
  // ...or the name is introduced in the 'aliases' section
  require_setup_throws(output_params("io_dexpr_collision",gname,{"out := (Rair*qv)"},{"Rair := T_mid"}),both_msg);

  {
    // The constant alone is fine
    auto cst = output_params("io_dexpr_constant",gname,{"x := T_mid*gravit"});
    OutputManager om_cst;
    om_cst.initialize(comm,cst,t0,false);
    REQUIRE_NOTHROW (om_cst.setup(fm,gm->get_grid_names()));
    om_cst.finalize();
  }

  scorpio::finalize_subsystem();
}

TEST_CASE ("io_with_math_functions")
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

  // A pressure, constant in time, to compare pow(..) with the Exner diag
  const auto layout3d = grid->get_3d_scalar_layout(LEV);
  const int nlocal = grid->get_num_local_dofs();
  auto p_of = [](const int i, const int k) { return Real(50000 + 1000*(i+k)); };
  {
    Field p_mid (FieldIdentifier("p_mid",layout3d,Pa,gname));
    p_mid.allocate_view();
    auto p_h = p_mid.get_view<Real**,Host>();
    for (int i=0; i<nlocal; ++i) {
      for (int k=0; k<nlevs; ++k) {
        p_h(i,k) = p_of(i,k);
      }
    }
    p_mid.sync_to_dev();
    p_mid.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(p_mid);
  }

  // A wind, constant in time, to compare wind_speed with sqrt(U**2+V**2). Expressions cannot
  // index a component, so U and V are fields of their own, subfields of horiz_winds, as the
  // atmosphere driver makes them.
  auto u_of = [](const int i, const int k) { return Real(3 + i - 0.5*k); };
  auto v_of = [](const int i, const int k) { return Real(-4 + 0.25*i + k); };
  {
    const auto layout_uv = grid->get_3d_vector_layout(LEV,2,"cmp");
    Field winds (FieldIdentifier("horiz_winds",layout_uv,m/s,gname));
    winds.allocate_view();
    auto w_h = winds.get_view<Real***,Host>();
    for (int i=0; i<nlocal; ++i) {
      for (int k=0; k<nlevs; ++k) {
        w_h(i,0,k) = u_of(i,k);
        w_h(i,1,k) = v_of(i,k);
      }
    }
    winds.sync_to_dev();
    winds.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(winds);
    fm->add_field(winds.get_component(0).alias("U"));
    fm->add_field(winds.get_component(1).alias("V"));
  }

  // Rd/cp, written with enough digits to round trip
  using C = physics::Constants<Real>;
  const double kappa = static_cast<double>(C::RD.value*C::INV_CP.value);
  char kappa_str[64];
  std::snprintf(kappa_str,sizeof(kappa_str),"%.17g",kappa);

  const std::string prefix = "io_dexpr_math";
  const int dt = 1;
  const int nsteps = 2;

  // The output names, and what each one is
  const std::vector<std::string> names = {
    "lq", "eq", "sq", "aT", "lT", "rt_exp", "rt_sqrt", "lq_avg", "avg_lq",
    "T2", "T2p", "cbrt", "ex", "ex_hand", "ex_ref", "ws_hand", "ws_ref",
  };
  auto params = output_params(prefix,gname,
      {
        // qv>=1, so these are all in their domain. qv is in kg/kg, so it is
        // dimensionless, as log and exp require
        "lq := log(qv)",
        "eq := exp(-qv)",
        "sq := sqrt(qv)",
        // T_mid takes both signs
        "aT := abs(T_mid)",
        // T_mid has units and invalid entries; both are dealt with by keywords
        "lT := log(T_mid,invalid='mask',units='ignore')",
        // Functions compose with each other and with operators: round trips
        "rt_exp := exp(log(qv))",
        "rt_sqrt := sqrt(qv)*sqrt(qv)",
        // Powers: the operator and the function are the same diag. A cube root
        // is the exponent as a fraction.
        "T2 := T_mid**2",
        "T2p := pow(T_mid,2)",
        "cbrt := qv**(1/3)",
        // Exner, three ways: as a power, as exp(kappa*log(..)), and as the diag itself
        std::string("ex := pow(p_mid/P0,") + kappa_str + ")",
        std::string("ex_hand := exp(") + kappa_str + "*log(p_mid/P0))",
        "ex_ref := Exner",
        // The wind speed diag, and the same by hand: units (m/s)^2 and back
        "ws_hand := sqrt(U**2+V**2)",
        "ws_ref := wind_speed",
        // Functions and methods compose, either way around
        "lq_avg := log(qv).mean('lev')",
        "avg_lq := log(qv.mean('lev'))",
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
    fm->get_field("p_mid").get_header().get_tracking().update_time_stamp(t);
    for (const auto& name : {"horiz_winds","U","V"}) {
      fm->get_field(name).get_header().get_tracking().update_time_stamp(t);
    }

    om.run(t);
  }
  om.finalize();
  }

  const auto filename = prefix + ".INSTANT.nsteps_x1.np" +
                        std::to_string(comm.size()) + "." + (t0+dt).to_string() + ".nc";
  scorpio::register_file(filename,scorpio::Read);
  for (const auto& name : names) {
    REQUIRE (scorpio::has_var(filename,name));
  }
  scorpio::release_file(filename);

  // NOTE: scoped, so the fields are gone before Kokkos is finalized
  {
  const auto layout2d = grid->get_2d_scalar_layout();
  std::map<std::string,Field> f;
  std::vector<Field> to_read;
  for (const auto& name : names) {
    const bool is2d = name=="lq_avg" or name=="avg_lq";
    f[name] = Field(FieldIdentifier(name,is2d ? layout2d : layout3d,none,gname));
    f[name].allocate_view();
    f[name].get_header().get_tracking().update_time_stamp(t0);
    to_read.push_back(f[name]);
  }

  FieldReader reader;
  reader.set_file_specs(filename);
  reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
  reader.set_fields(to_read);

  constexpr auto fv = constants::fill_value<Real>;
  for (int n=1; n<=nsteps; ++n) {
    reader.read(n-1);
    for (auto& it : f) {
      it.second.sync_to_host();
    }
    auto v2 = [&](const std::string& name) { return f.at(name).get_view<const Real**,Host>(); };
    auto v1 = [&](const std::string& name) { return f.at(name).get_view<const Real*,Host>(); };

    for (int i=0; i<nlocal; ++i) {
      Real qv_col_sum = 0;
      Real lq_col_sum = 0;
      for (int k=0; k<nlevs; ++k) {
        const Real qv_v = n + i*(k+1);
        const Real T_v  = n - i*(k+1);
        qv_col_sum += qv_v;
        lq_col_sum += std::log(qv_v);

        REQUIRE (v2("lq")(i,k)==Approx(std::log(qv_v)));
        REQUIRE (v2("eq")(i,k)==Approx(std::exp(-qv_v)));
        REQUIRE (v2("sq")(i,k)==Approx(std::sqrt(qv_v)));
        REQUIRE (v2("aT")(i,k)==std::abs(T_v));
        if (T_v>0) {
          REQUIRE (v2("lT")(i,k)==Approx(std::log(T_v)));
        } else {
          REQUIRE (v2("lT")(i,k)==fv);
        }
        REQUIRE (v2("rt_exp")(i,k)==Approx(qv_v));
        REQUIRE (v2("rt_sqrt")(i,k)==Approx(qv_v));
        REQUIRE (v2("T2")(i,k)==Approx(T_v*T_v));
        REQUIRE (v2("T2p")(i,k)==v2("T2")(i,k));
        REQUIRE (v2("cbrt")(i,k)==Approx(std::cbrt(qv_v)));
        // The three Exners agree, and match the formula
        REQUIRE (v2("ex")(i,k)==Approx(v2("ex_ref")(i,k)));
        REQUIRE (v2("ex_hand")(i,k)==Approx(v2("ex_ref")(i,k)));
        REQUIRE (v2("ex")(i,k)==Approx(std::pow(p_of(i,k)/C::P0.value,kappa)));
        // The hand-rolled wind speed matches the diag, and the formula
        REQUIRE (v2("ws_hand")(i,k)==Approx(v2("ws_ref")(i,k)));
        REQUIRE (v2("ws_hand")(i,k)==Approx(std::sqrt(u_of(i,k)*u_of(i,k) + v_of(i,k)*v_of(i,k))));
      }
      REQUIRE (v1("lq_avg")(i)==Approx(lq_col_sum/nlevs));
      REQUIRE (v1("avg_lq")(i)==Approx(std::log(qv_col_sum/nlevs)));
    }
  }
  reader.clean_up();
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
