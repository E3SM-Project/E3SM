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

#include <fstream>
#include <string>
#include <vector>

/*
 * Fields interpolated to the same pressure level share a valid-sample count,
 * also when aliased, except for a field whose input is masked or may be filled.
 * Midpoint and interface fields do not share, since their ranges differ.
 * An alias of a filled model field uses the count of that field.
 * Model fields have a count of their own, also if they share a mask, since
 * their fill pattern need not follow it.
 * A vector field at a level shares the level's count name; the count adds the
 * vector dim to the name, as for any count.
 * In MAX/MIN streams, an entry that is never valid is fill_value, also when
 * its count is shared with a field that is valid there.
 */

namespace scream {

namespace {

constexpr int ncols  = 4;
constexpr int nlevs  = 4;
constexpr int nsteps = 4;

constexpr Real FILL = constants::fill_value<Real>;

// 200/400/600/800 hPa, so the requested 500 hPa is strictly interior and is
// bracketed by levels 1 and 2 in every column. The interfaces reach 900 hPa,
// so 850 hPa is in range for an interface field, but not for a midpoint one.
Real p_mid_at (const int k) { return 20000.0 + 20000.0*k; }
Real p_int_at (const int k) { return 10000.0 + 20000.0*k; }

// X is the step number, constant in the vertical, so interpolation is exact.
// C is never positive, so X.where(C>0) is fill_value everywhere.
void calc_fields (Field& X, Field& C, const int n)
{
  auto X_h = X.get_view<Real**,Host>();
  auto C_h = C.get_view<Real**,Host>();
  const int nl = X.get_header().get_identifier().get_layout().dim(0);
  for (int i=0; i<nl; ++i) {
    for (int k=0; k<nlevs; ++k) {
      X_h(i,k) = n;
      C_h(i,k) = -1.0;
    }
  }
  X.sync_to_dev();
  C.sync_to_dev();
}

// The step number, except in the first local column, which is filled
void calc_F (Field& F, const int n)
{
  auto F_h = F.get_view<Real**,Host>();
  const int nl = F.get_header().get_identifier().get_layout().dim(0);
  for (int i=0; i<nl; ++i) {
    for (int k=0; k<nlevs; ++k) {
      F_h(i,k) = i==0 ? FILL : Real(n);
    }
  }
  F.sync_to_dev();
}

std::shared_ptr<FieldManager>
create_test_fm (const std::shared_ptr<const GridsManager>& gm,
                const std::shared_ptr<const AbstractGrid>& grid,
                const util::TimeStamp& t0)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto layout_mid = grid->get_3d_scalar_layout(LEV);
  const auto layout_int = grid->get_3d_scalar_layout(ILEV);
  auto fm = std::make_shared<FieldManager>(gm);

  Field X(FieldIdentifier("X",layout_mid,K,grid->name()));
  Field C(FieldIdentifier("C",layout_mid,K,grid->name()));
  // FieldAtPressureLevel lists both p_mid and p_int as inputs, since it does
  // not know which it needs until it sees the field's layout. Provide both.
  Field pm(FieldIdentifier("p_mid",layout_mid,Pa,grid->name()));
  Field pi(FieldIdentifier("p_int",layout_int,Pa,grid->name()));
  for (auto* f : {&X,&C,&pm,&pi}) {
    f->allocate_view();
    f->get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(*f);
  }

  auto pm_h = pm.get_view<Real**,Host>();
  auto pi_h = pi.get_view<Real**,Host>();
  const int nl = grid->get_num_local_dofs();
  for (int i=0; i<nl; ++i) {
    for (int k=0; k<nlevs; ++k)   pm_h(i,k) = p_mid_at(k);
    for (int k=0; k<nlevs+1; ++k) pi_h(i,k) = p_int_at(k);
  }
  pm.sync_to_dev();
  pi.sync_to_dev();

  // A field on interfaces (set to the step number in the time loop)
  Field Xi(FieldIdentifier("Xi",layout_int,K,grid->name()));
  Xi.allocate_view();
  Xi.get_header().get_tracking().update_time_stamp(t0);
  Xi.deep_copy(0);
  fm->add_field(Xi);

  // A vector field (set to the step number in the time loop)
  Field W(FieldIdentifier("W",grid->get_3d_vector_layout(LEV,2),K,grid->name()));
  W.allocate_view();
  W.get_header().get_tracking().update_time_stamp(t0);
  W.deep_copy(0);
  fm->add_field(W);

  // A model field that is filled in the first local column
  Field M(FieldIdentifier("m2d",grid->get_2d_scalar_layout(),K,grid->name()));
  M.allocate_view();
  M.get_header().set_may_be_filled(true);
  M.get_header().get_tracking().update_time_stamp(t0);
  M.deep_copy(7.0);
  M.sync_to_host();
  M.get_view<Real*,Host>()(0) = FILL;
  M.sync_to_dev();
  fm->add_field(M);

  // A model field that may be filled, with no mask: filled in the first local column
  Field F(FieldIdentifier("F",layout_mid,K,grid->name()));
  F.allocate_view();
  F.get_header().set_may_be_filled(true);
  F.get_header().get_tracking().update_time_stamp(t0);
  calc_F(F,0);
  fm->add_field(F);

  // Two model fields with the same mask, valid everywhere, but S2 is filled
  // everywhere (like a COSP output, which is refreshed less often than its mask)
  Field shared_mask(FieldIdentifier("shared_mask",grid->get_2d_scalar_layout(),none,grid->name(),DataType::IntType));
  shared_mask.allocate_view();
  shared_mask.deep_copy(1);
  for (std::string n : {"S1","S2"}) {
    Field S(FieldIdentifier(n,grid->get_2d_scalar_layout(),K,grid->name()));
    S.allocate_view();
    S.set_valid_mask(shared_mask);
    S.get_header().set_may_be_filled(true);
    S.get_header().get_tracking().update_time_stamp(t0);
    S.deep_copy(n=="S1" ? Real(3.0) : FILL);
    fm->add_field(S);
  }

  calc_fields(X,C,0);
  return fm;
}

// What a stream of the given type writes for the step numbers 1..nsteps
Real expected_value (const std::string& avg_type)
{
  if (avg_type=="MAX") return nsteps;
  if (avg_type=="MIN") return 1;
  Real sum = 0;
  for (int n=1; n<=nsteps; ++n) sum += n;
  return sum/nsteps;
}

void run_and_check (const std::string& avg_type,
                    const std::shared_ptr<const GridsManager>& gm,
                    const ekat::Comm& comm)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  auto grid = gm->get_grid("point_grid");
  const auto gname = grid->name();

  util::TimeStamp t0({2023,1,1},{0,0,0});
  auto fm = create_test_fm(gm,grid,t0);

  const std::string prefix = "io_dexpr_plev_avg_cnt";
  const int dt = 1;

  ekat::ParameterList params;
  params.set<std::string>("filename_prefix",prefix);
  params.set<std::string>("averaging_type",avg_type);
  params.set<std::string>("floating_point_precision","real");
  auto& f_pl = params.sublist("fields").sublist(gname);
  // The masked source is an intermediate (not written)
  f_pl.set<std::vector<std::string>>("aliases",{
      "Xmasked := X.where(C>0)",
  });
  f_pl.set<std::vector<std::string>>("field_names",{
      "X_at_500hPa",
      "C500 := C_at_500hPa",
      "Xmasked_at_500hPa",
      // NOTE: the alias sorts before its target, and is processed first
      "AM := m2d",
      "X_at_850hPa",
      "Xi_at_850hPa",
      "F_at_500hPa",
      "S1",
      "S2",
      "W_at_500hPa",
  });
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
      auto X = fm->get_field("X");
      auto C = fm->get_field("C");
      calc_fields(X,C,n);
      auto Xi = fm->get_field("Xi");
      Xi.deep_copy(n);
      auto F = fm->get_field("F");
      calc_F(F,n);
      auto W = fm->get_field("W");
      W.deep_copy(n);
      for (auto f : {X,C,Xi,F,W,fm->get_field("p_mid"),fm->get_field("p_int"),fm->get_field("m2d"),
                     fm->get_field("S1"),fm->get_field("S2")}) {
        f.get_header().get_tracking().update_time_stamp(t);
      }
      om.run(t);
    }
    om.finalize();
  }

  const auto filename = prefix + "." + avg_type + ".nsteps_x" + std::to_string(nsteps) +
                        ".np" + std::to_string(comm.size()) + "." +
                        t0.to_string() + ".nc";
  std::ifstream file_check(filename);
  REQUIRE (file_check.good());
  file_check.close();

  scorpio::register_file(filename,scorpio::Read);
  REQUIRE (scorpio::has_var(filename,"X_at_500hPa"));
  REQUIRE (scorpio::has_var(filename,"C500"));
  REQUIRE (scorpio::has_var(filename,"Xmasked_at_500hPa"));
  // The unmasked fields share the count keyed on the level, alias included.
  // CHECK, so a failure still reports the values below
  CHECK (scorpio::has_var(filename,"avg_count_500hPa_ncol"));
  CHECK_FALSE (scorpio::has_var(filename,"avg_count_X_at_500hPa_ncol"));
  CHECK_FALSE (scorpio::has_var(filename,"avg_count_C500_ncol"));
  CHECK_FALSE (scorpio::has_var(filename,"avg_count_ncol"));
  // ... while the field with a masked input has its own
  CHECK (scorpio::has_var(filename,"avg_count_Xmasked_at_500hPa_ncol"));
  // ... and so does the alias of the filled model field
  CHECK (scorpio::has_var(filename,"avg_count_AM_ncol"));
  // Midpoint and interface fields at the same level do not share
  CHECK (scorpio::has_var(filename,"avg_count_850hPa_ncol"));
  CHECK (scorpio::has_var(filename,"avg_count_850hPa_int_ncol"));
  // A field whose input may be filled has its own count, also without a mask
  CHECK (scorpio::has_var(filename,"avg_count_F_at_500hPa_ncol"));
  // Model fields with the same mask do not share
  CHECK (scorpio::has_var(filename,"avg_count_S1_ncol"));
  CHECK (scorpio::has_var(filename,"avg_count_S2_ncol"));
  CHECK_FALSE (scorpio::has_var(filename,"avg_count_shared_mask_ncol"));
  // The vector field at 500 hPa: the count is named after the level and the layout
  CHECK (scorpio::has_var(filename,"avg_count_500hPa_ncol_dim2"));
  CHECK_FALSE (scorpio::has_var(filename,"avg_count_500hPa_dim2_ncol_dim2"));
  scorpio::release_file(filename);

  {
    const auto layout2d = grid->get_2d_scalar_layout();
    Field xu(FieldIdentifier("X_at_500hPa",layout2d,K,gname));
    Field xm(FieldIdentifier("Xmasked_at_500hPa",layout2d,K,gname));
    Field c5(FieldIdentifier("C500",layout2d,K,gname));
    Field am(FieldIdentifier("AM",layout2d,K,gname));
    Field x8(FieldIdentifier("X_at_850hPa",layout2d,K,gname));
    Field xi8(FieldIdentifier("Xi_at_850hPa",layout2d,K,gname));
    Field f5(FieldIdentifier("F_at_500hPa",layout2d,K,gname));
    Field s1(FieldIdentifier("S1",layout2d,K,gname));
    Field s2(FieldIdentifier("S2",layout2d,K,gname));
    for (auto* f : {&xu,&xm,&c5,&am,&x8,&xi8,&f5,&s1,&s2}) {
      f->allocate_view();
      f->get_header().get_tracking().update_time_stamp(t0);
    }

    FieldReader reader;
    reader.set_file_specs(filename);
    reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
    reader.set_fields({xu,xm,c5,am,x8,xi8,f5,s1,s2});
    reader.read(0);
    for (auto* f : {&xu,&xm,&c5,&am,&x8,&xi8,&f5,&s1,&s2}) {
      f->sync_to_host();
    }
    auto xu_h = xu.get_view<const Real*,Host>();
    auto xm_h = xm.get_view<const Real*,Host>();
    auto c5_h = c5.get_view<const Real*,Host>();
    auto am_h = am.get_view<const Real*,Host>();
    auto x8_h = x8.get_view<const Real*,Host>();
    auto xi8_h = xi8.get_view<const Real*,Host>();
    auto f5_h = f5.get_view<const Real*,Host>();
    auto s1_h = s1.get_view<const Real*,Host>();
    auto s2_h = s2.get_view<const Real*,Host>();

    const Real expected = expected_value(avg_type);

    const int nlocal = grid->get_num_local_dofs();
    for (int i=0; i<nlocal; ++i) {
      // The unmasked field is valid, so 500hPa is in range
      REQUIRE (xu_h(i)!=FILL);
      REQUIRE (xu_h(i)==Approx(expected));
      // (with a count shared with AM, column 0 would be FILL)
      REQUIRE (c5_h(i)==Approx(-1.0));
      // The filled column has no valid samples
      REQUIRE (am_h(i)==(i==0 ? FILL : Real(7.0)));
      if (i==0) {
        REQUIRE (f5_h(i)==FILL);
      } else {
        REQUIRE (f5_h(i)==Approx(expected));
      }

      // 850 hPa is below the lowest midpoint, but above the lowest interface
      REQUIRE (x8_h(i)==FILL);
      REQUIRE (xi8_h(i)==Approx(expected));

      // No valid samples: must be fill_value (it was 0 with a shared count)
      REQUIRE (xm_h(i)==FILL);

      // S2 is never valid (with the count of S1, its average would be 0)
      REQUIRE (s1_h(i)==Approx(3.0));
      REQUIRE (s2_h(i)==FILL);
    }
  }
}

void run_vremap_and_check (const std::shared_ptr<const GridsManager>& gm,
                           const ekat::Comm& comm)
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  auto grid = gm->get_grid("point_grid");
  const auto gname = grid->name();

  // 850 hPa is below the lowest midpoint, but above the lowest interface
  const std::vector<Real> p_tgt = {50000.0, 85000.0};
  const std::string remap_file = "io_dexpr_plev_avg_cnt_vremap_np" + std::to_string(comm.size()) + ".nc";
  scorpio::register_file(remap_file,scorpio::FileMode::Write);
  const int nlevs_tgt = p_tgt.size();
  scorpio::define_dim(remap_file,"lev",nlevs_tgt);
  scorpio::define_var(remap_file,"p_levs",{"lev"},"real");
  scorpio::enddef(remap_file);
  scorpio::write_var(remap_file,"p_levs",p_tgt.data());
  scorpio::release_file(remap_file);

  // After the remap X and Xi have the same layout, hence the same count, which
  // is 0 at 850 hPa, where X is masked and Xi is not
  for (const std::string avg_type : {"MAX","MIN"}) {
    util::TimeStamp t0({2023,1,1},{0,0,0});
    auto fm = create_test_fm(gm,grid,t0);

    const std::string prefix = "io_dexpr_plev_avg_cnt_vremap";
    const int dt = 1;

    ekat::ParameterList params;
    params.set<std::string>("filename_prefix",prefix);
    params.set<std::string>("averaging_type",avg_type);
    params.set<std::string>("floating_point_precision","real");
    params.set<std::string>("vertical_remap_file",remap_file);
    params.sublist("fields").sublist(gname).set<std::vector<std::string>>("field_names",{"X","Xi"});
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
        auto X = fm->get_field("X");
        auto C = fm->get_field("C");
        calc_fields(X,C,n);
        auto Xi = fm->get_field("Xi");
        Xi.deep_copy(n);
        for (auto f : {X,C,Xi,fm->get_field("p_mid"),fm->get_field("p_int")}) {
          f.get_header().get_tracking().update_time_stamp(t);
        }
        om.run(t);
      }
      om.finalize();
    }

    const auto filename = prefix + "." + avg_type + ".nsteps_x" + std::to_string(nsteps) +
                          ".np" + std::to_string(comm.size()) + "." +
                          t0.to_string() + ".nc";
    scorpio::register_file(filename,scorpio::Read);
    REQUIRE (scorpio::has_var(filename,"avg_count_ncol_lev"));
    scorpio::release_file(filename);

    // Read back on a grid with the target levels
    auto gm_tgt = create_mesh_free_grids_manager(comm,0,0,nlevs_tgt,ncols);
    gm_tgt->build_grids();
    auto grid_tgt = gm_tgt->get_grid("point_grid");
    const auto layout = grid_tgt->get_3d_scalar_layout(LEV);
    Field x(FieldIdentifier("X",layout,K,grid_tgt->name()));
    Field xi(FieldIdentifier("Xi",layout,K,grid_tgt->name()));
    for (auto* f : {&x,&xi}) {
      f->allocate_view();
      f->get_header().get_tracking().update_time_stamp(t0);
    }
    read_fields(filename,{x,xi},grid_tgt->get_partitioned_dim_gids(),comm);
    x.sync_to_host();
    xi.sync_to_host();
    auto x_h  = x.get_view<const Real**,Host>();
    auto xi_h = xi.get_view<const Real**,Host>();

    const Real expected = expected_value(avg_type);
    for (int i=0; i<grid_tgt->get_num_local_dofs(); ++i) {
      REQUIRE (x_h(i,0)==Approx(expected));
      REQUIRE (x_h(i,1)==FILL);
      REQUIRE (xi_h(i,0)==Approx(expected));
      REQUIRE (xi_h(i,1)==Approx(expected));
    }
  }
}

} // anonymous namespace

TEST_CASE ("plev_diags_share_an_avg_count_unless_input_is_masked")
{
  register_diagnostics();

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ncols);
  gm->build_grids();

  for (const std::string avg_type : {"AVERAGE","MAX","MIN"}) {
    run_and_check(avg_type,gm,comm);
  }
  run_vremap_and_check(gm,comm);

  scorpio::finalize_subsystem();
}

} // namespace scream
