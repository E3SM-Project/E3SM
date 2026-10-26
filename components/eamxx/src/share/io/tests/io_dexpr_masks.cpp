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

#include <algorithm>
#include <fstream>
#include <string>
#include <vector>

/*
 * Conditional sampling in AVERAGE/MAX/MIN streams, against hand-computed values:
 *  - only the valid steps count, and an entry that was never valid is fill_value
 *  - each masked field gets its own valid-sample count, also if it is an alias
 *  - X.where(A).where(B) is the intersection of the two conditions
 */

namespace scream {

namespace {

constexpr int ncols = 4;
constexpr int nlevs = 3;
constexpr int nsteps = 4;

constexpr Real FILL = constants::fill_value<Real>;

// X is the step number. The two conditions differ, so a1 and a5 have different masks:
//   C1>0 on step n iff (n+k) is even
//   C2>0 on step n iff n>=2, except at level 0, where it never holds
void calc_fields (Field& X, Field& C1, Field& C2, const int n)
{
  auto X_h  = X .get_view<Real**,Host>();
  auto C1_h = C1.get_view<Real**,Host>();
  auto C2_h = C2.get_view<Real**,Host>();
  const int nl = X.get_header().get_identifier().get_layout().dim(0);
  for (int i=0; i<nl; ++i) {
    for (int k=0; k<nlevs; ++k) {
      X_h (i,k) = n;
      C1_h(i,k) = ((n+k)%2==0) ? 1.0 : -1.0;
      C2_h(i,k) = (k>0 and n>=2) ? 1.0 : -1.0;
    }
  }
  X.sync_to_dev();
  C1.sync_to_dev();
  C2.sync_to_dev();
}

// The same predicates, for the expected values
bool c1_holds (int n, int k) { return (n+k)%2==0; }
bool c2_holds (int n, int k) { return k>0 and n>=2; }

std::shared_ptr<FieldManager>
create_test_fm (const std::shared_ptr<const GridsManager>& gm,
                const std::shared_ptr<const AbstractGrid>& grid,
                const util::TimeStamp& t0)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto layout3d = grid->get_3d_scalar_layout(LEV);
  auto fm = std::make_shared<FieldManager>(gm);

  Field X(FieldIdentifier("X" ,layout3d,K,grid->name()));
  Field C1(FieldIdentifier("C1",layout3d,K,grid->name()));
  Field C2(FieldIdentifier("C2",layout3d,K,grid->name()));
  for (auto* f : {&X,&C1,&C2}) {
    f->allocate_view();
    f->get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(*f);
  }
  calc_fields(X,C1,C2,0);
  return fm;
}

// Run one stream with the given averaging type, and check what it wrote
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

  const std::string prefix = "io_dexpr_masks";
  const int dt = 1;

  ekat::ParameterList params;
  params.set<std::string>("filename_prefix",prefix);
  params.set<std::string>("averaging_type",avg_type);
  params.set<std::string>("floating_point_precision","real");
  auto& f_pl = params.sublist("fields").sublist(gname);
  f_pl.set<std::vector<std::string>>("field_names",{
      "a1 := X.where(C1>0)",
      "a5 := X.where(C1>0).where(C2>0)",
      // An alias of an alias
      "b1 := a1",
      // The vertical sum reads the mask of its input, not the fill values
      "s5 := X.where(C1>0).where(C2>0).sum('lev')",
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
      auto X  = fm->get_field("X");
      auto C1 = fm->get_field("C1");
      auto C2 = fm->get_field("C2");
      calc_fields(X,C1,C2,n);
      for (auto f : {X,C1,C2}) {
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

  // Each masked field gets its own count, and no count is shared by layout
  scorpio::register_file(filename,scorpio::Read);
  for (const std::string n : {"a1","a5","b1"}) {
    REQUIRE (scorpio::has_var(filename,n));
    REQUIRE (scorpio::has_var(filename,"avg_count_" + n + "_ncol_lev"));
  }
  REQUIRE (scorpio::has_var(filename,"s5"));
  REQUIRE (scorpio::has_var(filename,"avg_count_s5_ncol"));
  REQUIRE_FALSE (scorpio::has_var(filename,"avg_count_ncol_lev"));
  REQUIRE_FALSE (scorpio::has_var(filename,"avg_count_ncol"));
  // Fields that may hold fill values declare it
  CHECK (scorpio::has_attribute(filename,"a1","_FillValue"));
  CHECK (scorpio::has_attribute(filename,"s5","_FillValue"));
  scorpio::release_file(filename);

  const auto layout3d = grid->get_3d_scalar_layout(LEV);
  const auto layout2d = grid->get_2d_scalar_layout();
  Field a1(FieldIdentifier("a1",layout3d,K,gname));
  Field a5(FieldIdentifier("a5",layout3d,K,gname));
  Field b1(FieldIdentifier("b1",layout3d,K,gname));
  Field s5(FieldIdentifier("s5",layout2d,K,gname));
  for (auto* f : {&a1,&a5,&b1,&s5}) {
    f->allocate_view();
    f->get_header().get_tracking().update_time_stamp(t0);
  }

  FieldReader reader;
  reader.set_file_specs(filename);
  reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
  reader.set_fields({a1,a5,b1,s5});
  reader.read(0);
  for (auto* f : {&a1,&a5,&b1,&s5}) {
    f->sync_to_host();
  }
  auto a1_h = a1.get_view<const Real**,Host>();
  auto a5_h = a5.get_view<const Real**,Host>();
  auto b1_h = b1.get_view<const Real**,Host>();
  auto s5_h = s5.get_view<const Real*,Host>();

  // Accumulate the valid samples like the stream does
  struct Acc {
    std::string type;
    Real val = 0; int cnt = 0;
    void add (const Real v) {
      val = cnt==0 ? v : (type=="MAX" ? std::max(val,v) : type=="MIN" ? std::min(val,v) : val+v);
      ++cnt;
    }
    Real get () const { return cnt==0 ? FILL : (type=="AVERAGE" ? val/cnt : val); }
  };

  const int nlocal = grid->get_num_local_dofs();
  for (int i=0; i<nlocal; ++i) {
    // s5: at each step, the sum over the levels where both C1 and C2 hold (if any)
    Acc acc_s5{avg_type};
    for (int n=1; n<=nsteps; ++n) {
      int nvalid = 0;
      for (int k=0; k<nlevs; ++k) {
        if (c1_holds(n,k) and c2_holds(n,k)) ++nvalid;
      }
      if (nvalid>0) acc_s5.add(n*nvalid);
    }
    REQUIRE (acc_s5.cnt>0);
    REQUIRE (s5_h(i)==Approx(acc_s5.get()));

    for (int k=0; k<nlevs; ++k) {
      // a1: the steps where C1 held; a5: the steps where both C1 and C2 held
      Acc acc1{avg_type}, acc5{avg_type};
      for (int n=1; n<=nsteps; ++n) {
        if (c1_holds(n,k))                   acc1.add(n);
        if (c1_holds(n,k) and c2_holds(n,k)) acc5.add(n);
      }

      // Never valid (e.g., a5 at level 0) means fill_value, whatever the averaging type
      if (acc1.cnt==0) {
        REQUIRE (a1_h(i,k)==FILL);
      } else {
        REQUIRE (a1_h(i,k)==Approx(acc1.get()));
      }
      if (acc5.cnt==0) {
        REQUIRE (a5_h(i,k)==FILL);
      } else {
        REQUIRE (a5_h(i,k)==Approx(acc5.get()));
      }
      REQUIRE (b1_h(i,k)==a1_h(i,k));
    }
  }
}

} // anonymous namespace

TEST_CASE ("masked_conditional_sampling")
{
  register_diagnostics();

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ncols);
  gm->build_grids();

  for (const std::string avg_type : {"AVERAGE","MAX","MIN"}) {
    run_and_check(avg_type,gm,comm);
  }

  // Aliases that form a cycle are an error
  {
    auto grid = gm->get_grid("point_grid");
    util::TimeStamp t0({2023,1,1},{0,0,0});
    auto fm = create_test_fm(gm,grid,t0);

    ekat::ParameterList params;
    params.set<std::string>("filename_prefix","io_dexpr_masks_cycle");
    params.set<std::string>("averaging_type","AVERAGE");
    params.sublist("fields").sublist(grid->name()).set<std::vector<std::string>>("field_names",{
        "u := v",
        "v := u",
    });
    auto& ctrl_pl = params.sublist("output_control");
    ctrl_pl.set<std::string>("frequency_units","nsteps");
    ctrl_pl.set<int>("frequency",nsteps);

    OutputManager om;
    om.initialize(comm,params,t0,false);
    REQUIRE_THROWS (om.setup(fm,gm->get_grid_names()));
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
