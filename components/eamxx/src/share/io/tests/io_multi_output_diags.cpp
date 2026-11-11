#include <catch2/catch.hpp>

#include "share/diagnostics/abstract_diagnostic.hpp"
#include "share/diagnostics/register_diagnostics.hpp"
#include "share/data_managers/field_manager.hpp"
#include "share/data_managers/mesh_free_grids_manager.hpp"
#include "share/field/field.hpp"
#include "share/field/field_reader.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/io/scorpio_output.hpp"
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
 * A diag with several outputs, through output streams:
 *  - requests of any of its outputs, in one or more streams, share one instance,
 *    which is evaluated once per step
 *  - only the requested outputs reach the stream (and its avg counts)
 *  - outputs that share a mask share an avg count, whatever the request order
 *  - destroying the streams clears all its output names from the diag repo
 */

namespace scream {

namespace {

constexpr int ncols = 4;
constexpr int nlevs = 3;
constexpr int nsteps = 4;

constexpr Real FILL = constants::fill_value<Real>;

// On step n, at column g (global id) and level k:
//  - a = n+g+k
//  - b = 1 where the masked outputs are valid, else -1. The last column is never valid.
Real a_val (const int n, const int g, const int k) { return n+g+k; }
bool valid (const int n, const int g) { return g<ncols-1 and (n+g)%2==0; }

// Expected outputs of MOTest
Real prof_val (const int n, const int g, const int k) { return a_val(n,g,k) + (valid(n,g) ? 1 : -1); }
Real col_val  (const int n, const int g) {
  Real s = 0;
  for (int k=0; k<nlevs; ++k) {
    s += a_val(n,g,k);
  }
  return s;
}

// A diag with three outputs on different layouts, one of them masked:
//  - MOT_prof   (COL,LEV): a+b
//  - MOT_col    (COL):     the vertical sum of a
//  - MOT_masked (COL):     a at the first level, valid where b>0
class MOTest : public AbstractDiagnostic
{
public:
  MOTest (const ekat::Comm& comm, const ekat::ParameterList& params,
          const std::shared_ptr<const AbstractGrid>& grid)
   : AbstractDiagnostic(comm,params,grid)
  {
    using namespace ShortFieldTagsNames;
    using namespace ekat::units;
    m_field_in_names = {"a","b"};

    const auto& gn = m_grid->name();
    m_prof   = Field(FieldIdentifier("MOT_prof",  m_grid->get_3d_scalar_layout(LEV),K,gn),true);
    m_col    = Field(FieldIdentifier("MOT_col",   m_grid->get_2d_scalar_layout(),K,gn),true);
    m_masked = Field(FieldIdentifier("MOT_masked",m_grid->get_2d_scalar_layout(),K,gn),true);
    m_masked.create_valid_mask();
    m_masked.get_header().set_may_be_filled(true);

    add_output(m_prof);
    add_output(m_col);
    add_output(m_masked);

    ++num_instances;
  }

  std::string name () const override { return "MOTest"; }

  // Over all instances
  static int num_instances;
  static int num_evaluations;

protected:
  void compute_impl () override {
    const auto& a = m_fields_in.at("a");
    const auto& b = m_fields_in.at("b");
    a.sync_to_host();
    b.sync_to_host();
    auto a_h    = a.get_view<const Real**,Host>();
    auto b_h    = b.get_view<const Real*,Host>();
    auto prof_h = m_prof.get_view<Real**,Host>();
    auto col_h  = m_col.get_view<Real*,Host>();
    auto msk_h  = m_masked.get_view<Real*,Host>();
    auto mask_h = m_masked.get_valid_mask().get_view<int*,Host>();
    for (int icol=0; icol<m_grid->get_num_local_dofs(); ++icol) {
      col_h(icol) = 0;
      for (int ilev=0; ilev<nlevs; ++ilev) {
        prof_h(icol,ilev) = a_h(icol,ilev) + b_h(icol);
        col_h(icol) += a_h(icol,ilev);
      }
      mask_h(icol) = b_h(icol)>0 ? 1 : 0;
      msk_h(icol)  = b_h(icol)>0 ? a_h(icol,0) : FILL;
    }
    m_prof.sync_to_dev();
    m_col.sync_to_dev();
    m_masked.sync_to_dev();
    m_masked.get_valid_mask().sync_to_dev();
    ++num_evaluations;
  }

  Field m_prof, m_col, m_masked;
};
int MOTest::num_instances = 0;
int MOTest::num_evaluations = 0;

// A diag with two outputs sharing one mask, which has a custom name:
//  - MOS_a (COL): a at the first level, valid where b>0
//  - MOS_z (COL): a at the last level, valid where b>0
class MOShared : public AbstractDiagnostic
{
public:
  MOShared (const ekat::Comm& comm, const ekat::ParameterList& params,
            const std::shared_ptr<const AbstractGrid>& grid)
   : AbstractDiagnostic(comm,params,grid)
  {
    using namespace ekat::units;
    m_field_in_names = {"a","b"};

    const auto& gn = m_grid->name();
    m_a = Field(FieldIdentifier("MOS_a",m_grid->get_2d_scalar_layout(),K,gn),true);
    m_z = Field(FieldIdentifier("MOS_z",m_grid->get_2d_scalar_layout(),K,gn),true);
    m_a.create_valid_mask("mo_region_mask");
    m_z.set_valid_mask(m_a.get_valid_mask());
    for (auto* f : {&m_a,&m_z}) {
      f->get_header().set_may_be_filled(true);
      add_output(*f);
    }
  }

  std::string name () const override { return "MOShared"; }

protected:
  void compute_impl () override {
    const auto& a = m_fields_in.at("a");
    const auto& b = m_fields_in.at("b");
    a.sync_to_host();
    b.sync_to_host();
    auto a_h    = a.get_view<const Real**,Host>();
    auto b_h    = b.get_view<const Real*,Host>();
    auto fa_h   = m_a.get_view<Real*,Host>();
    auto fz_h   = m_z.get_view<Real*,Host>();
    auto mask_h = m_a.get_valid_mask().get_view<int*,Host>();
    for (int icol=0; icol<m_grid->get_num_local_dofs(); ++icol) {
      mask_h(icol) = b_h(icol)>0 ? 1 : 0;
      fa_h(icol) = b_h(icol)>0 ? a_h(icol,0) : FILL;
      fz_h(icol) = b_h(icol)>0 ? a_h(icol,nlevs-1) : FILL;
    }
    m_a.sync_to_dev();
    m_z.sync_to_dev();
    m_a.get_valid_mask().sync_to_dev();
  }

  Field m_a, m_z;
};

// A diag with an output of the same name as an output of MOTest. It is
// never evaluated: requesting it next to MOTest must fail at setup.
class MOClash : public AbstractDiagnostic
{
public:
  MOClash (const ekat::Comm& comm, const ekat::ParameterList& params,
           const std::shared_ptr<const AbstractGrid>& grid)
   : AbstractDiagnostic(comm,params,grid)
  {
    using namespace ekat::units;
    m_field_in_names = {"a"};

    const auto& gn = m_grid->name();
    add_output(Field(FieldIdentifier("MOX_a",     m_grid->get_2d_scalar_layout(),K,gn),true));
    add_output(Field(FieldIdentifier("MOT_masked",m_grid->get_2d_scalar_layout(),K,gn),true));
  }

  std::string name () const override { return "MOClash"; }

protected:
  void compute_impl () override {}
};

// Exposes the diags of a stream, and the diag repo shared by all streams
class OutputTester : public AtmosphereOutput
{
public:
  OutputTester(const ekat::Comm& comm, const ekat::ParameterList& params,
               const std::shared_ptr<const fm_type>& field_mgr, const std::string& grid_name)
    : AtmosphereOutput (comm,params,field_mgr,grid_name)
  { /* Nothing to do here */ }

  std::list<diag_ptr_type> get_diags () const { return m_diagnostics; }

  static const strmap_t<diag_ptr_type>& get_repo () { return m_diag_repo; }
};

void register_test_diags ()
{
  // For the expressions
  register_diagnostics();

  auto& factory = DiagnosticFactory::instance();
  for (const std::string n : {"MOT_prof","MOT_col","MOT_masked"}) {
    factory.register_product(n,&create_diagnostic<MOTest>);
  }
  for (const std::string n : {"MOS_a","MOS_z"}) {
    factory.register_product(n,&create_diagnostic<MOShared>);
  }
  // A name that creates MOTest, but is not one of its outputs
  factory.register_product("MOT_other",&create_diagnostic<MOTest>);
  // A diag with an output named like one of MOTest
  factory.register_product("MOX_a",&create_diagnostic<MOClash>);
}

void calc_fields (Field& a, Field& b, const AbstractGrid& grid, const int n)
{
  using gid_type = AbstractGrid::gid_type;
  auto gids_h = grid.get_dofs_gids().get_view<const gid_type*,Host>();
  auto a_h = a.get_view<Real**,Host>();
  auto b_h = b.get_view<Real*,Host>();
  for (int i=0; i<grid.get_num_local_dofs(); ++i) {
    const int g = gids_h(i);
    b_h(i) = valid(n,g) ? 1 : -1;
    for (int k=0; k<nlevs; ++k) {
      a_h(i,k) = a_val(n,g,k);
    }
  }
  a.sync_to_dev();
  b.sync_to_dev();
}

std::shared_ptr<FieldManager>
create_test_fm (const std::shared_ptr<const GridsManager>& gm,
                const util::TimeStamp& t0,
                const std::vector<std::string>& extra_fields = {})
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  auto grid = gm->get_grid("point_grid");
  auto fm = std::make_shared<FieldManager>(gm);

  Field a(FieldIdentifier("a",grid->get_3d_scalar_layout(LEV),K,grid->name()));
  Field b(FieldIdentifier("b",grid->get_2d_scalar_layout(),K,grid->name()));
  std::vector<Field> fields = {a,b};
  for (const auto& n : extra_fields) {
    fields.emplace_back(FieldIdentifier(n,grid->get_2d_scalar_layout(),K,grid->name()));
  }
  for (auto& f : fields) {
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }
  calc_fields(a,b,*grid,0);
  return fm;
}

ekat::ParameterList
stream_params (const std::string& prefix, const std::string& avg_type,
               const std::vector<std::string>& names, const std::string& grid_name,
               const int freq)
{
  ekat::ParameterList params;
  params.set<std::string>("filename_prefix",prefix);
  params.set<std::string>("averaging_type",avg_type);
  params.set<std::string>("floating_point_precision","real");
  params.sublist("fields").sublist(grid_name).set<std::vector<std::string>>("field_names",names);
  auto& ctrl_pl = params.sublist("output_control");
  ctrl_pl.set<std::string>("frequency_units","nsteps");
  ctrl_pl.set<int>("frequency",freq);
  ctrl_pl.set<bool>("save_grid_data",false);
  return params;
}

// Run one stream for nsteps steps, and return the name of the file it wrote
std::string run_stream (const ekat::ParameterList& params,
                        const std::shared_ptr<const GridsManager>& gm,
                        const ekat::Comm& comm)
{
  auto grid = gm->get_grid("point_grid");
  util::TimeStamp t0({2023,1,1},{0,0,0});
  auto fm = create_test_fm(gm,t0);

  const int dt = 1;
  auto t = t0;
  OutputManager om;
  om.initialize(comm,params,t0,false);
  om.setup(fm,gm->get_grid_names());
  for (int n=1; n<=nsteps; ++n) {
    om.init_timestep(t,dt);
    t += dt;
    auto a = fm->get_field("a");
    auto b = fm->get_field("b");
    calc_fields(a,b,*grid,n);
    for (auto f : {a,b}) {
      f.get_header().get_tracking().update_time_stamp(t);
    }
    om.run(t);
  }
  om.finalize();

  const auto& ctrl = params.sublist("output_control");
  const auto filename = params.get<std::string>("filename_prefix") + "." +
                        params.get<std::string>("averaging_type") +
                        ".nsteps_x" + std::to_string(ctrl.get<int>("frequency")) +
                        ".np" + std::to_string(comm.size()) + "." +
                        t0.to_string() + ".nc";
  std::ifstream file_check(filename);
  REQUIRE (file_check.good());
  return filename;
}

// Accumulates samples like a stream of the given averaging type does
struct Acc {
  std::string type;
  Real val = 0; int cnt = 0;
  void add (const Real v) {
    val = cnt==0 ? v : (type=="MAX" ? std::max(val,v) : type=="MIN" ? std::min(val,v) : val+v);
    ++cnt;
  }
  Real get () const { return cnt==0 ? FILL : (type=="AVERAGE" ? val/cnt : val); }
};

} // anonymous namespace

TEST_CASE ("multi_output_instant")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);
  register_test_diags();

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ncols);
  gm->build_grids();
  auto grid = gm->get_grid("point_grid");
  const auto gname = grid->name();

  // All outputs, an expression of one, and an alias of another
  MOTest::num_instances = MOTest::num_evaluations = 0;
  auto params = stream_params("io_multi_output_diags","INSTANT",
      {"MOT_prof","MOT_col","MOT_masked","d := MOT_col*2.0","e := MOT_prof"},gname,1);
  const auto filename = run_stream(params,gm,comm);

  // One instance, evaluated once per output step (t0 included)
  REQUIRE (MOTest::num_instances==1);
  REQUIRE (MOTest::num_evaluations==nsteps+1);

  // The stream was destroyed: none of the outputs is left in the repo
  for (const std::string n : {"MOT_prof","MOT_col","MOT_masked"}) {
    REQUIRE (OutputTester::get_repo().count(n)==0);
  }

  const auto l2d = grid->get_2d_scalar_layout();
  const auto l3d = grid->get_3d_scalar_layout(LEV);
  Field prof(FieldIdentifier("MOT_prof",l3d,K,gname));
  Field col (FieldIdentifier("MOT_col",l2d,K,gname));
  Field msk (FieldIdentifier("MOT_masked",l2d,K,gname));
  Field d   (FieldIdentifier("d",l2d,K,gname));
  Field e   (FieldIdentifier("e",l3d,K,gname));
  for (auto* f : {&prof,&col,&msk,&d,&e}) {
    f->allocate_view();
  }

  bool saw_valid = false, saw_invalid = false;
  {
    FieldReader reader;
    reader.set_file_specs(filename);
    reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
    reader.set_fields({prof,col,msk,d,e});

    using gid_type = AbstractGrid::gid_type;
    auto gids_h = grid->get_dofs_gids().get_view<const gid_type*,Host>();
    for (int n=0; n<=nsteps; ++n) {
      reader.read(n);
      for (auto* f : {&prof,&col,&msk,&d,&e}) {
        f->sync_to_host();
      }
      auto prof_h = prof.get_view<const Real**,Host>();
      auto col_h  = col.get_view<const Real*,Host>();
      auto msk_h  = msk.get_view<const Real*,Host>();
      auto d_h    = d.get_view<const Real*,Host>();
      auto e_h    = e.get_view<const Real**,Host>();
      for (int i=0; i<grid->get_num_local_dofs(); ++i) {
        const int g = gids_h(i);
        for (int k=0; k<nlevs; ++k) {
          REQUIRE (prof_h(i,k)==Approx(prof_val(n,g,k)));
          REQUIRE (e_h(i,k)==prof_h(i,k));
        }
        REQUIRE (col_h(i)==Approx(col_val(n,g)));
        REQUIRE (d_h(i)==Approx(2*col_val(n,g)));
        if (valid(n,g)) {
          REQUIRE (msk_h(i)==Approx(a_val(n,g,0)));
          saw_valid = true;
        } else {
          REQUIRE (msk_h(i)==FILL);
          saw_invalid = true;
        }
      }
    }
  }
  // Make sure both branches were exercised
  if (comm.size()==1) {
    REQUIRE (saw_valid);
    REQUIRE (saw_invalid);
  }

  scorpio::finalize_subsystem();
}

TEST_CASE ("multi_output_averaged")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);
  register_test_diags();

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ncols);
  gm->build_grids();
  auto grid = gm->get_grid("point_grid");
  const auto gname = grid->name();

  using gid_type = AbstractGrid::gid_type;
  auto gids_h = grid->get_dofs_gids().get_view<const gid_type*,Host>();
  const int nlocal = grid->get_num_local_dofs();

  for (const std::string avg_type : {"AVERAGE","MAX","MIN"}) {
    // All outputs
    auto params = stream_params("io_multi_output_diags_all",avg_type,
        {"MOT_prof","MOT_col","MOT_masked"},gname,nsteps);
    const auto filename = run_stream(params,gm,comm);

    // The masked output gets its own count, and is filled where never valid
    scorpio::register_file(filename,scorpio::Read);
    REQUIRE (scorpio::has_var(filename,"avg_count_MOT_masked_ncol"));
    CHECK (scorpio::has_attribute(filename,"MOT_masked","_FillValue"));
    REQUIRE_FALSE (scorpio::has_var(filename,"avg_count_MOT_col_ncol"));
    REQUIRE_FALSE (scorpio::has_var(filename,"avg_count_MOT_prof_ncol_lev"));
    scorpio::release_file(filename);

    Field prof(FieldIdentifier("MOT_prof",grid->get_3d_scalar_layout(LEV),K,gname));
    Field col (FieldIdentifier("MOT_col",grid->get_2d_scalar_layout(),K,gname));
    Field msk (FieldIdentifier("MOT_masked",grid->get_2d_scalar_layout(),K,gname));
    for (auto* f : {&prof,&col,&msk}) {
      f->allocate_view();
    }
    FieldReader reader;
    reader.set_file_specs(filename);
    reader.set_dim_decomp(grid->get_partitioned_dim_gids(),comm);
    reader.set_fields({prof,col,msk});
    reader.read(0);
    for (auto* f : {&prof,&col,&msk}) {
      f->sync_to_host();
    }
    auto prof_h = prof.get_view<const Real**,Host>();
    auto col_h  = col.get_view<const Real*,Host>();
    auto msk_h  = msk.get_view<const Real*,Host>();

    bool saw_partial = false, saw_empty = false;
    for (int i=0; i<nlocal; ++i) {
      const int g = gids_h(i);
      Acc acc_col{avg_type}, acc_msk{avg_type};
      for (int n=1; n<=nsteps; ++n) {
        acc_col.add(col_val(n,g));
        if (valid(n,g)) {
          acc_msk.add(a_val(n,g,0));
        }
      }
      REQUIRE (col_h(i)==Approx(acc_col.get()));
      if (acc_msk.cnt==0) {
        REQUIRE (msk_h(i)==FILL);
        saw_empty = true;
      } else {
        REQUIRE (msk_h(i)==Approx(acc_msk.get()));
        REQUIRE (acc_msk.cnt<nsteps);
        saw_partial = true;
      }
      for (int k=0; k<nlevs; ++k) {
        Acc acc_prof{avg_type};
        for (int n=1; n<=nsteps; ++n) {
          acc_prof.add(prof_val(n,g,k));
        }
        REQUIRE (prof_h(i,k)==Approx(acc_prof.get()));
      }
    }
    // Make sure both branches were exercised
    if (comm.size()==1) {
      REQUIRE (saw_partial);
      REQUIRE (saw_empty);
    }

    // Only one output: the masked one was not requested, so there is no count at all
    auto params_col = stream_params("io_multi_output_diags_col",avg_type,
        {"MOT_col"},gname,nsteps);
    const auto filename_col = run_stream(params_col,gm,comm);
    scorpio::register_file(filename_col,scorpio::Read);
    REQUIRE (scorpio::has_var(filename_col,"MOT_col"));
    REQUIRE_FALSE (scorpio::has_var(filename_col,"MOT_masked"));
    REQUIRE_FALSE (scorpio::has_var(filename_col,"avg_count_ncol"));
    REQUIRE_FALSE (scorpio::has_var(filename_col,"avg_count_MOT_masked_ncol"));
    scorpio::release_file(filename_col);

    // Outputs sharing a mask share a count. Each must take it from the mask:
    // the second one requested must not get a count of its own
    auto params_shared = stream_params("io_multi_output_diags_shared",avg_type,
        {"MOS_a","MOS_z"},gname,nsteps);
    const auto filename_shared = run_stream(params_shared,gm,comm);
    scorpio::register_file(filename_shared,scorpio::Read);
    REQUIRE (scorpio::has_var(filename_shared,"avg_count_mo_region_ncol"));
    REQUIRE_FALSE (scorpio::has_var(filename_shared,"avg_count_MOS_a_ncol"));
    REQUIRE_FALSE (scorpio::has_var(filename_shared,"avg_count_MOS_z_ncol"));
    REQUIRE_FALSE (scorpio::has_var(filename_shared,"avg_count_ncol"));
    scorpio::release_file(filename_shared);
  }

  scorpio::finalize_subsystem();
}

TEST_CASE ("multi_output_sharing")
{
  // NOTE: like diags_sharing_check in io_diags.cpp, this test does no scorpio call:
  //       AtmosphereOutput objects do not open files until run() writes
  ekat::Comm comm(MPI_COMM_WORLD);
  register_test_diags();

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ncols);
  gm->build_grids();
  auto grid = gm->get_grid("point_grid");
  const auto gname = grid->name();

  auto time = util::TimeStamp({2023,1,1},{0,0,0});
  auto fm = create_test_fm(gm,time);

  auto params_for = [&](const std::vector<std::string>& names) {
    ekat::ParameterList params;
    params.set<std::string>("averaging_type","average");
    params.sublist("fields").sublist(gname).set<std::vector<std::string>>("field_names",names);
    return params;
  };
  const auto& repo = OutputTester::get_repo();
  auto in_repo = [&](const std::string& n) { return repo.count(n)==1; };

  SECTION ("across_streams") {
    MOTest::num_instances = MOTest::num_evaluations = 0;
    {
      // Different outputs in different streams: one instance, evaluated once
      OutputTester out1(comm,params_for({"MOT_col","MOT_prof"}),fm,gname);
      OutputTester out2(comm,params_for({"MOT_masked"}),fm,gname);

      auto diags1 = out1.get_diags();
      auto diags2 = out2.get_diags();
      REQUIRE (diags1.size()==1);
      REQUIRE (diags2.size()==1);
      REQUIRE (diags1.front()==diags2.front());
      REQUIRE (MOTest::num_instances==1);

      out1.init_timestep(time);
      out2.init_timestep(time);
      time += 3600;
      fm->get_field("a").get_header().get_tracking().update_time_stamp(time);
      fm->get_field("b").get_header().get_tracking().update_time_stamp(time);
      out1.run("UNUSED", time, false, false, 0, false);
      out2.run("UNUSED", time, false, false, 0, false);
      REQUIRE (MOTest::num_evaluations==1);

      // All outputs are stored under their name, also the one no stream requested
      for (const std::string n : {"MOT_prof","MOT_col","MOT_masked"}) {
        REQUIRE (in_repo(n));
        REQUIRE (repo.at(n)==diags1.front());
      }
    }
    // The streams are gone, and so are all the names of their diags
    for (const std::string n : {"MOT_prof","MOT_col","MOT_masked"}) {
      REQUIRE_FALSE (in_repo(n));
    }

    // A new stream builds a new instance
    OutputTester out3(comm,params_for({"MOT_masked"}),fm,gname);
    REQUIRE (MOTest::num_instances==2);
    REQUIRE (in_repo("MOT_prof"));
  }

  SECTION ("unrequested_output_in_repo") {
    // One stream requests one output: the others are stored too, but not in its FM
    MOTest::num_instances = 0;
    OutputTester out1(comm,params_for({"MOT_col"}),fm,gname);
    REQUIRE (out1.get_diags().size()==1);
    REQUIRE (in_repo("MOT_masked"));

    // ...so another stream requesting another output reuses the instance
    OutputTester out2(comm,params_for({"MOT_masked","MOT_prof"}),fm,gname);
    REQUIRE (MOTest::num_instances==1);
    REQUIRE (out2.get_diags().size()==1);
    REQUIRE (out1.get_diags().front()==out2.get_diags().front());
  }

  SECTION ("alias_named_like_an_output") {
    // An alias is local to its stream, and may have the name of an output the
    // stream does not request. Names are processed in alphabetical order, so the
    // alias MOT_col is processed before the diag, and the alias MOT_prof after it
    MOTest::num_instances = 0;
    for (const auto& names : {std::vector<std::string>{"MOT_col := a","MOT_prof"},
                              std::vector<std::string>{"MOT_prof := a","MOT_col"}}) {
      OutputTester out(comm,params_for(names),fm,gname);
      REQUIRE (out.get_diags().size()==1);
      REQUIRE (out.get_diags().front()->name()=="MOTest");
    }
    REQUIRE (MOTest::num_instances==2);
  }

  SECTION ("errors") {
    // A name that resolves to MOTest, but is not one of its outputs
    REQUIRE_THROWS_WITH (OutputTester(comm,params_for({"MOT_other"}),fm,gname),
        Catch::Contains("not an output") and Catch::Contains("MOT_masked"));
    REQUIRE_FALSE (in_repo("MOT_other"));
    REQUIRE_FALSE (in_repo("MOT_prof"));

    // An output with the name of a model field would shadow it
    auto fm_clash = create_test_fm(gm,time,{"MOT_col"});
    REQUIRE_THROWS_WITH (OutputTester(comm,params_for({"MOT_prof"}),fm_clash,gname),
        Catch::Contains("has the name of a model field") and Catch::Contains("field name: MOT_col"));
    REQUIRE_FALSE (in_repo("MOT_prof"));
    REQUIRE_FALSE (in_repo("MOT_masked"));

    // Two diags with an output of the same name. MOT_col comes first, so
    // MOTest is stored under MOT_masked before MOClash is created
    REQUIRE_THROWS_WITH (OutputTester(comm,params_for({"MOT_col","MOX_a"}),fm,gname),
        Catch::Contains("Two diagnostics produce a field with the same name") and
        Catch::Contains("field name: MOT_masked") and
        Catch::Contains("MOTest") and Catch::Contains("MOClash"));
    for (const std::string n : {"MOT_prof","MOT_col","MOT_masked","MOX_a"}) {
      REQUIRE_FALSE (in_repo(n));
    }

    // A stream that fails after storing the outputs (an expression to write
    // needs a name) leaves none of them in the repo
    REQUIRE_THROWS_WITH (OutputTester(comm,params_for({"MOT_prof","b*2.0"}),fm,gname),
        Catch::Contains("An expression must be given an output name"));
    for (const std::string n : {"MOT_prof","MOT_col","MOT_masked"}) {
      REQUIRE_FALSE (in_repo(n));
    }
  }
}

TEST_CASE ("number_paths")
{
  // NumberPath computes all three number paths: requesting two gives one diag
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);
  register_diagnostics();

  auto gm = create_mesh_free_grids_manager(comm,0,0,nlevs,ncols);
  gm->build_grids();
  auto grid = gm->get_grid("point_grid");
  const auto gname = grid->name();

  auto t0 = util::TimeStamp({2023,1,1},{0,0,0});
  auto fm = std::make_shared<FieldManager>(gm);
  for (const std::string n : {"pseudo_density","qc","nc","qi","ni","qr","nr"}) {
    Field f(FieldIdentifier(n,grid->get_3d_scalar_layout(LEV),ekat::units::none,gname));
    f.allocate_view();
    f.get_header().get_tracking().update_time_stamp(t0);
    fm->add_field(f);
  }

  ekat::ParameterList params;
  params.set<std::string>("averaging_type","average");
  params.sublist("fields").sublist(gname).set<std::vector<std::string>>("field_names",
      {"LiqNumberPath","RainNumberPath"});
  OutputTester out(comm,params,fm,gname);

  auto diags = out.get_diags();
  REQUIRE (diags.size()==1);
  REQUIRE (diags.front()->name()=="NumberPath");
  REQUIRE (OutputTester::get_repo().at("IceNumberPath")==diags.front());
}

} // namespace scream
