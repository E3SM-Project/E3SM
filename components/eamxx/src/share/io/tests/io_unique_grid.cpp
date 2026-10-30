#include <catch2/catch.hpp>

#include "share/io/eamxx_output_manager.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include "share/data_managers/field_manager.hpp"
#include "share/field/field_reader.hpp"
#include "share/field/field_utils.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/se_grid.hpp"

#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_comm.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_units.hpp>

#include <algorithm>
#include <numeric>

namespace {

using namespace scream;
using namespace ShortFieldTagsNames;

using gid_type = AbstractGrid::gid_type;

constexpr int num_elems = 3;
constexpr int num_gp    = 4;
constexpr int num_levs  = 2*SCREAM_PACK_SIZE+1;
constexpr int vec_dim   = 2;

/*
 * In this test we do a SE->unique grid remap "on the fly" during the output phase
 * (via the output_data_layout option), for fields defined on a SE grid.
 * Then, we load the generated file, and compare the values against the ones
 * we expect. Since the values of the SE fields are a function of the unique gid
 * of each dof, we know what the unique grid fields should contain.
 */

// A grids manager where we can stuff grids
class TestGridsManager : public GridsManager {
public:
  std::string name () const override { return "Test grids manager"; }
  void build_grids () override {}
  void add (const std::shared_ptr<AbstractGrid>& g) { add_nonconst_grid(g); }
};

// The unique gid of the gp (ie,ip,jp) on this rank. Elements form a strip,
// where consecutive elements share a edge: ip=num_gp-1 on element ie is
// the same dof as ip=0 on element ie+1. Different ranks have disjoint dofs.
int num_unique_cols_per_rank () { return (num_elems*(num_gp-1)+1)*num_gp; }
gid_type unique_gid (const int rank, const int ie, const int ip, const int jp)
{
  return rank*num_unique_cols_per_rank() + (ie*(num_gp-1)+ip)*num_gp + jp;
}

// The value stored in a field at a given gid. The weights make sure that
// different components/levels contain different values
// (we assume icmp<10 and ilev<100)
Real value (const gid_type gid, const int icmp, const int ilev)
{
  // NOTE: these are integers (stored as Real), so that they are exact, no matter how the
  //       compiler contracts/vectorizes the arithmetic at different call sites
  return Real(gid)*1000 + icmp*100 + ilev;
}

std::shared_ptr<SEGrid> build_se_grid (const ekat::Comm& comm, Field& unique_gids)
{
  using namespace ShortFieldTagsNames;
  auto se_grid = std::make_shared<SEGrid>("se",num_elems,num_gp,num_levs,comm);

  auto dofs_gids   = se_grid->get_dofs_gids();
  unique_gids = Field(FieldIdentifier("unique_grid_gids",FieldLayout({CMP},{se_grid->get_num_local_dofs()}),
                                      ekat::units::none,se_grid->name(),DataType::IntType));
  unique_gids.allocate_view();
  auto elem_gids = se_grid->get_partitioned_dim_gids();
  auto lid2idx   = se_grid->get_lid_to_idx_map();

  auto h_dofs   = dofs_gids.get_view<gid_type*,Host>();
  auto h_unique   = unique_gids.get_view<gid_type*,Host>();
  auto h_elem = elem_gids.get_view<gid_type*,Host>();
  auto h_l2i  = lid2idx.get_view<int**,Host>();

  const int rank = comm.rank();
  for (int ie=0; ie<num_elems; ++ie) {
    h_elem(ie) = ie + num_elems*rank;
    for (int ip=0; ip<num_gp; ++ip) {
      for (int jp=0; jp<num_gp; ++jp) {
        const int idof = ie*num_gp*num_gp + ip*num_gp + jp;
        h_dofs(idof) = rank*num_elems*num_gp*num_gp + idof;
        h_unique(idof) = unique_gid(rank,ie,ip,jp);
        h_l2i(idof,0) = ie;
        h_l2i(idof,1) = ip;
        h_l2i(idof,2) = jp;
      }
    }
  }
  dofs_gids.sync_to_dev();
  unique_gids.sync_to_dev();
  elem_gids.sync_to_dev();
  lid2idx.sync_to_dev();

  return se_grid;
}

// A point grid with all the unique dofs of the SE grid
std::shared_ptr<PointGrid> build_unique_grid (const ekat::Comm& comm)
{
  auto unique_grid = std::make_shared<PointGrid>("unique",num_unique_cols_per_rank(),num_levs,comm);

  auto unique_gids = unique_grid->get_dofs_gids();
  auto h_gids  = unique_gids.get_view<gid_type*,Host>();
  std::iota(h_gids.data(),h_gids.data()+num_unique_cols_per_rank(),comm.rank()*num_unique_cols_per_rank());
  unique_gids.sync_to_dev();

  return unique_grid;
}

TEST_CASE("io_unique_grid")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  // Initialize the pio_subsystem for this test:
  scorpio::init_subsystem(comm);

  Field se_unique_gids;
  auto se_grid = build_se_grid(comm,se_unique_gids);
  auto unique_grid = build_unique_grid(comm);
  se_grid->set_unique_grid(unique_grid,se_unique_gids);

  auto gm = std::make_shared<TestGridsManager>();
  gm->add(se_grid);
  gm->add(unique_grid);

  // Create the SE fields
  auto nondim = ekat::units::none;
  const int ps = SCREAM_PACK_SIZE;
  FieldIdentifier fid_1 ("field_1",se_grid->get_3d_scalar_layout(LEV),nondim,se_grid->name());
  FieldIdentifier fid_2 ("field_2",se_grid->get_3d_vector_layout(LEV,vec_dim,"vec"),nondim,se_grid->name());
  FieldIdentifier fid_3 ("field_3",se_grid->get_2d_scalar_layout(),nondim,se_grid->name());

  auto fm = std::make_shared<FieldManager>(gm);
  fm->register_field(FieldRequest(fid_1,ps));
  fm->register_field(FieldRequest(fid_2,ps));
  fm->register_field(FieldRequest(fid_3));
  fm->registration_ends();

  util::TimeStamp t0({2000,1,1},{0,0,0});
  fm->init_fields_time_stamp(t0);

  // Fill the SE fields with values that depend on the unique gid of each dof
  const auto se_gids = se_grid->get_unique_grid_gids().get_view<const gid_type*,Host>();
  const auto l2i     = se_grid->get_lid_to_idx_map().get_view<const int**,Host>();
  {
    auto f1 = fm->get_field("field_1",se_grid->name());
    auto f2 = fm->get_field("field_2",se_grid->name());
    auto f3 = fm->get_field("field_3",se_grid->name());
    auto h1 = f1.get_view<Real****,Host>();
    auto h2 = f2.get_view<Real*****,Host>();
    auto h3 = f3.get_view<Real***,Host>();
    for (int idof=0; idof<se_grid->get_num_local_dofs(); ++idof) {
      const int ie = l2i(idof,0);
      const int ip = l2i(idof,1);
      const int jp = l2i(idof,2);
      const auto gid = se_gids(idof);
      h3(ie,ip,jp) = value(gid,0,0);
      for (int il=0; il<num_levs; ++il) {
        h1(ie,ip,jp,il) = value(gid,0,il);
        for (int ic=0; ic<vec_dim; ++ic) {
          h2(ie,ic,ip,jp,il) = value(gid,ic,il);
        }
      }
    }
    for (auto f : {f1,f2,f3}) {
      f.sync_to_dev();
    }
  }

  // Write all fields to file, remapping them to the cg grid
  std::vector<std::string> fnames = {"field_1", "field_2", "field_3"};
  ekat::ParameterList out_params;
  out_params.set<std::string>("averaging_type","instant");
  out_params.set<std::string>("filename_prefix","io_unique_grid");
  out_params.sublist("fields").sublist(se_grid->name()).set<std::vector<std::string>>("field_names",fnames);
  out_params.sublist("fields").sublist(se_grid->name()).set<std::string>("output_data_layout","unique");

  out_params.sublist("output_control").set<int>("frequency",1);
  out_params.sublist("output_control").set<std::string>("frequency_units","nsteps");
  out_params.set<std::string>("floating_point_precision","real");

  OutputManager output;
  output.initialize(comm, out_params, t0, false);
  output.setup (fm, {se_grid->name()});
  output.run(t0);
  output.finalize();

  // Next, load all fields from file into fields on the cg grid
  std::string filename = "io_unique_grid.INSTANT.nsteps_x1.np" + std::to_string(comm.size()) + "." + t0.to_string() + ".nc";
  filename.erase(std::remove(filename.begin(),filename.end(),':'),filename.end());

  auto make_unique_field = [&](const FieldIdentifier& fid_se) {
    const auto& l = fid_se.get_layout();
    FieldLayout cg_layout = unique_grid->equivalent_layout(l);
    Field f(FieldIdentifier(fid_se.name(),cg_layout,nondim,unique_grid->name()));
    f.allocate_view();
    return f;
  };
  auto r1 = make_unique_field(fid_1);
  auto r2 = make_unique_field(fid_2);
  auto r3 = make_unique_field(fid_3);
  read_fields (filename,{r1,r2,r3},unique_grid->get_partitioned_dim_gids(),comm);

  // Build the expected fields
  auto e1 = r1.clone("e1");
  auto e2 = r2.clone("e2");
  auto e3 = r3.clone("e3");
  {
    const auto unique_gids = unique_grid->get_dofs_gids().get_view<const gid_type*,Host>();
    auto h1 = e1.get_view<Real**,Host>();
    auto h2 = e2.get_view<Real***,Host>();
    auto h3 = e3.get_view<Real*,Host>();
    for (int icol=0; icol<unique_grid->get_num_local_dofs(); ++icol) {
      const auto gid = unique_gids(icol);
      h3(icol) = value(gid,0,0);
      for (int il=0; il<num_levs; ++il) {
        h1(icol,il) = value(gid,0,il);
        for (int ic=0; ic<vec_dim; ++ic) {
          h2(icol,ic,il) = value(gid,ic,il);
        }
      }
    }
    for (auto f : {e1,e2,e3}) {
      f.sync_to_dev();
    }
  }

  // NOTE: CHECK (rather than REQUIRE), so that all ranks reach the scorpio finalize call below
  CHECK (views_are_equal(r1,e1));
  CHECK (views_are_equal(r2,e2));
  CHECK (views_are_equal(r3,e3));

  // Cleanup everything
  scorpio::finalize_subsystem();
}

} // anonymous namespace
