#include <catch2/catch.hpp>

#include "share/remap/horizontal_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <numeric>

namespace scream {

namespace {

void root_print (const std::string& msg, const ekat::Comm& comm) {
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

// A coarsening map: every tgt dof averages two consecutive src dofs.
void create_remap_file (const std::string& filename, const int ngdofs_tgt)
{
  const int ngdofs_src = ngdofs_tgt + 1;
  const int nnz = 2*ngdofs_tgt;

  scorpio::register_file(filename, scorpio::FileMode::Write);

  scorpio::define_dim(filename,"n_a", ngdofs_src);
  scorpio::define_dim(filename,"n_b", ngdofs_tgt);
  scorpio::define_dim(filename,"n_s", nnz);

  scorpio::define_var(filename,"col",{"n_s"},"int");
  scorpio::define_var(filename,"row",{"n_s"},"int");
  scorpio::define_var(filename,"S"  ,{"n_s"},"double");

  scorpio::enddef(filename);

  std::vector<int> col(nnz), row(nnz);
  std::vector<double> S(nnz,0.5);
  const int gid_base = 1;
  for (int i=0; i<ngdofs_tgt; ++i) {
    row[2*i]   = gid_base + i;
    row[2*i+1] = gid_base + i;
    col[2*i]   = gid_base + i;
    col[2*i+1] = gid_base + i+1;
  }

  scorpio::write_var(filename,"row",row.data());
  scorpio::write_var(filename,"col",col.data());
  scorpio::write_var(filename,"S",    S.data());

  scorpio::release_file(filename);
}

std::shared_ptr<AbstractGrid>
build_src_grid (const ekat::Comm& comm, const int ngdofs)
{
  using gid_type = AbstractGrid::gid_type;
  const int nlevs = 20;

  int nldofs   = ngdofs / comm.size();
  int remainder = ngdofs % comm.size();
  int offset   = nldofs * comm.rank() + std::min(comm.rank(),remainder);
  if (comm.rank()<remainder) {
    ++nldofs;
  }

  auto grid = std::make_shared<PointGrid>("src",nldofs,nlevs,comm);
  auto dofs_h = grid->get_dofs_gids().get_view<gid_type*,Host>();
  std::iota(dofs_h.data(),dofs_h.data()+nldofs,offset+1);
  grid->get_dofs_gids().sync_to_dev();

  return grid;
}

// Remap f=gid: every tgt dof must get the average of its two src dofs
void check_remap (const std::shared_ptr<HorizontalRemapper>& remap)
{
  using gid_type = AbstractGrid::gid_type;

  auto src_grid = remap->get_src_grid();
  auto tgt_grid = remap->get_tgt_grid();

  FieldIdentifier fid("f",src_grid->get_2d_scalar_layout(),ekat::units::none,src_grid->name());
  Field src(fid);
  src.allocate_view();
  auto src_h = src.get_view<Real*,Host>();
  auto src_gids_h = src_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  for (int i=0; i<src_grid->get_num_local_dofs(); ++i) {
    src_h(i) = src_gids_h(i);
  }
  src.sync_to_dev();

  auto tgt = remap->register_field_from_src(src);
  remap->registration_ends();
  remap->remap_fwd();

  tgt.sync_to_host();
  auto tgt_h = tgt.get_view<const Real*,Host>();
  auto tgt_gids_h = tgt_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  for (int i=0; i<tgt_grid->get_num_local_dofs(); ++i) {
    REQUIRE (tgt_h(i)==tgt_gids_h(i)+0.5);
  }
}

// Mimic VerticalRemapper::create_tgt_grid: same GIDs, different object
std::shared_ptr<AbstractGrid>
clone_as_vremap_tgt (const std::shared_ptr<AbstractGrid>& src_grid)
{
  auto g = src_grid->clone(src_grid->name()+"_vremap_tgt",true);
  g->reset_vertical_configuration(10, AbstractGrid::VKind::Pressure);
  return g;
}

} // anonymous namespace

// Two remappers created from the same map file share the remap data. Both must
// map src->tgt, regardless of which is created first, and even if the grid of the
// second is a clone of the src grid.
// NOTE: both remappers must be alive at the same time, or the data is not shared.
TEST_CASE("horiz_remap_shared_map_file")
{
  ekat::Comm comm(MPI_COMM_WORLD);

  root_print ("\n +-------------------------------------------+\n",comm);
  root_print (" |   Horiz remap: two remappers, one map     |\n",comm);
  root_print (" +-------------------------------------------+\n\n",comm);

  // NOTE: the TEST_CASE body runs once per SECTION
  if (not scorpio::is_subsystem_inited()) {
    scorpio::init_subsystem(comm);
  }

  const int nldofs_tgt = 3;
  const int ngdofs_tgt = nldofs_tgt*comm.size();
  const int ngdofs_src = ngdofs_tgt+1;

  std::string filename = "hr_shared_map_tests." + std::to_string(comm.size()) + ".nc";
  create_remap_file(filename, ngdofs_tgt);

  auto src_grid = build_src_grid(comm, ngdofs_src);

  SECTION ("plain first, clone second") {
    auto remap_plain = std::make_shared<HorizontalRemapper>(src_grid,filename);
    REQUIRE (remap_plain->get_tgt_grid()->get_num_global_dofs()==ngdofs_tgt);

    // Keep alive, so the remap data is shared
    auto vclone = clone_as_vremap_tgt(src_grid);
    auto remap_clone = std::make_shared<HorizontalRemapper>(vclone,filename);

    REQUIRE (remap_clone->get_src_grid()->get_num_global_dofs()==ngdofs_src);
    REQUIRE (remap_clone->get_tgt_grid()->get_num_global_dofs()==ngdofs_tgt);
    // The tgt grid has the vertical configuration of the grid the remapper was built from
    REQUIRE (remap_plain->get_tgt_grid()->get_num_vertical_levels()==20);
    REQUIRE (remap_clone->get_tgt_grid()->get_num_vertical_levels()==10);

    check_remap(remap_plain);
    check_remap(remap_clone);
  }

  SECTION ("clone first, plain second") {
    auto vclone = clone_as_vremap_tgt(src_grid);
    auto remap_clone = std::make_shared<HorizontalRemapper>(vclone,filename);
    REQUIRE (remap_clone->get_tgt_grid()->get_num_global_dofs()==ngdofs_tgt);

    auto remap_plain = std::make_shared<HorizontalRemapper>(src_grid,filename);

    REQUIRE (remap_plain->get_src_grid()->get_num_global_dofs()==ngdofs_src);
    REQUIRE (remap_plain->get_tgt_grid()->get_num_global_dofs()==ngdofs_tgt);

    check_remap(remap_clone);
    check_remap(remap_plain);
  }

  if (scorpio::is_subsystem_inited()) {
    scorpio::finalize_subsystem();
  }
}

} // namespace scream
