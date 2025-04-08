#include <catch2/catch.hpp>

#include "share/io/scorpio_scm_input.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

#include "share/grid/point_grid.hpp"
#include "share/field/field.hpp"
#include "share/field/field_utils.hpp"

#include "share/util/eamxx_setup_random_test.hpp"

namespace scream {

// Returns fields after initialization
void write (const int seed, const ekat::Comm& comm)
{
  using ekat::units::Units;
  using RPDF = std::uniform_real_distribution<Real>;
  using IPDF = std::uniform_int_distribution<int>;
  using Engine = std::mt19937_64;

  Engine engine(seed);

  RPDF lat_pdf(-90.0,90.0);
  RPDF lon_pdf(-180.0,180.0);

  int nlevs = IPDF(5,8)(engine);
  int ncols = IPDF(4,5)(engine);

  comm.broadcast(&ncols,1,comm.root_rank());
  comm.broadcast(&nlevs,1,comm.root_rank());

  // Create grid
  auto grid = create_point_grid("test",ncols,nlevs,comm);

  // Create lat/lon grid
  Units deg (Units::nondimensional(),"deg");
  auto lat = grid->create_geometry_data("lat",grid->get_2d_scalar_layout(),deg);
  auto lon = grid->create_geometry_data("lon",grid->get_2d_scalar_layout(),deg);
  randomize(lat,engine,lat_pdf);
  randomize(lon,engine,lon_pdf);

  // Create variable data
  FieldIdentifier fid("var",grid->get_3d_scalar_layout(true),Units::nondimensional(),"");
  Field var(fid);
  var.allocate_view();
  randomize(var,engine,RPDF(-1.0,1.0));

  // Create file
  auto filename = "io_scm_np" + std::to_string(comm.size()) + ".nc";
  scorpio::register_file(filename,scorpio::Write,scorpio::DefaultIOType);
  
  scorpio::define_dim(filename,"ncol",ncols);
  scorpio::define_dim(filename,"lev",nlevs);

  scorpio::define_var(filename,"lat",{"ncol"},"real");
  scorpio::define_var(filename,"lon",{"ncol"},"real");

  scorpio::define_var(filename,"var",{"ncol","lev"},"real");

  auto my_col_gids = grid->get_partitioned_dim_gids().get_view<const AbstractGrid::gid_type*,Host>();
  std::vector<scorpio::offset_t> my_offsets(my_col_gids.size());
  for (size_t i=0; i<my_offsets.size(); ++i) {
    my_offsets[i] = my_col_gids[i] - grid->get_global_min_partitioned_dim_gid ();
  }
  scorpio::set_dim_decomp(filename,"ncol",my_offsets);
  scorpio::enddef(filename);

  // Write to file
  scorpio::write_var(filename,"lat",lat.get_internal_view_data<Real,Host>());
  scorpio::write_var(filename,"lon",lon.get_internal_view_data<Real,Host>());
  scorpio::write_var(filename,"var",var.get_internal_view_data<Real,Host>());

  scorpio::release_file(filename);
}

void read (const int seed, const ekat::Comm& comm)
{
  using ekat::units::Units;
  using IPDF = std::uniform_int_distribution<int>;
  using Engine = std::mt19937_64;

  // Pick a col index, and find its lat/lon from file
  auto filename = "io_scm_np" + std::to_string(comm.size()) + ".nc";
  scorpio::register_file(filename,scorpio::Read,scorpio::DefaultIOType);

  int ncols = scorpio::get_dimlen(filename,"ncol");
  int nlevs = scorpio::get_dimlen(filename,"lev");

  Engine engine(seed);

  // Read lat/lon/var on ALL ranks. Rank 0 decides what's the tgt lat/lon
  std::vector<Real> lat(ncols), lon(ncols), var(ncols*nlevs);
  scorpio::read_var(filename,"lat",lat.data());
  scorpio::read_var(filename,"lon",lon.data());
  scorpio::read_var(filename,"var",var.data());

  auto tgt_col = IPDF(0,ncols-1)(engine);
  comm.broadcast(&tgt_col,1,comm.root_rank());

  auto tgt_lat = lat[tgt_col];
  auto tgt_lon = lon[tgt_col];

  // Create single-col grid
  auto grid = std::make_shared<PointGrid>("scm_grid",1,nlevs,comm);
  auto dofs_gids = grid->get_dofs_gids();
  dofs_gids.deep_copy(0);

  // Create field to read
  FieldIdentifier fid("var",grid->get_3d_scalar_layout(true),Units::nondimensional(),"");
  Field var_f(fid);
  var_f.allocate_view();

  // Read field
  SCMInput reader(filename,tgt_lat,tgt_lon,{var_f},comm);
  reader.read_variables();

  // Check
  auto var_h = var_f.get_view<const Real**,Host>();
  for (int ilev=0; ilev<nlevs; ++ilev) {
    REQUIRE (var_h(0,ilev)==var[tgt_col*nlevs+ilev]);
  }

  scorpio::release_file(filename);
}

TEST_CASE ("scm_io") {
  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto seed = get_random_test_seed(&comm);

  write(seed,comm);
  read(seed,comm);

  scorpio::finalize_subsystem();
}

} // anonymous namespace
