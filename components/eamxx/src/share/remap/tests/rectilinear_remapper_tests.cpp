#include <catch2/catch.hpp>

#include "share/remap/horizontal_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_reader.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_test_utils.hpp>

// Tests remap to a generic 2d rectilinear (HRRR-like) grid. Unlike a lat-lon grid,
// the nx*ny tgt points do NOT have separable lat/lon (so nx!=nlon and ny!=nlat).

namespace scream {

namespace {

// Smooth function on the sphere, which the barycentric weights in the map file
// reproduce with an error of a few 1e-3 (the src grid is ~7 degrees wide)
Real smooth_fn (const Real lat, const Real lon) {
  constexpr Real pi = 3.14159265358979323846;
  return std::cos(lat*pi/180)*std::cos(lon*pi/180);
}

} // anonymous namespace

TEST_CASE ("rectilinear_remapper")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using gid_type = AbstractGrid::gid_type;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto& ts = ekat::TestSession::get();
  const std::string map_file = ts.params.at("map-file");
  const int nx = std::stoi(ts.params.at("nx"));
  const int ny = std::stoi(ts.params.at("ny"));

  const int nlevs = 3;

  // Create src grid, and read lat/lon from map file
  const int ncols_src = scorpio::get_dimlen(map_file,"n_a");
  const int ncols_tgt = scorpio::get_dimlen(map_file,"n_b");
  REQUIRE (ncols_tgt==nx*ny);
  auto src_grid = create_point_grid("ne4pg2",ncols_src,nlevs,comm,1);
  {
    auto deg = none.rename("deg");
    const auto& l2d = src_grid->get_2d_scalar_layout();
    auto lat = src_grid->create_geometry_data("lat",l2d,deg);
    auto lon = src_grid->create_geometry_data("lon",l2d,deg);
    std::map<FieldTag,std::string> io_dims = { {COL,"n_a"} };
    auto io_gids = src_grid->get_partitioned_dim_gids().alias("gids",io_dims);
    auto io_lat  = lat.alias("yc_a",io_dims);
    auto io_lon  = lon.alias("xc_a",io_dims);
    read_fields(map_file,{io_lat,io_lon},io_gids,comm);
  }
  scorpio::clear_unused_decomps();

  // Create src fields: a 2d one, and a 3d one
  auto create_src_field = [&](const std::string& name, const FieldLayout& fl) {
    Field f(FieldIdentifier(name,fl,none,src_grid->name()));
    f.allocate_view();
    return f;
  };
  auto f2d = create_src_field("f2d",src_grid->get_2d_scalar_layout());
  auto f3d = create_src_field("f3d",src_grid->get_3d_scalar_layout(LEV));
  {
    const int nl = src_grid->get_num_local_dofs();
    auto lat = src_grid->get_geometry_data("lat").get_view<const Real*,Host>();
    auto lon = src_grid->get_geometry_data("lon").get_view<const Real*,Host>();
    auto f2d_h = f2d.get_view<Real*,Host>();
    auto f3d_h = f3d.get_view<Real**,Host>();
    for (int i=0; i<nl; ++i) {
      f2d_h(i) = smooth_fn(lat(i),lon(i));
      for (int k=0; k<nlevs; ++k) {
        f3d_h(i,k) = (k+1)*f2d_h(i);
      }
    }
    f2d.sync_to_dev();
    f3d.sync_to_dev();
  }

  SECTION ("bad_sizes") {
    // nx*ny does not match the tgt grid
    REQUIRE_THROWS (HorizontalRemapper(src_grid,map_file,false,{nx,ny+1}));
    // Invalid sizes
    REQUIRE_THROWS (HorizontalRemapper(src_grid,map_file,false,{nx}));
    REQUIRE_THROWS (HorizontalRemapper(src_grid,map_file,false,{-nx,-ny}));
  }

  SECTION ("inconsistent_sizes_across_remappers") {
    // Data is shared across remappers with the same map file. Asking for different sizes is an error
    HorizontalRemapper r1(src_grid,map_file,false,{nx,ny});
    REQUIRE_THROWS (HorizontalRemapper(src_grid,map_file,false,{ny,nx}));
    REQUIRE_THROWS (HorizontalRemapper(src_grid,map_file,false));
  }

  SECTION ("remap") {
    HorizontalRemapper remap(src_grid,map_file,false,{nx,ny});
    auto tgt_grid = remap.get_tgt_grid();

    // The tgt grid must have the geo data to output with (y,x) layout, and NOT the one for lat-lon output
    REQUIRE (tgt_grid->get_num_global_dofs()==nx*ny);
    REQUIRE (tgt_grid->has_geometry_data("x_idx"));
    REQUIRE (tgt_grid->has_geometry_data("y_idx"));
    REQUIRE (not tgt_grid->has_geometry_data("lat_idx"));
    REQUIRE (not tgt_grid->has_geometry_data("lon_idx"));

    // lat/lon remain regular col-dependent fields (not 1d coordinate arrays)
    REQUIRE (tgt_grid->has_geometry_data("lat"));
    REQUIRE (tgt_grid->has_geometry_data("lon"));
    REQUIRE (tgt_grid->get_geometry_data("lat").get_header().get_identifier().get_layout()==tgt_grid->get_2d_scalar_layout());
    REQUIRE (tgt_grid->get_geometry_data("lon").get_header().get_identifier().get_layout()==tgt_grid->get_2d_scalar_layout());

    auto x_idx = tgt_grid->get_geometry_data("x_idx");
    auto y_idx = tgt_grid->get_geometry_data("y_idx");
    REQUIRE (x_idx.get_header().get_extra_data<int>("rectilinear_extent")==nx);
    REQUIRE (y_idx.get_header().get_extra_data<int>("rectilinear_extent")==ny);

    // Check x_idx/y_idx: x is the fastest varying index
    const int nl_tgt = tgt_grid->get_num_local_dofs();
    auto gids = tgt_grid->get_dofs_gids().get_view<const gid_type*,Host>();
    auto xi = x_idx.get_view<const int*,Host>();
    auto yi = y_idx.get_view<const int*,Host>();
    const auto min_gid = tgt_grid->get_global_min_dof_gid();
    for (int i=0; i<nl_tgt; ++i) {
      REQUIRE (xi(i)>=0);
      REQUIRE (xi(i)<nx);
      REQUIRE (yi(i)>=0);
      REQUIRE (yi(i)<ny);
      REQUIRE (yi(i)*nx+xi(i)==gids(i)-min_gid);
    }

    // Remap, and compare with exact smooth function, evaluated at tgt lat/lon
    auto t2d = remap.register_field_from_src(f2d);
    auto t3d = remap.register_field_from_src(f3d);
    remap.registration_ends();
    remap.remap_fwd();
    t2d.sync_to_host();
    t3d.sync_to_host();

    auto lat = tgt_grid->get_geometry_data("lat").get_view<const Real*,Host>();
    auto lon = tgt_grid->get_geometry_data("lon").get_view<const Real*,Host>();
    auto t2d_h = t2d.get_view<const Real*,Host>();
    auto t3d_h = t3d.get_view<const Real**,Host>();
    constexpr Real tol = 5e-3;
    for (int i=0; i<nl_tgt; ++i) {
      const Real exact = smooth_fn(lat(i),lon(i));
      REQUIRE (std::abs(t2d_h(i)-exact)<tol);
      for (int k=0; k<nlevs; ++k) {
        REQUIRE (std::abs(t3d_h(i,k)-(k+1)*exact)<(k+1)*tol);
      }
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
