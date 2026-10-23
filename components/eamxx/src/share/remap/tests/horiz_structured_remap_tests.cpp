#include <catch2/catch.hpp>

#include "share/remap/horizontal_remapper.hpp"
#include "share/grid/point_grid.hpp"
#include "share/field/field_reader.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_test_utils.hpp>

// Tests remap to structured grids, which are autodetected from the map file (dst_grid_rank=2).
// A structured grid is either
//  - lat-lon: lat (lon) only depends on y (x). Here lat/lon are 1d coordinate arrays.
//  - curvilinear (e.g., HRRR-like): lat and lon depend on both x and y. Here lat/lon are col-dependent fields.
// The curvilinear test map has lat/lon values that are all different, so that nx!=nlon and ny!=nlat.

namespace scream {

namespace {

// Some smooth function on the sphere, to use as src field
Real smooth_fn (const Real lat, const Real lon) {
  constexpr Real pi = 3.14159265358979323846;
  return std::cos(lat*pi/180)*std::cos(lon*pi/180);
}

// Get [ny,nx] of a structured grid, as stored by the remapper in the geo data structured_dims
struct StructuredDims { int ny, nx; };
StructuredDims get_structured_dims (const AbstractGrid& grid)
{
  auto h = grid.get_geometry_data("structured_dims").get_view<const int*,Host>();
  return {h(0),h(1)};
}

// Read the whole var from the file (on all ranks). The var must not be decomposed.
template<typename T>
std::vector<T> read_full_var (const std::string& filename, const std::string& varname, const int n)
{
  const bool was_open = scorpio::is_file_open(filename);
  if (not was_open) {
    scorpio::register_file(filename,scorpio::Read);
  }
  std::vector<T> v(n);
  scorpio::read_var(filename,varname,v.data());
  if (not was_open) {
    scorpio::release_file(filename);
  }
  return v;
}

// Create a map file where the tgt grid is nx*ny points, and the src grid has one more point.
// The first nx*ny src points are mapped to the tgt points with weight 1. The tgt points
// are lat-lon, with lat_j=-60+20*j and lon_i=10+30*i, but they are stored with either
//  - lon fastest (lat-major): gid = j*nx + i (standard, a lat-lon grid)
//  - lat fastest (lon-major): gid = i*ny + j (this breaks the assumption that x is the fastest index)
void create_latlon_map (const std::string& filename, const int nx, const int ny, const bool lat_major)
{
  const int n_b = nx*ny;
  const int n_a = n_b+1;

  scorpio::register_file(filename, scorpio::FileMode::Write);
  scorpio::define_dim(filename, "n_a", n_a);
  scorpio::define_dim(filename, "n_b", n_b);
  scorpio::define_dim(filename, "n_s", n_b);
  scorpio::define_dim(filename, "dst_grid_rank", 2);

  scorpio::define_var(filename, "dst_grid_dims", {"dst_grid_rank"}, "int");
  scorpio::define_var(filename, "xc_b",   {"n_b"}, "double");
  scorpio::define_var(filename, "yc_b",   {"n_b"}, "double");
  scorpio::define_var(filename, "area_b", {"n_b"}, "double");
  scorpio::define_var(filename, "col", {"n_s"}, "int");
  scorpio::define_var(filename, "row", {"n_s"}, "int");
  scorpio::define_var(filename, "S",   {"n_s"}, "double");
  scorpio::enddef(filename);

  std::vector<double> lat(n_b), lon(n_b), area(n_b,1.0), S(n_b,1.0);
  std::vector<int> row(n_b), col(n_b);
  for (int j=0; j<ny; ++j) {
    for (int i=0; i<nx; ++i) {
      const int gid = lat_major ? j*nx + i : i*ny + j;
      lat[gid] = -60 + 20*j;
      lon[gid] =  10 + 30*i;
    }
  }
  for (int i=0; i<n_b; ++i) {
    row[i] = col[i] = i+1;
  }
  int dims[2] = {nx,ny};
  scorpio::write_var(filename,"dst_grid_dims",dims);
  scorpio::write_var(filename,"xc_b",lon.data());
  scorpio::write_var(filename,"yc_b",lat.data());
  scorpio::write_var(filename,"area_b",area.data());
  scorpio::write_var(filename,"row",row.data());
  scorpio::write_var(filename,"col",col.data());
  scorpio::write_var(filename,"S",S.data());
  scorpio::release_file(filename);
}

} // anonymous namespace

TEST_CASE ("horiz_structured_remap")
{
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;
  using gid_type = AbstractGrid::gid_type;

  ekat::Comm comm(MPI_COMM_WORLD);
  scorpio::init_subsystem(comm);

  auto& ts = ekat::TestSession::get();
  const std::string curv_map_file = ts.params.at("curvilinear-map-file");
  const std::string ll_map_file = ts.params.at("latlon-map-file");

  // The size of the tgt grid is in the map file
  const auto dst_grid_dims = read_full_var<int>(curv_map_file,"dst_grid_dims",2);
  const int nx = dst_grid_dims[0];
  const int ny = dst_grid_dims[1];

  const int nlevs = 3;

  // Create src grid, and read lat/lon from map file
  const int ncols_src = scorpio::get_dimlen(curv_map_file,"n_a");
  const int ncols_tgt = scorpio::get_dimlen(curv_map_file,"n_b");
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
    read_fields(curv_map_file,{io_lat,io_lon},io_gids,comm);
  }
  scorpio::clear_unused_decomps();

  // Read the map file content on all ranks, to compute a reference solution independently of
  // how the remapper works. This does not depend on the lat/lon values of the tgt grid points,
  // so it works for any weights.
  const int nnz = scorpio::get_dimlen(curv_map_file,"n_s");
  const auto map_row = read_full_var<int>(curv_map_file,"row",nnz);
  const auto map_col = read_full_var<int>(curv_map_file,"col",nnz);
  const auto map_S   = read_full_var<double>(curv_map_file,"S",nnz);
  const auto map_xc_a = read_full_var<double>(curv_map_file,"xc_a",ncols_src);
  const auto map_yc_a = read_full_var<double>(curv_map_file,"yc_a",ncols_src);
  const auto map_xc_b = read_full_var<double>(curv_map_file,"xc_b",ncols_tgt);
  const auto map_yc_b = read_full_var<double>(curv_map_file,"yc_b",ncols_tgt);
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

  // Check that gid=iy*nx+ix for the (y,x) dims of the grid
  auto check_ordering = [&](const AbstractGrid& tgt_grid, const int nx_, const int ny_) {
    // Nothing to check at the local level, other than gids being in the expected range
    const auto min_gid = tgt_grid.get_global_min_dof_gid();
    auto gids = tgt_grid.get_dofs_gids().get_view<const gid_type*,Host>();
    for (int i=0; i<tgt_grid.get_num_local_dofs(); ++i) {
      REQUIRE (gids(i)-min_gid>=0);
      REQUIRE (gids(i)-min_gid<nx_*ny_);
    }
  };

  SECTION ("curvilinear") {
    // The map file has dst_grid_rank=2, but lat/lon cannot be factored: lat and lon depend on both x and y.
    HorizontalRemapper remap(src_grid,curv_map_file);
    auto tgt_grid = remap.get_tgt_grid();

    REQUIRE (tgt_grid->get_num_global_dofs()==nx*ny);
    REQUIRE (tgt_grid->has_geometry_data("structured_dims"));
    const auto dims = get_structured_dims(*tgt_grid);
    REQUIRE (dims.ny==ny);
    REQUIRE (dims.nx==nx);
    check_ordering(*tgt_grid,nx,ny);

    // lat/lon remain regular fields over the columns
    for (const auto& n : {"lat","lon"}) {
      REQUIRE (tgt_grid->has_geometry_data(n));
      const auto& fid = tgt_grid->get_geometry_data(n).get_header().get_identifier();
      REQUIRE (fid.get_layout()==tgt_grid->get_2d_scalar_layout());
    }

    // Remap, and compare with the reference solution, computed from the map file content
    auto t2d = remap.register_field_from_src(f2d);
    auto t3d = remap.register_field_from_src(f3d);
    remap.registration_ends();
    remap.remap_fwd();
    t2d.sync_to_host();
    t3d.sync_to_host();

    std::vector<Real> ref(ncols_tgt,0);
    for (int s=0; s<nnz; ++s) {
      // Map files are 1-based
      ref[map_row[s]-1] += map_S[s]*smooth_fn(map_yc_a[map_col[s]-1],map_xc_a[map_col[s]-1]);
    }

    const auto min_gid = tgt_grid->get_global_min_dof_gid();
    auto gids = tgt_grid->get_dofs_gids().get_view<const gid_type*,Host>();
    auto lat = tgt_grid->get_geometry_data("lat").get_view<const Real*,Host>();
    auto lon = tgt_grid->get_geometry_data("lon").get_view<const Real*,Host>();
    auto t2d_h = t2d.get_view<const Real*,Host>();
    auto t3d_h = t3d.get_view<const Real**,Host>();
#ifdef SCREAM_DOUBLE_PRECISION
    constexpr Real tol = 1e-12;
#else
    constexpr Real tol = 2e-5;
#endif
    for (int i=0; i<tgt_grid->get_num_local_dofs(); ++i) {
      const int idx = gids(i) - min_gid;
      // The geo data must follow the same ordering as the map file
      REQUIRE (std::abs(lat(i)-map_yc_b[idx])<tol);
      REQUIRE (std::abs(lon(i)-map_xc_b[idx])<tol);

      REQUIRE (std::abs(t2d_h(i)-ref[idx])<tol);
      for (int k=0; k<nlevs; ++k) {
        REQUIRE (std::abs(t3d_h(i,k)-(k+1)*ref[idx])<(k+1)*tol);
      }
    }
  }

  SECTION ("latlon") {
    // lat (lon) only depends on y (x): lat/lon are 1d coordinate arrays
    const int nlon = 20;
    const int nlat = 10;
    REQUIRE (scorpio::get_dimlen(ll_map_file,"n_a")==ncols_src);
    REQUIRE (scorpio::get_dimlen(ll_map_file,"n_b")==nlon*nlat);

    HorizontalRemapper remap(src_grid,ll_map_file);
    auto tgt_grid = remap.get_tgt_grid();
    REQUIRE (tgt_grid->get_num_global_dofs()==nlon*nlat);
    REQUIRE (tgt_grid->has_geometry_data("structured_dims"));
    const auto dims = get_structured_dims(*tgt_grid);
    REQUIRE (dims.ny==nlat);
    REQUIRE (dims.nx==nlon);
    check_ordering(*tgt_grid,nlon,nlat);

    // lat/lon are 1d coordinate arrays, following the order of the points in the map file
    auto lat_f = tgt_grid->get_geometry_data("lat");
    auto lon_f = tgt_grid->get_geometry_data("lon");
    REQUIRE (lat_f.get_header().get_identifier().get_layout().size()==nlat);
    REQUIRE (lon_f.get_header().get_identifier().get_layout().size()==nlon);
    auto lat = lat_f.get_view<const Real*,Host>();
    auto lon = lon_f.get_view<const Real*,Host>();
    REQUIRE (std::abs(lat(0)+78.1216682)<1e-4);
    REQUIRE (std::abs(lon(0)-9)<1e-4);
    for (int i=1; i<nlat; ++i) {
      REQUIRE (lat(i)>lat(i-1));
    }
    for (int i=1; i<nlon; ++i) {
      REQUIRE (lon(i)>lon(i-1));
    }
  }

  SECTION ("latlon_ordering") {
    // A lat-lon grid stored with lat as the fastest varying index cannot be seen as an (nx,ny) lat-lon
    // grid. We still get a valid grid: a curvilinear one, with lat/lon being col-dependent fields.
    const int nlon = 4;
    const int nlat = 3;
    for (bool lat_major : {true,false}) {
      const std::string fname = std::string("latlon_map_") + (lat_major ? "lat_major" : "lon_major") + ".nc";
      create_latlon_map(fname,nlon,nlat,lat_major);

      auto grid = create_point_grid("src",nlon*nlat+1,nlevs,comm,1);
      HorizontalRemapper remap(grid,fname);
      auto tgt_grid = remap.get_tgt_grid();

      REQUIRE (tgt_grid->has_geometry_data("structured_dims"));
      const auto dims = get_structured_dims(*tgt_grid);
      // The dims have the same sizes for lat-lon and curvilinear grids
      REQUIRE (dims.ny==nlat);
      REQUIRE (dims.nx==nlon);
      check_ordering(*tgt_grid,nlon,nlat);
      // What differs is the layout of lat/lon: 1d coordinate arrays vs col-dependent fields
      if (lat_major) {
        REQUIRE (tgt_grid->get_geometry_data("lat").get_header().get_identifier().get_layout().size()==nlat);
        REQUIRE (tgt_grid->get_geometry_data("lon").get_header().get_identifier().get_layout().size()==nlon);
      } else {
        REQUIRE (tgt_grid->get_geometry_data("lat").get_header().get_identifier().get_layout()==tgt_grid->get_2d_scalar_layout());
        REQUIRE (tgt_grid->get_geometry_data("lon").get_header().get_identifier().get_layout()==tgt_grid->get_2d_scalar_layout());
      }
    }
  }

  scorpio::finalize_subsystem();
}

} // namespace scream
