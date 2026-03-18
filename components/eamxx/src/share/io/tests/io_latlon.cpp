#include <catch2/catch.hpp>
#include <memory>

#include "share/io/scorpio_input.hpp"
#include "share/io/eamxx_output_manager.hpp"
#include "share/grid/point_grid.hpp"

#include <ekat_test_utils.hpp>

namespace scream {

// Create a point grid with ncols equal to the src or tgt grid of the map file
// Also read in lat/lon from map file
inline std::shared_ptr<AbstractGrid>
create_grid(const ekat::Comm& comm, const std::string& map_file, const std::string& name, bool src)
{
  // Note: map files are 1-based, so create pt grid with 1-based gids
  auto nondim = ekat::units::Units::nondimensional();
  ekat::units::Units deg(nondim,"degrees");
  std::string suffix = src ? "_a" : "_b";
  int ncols = scorpio::get_dimlen(map_file,"n"+suffix);
  auto grid = create_point_grid(name,ncols,1,comm,1);

  auto lat = grid->create_geometry_data("lat",grid->get_2d_scalar_layout(),deg);
  auto lon = grid->create_geometry_data("lon",grid->get_2d_scalar_layout(),deg);
  auto io_grid = grid->clone(grid->name(),true);
  io_grid->reset_field_tag_name(FieldTag::Column,"n"+suffix);
  AtmosphereInput reader(map_file,io_grid,std::vector<Field>{lat.alias("yc"+suffix),lon.alias("xc"+suffix)});
  reader.read_variables();
  return grid;
}

// Create a fm with a copy of lat/lon from the grid in it
std::shared_ptr<FieldManager> create_fm(const std::shared_ptr<const AbstractGrid>& grid)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;
  using FID = FieldIdentifier;

  // Create a fm
  auto fm = std::make_shared<FieldManager>(grid,RepoState::Closed);

  Field f1(FID("f_lat", grid->get_2d_scalar_layout(),m,grid->name()),true);
  Field f2(FID("f_lon", grid->get_2d_scalar_layout(),m,grid->name()),true);
  f1.deep_copy(grid->get_geometry_data("lat"));
  f2.deep_copy(grid->get_geometry_data("lon"));

  fm->add_field(f1);
  fm->add_field(f2);
  return fm;
}

void print (const std::string& msg, const ekat::Comm& comm) {
  if (comm.am_i_root()) {
    printf("%s",msg.c_str());
  }
}

TEST_CASE("io_latlon")
{
  using strvec_t = std::vector<std::string>;

  // Setup the global structure
  ekat::Comm comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.
  scorpio::init_subsystem(comm);

  print (" -> Test Setup ...\n",comm);

  auto& ts = ekat::TestSession::get();
  const std::string& map_file = ts.params.at("map-file");

  // Construct a timestamp
  util::TimeStamp t0 ({2000,1,1},{0,0,0});

  auto src_grid = create_grid(comm,map_file,"ne4pg2",true);
  auto tgt_grid = create_grid(comm,map_file,"latlon",false);

  // Later, we'll open the map file again to read lat/lon in the remapper,
  // but we may have a different grid distribution, so we need to ensure that
  // the var decomp will be set from scratch.
  scorpio::clear_unused_decomps();

  auto fm = create_fm(src_grid);
  fm->init_fields_time_stamp(t0);
  print (" -> Test Setup ... done\n",comm);

  // Setup remapped output streams and run them
  print (" -> Write remapped output ... \n",comm);
  OutputManager om;

  ekat::ParameterList params;
  params.set<std::string>("filename_prefix","io_latlon");
  params.set<std::string>("averaging_type","instant");
  params.set<std::string>("file_max_storage_type","one_year");
  auto& oc = params.sublist("output_control");
  oc.set<int>("frequency",1);
  oc.set<std::string>("frequency_units","nsteps");
  params.set<strvec_t>("field_names",{"f_lat","f_lon"});
  params.set<std::string>("horiz_remap_file",map_file);

  om.initialize(comm,params,t0,false);
  om.setup(fm,{src_grid->name()});
  om.finalize();
  print (" -> Write remapped putput ... done!\n",comm);
  scorpio::finalize_subsystem();

}

} // namespace scream
