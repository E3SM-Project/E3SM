#include <catch2/catch.hpp>

#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "share/grid/mesh_free_grids_manager.hpp"

#include "share/field/field_identifier.hpp"
#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_utils.hpp"

#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"

#include "ekat/ekat_pack.hpp"
#include "ekat/util/ekat_units.hpp"
#include "ekat/io/ekat_yaml.hpp"
#include "ekat/ekat_parameter_list.hpp"

namespace {

using namespace scream;
using namespace ekat::units;
const int packsize = SCREAM_PACK_SIZE;
using Pack         = ekat::Pack<Real,packsize>;

std::shared_ptr<FieldManager>
get_test_fm(const std::shared_ptr<const AbstractGrid>& grid,
            const util::TimeStamp& t0, const bool do_randomize);

std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const int num_my_elems, const int np, const int num_levs);

ekat::ParameterList get_in_params(const ekat::Comm& comm,
                                  const util::TimeStamp& t0);

TEST_CASE("se_grid_io")
{
  using strvec_t = std::vector<std::string>;
  ekat::Comm io_comm(MPI_COMM_WORLD);

  int num_my_elems = 2;
  int np = 4;
  int num_levs = 2 + SCREAM_PACK_SIZE;
  int dt = 10;

  // Initialize the pio_subsystem for this test:
  scorpio::init_subsystem(io_comm);

  // First set up a field manager and grids manager to interact with the output functions
  auto gm = get_test_gm(io_comm,num_my_elems,np,num_levs);
  auto grid = gm->get_grid("SE Grid");

  // Construct a timestamp
  util::TimeStamp t0 ({2000,1,1},{0,0,0});

  auto fm0 = get_test_fm(grid,t0,true);
  ekat::ParameterList params;
  params.set<std::string>("filename_prefix","io_se_grid");
  params.set<std::string>("Averaging Type","Instant");
  params.set<int>("Max Snapshots Per File",1);
  params.set<strvec_t>("Field Names",{"field_1","field_2","field_3","field_packed"});
  params.set<std::string>("Floating Point Precision","real");
  auto& ctl_pl = params.sublist("output_control");
  ctl_pl.set("Frequency",1);
  ctl_pl.set<std::string>("frequency_units","nsteps");

  OutputManager om;
  om.initialize(io_comm,params,t0,false);
  om.setup(fm0,gm);
  om.init_timestep(t0,dt);
  om.run(t0+dt);
  om.finalize();

  // Get a fresh new field manager, and set fields to NaN
  auto fm1 = get_test_fm(grid,t0,false);
  const auto fnames = {"field_1", "field_2", "field_3", "field_packed"};
  for (const auto& fname : fnames) {
    auto f = fm1->get_field(fname);
    f.deep_copy(ekat::ScalarTraits<Real>::invalid());
  }

  // Check fields were written correctly
  auto in_params = get_in_params(io_comm,t0);
  AtmosphereInput ins_input(in_params,fm1);
  ins_input.read_variables();

  for (const auto& fname : fnames) {
    auto f0 = fm0->get_field(fname);
    auto f1 = fm1->get_field(fname);
    REQUIRE (views_are_equal(f0,f1));
  }
  ins_input.finalize();

  // All Done
  scorpio::finalize_subsystem();
}

/*===================================================================================================*/
std::shared_ptr<FieldManager>
get_test_fm(const std::shared_ptr<const AbstractGrid>& grid,
            const util::TimeStamp& t0, const bool do_randomize)
{
  using namespace ShortFieldTagsNames;
  using FL = FieldLayout;
  using FR = FieldRequest;

  const auto& comm = grid->get_comm();

  // Create a fm
  auto fm = std::make_shared<FieldManager>(grid);

  // Create some fields for this fm
  const int nlevs = grid->get_num_vertical_levels();
  const std::string& gn = grid->name();

  FieldIdentifier fid1("field_1",grid->get_2d_scalar_layout(),kg,gn);
  FieldIdentifier fid2("field_2",FL{{LEV},{nlevs}},kg,gn);
  FieldIdentifier fid3("field_3",grid->get_3d_scalar_layout(true),kg/m,gn);
  FieldIdentifier fid4("field_packed",grid->get_3d_scalar_layout(true),kg/m,gn);

  // Register fields with fm
  fm->registration_begins();
  fm->register_field(FR{fid1});
  fm->register_field(FR{fid2});
  fm->register_field(FR{fid3});
  fm->register_field(FR{fid4,Pack::n}); // Register field as packed
  fm->registration_ends();

  // Randomize fields
  const auto fnames = {"field_1", "field_2", "field_3", "field_packed"};
  if (do_randomize) {
    auto engine = setup_random_test (&comm);
    using RPDF = std::uniform_real_distribution<Real>;
    RPDF pdf(0.01,0.99);

    for (const auto& fname : fnames) {
      auto f = fm->get_field(fname);
      randomize(f,engine,pdf);
      f.get_header().get_tracking().update_time_stamp(t0);
    }
  } else {
    for (const auto& fname : fnames) {
      auto f = fm->get_field(fname);
      f.deep_copy(-1);
      f.get_header().get_tracking().update_time_stamp(t0);
    }
  }

  // field_2 is not partitioned, so let's sync it across ranks
  auto f2 = fm->get_field("field_2");
  auto v2 = f2.get_view<Real*>();
  comm.all_reduce(v2.data(),nlevs,MPI_MAX);

  return fm;
}
/*==========================================================================================================*/
std::shared_ptr<GridsManager>
get_test_gm(const ekat::Comm& io_comm, const int num_my_elems, const int np, const int num_levs)
{
  auto gm = create_mesh_free_grids_manager(io_comm,num_my_elems,np,num_levs,0);
  gm->build_grids();

  return gm;
}
/*==================================================================================================*/
ekat::ParameterList get_in_params(const ekat::Comm& comm,
                                  const util::TimeStamp& t0)
{
  using vos_type = std::vector<std::string>;
  ekat::ParameterList in_params("Input Parameters");

  std::string filename = "io_se_grid.INSTANT.nsteps_x1.np"
                       + std::to_string(comm.size())
                       + "." + t0.to_string() + ".nc";

  in_params.set<std::string>("Filename",filename);
  in_params.set<vos_type>("Field Names",{"field_1", "field_2", "field_3", "field_packed"});
  in_params.set<std::string>("Floating Point Precision","real");
  return in_params;
}

} // anonymous namespace
