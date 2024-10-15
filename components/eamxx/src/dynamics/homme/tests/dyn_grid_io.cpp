#include <catch2/catch.hpp>

#include "TimeLevel.hpp"

#include "dynamics/homme/homme_grids_manager.hpp"
#include "dynamics/homme/homme_dimensions.hpp"
#include "dynamics/homme/interface/scream_homme_interface.hpp"

#include "share/io/scorpio_input.hpp"
#include "share/io/scream_output_manager.hpp"

#include "share/field/field_utils.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/util/scream_setup_random_test.hpp"
#include "share/util/scream_time_stamp.hpp"

#include "ekat/util/ekat_units.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"

extern "C" {
// These are specific C/F calls for these tests (i.e., not part of scream_homme_interface.hpp)
void init_test_params_f90 ();
void cleanup_test_f90 ();
}

namespace {

/*
 * In this test we do a dyn->physGLL remap "on the fly" during the output phase.
 * Then, we load the generated file on a separate FieldManager, and compare the
 * values against the ones we would get manually running the d2p remapper.
 */

TEST_CASE("dyn_grid_io")
{
  using namespace scream;
  using namespace ShortFieldTagsNames;

  ekat::Comm comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.

  // Initialize the pio_subsystem for this test:
  scorpio::init_subsystem(comm);

  // Init homme context
  if (!is_parallel_inited_f90()) {
    auto comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);
    init_parallel_f90(comm_f);
  }
  init_test_params_f90 ();

  // Set parameters
  constexpr int ne = 2;
  set_homme_param("ne",ne);

  // Create the grids
  ekat::ParameterList params;
  params.set<std::string>("physics_grid_type","GLL");
  params.set<std::string>("vertical_coordinate_filename","NONE");
  auto gm = std::make_shared<HommeGridsManager>(comm,params);
  gm->build_grids();

  auto dyn_grid  = gm->get_grid("Dynamics");
  auto phys_grid = gm->get_grid("Physics GLL");

  // Local counters
  EKAT_REQUIRE_MSG(phys_grid->get_num_local_dofs()>0, "Internal test error! Fix dyn_grid_io, please.\n");
  EKAT_REQUIRE_MSG(get_num_local_elems_f90()>0, "Internal test error! Fix dyn_grid_io, please.\n");

  // Create physics and dynamics Fields
  auto dyn_scalar3d_mid = dyn_grid->get_3d_scalar_layout(true);
  auto dyn_vector3d_mid = dyn_grid->get_3d_vector_layout(true,2);
  auto dyn_scalar2d     = dyn_grid->get_2d_scalar_layout();

  auto phys_scalar3d_mid = phys_grid->get_3d_scalar_layout(true);
  auto phys_vector3d_mid = phys_grid->get_3d_vector_layout(true,2);
  auto phys_scalar2d     = phys_grid->get_2d_scalar_layout();

  auto nondim = ekat::units::Units::nondimensional();
  FieldIdentifier fid_dyn_1 ("field_1",dyn_scalar3d_mid,nondim,dyn_grid->name());
  FieldIdentifier fid_dyn_2 ("field_2",dyn_vector3d_mid,nondim,dyn_grid->name());
  FieldIdentifier fid_dyn_3 ("field_3",dyn_scalar2d    ,nondim,dyn_grid->name());

  FieldIdentifier fid_phys_1 ("field_1",phys_scalar3d_mid,nondim,phys_grid->name());
  FieldIdentifier fid_phys_2 ("field_2",phys_vector3d_mid,nondim,phys_grid->name());
  FieldIdentifier fid_phys_3 ("field_3",phys_scalar2d    ,nondim,phys_grid->name());

  // The starting FM
  auto fm_dyn = std::make_shared<FieldManager> (dyn_grid);

  // The FM we will manually remap onto
  auto fm_ctrl= std::make_shared<FieldManager> (phys_grid);

  // The FM we will read into, and compare against the previous
  auto fm_phys= std::make_shared<FieldManager> (phys_grid);

  fm_dyn->registration_begins();
  fm_phys->registration_begins();
  fm_ctrl->registration_begins();

  const int ps = HOMMEXX_PACK_SIZE;
  util::TimeStamp t0({2000,1,1},{0,0,0});

  fm_dyn->register_field(FieldRequest(fid_dyn_1,ps));
  fm_dyn->register_field(FieldRequest(fid_dyn_2,ps));
  fm_dyn->register_field(FieldRequest(fid_dyn_3));

  fm_phys->register_field(FieldRequest(fid_phys_1,ps));
  fm_phys->register_field(FieldRequest(fid_phys_2,ps));
  fm_phys->register_field(FieldRequest(fid_phys_3));

  fm_ctrl->register_field(FieldRequest(fid_phys_1,ps));
  fm_ctrl->register_field(FieldRequest(fid_phys_2,ps));
  fm_ctrl->register_field(FieldRequest(fid_phys_3));

  fm_dyn->registration_ends();
  fm_phys->registration_ends();
  fm_ctrl->registration_ends();
  fm_dyn->init_fields_time_stamp(t0);
  fm_phys->init_fields_time_stamp(t0);
  fm_ctrl->init_fields_time_stamp(t0);

  std::vector<std::string> fnames = {"field_1", "field_2", "field_3"};

  // Randomize dyn fields, then remap to ctrl fields
  std::uniform_real_distribution<Real> pdf(0.01,100.0);
  auto engine = setup_random_test(&comm);
  auto dyn2ctrl = gm->create_remapper(dyn_grid,phys_grid);
  dyn2ctrl->registration_begins();
  for (const auto& fn : fnames) {
    auto fd = fm_dyn->get_field(fn);
    auto fc = fm_ctrl->get_field(fn);
    dyn2ctrl->register_field(fd,fc);
    randomize(fd,engine,pdf);

    // Init phys field to something obviously wrong
    auto fp = fm_phys->get_field(fn);
    fp.deep_copy(-1.0);
  }
  dyn2ctrl->registration_ends();
  dyn2ctrl->remap(true);

  // Now try to write all fields to file from the dyn grid fm
  // Note: add MPI ranks to filename, to allow MPI tests to run in parallel
  ekat::ParameterList out_params;
  out_params.set<std::string>("Averaging Type","Instant");
  out_params.set<std::string>("filename_prefix","dyn_grid_io");
  out_params.sublist("Fields").sublist("Dynamics").set<std::vector<std::string>>("Field Names",fnames);
  out_params.sublist("Fields").sublist("Dynamics").set<std::string>("IO Grid Name","Physics GLL");

  out_params.sublist("output_control").set<int>("Frequency",1);
  out_params.sublist("output_control").set<std::string>("frequency_units","nsteps");
  out_params.set<std::string>("Floating Point Precision","real");

  OutputManager output;
  output.initialize(comm, out_params, t0, false);
  output.setup (fm_dyn, gm);
  output.run(t0);
  output.finalize();

  // Next, let's load all fields from file directly into the dyn grid fm
  std::string filename = "dyn_grid_io.INSTANT.nsteps_x1.np" + std::to_string(comm.size()) + "." + t0.to_string() + ".nc";
  filename.erase(std::remove(filename.begin(),filename.end(),':'),filename.end());

  ekat::ParameterList in_params;
  in_params.set<std::string>("Filename",filename);
  in_params.set<std::vector<std::string>>("Field Names",fnames);
  AtmosphereInput input (in_params,fm_phys);
  input.read_variables();
  input.finalize();

  // Compare against ctrl fields
  for (const auto& fn : fnames) {
    auto fp = fm_phys->get_field(fn);
    auto fc = fm_ctrl->get_field(fn);
    REQUIRE(views_are_equal(fp,fc));
  }

  // Cleanup everything
  scorpio::finalize_subsystem();
  Homme::Context::finalize_singleton();
  cleanup_test_f90();
}

} // anonymous namespace
