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

TEST_CASE("dyn_grid_io")
{
  using namespace scream;
  using namespace ShortFieldTagsNames;
  constexpr int np  = HOMMEXX_NP;
  constexpr int nlev = HOMMEXX_NUM_PHYSICAL_LEV;

  ekat::Comm comm(MPI_COMM_WORLD);  // MPI communicator group used for I/O set as ekat object.

  // Initialize the pio_subsystem for this test:
  MPI_Fint fcomm = MPI_Comm_c2f(comm.mpi_comm());
  scorpio::eam_init_pio_subsystem(fcomm);

  // Init homme context
  if (!is_parallel_inited_f90()) {
    auto comm_f = MPI_Comm_c2f(MPI_COMM_WORLD);
    init_parallel_f90(comm_f);
  }
  init_test_params_f90 ();

  // The homme context
  auto& c = Homme::Context::singleton();

  // The TimeLevel structure is needed by the PD remapper
  // Note: we don't remap states, so doesn't matter what values we pick
  c.create_if_not_there<Homme::TimeLevel>();
  auto& tl = c.get<Homme::TimeLevel>();
  tl.np1 = 0;
  tl.nm1 = 1;
  tl.n0  = 2;

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
  const int nelems = get_num_local_elems_f90();
  const int ncols  = phys_grid->get_num_local_dofs();
  EKAT_REQUIRE_MSG(ncols>0, "Internal test error! Fix dyn_grid_io, please.\n");
  EKAT_REQUIRE_MSG(nelems>0, "Internal test error! Fix dyn_grid_io, please.\n");

  // Create physics and dynamics Fields
  FieldLayout layout_dyn_1({EL,GP,GP,LEV},{nelems,np,np,nlev});
  FieldLayout layout_dyn_2({EL,CMP,GP,GP,LEV},{nelems,2,np,np,nlev});
  FieldLayout layout_dyn_3({EL,GP,GP},{nelems,np,np});

  FieldLayout layout_phys_1({COL,LEV},{ncols,nlev});
  FieldLayout layout_phys_2({COL,CMP,LEV},{ncols,2,nlev});
  FieldLayout layout_phys_3({COL},{ncols});

  auto nondim = ekat::units::Units::nondimensional();
  FieldIdentifier fid_dyn_1 ("field_1",layout_dyn_1,nondim,dyn_grid->name());
  FieldIdentifier fid_dyn_2 ("field_2",layout_dyn_2,nondim,dyn_grid->name());
  FieldIdentifier fid_dyn_3 ("field_3",layout_dyn_3,nondim,dyn_grid->name());

  FieldIdentifier fid_phys_1 ("field_1",layout_phys_1,nondim,phys_grid->name());
  FieldIdentifier fid_phys_2 ("field_2",layout_phys_2,nondim,phys_grid->name());
  FieldIdentifier fid_phys_3 ("field_3",layout_phys_3,nondim,phys_grid->name());

  auto fm_dyn = std::make_shared<FieldManager> (dyn_grid);
  auto fm_phys= std::make_shared<FieldManager> (phys_grid);
  auto fm_ctrl= std::make_shared<FieldManager> (phys_grid);

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

  // Randomize control fields, then remap to dyn fields
  std::uniform_real_distribution<Real> pdf(0.01,100.0);
  auto engine = setup_random_test(&comm);
  auto dyn2phys = gm->create_remapper(dyn_grid,phys_grid);
  auto dyn2ctrl = gm->create_remapper(dyn_grid,phys_grid);
  dyn2phys->registration_begins();
  dyn2ctrl->registration_begins();
  for (const auto& fn : fnames) {
    auto fd = fm_dyn->get_field(fn);
    auto fp = fm_phys->get_field(fn);
    auto fc = fm_ctrl->get_field(fn);
    dyn2ctrl->register_field(fd,fc);
    dyn2phys->register_field(fd,fp);
    randomize(fc,engine,pdf);

  }
  dyn2phys->registration_ends();
  dyn2ctrl->registration_ends();
  dyn2ctrl->remap(false); // Remap bwd to get data from ctrl to dyn

  // Now try to write all fields to file from the dyn grid fm
  // IMPORTANT! Make the file name dependent on comm size, so all tests
  //            can run in parallel without race conditions on the nc file.
  ekat::ParameterList io_params;
  io_params.set<int>("Max Snapshots Per File",1);
  io_params.set<std::string>("Averaging Type","Instant");
  io_params.set<std::vector<std::string>>("Grids",{"Dynamics"});
  io_params.set<std::string>("Casename","dyn_grid_io_np" + std::to_string(comm.size()));
  io_params.sublist("Fields").sublist("Dynamics").set<std::vector<std::string>>("Field Names",fnames);
  io_params.sublist("Fields").sublist("Dynamics").set<std::string>("IO Grid Name","Physics GLL");

  io_params.sublist("output_control").set<int>("Frequency",1);
  io_params.sublist("output_control").set<std::string>("frequency_units","nsteps");
  io_params.set<std::string>("Floating Point Precision","real");

  OutputManager output;
  // AtmosphereOutput output(comm,io_params,fm_dyn,gm);
  output.setup (comm, io_params, fm_dyn, gm, t0, t0, false);
  output.run(t0);
  output.finalize();

  // Clear the content of the dyn fields, to avoid seeing the same numbers
  // only b/c nothing was in fact read.
  for (const auto& fn : fnames) {
    auto f = fm_dyn->get_field(fn);
    f.deep_copy(-1.0);
  }

  // Next, let's load all fields from file directly into the dyn grid fm
  std::string filename = "dyn_grid_io_np" + std::to_string(comm.size())
                       + ".INSTANT.nsteps_x1." + t0.to_string() + ".nc";
  filename.erase(std::remove(filename.begin(),filename.end(),':'),filename.end());

  io_params.set<std::string>("Filename",filename);
  io_params.set<std::string>("Grid",dyn_grid->name());
  AtmosphereInput input (io_params,fm_dyn, gm);
  input.read_variables();
  input.finalize();

  // Remap dyn->phys, and compare against ctrl
  dyn2phys->remap(true);
  for (const auto& fn : fnames) {
    auto fp = fm_phys->get_field(fn);
    auto fc = fm_ctrl->get_field(fn);
    REQUIRE(views_are_equal(fp,fc));
  }

  // Cleanup everything
  scorpio::eam_pio_finalize();
  Homme::Context::finalize_singleton();
  cleanup_test_f90();
}

} // anonymous namespace
