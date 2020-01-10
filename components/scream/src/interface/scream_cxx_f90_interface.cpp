#include <catch2/catch.hpp>
#include "interface/ScreamContext.hpp"
#include "control/atmosphere_driver.hpp"
#include "share/mpi/scream_comm.hpp"
#include "share/parameter_list.hpp"
#include "share/scream_session.hpp"

#include "share/grid/user_provided_grids_manager.hpp"
#include "control/tests/dummy_grid.hpp"
#include "control/tests/dummy_atm_proc.hpp"

#include <numeric>

extern "C"
{

/*===============================================================================================*/
void scream_init (const MPI_Fint& f_comm, const int& start_ymd, const int& start_tod) {

  // First of all, initialize the scream session
  scream::initialize_scream_session();

  // Get the context
  auto& c = scream::ScreamContext::singleton();

  using namespace scream;
  using device_type = scream::control::AtmosphereDriver::device_type;

  constexpr int num_cols   = 2;
  constexpr int vec_length = 4;

  // Create the C MPI_Comm from the Fortran one
  MPI_Comm mpi_comm_c = MPI_Comm_f2c(f_comm);
  auto& comm = c.create<scream::Comm>(mpi_comm_c);

  // Create a parameter list for inputs
  scream::ParameterList ad_params("Atmosphere Driver");
  auto& proc_params = ad_params.sublist("Atmosphere Processes");

  proc_params.set("Number of Entries",2);
  proc_params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = proc_params.sublist("Process 0");
  p0.set<std::string>("Process Name", "Physics_fwd");
  p0.set<int>("Number of vector components",vec_length);

  auto& p1 = proc_params.sublist("Process 1");
  p1.set<std::string>("Process Name", "Physics_bwd");
  p1.set<int>("Number of vector components",vec_length);

  auto& gm_params = ad_params.sublist("Grids Manager");
  gm_params.set<std::string>("Reference Grid","Physics_fwd");
  gm_params.set<std::string>("Type","User Provided");

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("physics_fwd",&create_atmosphere_process<DummyProcess<device_type,2,true>>);
  proc_factory.register_product("physics_bwd",&create_atmosphere_process<DummyProcess<device_type,4,false>>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("User Provided",create_user_provided_grids_manager);

  // Set the dummy grid in the UserProvidedGridManager
  // Recall that this class stores *static* members, so whatever
  // we set here, will be reflected in the GM built by the factory.
  auto& upgm = c.create<UserProvidedGridsManager>(); // upgm;
  auto dummy_grid_fwd = std::make_shared<DummyPhysicsGrid>(num_cols,true);
  auto dummy_grid_bwd = std::make_shared<DummyPhysicsGrid>(num_cols,false);

  upgm.set_grid(dummy_grid_fwd);
  upgm.set_grid(dummy_grid_bwd);
  upgm.set_reference_grid("Physics_fwd");
  using remapper_type = DummyPhysicsGridRemapper<Real,device_type>;
  upgm.set_remapper(std::make_shared<remapper_type>(dummy_grid_fwd,dummy_grid_bwd));
  upgm.set_remapper(std::make_shared<remapper_type>(dummy_grid_bwd,dummy_grid_fwd));

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  // Create the bare ad, then init it
  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  auto& ad = c.create<scream::control::AtmosphereDriver>();

  // Init and run (to finalize, wait till checks are completed,
  // or you'll clear the field repo!)
  util::TimeStamp init_time(2019,0,0);
  const Real dt = 3500.0; // This should trigger an adjustment of the last time step
  ad.initialize(atm_comm,ad_params,init_time);

  (void) start_ymd;
  (void) start_tod;
  (void) f_comm;
}
/*===============================================================================================*/
void scream_run (const double& dt) {
  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  using namespace scream;

  constexpr int num_cols   = 2;
  constexpr int vec_length = 4;
  // Get the context
  auto& c = scream::ScreamContext::singleton();

  // Get the AD, and run it
  auto& ad = c.getNonConst<scream::control::AtmosphereDriver>();

  const auto& init_time = ad.get_atm_time_stamp();
  util::TimeStamp end_time (2019,1,0);
  const Real dt_tmp = 3500.0; // This should trigger an adjustment of the last time step

  // Fill the field with initial guess
  std::vector<FieldTag> tags = {FieldTag::Column,FieldTag::Component};
  std::vector<int> dims = {num_cols, vec_length};
  FieldLayout layout (tags,dims);
  FieldIdentifier fid("field_0",layout,units::m,"Physics_fwd");
  const auto& repo = ad.get_field_repo();
  const auto& field = repo.get_field(fid);
  auto d_view = field.get_view();
  auto h_view = Kokkos::create_mirror_view(d_view);
  const int size = h_view.size();
  decltype(h_view) answer("",size);

  Kokkos::deep_copy(d_view,0.0);
  Kokkos::deep_copy(h_view,0.0);
#ifdef SCREAM_DEBUG
  Kokkos::deep_copy(ad.get_bkp_field_repo().get_field(fid).get_view(),d_view);
#endif

  // Run
  int iter = 0;
  for (auto time=init_time; time<end_time; time+=dt_tmp, ++iter) {
    const auto& ts = ad.get_atm_time_stamp();
    std::cout << " -------------------------------------------------------\n";
    std::cout << "   Start of atmosphere time step\n";
    std::cout << "    - current time: " << ts.to_string() << "\n";
    Real dt_adj = dt_tmp;
    bool adjusted = false;
    if (end_time<(ts+dt_tmp)) {
      dt_adj = (end_time-ts).get_seconds();
      adjusted = true;
    }
    std::cout << "    - time step: " << dt_adj << (adjusted ? "s (adjusted to hit final time)\n" : "s\n");
    std::cout << "    - end time: " << (ts+dt_adj).to_string() << "\n";
    ad.run(dt_adj);

    // Prepare the answer
    // Every atm proc does out(:) = sin(in(:)). There are 2 atm process, with
    // input and output swapped, that do this update, so at every
    // time step we should get f(i) = sin(sin(f(i))).
    int rem = iter % 4;
    for (int i=0; i<size; ++i) {
      switch (rem) {
        // Do update twice, since both atm procs do the same update
        case 0: answer[i] += 2.0; answer[i] += 2.0; break;
        case 1: answer[i] *= 2.0; answer[i] *= 2.0; break;
        case 2: answer[i] -= 2.0; answer[i] -= 2.0; break;
        case 3: answer[i] /= 2.0; answer[i] /= 2.0; break;
      }
    }
  }

  // Check the answer
  Kokkos::deep_copy(h_view,d_view);
  for (int i=0; i<h_view.extent_int(0); ++i) {
    if (h_view(i) != answer[i]) {
      std::cout << "Error in ping-pong test for test: " <<  i << "\n";
    };
  }

  (void) dt;
}

/*===============================================================================================*/
void scream_finalize (/* args ? */) {
  // TODO: uncomment once you have valid inputs. I fear AD may crash with no inputs.
  // Get the context
  auto& c = scream::ScreamContext::singleton();

  // Get the AD, and finalize it
  auto& ad = c.getNonConst<scream::control::AtmosphereDriver>();
  auto& upgm = c.getNonConst<scream::UserProvidedGridsManager>();
  ad.finalize();
  upgm.clean_up();
}

} // extern "C"
