#include <catch2/catch.hpp>
#include "share/atmosphere_process.hpp"
#include "share/scream_pack.hpp"
#include "share/grid/user_provided_grids_manager.hpp"
#include "share/grid/default_grid.hpp"
#include "control/atmosphere_driver.hpp"

namespace scream {

// === A dummy atm process, on Physics grid === //

template<typename DeviceType, int PackSize>
class DummyProcess : public scream::AtmosphereProcess {
public:
  using device_type = DeviceType;

  DummyProcess (const Comm& comm, const ParameterList& params)
   : m_input(FieldIdentifier("INVALID",{FieldTag::Invalid}))
   , m_output(FieldIdentifier("INVALID",{FieldTag::Invalid}))
   , m_comm(comm)
  {
    m_params = params;
    m_id = comm.rank();
  }

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  std::set<std::string> get_required_grids () const {
    // TODO: define what grid the coupling runs on. Check with MOAB folks.
    static std::set<std::string> s;
    s.insert("Physics");
    return s;
  }

  // Return some sort of name, linked to PType
  std::string name () const { return "Dummy process"; }

  // The communicator associated with this atm process
  const Comm& get_comm () const { return m_comm; }

  void set_grid (const std::shared_ptr<const GridsManager> grids_manager) {
    m_grid = grids_manager->get_grid("Physics");

    auto num_cols = m_grid->num_dofs();

    std::vector<FieldTag> tags = {FieldTag::Column,FieldTag::Component};
    std::vector<int> dims = {num_cols, m_params.get<int>("Number of vector components")};
    FieldLayout layout (tags,dims);

    std::string in_name = "field_";
    std::string out_name = "field_";
    auto size = m_comm.size();
    if (size==1) {
      in_name  += "0";
      out_name += "1";
    } else {
      in_name  += std::to_string(m_id);
      out_name += std::to_string( (m_id + size - 1) % size );
    }

    m_input_fids.emplace(in_name,layout,"Physics");
    m_output_fids.emplace(out_name,layout,"Physics");
  }

  void initialize (const util::TimeStamp& t0) {
    m_time_stamp = t0;
  }

  void run (const double dt) {
    auto in = m_input.get_view();
    auto out = m_output.get_view();
    auto id = m_id;
    Kokkos::parallel_for(Kokkos::RangePolicy<>(0,16),
      KOKKOS_LAMBDA(const int i) {
        out(i) = sin(in(i)+id);
    });
    Kokkos::fence();

    m_time_stamp += dt;
    m_output.get_header().get_tracking().update_time_stamp(m_time_stamp);
  }

  // Clean up
  void finalize ( ) {}

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const {
    using pack_type = pack::Pack<Real,PackSize>;
    field_repo.template register_field<pack_type>(*m_input_fids.begin());
    field_repo.template register_field<pack_type>(*m_output_fids.begin());
  }

  // Providing a list of required and computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_input_fids; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_output_fids; }

protected:

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real, device_type>& f) {
    m_input = f;
  }
  void set_computed_field_impl (const Field<      Real, device_type>& f) {
    m_output = f;
  }

  util::TimeStamp           m_time_stamp;

  std::set<FieldIdentifier> m_input_fids;
  std::set<FieldIdentifier> m_output_fids;

  Field<const Real,device_type> m_input;
  Field<Real,device_type>       m_output;

  std::shared_ptr<AbstractGrid> m_grid;

  ParameterList m_params;
  int     m_id;

  Comm    m_comm;
};

// === A dummy physics grids for this test === //

class DummyPhysicsGrid : public DefaultGrid<GridType::Physics>
{
public:
  DummyPhysicsGrid (const int num_cols)
   : DefaultGrid<GridType::Physics>("Physics")
  {
    m_num_dofs = num_cols;
  }
  ~DummyPhysicsGrid () = default;

protected:
};

TEST_CASE("ping-pong", "") {
  using namespace scream;
  using namespace scream::control;

  using device_type = AtmosphereDriver::device_type;

  constexpr int num_cols  = 32;

  // Create a parameter list for inputs
  ParameterList ad_params("Atmosphere Driver");
  auto& proc_params = ad_params.sublist("Atmosphere Processes");

  proc_params.set("Number of Entries",2);
  proc_params.set<std::string>("Schedule Type","Sequential");

  auto& p0 = proc_params.sublist("Process 0");
  p0.set<std::string>("Process Name", "Dummy");
  p0.set<int>("Number of vector components",2);

  auto& p1 = proc_params.sublist("Process 1");
  p1.set<std::string>("Process Name", "Dummy");
  p1.set<int>("Number of vector components",2);

  auto& gm_params = ad_params.sublist("Grids Manager");
  gm_params.set<std::string>("Type","User Provided");

  // Need to register products in the factory *before* we create any AtmosphereProcessGroup,
  // which rely on factory for process creation. The initialize method of the AD does that.
  // While we're at it, check that the case insensitive key of the factory works.
  auto& proc_factory = AtmosphereProcessFactory::instance();
  proc_factory.register_product("duMmy",&create_atmosphere_process<DummyProcess<device_type,SCREAM_PACK_SIZE>>);

  // Need to register grids managers before we create the driver
  auto& gm_factory = GridsManagerFactory::instance();
  gm_factory.register_product("User Provided",create_user_provided_grids_manager);

  // Set the dummy grid in the UserProvidedGridManager
  // Recall that this class stores *static* members, so whatever
  // we set here, will be reflected in the GM built by the factory.
  UserProvidedGridsManager upgm;
  upgm.set_grid(std::make_shared<DummyPhysicsGrid>(num_cols));
  upgm.set_reference_grid("Physics");

  // Create a comm
  Comm atm_comm (MPI_COMM_WORLD);

  // Create the driver
  AtmosphereDriver ad;

  // Init and run (do not finalize, or you'll clear the field repo!)
  util::TimeStamp init_time(2019,0,0);
  util::TimeStamp end_time (2019,1,0);
  const double dt = 1200.0;
  ad.initialize(atm_comm,ad_params,init_time);
  for (auto time=init_time+dt; time<end_time; time+=dt) {
    ad.run(dt);
  }

  // Every atm proc does out(:) = sin(in(:)+rank)
  Real answer = 0;
  for (auto time=init_time+dt; time<end_time; time+=dt) {
    for (int pid=0; pid<atm_comm.size(); ++pid) {
      answer = std::sin(answer+pid);
    }
  }

  // Get the field repo, and check the answer
  const auto& repo = ad.get_field_repo();

  std::vector<FieldTag> tags = {FieldTag::Column,FieldTag::Component};
  std::vector<int> dims = {num_cols, 2};
  FieldLayout layout (tags,dims);
  FieldIdentifier final_fid("field_" + std::to_string(atm_comm.size()-1),layout,"Physics");
  const auto& final_field = repo.get_field(final_fid);

  auto h_view = Kokkos::create_mirror_view(final_field.get_view());
  for (int i=0; i<h_view.extent_int(0); ++i) {
    REQUIRE (h_view(i) == answer);
  }

  // Finalize 
  ad.finalize();
}

} // empty namespace
