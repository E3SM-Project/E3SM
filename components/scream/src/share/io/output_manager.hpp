#ifndef SCREAM_OUTPUT_MANAGER_HPP
#define SCREAM_OUTPUT_MANAGER_HPP

#include "share/io/scorpio.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "share/field/field_repository.hpp"
#include "share/grid/grids_manager.hpp"

namespace scream
{

class OutputManager
{
public:
  using output_type = AtmosphereOutput;
  using device_type = DefaultDevice; 

  OutputManager () = default;
  virtual ~OutputManager () = default;

  void init();
  void run(util::TimeStamp& current_ts);
  void run(Real current_ts);
  void finalize();

  void set_comm(const ekat::Comm& comm) { atm_comm = comm; }
  void set_params(const ekat::ParameterList& params) { m_params=params; param_set=true;}
  void set_grids(const std::shared_ptr<const GridsManager>& gm) { m_grids_manager = gm; gm_set=true;}
  void set_repo(const std::shared_ptr<const FieldRepository<Real,device_type>>& repo) { m_device_field_repo = repo; repo_set=true;}

protected:
  std::vector<output_type>             m_output_streams;
  ekat::Comm                           atm_comm;
  ekat::Comm                           pio_comm; 
  ekat::ParameterList                  m_params;
  std::shared_ptr<const FieldRepository<Real,device_type>> m_device_field_repo;
  std::shared_ptr<const GridsManager>  m_grids_manager;

  bool                                 param_set = false;
  bool                                 gm_set    = false;
  bool                                 repo_set  = false;
  bool                                 no_output = true;  //TODO, when restarts are active then this paradigm won't work.  We will still need an output stream for restarts.  Possible solution, have restart information be a part of YAML control entries for output.  Then this will still work.

}; // class OutputManager

inline void OutputManager::init()
{
  using namespace scorpio;
  // Make sure params, repo and grids manager are already set
  EKAT_REQUIRE_MSG(param_set and gm_set and repo_set,"Error! Output manager requires a parameter list, grids manager and field repo to be set before initialization");
  // First parse the param list for output.
  // Starting with the stride.  NOTE: Only a stride of 1 is supported right now.  TODO: Allow for stride of any value.
  Int stride = m_params.get<Int>("PIO Stride");
  EKAT_REQUIRE_MSG(stride==1,"Error! Output only supports a PIO Stride of 1");  // TODO: Delete when no longer valid.
  Int comm_color = atm_comm.rank() % stride;
  auto pio_comm = atm_comm;  // TODO, EKAT should have a comm_split option (.split(comm_color));
  // PIO requires a subsystem to begin.  TODO, the component coupler actually inits the subsystem for the ATM,
  //                                     When the surface coupling is complete we can pass compid from the coupler
  //                                     to the manager and switch the third arguement below from true to false.
  int compid=0;  // For CIME based builds this will be the integer ID assigned to the atm by the component coupler.  For testing we simply set to 0
  MPI_Fint fcomm = MPI_Comm_c2f(pio_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
  eam_init_pio_subsystem(fcomm,compid,true);   // Gather the initial PIO subsystem data creater by component coupler
  
  // Construct and store an output stream instance for each output request.
  auto& list_of_files = m_params.get<std::vector<std::string>>("Output YAML Files");    // First grab the list of Output files from the control YAML
  for (auto& it : list_of_files)
  {
    ekat::ParameterList out_params(it);
    parse_yaml_file(it,out_params);
    output_type output_instance(pio_comm,out_params);
    m_output_streams.push_back(output_instance);
  }
  // Initialize each instance.
  for (auto& it : m_output_streams)
  {
    it.init(*m_device_field_repo, *m_grids_manager);
  }
  // Set no_output to false so that output will be created for the remainder of simulation
  no_output = false;
}
// Overload run to allow for passing a TimeStamp or just directly a Real
inline void OutputManager::run(util::TimeStamp& current_ts)
{
  if (no_output) { return; }
  for (auto& it : m_output_streams)
  {
    it.run(*m_device_field_repo, *m_grids_manager, current_ts.get_seconds());
  }
}
inline void OutputManager::run(Real current_ts)
{
  if (no_output) { return; }
  for (auto& it : m_output_streams)
  {
    it.run(*m_device_field_repo, *m_grids_manager, current_ts);
  }
}


inline void OutputManager::finalize()
{
  if (no_output) { return; }
  using namespace scorpio;
  for (auto& it : m_output_streams)
  {
    it.finalize();
  }
  // Finalize PIO overall
  eam_pio_finalize();
}
} // namespace scream
#endif // SCREAM_OUTPUT_MANAGER_HPP
