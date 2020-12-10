#ifndef SCREAM_OUTPUT_MANAGER_HPP
#define SCREAM_OUTPUT_MANAGER_HPP

#include "share/io/scorpio_output.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "share/field/field_repository.hpp"
#include "share/grid/grids_manager.hpp"

namespace scream
{

/*
 * The OutputManager Class handles all individual output streams.
 * Rather than have the SCREAM-AD or unit tests call the scorpio_output
 * class objects directly, the OutputManager stores these in a vector,
 * m_output_streams.
 *
 * Similar to a typical atmospheric process object, the OutputManager
 * has an init, run and finalize routine which is called during the AD
 * during those respective steps.
 *
 * PROPER USAGE:
 * The output manager requires a communication group, a parameter list of control
 * variables, a grid manager and a field repository.
 * Each of these four things are set using the setter function 'set_X' where
 *   X = comm, for the EKAT comm group.
 *     = params, for the parameter list.  In typical SCREAM runs this parameter
 *       list is a sublist of the scream control yaml file.
 *     = grids, for the grids mananger
 *     = repo, for the field repository.
 * The setup of the output manager in the SCREAM-AD is one of the last steps to
 * ensure that all four of the above objects have already been constructed.
 * see /control/atmospheric_driver.cpp for an example.
 *
 * For UNIT TESTS:
 * The output manager does require a comm group, list of parameters, grids manager
 * and field repository to work.  If output is desired in a unit test then these
 * must be established.  There are examples in /src/share/io/tests of how to 
 * establish a simple grids manager and field repository.  As well as how to
 * locally create a parameter list.
 *
 * Adding output streams mid-simulation:
 * It is possible to add an output stream after init has been called by calling
 * the internal function 'new_output' which takes an EKAT parameter list as input.
 * See comments in new_output below for more details.
 *
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 */
class OutputManager
{
public:
  using output_type = AtmosphereOutput;

  // Destructor
  virtual ~OutputManager () = default;
  // Constructors
  OutputManager() = default;

  OutputManager( const ekat::Comm& comm, const ekat::ParameterList& params,
                 const std::shared_ptr<const FieldRepository<Real>>& repo,
                 const std::shared_ptr<const GridsManager>& gm,
                 const bool runtype_restart)
  {
    atm_comm            = comm;
    m_params            = params;
    param_set           = true;
    m_device_field_repo = repo;
    repo_set            = true;
    m_grids_manager     = gm;
    gm_set              = true;
    m_runtype_restart   = runtype_restart;
  }

  OutputManager( const ekat::Comm& comm, const ekat::ParameterList& params,
                 const std::shared_ptr<const FieldRepository<Real>>& repo,
                 const std::shared_ptr<const GridsManager>& gm)
  {
    atm_comm            = comm;
    m_params            = params;
    param_set           = true;
    m_device_field_repo = repo;
    repo_set            = true;
    m_grids_manager     = gm;
    gm_set              = true;
  }

  void init();
  void run(util::TimeStamp& current_ts);
  void run(Real current_ts);
  void finalize();

/* Short function to add a new output stream to the output manager.  By making this an independent
 * function it is possible to add an output stream after initialization has been called.           */
  void new_output(const ekat::ParameterList& params);
  void new_output(const ekat::ParameterList& params,const bool runtype_restart);

  void set_comm(const ekat::Comm& comm) { atm_comm = comm; }
  void set_params(const ekat::ParameterList& params) { m_params=params; param_set=true;}
  void set_grids(const std::shared_ptr<const GridsManager>& gm) { m_grids_manager = gm; gm_set=true;}
  void set_repo(const std::shared_ptr<const FieldRepository<Real>>& repo) { m_device_field_repo = repo; repo_set=true;}
  void set_runtype_restart(const bool bval) { m_runtype_restart = bval; }
  void make_restart_param_list(ekat::ParameterList& params);

protected:
  std::vector<std::shared_ptr<output_type>>    m_output_streams;
  ekat::Comm                                   atm_comm;
  ekat::Comm                                   pio_comm; 
  ekat::ParameterList                          m_params;
  std::shared_ptr<const FieldRepository<Real>> m_device_field_repo;
  std::shared_ptr<const GridsManager>          m_grids_manager;

  bool                                 param_set = false;
  bool                                 gm_set    = false;
  bool                                 repo_set  = false;
  bool                        m_runtype_restart  = false;

}; // class OutputManager
/*===============================================================================================*/
/* Short function to add a new output stream to the output manager.  By making this an independent
 * function it is possible to add an output stream after initialization has been called.
 * Note: new_output requires a parameter list with all of the output control information for this
 * stream.  See scorpio_output.hpp for more information on what the parameter list needs.         */
inline void OutputManager::new_output(const ekat::ParameterList& params)
{
  auto output_instance = std::make_shared<output_type>(pio_comm,params,m_device_field_repo,m_grids_manager,m_runtype_restart);
  output_instance->init();
  m_output_streams.push_back(output_instance);
}
/* --------------------------------------------------------------------- */
inline void OutputManager::new_output(const ekat::ParameterList& params, const bool runtype_restart)
{
  auto output_instance = std::make_shared<output_type>(pio_comm,params,m_device_field_repo,m_grids_manager,runtype_restart);
  output_instance->init();
  m_output_streams.push_back(output_instance);
}
/*===============================================================================================*/
/*
 * init():
 *   1. Parses the parameter list that contains all output manager control information.
 *   2. Creates a restart file output stream, if defined in the parameter list.
 *   3. Cycles through all of the output streams listed in the parameter list.
 *      See scorpio_output.hpp for more information on output streams.
 */
inline void OutputManager::init()
{
  using namespace scorpio;
  if (!param_set and !gm_set and !repo_set) { return; }
  // Make sure params, repo and grids manager are already set
  EKAT_REQUIRE_MSG(param_set,"Error! Output manager requires a parameter list to be set before initialization");
  EKAT_REQUIRE_MSG(gm_set,   "Error! Output manager requires a grids manager to be set before initialization");
  EKAT_REQUIRE_MSG(repo_set, "Error! Output manager requires a field repo to be set before initialization");
  // First parse the param list for output.
  // Starting with the stride.  NOTE: Only a stride of 1 is supported right now.  TODO: Allow for stride of any value.
  Int stride = m_params.get<Int>("PIO Stride",1);
  EKAT_REQUIRE_MSG(stride==1,"Error! Output only supports a PIO Stride of 1");  // TODO: Delete when no longer valid.
  // TODO: Implement stride>1 for PIO.  The following commented lines spec out how this would look.
  // Int comm_color = atm_comm.rank() % stride;
  pio_comm = atm_comm;  // TODO, EKAT should have a comm_split option (.split(comm_color));
  // PIO requires a subsystem to begin.  TODO, the component coupler actually inits the subsystem for the ATM,
  // int compid=0;  // For CIME based builds this will be the integer ID assigned to the atm by the component coupler.  For testing we simply set to 0
//  MPI_Fint fcomm = MPI_Comm_c2f(pio_comm.mpi_comm());  // MPI communicator group used for I/O.  In our simple test we use MPI_COMM_WORLD, however a subset could be used.
//  eam_init_pio_subsystem(fcomm);   // Gather the initial PIO subsystem data creater by component coupler
  
  // Check for restart output:
  // Restarts are a special kind of output that have unique characteristics.  The following constructs the EKAT
  // parameter list to control the restart output stream.  
  // The user can trigger restart output at runtime by adding the node "Restart Control" to the SCREAM control YAML
  // file.  A typical restart control entry may look like:
  // ----
  // Restart Contol
  //   FREQUENCY
  //     OUT_N: INT
  //     OUT_OPTION: STRING
  // ----
  // where OUT_OPTION is the units of the output frequency (ex. Steps, Months, etc.) and OUT_N>0 is the actual frequency,
  // (ex. 1 would be every OUT_OPTION, 2 would be every other OUT_OPTION).
  Int restart_hist_N = 0;
  std::string restart_hist_OPTION = "None";
  if (m_params.isSublist("Restart Control"))
  {
    auto& out_params = m_params.sublist("Restart Control");
    make_restart_param_list(out_params);
    new_output(out_params,false);
    // Gather restart frequency info for restart history files
    auto& freq_params = out_params.sublist("FREQUENCY");
    restart_hist_N = freq_params.get<Int>("OUT_N");
    restart_hist_OPTION = freq_params.get<std::string>("OUT_OPTION");
  }
  // Construct and store an output stream instance for each output request.
  // Typical output is controlled by parameter list which is stored in individual YAML files.
  // A list of those files is passed to the output manager via the parameter list using the
  // key "Output YAML Files".  See srd/share/io/tests for examples.
  auto& list_of_files = m_params.get<std::vector<std::string>>("Output YAML Files");    // First grab the list of Output files from the control YAML
  for (auto& it : list_of_files)
  {
    ekat::ParameterList out_params(it);
    parse_yaml_file(it,out_params);
    out_params.set<Int>("restart_hist_N",restart_hist_N);
    out_params.set<std::string>("restart_hist_OPTION",restart_hist_OPTION);
    new_output(out_params);
  }
}
/*===============================================================================================*/
/*
 * run(time):
 *   1. Run is overloaded to allow for the passage of a timestamp (default) or time
 *   as a float.
 *   2. This routine simply calls 'run' for each output stream stored in the output
 *      manager.
 */
// Overload run to allow for passing a TimeStamp or just directly a Real
inline void OutputManager::run(util::TimeStamp& current_ts)
{
  for (auto& it : m_output_streams)
  {
    it->run(current_ts);
  }
}
/*-----------------------------------------------------------------------------------------------*/
inline void OutputManager::run(Real current_ts)
{
  for (auto& it : m_output_streams)
  {
    it->run(current_ts);
  }
}
/*===============================================================================================*/
/*
 * finalize():
 *   1. This routine calls 'finalize' for each output stream stored in the output
 *      manager.
 */
inline void OutputManager::finalize()
{
  for (auto& it : m_output_streams)
  {
    it->finalize();
    it = nullptr;
  }
}
/*===============================================================================================*/
inline void OutputManager::make_restart_param_list(ekat::ParameterList& params)
{
/*  
 * Given the unique nature of restart files, this routine sets up the specific requirements.
 */
  params.set<std::string>("FILENAME", "scorpio_restart_test");  //TODO change this to some sort of case_name, TODO: control so suffix is .r.nc
  // Parse the parameters that controls this output instance.
  params.set<std::string>("AVERAGING TYPE","Instant");  // Restart themselves are instant snapshots.  Restart histories are a different file and are controlled directly in scorpio_output.hpp
  params.set<std::string>("GRID","Physics"); // Restarts are only supported on the Physics grid for now.
  params.set<bool>("RESTART FILE",true);
  auto& freq_params = params.sublist("FREQUENCY");
  freq_params.get<Int>("OUT_MAX_STEPS",1);
  // Restart files don't also need a restart history file
  params.set<Int>("restart_hist_N",0);
  params.set<std::string>("restart_hist_OPTION","NONE");
  
  // set fields for restart
  // If a developer wants a field to be stored in the restart file than they must add the "restart" group tag
  // to that field when registering the field with the field repository.
  auto& restart_fields = m_device_field_repo->get_field_groups_names().at("restart");
  auto& field_params = params.sublist("FIELDS");
  field_params.set<Int>("Number of Fields",restart_fields.size());
  Int it_cnt = 0;
  for (auto it : restart_fields)
  {
    ++it_cnt;
    std::string f_it = "field ";
    f_it+=std::to_string(it_cnt);
    field_params.set<std::string>(f_it,it);
  }
}
/*===============================================================================================*/
} // namespace scream
#endif // SCREAM_OUTPUT_MANAGER_HPP
