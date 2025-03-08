#ifndef SCREAM_OUTPUT_MANAGER_HPP
#define SCREAM_OUTPUT_MANAGER_HPP

#include "share/io/scorpio_output.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/io/eamxx_io_file_specs.hpp"
#include "share/io/eamxx_io_control.hpp"

#include "share/field/field_manager.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include "ekat/logging/ekat_logger.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

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
 * variables, a grid manager and a field manager.
 * Each of these four things are set using the setter function 'set_X' where
 *   X = comm, for the EKAT comm group.
 *     = params, for the parameter list.  In typical SCREAM runs this parameter
 *       list is a sublist of the scream control yaml file.
 *     = grids, for the grids mananger
 *     = fm, for the field manager.
 * The setup of the output manager in the SCREAM-AD is one of the last steps to
 * ensure that all four of the above objects have already been constructed.
 * see /control/atmospheric_driver.cpp for an example.
 *
 * For UNIT TESTS:
 * The output manager does require a comm group, list of parameters, grids manager
 * and field manager to work.  If output is desired in a unit test then these
 * must be established.  There are examples in /src/share/io/tests of how to
 * establish a simple grids manager and field manager.  As well as how to
 * locally create a parameter list.
 *
 * Adding output streams mid-simulation:
 * TODO - This doesn't actually exist
 * It is possible to add an output stream after init has been called by calling
 * the internal function 'add_output_stream' which takes an EKAT parameter list as input.
 * See comments in add_output_stream below for more details.
 *
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 *  (2021-08-19) Luca Bertagna (SNL)
 */
class OutputManager
{
public:
  using fm_type = FieldManager;
  using gm_type = GridsManager;
  using globals_map_t = std::map<std::string,ekat::any>;

  // Constructor(s) & Destructor
  OutputManager() = default;
  virtual ~OutputManager ();

  // Initialize manager by storing class members that only depend on runtime parameters.
  // Inputs:
  //  - params: the parameter list with file/fields info, as well as method of output options
  //  - run_t0: the timestamp of the start of the current simulation
  //  - case_t0: the timestamp of the start of the overall simulation (precedes run_r0 for
  //             a restarted simulation. Restart logic is triggered *only* if case_t0<run_t0.
  //  - is_model_restart_output: whether this output stream is to write a model restart file
  void initialize (const ekat::Comm& io_comm, const ekat::ParameterList& params,
                   const util::TimeStamp& run_t0, const util::TimeStamp& case_t0,
                   const bool is_model_restart_output, const RunType run_type);

  // This overloads are to make certain unit tests easier
  void initialize (const ekat::Comm& io_comm, const ekat::ParameterList& params,
                   const util::TimeStamp& run_t0, const util::TimeStamp& case_t0,
                   const bool is_model_restart_output)
  {
    auto run_type = case_t0<run_t0 ? RunType::Restart : RunType::Initial;
    initialize(io_comm,params,run_t0,case_t0,is_model_restart_output,run_type);
  }
  void initialize (const ekat::Comm& io_comm, const ekat::ParameterList& params,
                   const util::TimeStamp& run_t0,
                   const bool is_model_restart_output)
  {
    initialize(io_comm, params, run_t0, run_t0, is_model_restart_output, RunType::Initial);
  }

  // Setup manager by creating the internal output streams using grids/field data
  // Inputs:
  //  - field_mgr/field_mgrs: field manager(s) storing fields to be outputed
  //  - grids_mgr: grid manager to create remapping from field managers grids onto the IO grid.
  //               This is needed, e.g., when outputing SEGrid fields without duplicating dofs.
  void setup (const std::shared_ptr<fm_type>& field_mgr,
              const std::shared_ptr<const gm_type>& grids_mgr);

  void setup (const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs,
              const std::shared_ptr<const gm_type>& grids_mgr);

  void set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& atm_logger) {
      m_atm_logger = atm_logger;
  }
  void add_global (const std::string& name, const ekat::any& global);

  void init_timestep (const util::TimeStamp& start_of_step, const Real dt);
  void run (const util::TimeStamp& current_ts);
  void finalize();

  long long res_dep_memory_footprint () const;

  bool is_restart () const { return m_is_model_restart_output; }

  // For debug and testing purposes
  const IOControl&   output_control    () const { return m_output_control;    }
  const IOFileSpecs& output_file_specs () const { return m_output_file_specs; }
protected:

  std::string compute_filename (const IOFileSpecs& file_specs,
                                const util::TimeStamp& timestamp) const;

  void set_file_header(const IOFileSpecs& file_specs);

  // Set internal class variables and processes the field_mgrs for restart fields
  // to add to the parameter list for a model restart managers.
  void setup_internals (const std::map<std::string,std::shared_ptr<fm_type>>& field_mgrs);

  void setup_file (      IOFileSpecs& filespecs,
                   const IOControl& control);

  // If a file can be closed (next snap won't fit) or needs flushing, do so
  void close_or_flush_if_needed (      IOFileSpecs& file_specs,
                                 const IOControl&   control) const;

  // Manage logging of info to atm.log
  void push_to_logger();

  using output_type     = AtmosphereOutput;
  using output_ptr_type = std::shared_ptr<output_type>;

  std::vector<output_ptr_type>   m_output_streams;
  std::vector<output_ptr_type>   m_geo_data_streams;

  globals_map_t                  m_globals;

  ekat::Comm                     m_io_comm;
  ekat::ParameterList            m_params;

  // The output filename root
  std::string       m_filename_prefix;

  std::vector<double> m_time_bnds;

  // How to combine multiple snapshots in the output: Instant, Max, Min, Average
  OutputAvgType     m_avg_type;

  // Whether this OutputManager handles a model restart file, or normal model output.
  bool m_is_model_restart_output;

  // Frequency of output and checkpointing
  // See eamxx_io_utils.hpp for details.
  IOControl m_output_control;
  IOControl m_checkpoint_control;

  // The specs (name, capacity, size) of output and checkpoint file.
  // See eamxx_io_utils.hpp for details.
  IOFileSpecs m_output_file_specs;
  IOFileSpecs m_checkpoint_file_specs;

  // Whether this run is the restart of a previous run, in which case
  // we might have to load an output checkpoint file (depending on avg type)
  RunType m_run_type;

  // Whether a restarted run can resume filling previous run output file (if not full)
  bool m_resume_output_file = false;

  // The initial time stamp of the simulation and run. For initial runs, they coincide,
  // but for restarted runs, run_t0>case_t0, with the former being the time at which the
  // restart happens, and the latter being the start time of the *original* run.
  util::TimeStamp   m_case_t0;
  util::TimeStamp   m_run_t0;

  // The logger to be used throughout the ATM to log message
  std::shared_ptr<ekat::logger::LoggerBase> m_atm_logger;

  // If true, we save grid data in output file
  bool m_save_grid_data;
};

} // namespace scream

#endif // SCREAM_OUTPUT_MANAGER_HPP
