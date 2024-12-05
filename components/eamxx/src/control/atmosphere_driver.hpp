#ifndef SCREAM_ATMOSPHERE_DRIVER_HPP
#define SCREAM_ATMOSPHERE_DRIVER_HPP

#include "control/surface_coupling_utils.hpp"
#include "share/iop/intensive_observation_period.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/util/scream_time_stamp.hpp"
#include "share/scream_types.hpp"
#include "share/io/scream_output_manager.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
#include "share/atm_process/SCDataManager.hpp"

#include "ekat/logging/ekat_logger.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <memory>

namespace scream {

// Forward declarations
class AtmosphereProcess;
class AtmosphereProcessGroup;

namespace control {

/*
 *  The driver for the atmosphere component.
 *
 *  This class (AD) is responsible to keep track of the different atmosphere
 *  subcomponents (ASC) (parametrizations and dynamics). The AD is responsible for ensuring
 *  that all the ASCs required by the current test case are correctly initialized,
 *  that they are called in the correct order, and that they do not leak memory at the end.
 *
 *  The AD is also responsible for handling the FieldManager(s) of the Atmosphere.
 *  It must keep track of the different fields, and of the ASCs that require
 *  each field, ensuring all dependencies are met (in a DAG fashion).
 */

class AtmosphereDriver
{
public:
  using field_mgr_type = FieldManager;
  using field_mgr_ptr  = std::shared_ptr<field_mgr_type>;

  AtmosphereDriver () = default;
  AtmosphereDriver (const ekat::Comm& atm_comm,
                    const ekat::ParameterList& params);

  // Must call finalize, so that, if AD is destroyed as part of uncaught
  // exception stack unwinding, we will still perform some cleanup ops,
  // among which, for instance, closing any open output file.
  ~AtmosphereDriver ();

  // ---- Begin initialization methods ---- //

  // Set comm for the whole atmosphere
  void set_comm (const ekat::Comm& atm_comm);

  // Return comm for the whole atmosphere
  const ekat::Comm& get_comm () const { return m_atm_comm; }

  // Set AD params
  void set_params (const ekat::ParameterList& params);

  // Init time stamps
  // run_type: -1: deduce from run/case t0, 0: initial, 1: restart
  void init_time_stamps (const util::TimeStamp& run_t0, const util::TimeStamp& case_t0, int run_type = -1);

  // Set AD params
  void init_scorpio (const int atm_id = 0);

  // Setup IntensiveObservationPeriod
  void setup_iop ();

  // Create atm processes, without initializing them
  void create_atm_processes ();

  // Create needed grids, based on processes needs.
  void create_grids ();

  // Create fields as requested by all processes
  void create_fields ();

  // Adds cpl import/export information to SCDataManager.
  void setup_surface_coupling_data_manager(SurfaceCouplingTransferType transfer_type,
                                           const int num_cpl_fields, const int num_scream_fields,
                                           const int field_size, Real* data_ptr,
#ifdef HAVE_MOAB
                                           Real* data_ptr_moab,
#endif
                                           char* names_ptr, int* cpl_indices_ptr, int* vec_comps_ptr,
                                           Real* constant_multiple_ptr, bool* do_transfer_during_init_ptr);

  // Find surface coupling processes and have
  // them setup internal SurfaceCoupling data.
  void setup_surface_coupling_processes() const;

  // Zero out precipitation flux
  void reset_accumulated_fields();

  // Create and add mass and energy conservation checks
  // and pass to m_atm_process_group.
  void setup_column_conservation_checks ();

  // If TMS process exists, creates link to SHOC for applying
  // tms' surface drag coefficient.
  void setup_shoc_tms_links();

  // Add column data to all pre/postcondition property checks
  // for use in output.
  void add_additional_column_data_to_property_checks ();

  void set_provenance_data (std::string caseid = "",
                            std::string rest_caseid = "",
                            std::string hostname = "",
                            std::string username = "",
                            std::string versionid = "");

  // Load initial conditions for atm inputs
  void initialize_fields ();

  // Create output managers
  void create_output_managers ();

  // Initialie I/O structures for output
  void initialize_output_managers ();

  // Call 'initialize' on all atm procs
  void initialize_atm_procs ();

  // ---- End of initialization methods ---- //

  // A wrapper of all of the above (except setting SurfaceCoupling),
  // which is handy for scream standalone tests.
  //  - atm_comm: the MPI comm containing all ranks assigned to the atmosphere
  //  - params: parameter list with all atm options (organized in sublists)
  //  - run_t0 : the time stamp where the run starts
  //  - case_t0: the time stamp where the original simulation started (for restarts)
  void initialize (const ekat::Comm& atm_comm,
                   const ekat::ParameterList& params,
                   const util::TimeStamp& run_t0,
                   const util::TimeStamp& case_t0);

  // Shortcut for tests not doing restart
  void initialize (const ekat::Comm& atm_comm,
                   const ekat::ParameterList& params,
                   const util::TimeStamp& t0) {
    initialize(atm_comm,params,t0,t0);
  }

  // The run method is responsible for advancing the atmosphere component by one atm time step
  // Inside here you should find calls to the run method of each subcomponent, including parameterizations
  // and dynamics (HOMME).
  // Note: dt is assumed to be in seconds
  void run (const int dt);

  // Clean up the driver (finalizes and cleans up all internals)
  // NOTE: if already finalized, this is a no-op
  void finalize ();

  field_mgr_ptr get_field_mgr (const std::string& grid_name) const;

  // Get atmosphere time stamp
  const util::TimeStamp& get_atm_time_stamp () const { return m_current_ts; }

  const std::shared_ptr<GridsManager>& get_grids_manager () const { return m_grids_manager; }

  const std::shared_ptr<ATMBufferManager>& get_memory_buffer() const { return m_memory_buffer; }

  const std::shared_ptr<AtmosphereProcessGroup>& get_atm_processes () const { return m_atm_process_group; }

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif
  void initialize_constant_field(const FieldIdentifier& fid, const ekat::ParameterList& ic_pl);
protected:

  void report_res_dep_memory_footprint () const;

  void create_logger ();
  void set_initial_conditions ();
  void restart_model ();

  // Read fields from a file when the names of the fields in
  // EAMxx do not match exactly with the .nc file. Example is
  // for topography data files, where GLL and PG2 grid have
  // different naming conventions for phis.
  void read_fields_from_file (const std::vector<std::string>& field_names_nc,
                              const std::vector<std::string>& field_names_eamxx,
                              const std::shared_ptr<const AbstractGrid>& grid,
                              const std::string& file_name,
                              const util::TimeStamp& t0);
  // Read fields from a file when the names of the fields in
  // EAMxx match with the .nc file.
  void read_fields_from_file (const std::vector<std::string>& field_names,
                              const std::shared_ptr<const AbstractGrid>& grid,
                              const std::string& file_name,
                              const util::TimeStamp& t0);
  void register_groups ();

  std::map<std::string,field_mgr_ptr>       m_field_mgrs;

  std::shared_ptr<AtmosphereProcessGroup>   m_atm_process_group;

  std::shared_ptr<GridsManager>             m_grids_manager;

  ekat::ParameterList                       m_atm_params;

  std::shared_ptr<OutputManager>            m_restart_output_manager;
  std::list<OutputManager>                  m_output_managers;

  std::shared_ptr<ATMBufferManager>         m_memory_buffer;
  std::shared_ptr<SCDataManager>            m_surface_coupling_import_data_manager;
  std::shared_ptr<SCDataManager>            m_surface_coupling_export_data_manager;

  std::shared_ptr<IntensiveObservationPeriod> m_iop;

  // This is the time stamp at the beginning of the time step.
  util::TimeStamp                           m_current_ts;

  // These are the time stamps of the beginning of this run and case
  // respectively. For initial runs, they are the same, but for
  // restarted runs, the latter is "older" than the former
  util::TimeStamp                           m_run_t0;
  util::TimeStamp                           m_case_t0;
  RunType                                   m_run_type;

  // This is the comm containing all (and only) the processes assigned to the atmosphere
  ekat::Comm                                m_atm_comm;

  // The logger to be used throughout the ATM to log message
  std::shared_ptr<ekat::logger::LoggerBase> m_atm_logger;

  // Some status flags, used to make sure we call the init functions in the right order
  static constexpr int s_comm_set       =    1;
  static constexpr int s_params_set     =    2;
  static constexpr int s_scorpio_inited =    4;
  static constexpr int s_procs_created  =    8;
  static constexpr int s_grids_created  =   16;
  static constexpr int s_fields_created =   32;
  static constexpr int s_sc_set         =   64;
  static constexpr int s_output_inited  =  128;
  static constexpr int s_fields_inited  =  256;
  static constexpr int s_procs_inited   =  512;
  static constexpr int s_ts_inited      = 1024;
  static constexpr int s_output_created = 2048;

  // Lazy version to ensure s_atm_inited & flag is true for every flag,
  // even if someone adds new flags later on
  static constexpr int s_atm_inited     =  ~0;

  // Utility function to check the ad status
  void check_ad_status (const int flag, const bool must_be_set = true);

  // Whether GPTL must be finalized by the AD (in certain standalone runs)
  bool m_gptl_externally_handled;

  // Current ad initialization status
  int m_ad_status = 0;

  // Current simulation casename
  std::string m_casename;
};

}  // namespace control
}  // namespace scream

#endif // SCREAM_ATMOSPHERE_DRIVER_HPP
