#ifndef SCREAM_SCORPIO_OUTPUT_HPP
#define SCREAM_SCORPIO_OUTPUT_HPP

#include <iostream>
#include <fstream>

#include "scream_config.h"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/scorpio_input.hpp"

#include "share/field/field_repository.hpp"
#include "share/field/field_header.hpp"
#include "share/field/field.hpp"
#include "share/field/field_identifier.hpp"

#include "share/grid/grids_manager.hpp"

/*  The AtmosphereOutput class handles an output stream in SCREAM.
 *  Typical usage is to register an AtmosphereOutput object with the OutputManager, see output_manager.hpp
 *
 *  Similar to the typical AtmosphereProcess, output streams have a init, run and finalize routines.
 *  These routines are called during the respective steps of the AD call, i.e. ad.init(), ad.run() and ad.finalize()/
 *
 *  Each AtmosphereOutput instance handles one output stream.
 *
 *  At construction time ALL output instances require
 *  1. an EKAT comm group and
 *  2. a EKAT parameter list
 *  3. a shared pointer to the field repository
 *  4. a shared pointer to the grids manager
 * 
 *  The EKAT parameter list must contain all of the top level output control information and follows the format:
 *  ------
 *  FILENAME: STRING
 *  AVERAGING TYPE: STRING
 *  GRID: STRING                   (optional)
 *  FREQUENCY:
 *    OUT_N: INT
 *    OUT_OPTION: STRING
 *    OUT_MAX_STEPS: INT
 *  FIELDS:
 *    Number of Fields: INT
 *    field 1: STRING
 *    ...
 *    field N: STRING
 *  restart_hist_N: INT            (optional)
 *  restart_hist_OPTION: STRING    (optional)
 *  RESTART FILE: BOOL             (optional)
 *  -----
 *  where,
 *  FILENAME is a string of the filename suffix.  TODO: change this to a casename associated with the whole run.
 *  AVERAGING TYPE is a string that describes which type of output, current options are:
 *    instant - no averaging, output each snap as is.
 *    average - average of the field over some interval described in frequency section.
 *    min     - minimum value of the field over some time interval
 *    max     - maximum value of the field over some time interval
 *  GRID is a string describing which grid to write on, currently the only option is the default Physics.
 *  FREQUENCY is a subsection that controls the frequency parameters of output.
 *    OUT_N is an integer of the frequency of output given the units from OUT_OPTION
 *    OUT_OPTION is a string for the units of output frequency, examples would be "Steps", "Months", "Years", "Hours", etc.
 *    OUT_MAX_STEPS is an integer of the maximum number of steps that can exist on a single file (controls files getting too big).
 *  FIELDS is a subsection that lists all the fields in this output stream.
 *    Number of Fields is an integer that specifies the number of fields in this output stream.
 *    field 1,...,field N is a list of each field in the output stream by name in field manager.
 *  restart_hist_N is an optional integer parameter that specifies the frequenct of restart history writes.
 *  restart_hist_OPTION is an optional string parameter for the units of restart history output.
 *  RESTART FILE is an optional boolean parameter that specifies if this output stream is a restart output, which is treated differently.
 *
 *  Usage of this class is to create an output file, write data to the file and close the file.
 *  This class keeps a running copy of data for all output fields locally to be used for the different averaging flags.
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 */

namespace scream
{

class AtmosphereOutput 
{
public:
  using dofs_list_type = AbstractGrid::dofs_list_type;
  using view_type_host = typename KokkosTypes<HostDevice>::view_1d<Real>;
  using input_type     = AtmosphereInput;

  virtual ~AtmosphereOutput () = default;

  // Constructor
  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params, 
                   const std::shared_ptr<const FieldRepository<Real>>& repo,
                   const std::shared_ptr<const GridsManager>& gm)
  {
    m_comm       = comm;
    m_params     = params;
    m_field_repo = repo;
    m_gm         = gm;
    m_read_restart_hist = false;
  }
  // Constructor
  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params, 
                   const std::shared_ptr<const FieldRepository<Real>>& repo,
                   const std::shared_ptr<const GridsManager>& gm,
                   const bool read_restart_hist)
  {
    m_comm       = comm;
    m_params     = params;
    m_field_repo = repo;
    m_gm         = gm;
    m_read_restart_hist = read_restart_hist;
  }

  // Main Functions
  void init();
  // Allow run to be called with a timestamp input (typical runs) and a floating point value (used for unit tests)
  void run(const Real time);
  void run(const util::TimeStamp& time);
  void finalize();

  // Helper Functions
  void check_status();
  std::map<std::string,Int> get_status() const { return m_status; }

protected:
  // Internal functions
  void register_dimensions(const std::string& name);
  void register_variables(const std::string& filename);
  void set_degrees_of_freedom(const std::string& filename);
  void register_views();
  void new_file(const std::string& filename);
  void run_impl(const Real time, const std::string& time_str);  // Actual run routine called by outward facing "run"
  void set_restart_hist_read( const bool bval ) { m_read_restart_hist = bval; }
  // Internal variables
  ekat::ParameterList                          m_params;
  ekat::Comm                                   m_comm;
  std::shared_ptr<const FieldRepository<Real>> m_field_repo;
  std::shared_ptr<const GridsManager>          m_gm;
  
  // Main output control data
  std::string m_casename;
  std::string m_avg_type;
  std::string m_grid_name;
  // Frequency of output control
  Int m_out_max_steps;
  Int m_out_frequency;
  std::string m_out_units;
  // How individual columns are distributed across MPI Ranks
  Int m_total_dofs;
  Int m_local_dofs;
  // Restart history control
  Int m_restart_hist_n;
  std::string m_restart_hist_option;
  // Internal maps to the output fields, how the columns are distributed, the file dimensions and the global ids.
  std::vector<std::string>               m_fields;
  std::map<std::string,Int>              m_dofs;
  std::map<std::string,Int>              m_dims;
  typename dofs_list_type::HostMirror    m_gids;
  // Local views of each field to be used for "averaging" output and writing to file.
  std::map<std::string,view_type_host>   m_view_local;

  // Manage when files are open and closed, and what type of file I am writing.
  bool m_is_init = false;
  bool m_is_restart_hist = false;  //TODO:  If instead we rely on the timestamp to determine how many steps are represented in the averaging value, or maybe filename, we won't need this.
  bool m_is_restart = false;
  bool m_read_restart_hist = false;

  // Helper map to monitor an output stream's status.  Most used fields are Snaps and Avg Count which are
  // used to monitor if a new file is needed and if output should be written according to the frequency settings.
  std::map<std::string,Int> m_status = {
                                  {"Init",          0},  // Records the number of files this output stream has managed
                                  {"Run",           0},  // Total number of times "Run" has been called
                                  {"Finalize",      0},  // Total number of times "Finalize" has been called (should always be 1)
                                  {"Snaps",         0},  // Total number of timesnaps saved to the currently open file.
                                  {"Avg Count",     0},  // Total number of timesnaps that have gone by since the last time output was written.
                                       }; 

}; // Class AtmosphereOutput

//// ====================== IMPLEMENTATION ===================== //
inline void AtmosphereOutput::check_status()
{
  printf("IO Status for Rank %5d, File - %.40s: (Init: %2d), (Run: %5d), (Finalize: %2d), (Avg. Count: %2d), (Snaps: %2d)\n",
                   m_comm.rank(),m_casename.c_str(),m_status["Init"],m_status["Run"],m_status["Finalize"],m_status["Avg Count"],m_status["Snaps"]);
} // check_status
} //namespace scream
#endif // SCREAM_SCORPIO_OUTPUT_HPP
