#ifndef SCREAM_SCORPIO_OUTPUT_HPP
#define SCREAM_SCORPIO_OUTPUT_HPP

#include "share/io/scream_scorpio_interface.hpp"
#include "scream_config.h"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"

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
 *  3. a shared pointer to the field manager
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
 *  FIELDS: [field_1, field_2, ..., field_N]
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
 *  FIELDS is a list of fields that need to be added to this output stream
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

// Forward declarations
template<typename T>
class FieldManager;

class FieldLayout;

class AbstractGrid;

namespace util { class TimeStamp; }

class AtmosphereOutput 
{
public:

  using KT = KokkosTypes<DefaultDevice>;
  template<int N>
  using view_Nd_host = typename KT::template view_ND<Real,N>::HostMirror;
  using view_1d_host = view_Nd_host<1>;

  virtual ~AtmosphereOutput () = default;

  // Constructor
  // Note on the last two inputs:
  //  - is_restart_run: if true, this is a restarted run. This is only importand
  //    if this output instance is *not* for writing model restart files,
  //    and *only* for particular choices of "Averaging Type".
  //  - is_model_restart_output: if true, this Output is for model restart files.
  //    In this case, we have to also create an "rpointer.atm" file (which
  //    contains metadata, and is expected by the component coupled)
  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params, 
                   const std::shared_ptr<const FieldManager<Real>>& field_mgr,
                   const bool is_restarted_run = false,
                   const bool is_model_restart_output = false);

  // Main Functions
  void init();
  void run (const util::TimeStamp& time);
  void finalize();

protected:

  // Internal functions
  void register_dimensions(const std::string& name);
  void register_variables(const std::string& filename);
  void set_degrees_of_freedom(const std::string& filename);
  std::vector<int> get_var_dof_offsets (const FieldLayout& layout);
  void register_views();
  void new_file(const std::string& filename);
  void run_impl(const Real time, const std::string& time_str);
  std::string compute_filename_root (const std::string& casename) const;
  void combine (const Real& new_val, Real& curr_val) const;

  // --- Internal variables --- //
  ekat::Comm                                  m_comm;
  ekat::ParameterList                         m_params;
  std::shared_ptr<const FieldManager<Real>>   m_field_mgr;
  std::shared_ptr<const AbstractGrid>         m_grid;
  
  // The output filename root
  std::string       m_casename;

  // How to combine multiple snapshots in the output: Instant, Max, Min, Average
  std::string       m_avg_type;

  // Frequency of output control
  int m_out_max_steps;
  int m_out_frequency;
  std::string m_out_units;

  // Internal maps to the output fields, how the columns are distributed, the file dimensions and the global ids.
  std::vector<std::string>            m_fields;
  std::map<std::string,FieldLayout>   m_layouts;
  std::map<std::string,int>           m_dofs;
  std::map<std::string,int>           m_dims;

  // Local views of each field to be used for "averaging" output and writing to file.
  std::map<std::string,view_1d_host>    m_host_views_1d;

  // Whether this Output object writes a model restart file, or normal model output.
  bool m_is_model_restart_output;

  // If this is normal output, whether data to restart the history is needed/generated.
  bool m_has_restart_data;
  int m_checkpoint_freq;

  // Whether this run is the restart of a previous run (in which case, we might load an output checkpoint)
  bool m_is_restarted_run;

  // Whether the output file is open.
  // Note: this is redundant, since it's equal to m_num_snapshots_in_file==0, but it makes code more readable.
  bool m_is_output_file_open = false;

  // When this equals m_out_frequency, it's time to write the output.
  int m_nsteps_since_last_output = 0;

  // When this equals m_checkpoint_freq, it's time to write a history restart file.
  int m_nsteps_since_last_checkpoint = 0;

  // When this equals m_out_max_steps, it's time to close the out file, and open a new one.
  int m_num_snapshots_in_file = 0;
};

inline void AtmosphereOutput::combine (const Real& new_val, Real& curr_val) const
{
  if (m_avg_type=="INSTANT" || m_nsteps_since_last_output == 1) {
    curr_val = new_val;
  } else {
    // Update local view given the averaging type.
    if (m_avg_type == "AVERAGE") {
      curr_val = (curr_val*(m_nsteps_since_last_output-1) + new_val)/(m_nsteps_since_last_output);
    } else if (m_avg_type == "MAX") {
      curr_val = std::max(curr_val,new_val);
    } else if (m_avg_type == "Min") {
      curr_val = std::min(curr_val,new_val);
    }
  }
}

} //namespace scream

#endif // SCREAM_SCORPIO_OUTPUT_HPP
