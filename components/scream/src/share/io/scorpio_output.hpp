#ifndef SCREAM_SCORPIO_OUTPUT_HPP
#define SCREAM_SCORPIO_OUTPUT_HPP

#include "share/io/scream_scorpio_interface.hpp"
#include "scream_config.h"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"

/*  The AtmosphereOutput class handles an output stream in SCREAM.
 *  Typical usage is to register an AtmosphereOutput object with the OutputManager (see scream_output_manager.hpp
 *
 *  Similar to other SCREAM classes, output streams have a init, run and finalize routines.
 *  These routines are called during the homonymous steps of the AD.
 *
 *  Each AtmosphereOutput instance handles one output stream. That means that each instance can
 *  only write to a single file, and can only write output from a single grid.
 *
 *  At construction time ALL output instances require at least an EKAT comm group, an
 *  EKAT parameter list, and a pointer to the field manager.
 * 
 *  The EKAT parameter list contains the following options to control output behavior
 *  ------
 *  Casename:                     STRING
 *  Averaging Type:               STRING
 *  Max Snapshots Per File:       INT                   (default: 1)
 *  Fields:                       ARRAY OF STRINGS
 *  Output:
      Frequency:                  INT
 *    Frequency Units:            STRING                (default: Steps)
 *  Checkpointing:
 *    Frequency:                  INT                   (default: 0)
 *    Frequency Units:            STRING                (default: ${Output->Frequency Units})
 *  Restart:
 *    Casename:                   STRING                (default: ${Casename})
 *    Perform Restart:            BOOL                  (default: true)
 *  -----
 *  The meaning of these parameters is the following:
 *  - Casename: the output filename root.
 *  - Averaging Type: a string that describes which type of output, current options are:
 *      instant - no averaging, output each snap as is.
 *      average - average of the field over some interval.
 *      min     - minimum value of the field over time interval.
 *      max     - maximum value of the field over time interval.
 *    Here, 'time interval' is described by ${Output Frequency} and ${Output Frequency Units}.
 *    E.g., with 'Output Frequency'=10 and 'Output Frequency Units'="Days", the time interval is 10 days.
 *  - Fields: a list of fields that need to be added to this output stream
 *  - Max Snapshots Per File: the maximum number of snapshots saved per file. After this many
 *  - Output: parameters for output control
 *    - Frequency: the frequency of output writes (in the units specified by ${Output Frequency Units})
 *    - Frequency Units: the units of output frequency (Steps, Months, Years, Hours, Days,...)
 *      snapshots have been written on a single nc file, the class will close the file, and open a new one
 *  - Checkpointing: parameters for checkpointing control
 *    - Frequency: the frequenct of checkpoints writes. This option is used/matters only if
 *      if Averaging Type is *not* Instant. A value of 0 is interpreted as 'no checkpointing'.
 *    - Frequency Units: the units of restart history output.
 *  - Restart: parameters for history restart
 *    - Casename: the history restart filename root.
 *    - Perform Restart: if this is a restarted run, and Averaging Type is not Instant, this flag
 *      determines whether we want to restart the output history or start from scrach. That is,
 *      you can set this to false to force a fresh new history, even in a restarted run.

 *  Note: you can specify lists (such as the 'Fields' list above) with either of the two syntaxes
 *    Fields: [field_name1, field_name2, ... , field_name_N]
 *    Fields:
 *      - field_name_1
 *      - field_name_2
 *        ...
 *      - field_name_N
 *
 *  Usage of this class is to create an output file, write data to the file and close the file.
 *  This class keeps a temp array for all output fields to be used to perform averaging.
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 *  (2021-08-19) Luca Bertagna (SNL)
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
  int m_out_frequency;
  std::string m_out_frequency_units;

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
  
  // Frequency of checkpoint writes
  int m_checkpoint_freq;
  std::string m_checkpoint_freq_units;

  // Whether this run is the restart of a previous run (in which case, we might load an output checkpoint)
  bool m_is_restarted_run;

  // Whether the output file is open.
  // Note: this is redundant, since it's equal to m_num_snapshots_in_file==0, but it makes code more readable.
  bool m_is_output_file_open = false;

  // When this equals m_out_frequency, it's time to write the output.
  int m_nsteps_since_last_output = 0;

  // When this equals m_checkpoint_freq, it's time to write a history restart file.
  int m_nsteps_since_last_checkpoint = 0;

  // To keep nc files small, we limit the number of snapshots in each nc file
  // When the number of snapshots in a file reaches m_out_max_steps, it's time
  // to close the out file, and open a new one.
  int m_max_snapshots_per_file;
  int m_num_snapshots_in_file = 0;
};

// ===================== IMPLEMENTATION ======================== //

// This helper function updates the current output val with a new one,
// according to the "averaging" type, and according to the number of
// model time steps since the last output step.
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
