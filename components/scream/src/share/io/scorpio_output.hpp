#ifndef SCREAM_SCORPIO_OUTPUT_HPP
#define SCREAM_SCORPIO_OUTPUT_HPP

#include "share/io/scream_scorpio_interface.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/util//scream_time_stamp.hpp"

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
 *  Fields:
 *     GRID_NAME_1:               ARRAY OF STRINGS
 *     GRID_NAME_2:               ARRAY OF STRINGS
 *     ...
 *     GRID_NAME_N:               ARRAY OF STRINGS
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
 *  - Grids: a list of grids that we are using to output to this file
 *  - Averaging Type: a string that describes which type of output, current options are:
 *      instant - no averaging, output each snap as is.
 *      average - average of the field over some interval.
 *      min     - minimum value of the field over time interval.
 *      max     - maximum value of the field over time interval.
 *    Here, 'time interval' is described by ${Output Frequency} and ${Output Frequency Units}.
 *    E.g., with 'Output Frequency'=10 and 'Output Frequency Units'="Days", the time interval is 10 days.
 *  - GRID_NAME_[1,...,N]: a list of fields that need to be added to the output stream for the grid
 *                         $GRID_NAME_[1,...,N]. Each grid name must appear in the 'Grids' list
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

 *  Note: you can specify lists (such as the 'Fields->GRID_NAME_1' list above) with either of the two syntaxes
 *    GRID_NAME_1: [field_name1, field_name2, ... , field_name_N]
 *    GRID_NAME_2:
 *      - field_name_1
 *      - field_name_2
 *        ...
 *      - field_name_N
 *
 *  Note: each instance of this class can only handle ONE grid, so if multiple grids are specified,
 &        you will need one instance per grid.
 *  Usage of this class is to create an output file, write data to the file and close the file.
 *  This class keeps a temp array for all output fields to be used to perform averaging.
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 *  (2021-08-19) Luca Bertagna (SNL)
 *  (2021-10-14) Luca Bertagna (SNL)
 */

namespace scream
{

class AtmosphereOutput 
{
public:
  using fm_type       = FieldManager<Real>;
  using grid_type     = AbstractGrid;
  using gm_type       = GridsManager;
  using remapper_type = AbstractRemapper<Real>;

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
                   const std::shared_ptr<const fm_type>& field_mgr,
                   const std::shared_ptr<const gm_type>& grids_mgr);

  // Main Functions
  void restart (const std::string& filename);
  void init();
  void setup_output_file (const std::string& filename);
  void run (const std::string& filename, const bool write, const int nsteps_since_last_output);
  void finalize() {}

protected:

  // Internal functions
  void set_params (const ekat::ParameterList& params, const std::string& grid_name);
  void set_field_manager (const std::shared_ptr<const fm_type>& field_mgr,
                          const std::shared_ptr<const gm_type>& grids_mgr);
  void set_grid (const std::shared_ptr<const AbstractGrid>& grid);
  void build_remapper (const std::shared_ptr<const gm_type>& grids_mgr);

  void register_dimensions(const std::string& name);
  void register_variables(const std::string& filename);
  void set_degrees_of_freedom(const std::string& filename);
  std::vector<int> get_var_dof_offsets (const FieldLayout& layout);
  void register_views();
  void combine (const Real& new_val, Real& curr_val, const int nsteps_since_last_output) const;

  // --- Internal variables --- //
  ekat::Comm                                  m_comm;
  std::shared_ptr<const FieldManager<Real>>   m_field_mgr;
  std::shared_ptr<const AbstractGrid>         m_grid;
  std::shared_ptr<remapper_type>              m_remapper;

  // How to combine multiple snapshots in the output: Instant, Max, Min, Average
  std::string       m_avg_type;

  // Internal maps to the output fields, how the columns are distributed, the file dimensions and the global ids.
  std::vector<std::string>            m_fields_names;
  std::map<std::string,FieldLayout>   m_layouts;
  std::map<std::string,int>           m_dofs;
  std::map<std::string,int>           m_dims;

  // Local views of each field to be used for "averaging" output and writing to file.
  std::map<std::string,view_1d_host>    m_host_views_1d;
};

// ===================== IMPLEMENTATION ======================== //

// This helper function updates the current output val with a new one,
// according to the "averaging" type, and according to the number of
// model time steps since the last output step.
inline void AtmosphereOutput::combine (const Real& new_val, Real& curr_val, const int nsteps_since_last_output) const
{
  if (m_avg_type=="INSTANT" || nsteps_since_last_output == 1) {
    curr_val = new_val;
  } else {
    // Update local view given the averaging type.
    if (m_avg_type == "AVERAGE") {
      curr_val = (curr_val*(nsteps_since_last_output-1) + new_val)/(nsteps_since_last_output);
    } else if (m_avg_type == "MAX") {
      curr_val = std::max(curr_val,new_val);
    } else if (m_avg_type == "Min") {
      curr_val = std::min(curr_val,new_val);
    }
  }
}

} //namespace scream

#endif // SCREAM_SCORPIO_OUTPUT_HPP
