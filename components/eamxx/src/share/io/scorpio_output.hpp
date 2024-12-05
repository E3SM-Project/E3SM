#ifndef SCREAM_SCORPIO_OUTPUT_HPP
#define SCREAM_SCORPIO_OUTPUT_HPP

#include "share/io/scream_scorpio_interface.hpp"
#include "share/io/scream_io_utils.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/util//scream_time_stamp.hpp"
#include "share/atm_process/atmosphere_diagnostic.hpp"

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
 *  filename_prefix:              STRING
 *  Averaging Type:               STRING
 *  Max Snapshots Per File:       INT                   (default: 1)
 *  Fields:
 *     GRID_NAME_1:
 *        Field Names:            ARRAY OF STRINGS
 *        IO Grid Name:           STRING                (optional)
 *     GRID_NAME_2:
 *        Field Names:            ARRAY OF STRINGS
 *        IO Grid Name:           STRING                (optional)
 *     ...
 *     GRID_NAME_N:
 *        Field Names:            ARRAY OF STRINGS
 *        IO Grid Name:           STRING                (optional)
 *  output_control:
 *    Frequency:                  INT
 *    frequency_units:            STRING                (default: nsteps)
 *  Restart:
 *    filename_prefix:            STRING                (default: ${filename_prefix})
 *    Perform Restart:            BOOL                  (default: true)
 *  -----
 *  The meaning of these parameters is the following:
 *  - filename_prefix: the output filename root.
 *  - Averaging Type: a string that describes which type of output, current options are:
 *      instant - no averaging, output each snap as is.
 *      average - average of the field over some interval.
 *      min     - minimum value of the field over time interval.
 *      max     - maximum value of the field over time interval.
 *    Here, 'time interval' is described by ${Output Frequency} and ${Output frequency_units}.
 *    E.g., with 'Output Frequency'=10 and 'Output frequency_units'="Days", the time interval is 10 days.
 *  - Fields: parameters specifying fields to output
 *     - GRID_NAME: parameters specifyign fields to output from grid $GRID_NAME
 *        - Field Names: names of fields defined on grid $grid_name that need to be outputed
 *        - IO Grid Name: if provided, remap fields to this grid before output (useful to remap
 *                        SEGrid fields to PointGrid fields on the fly, to save on output size)
 *  - Max Snapshots Per File: the maximum number of snapshots saved per file. After this many
 *    snapshots, the current files is closed and a new file created.
 *  - Output: parameters for output control
 *    - Frequency: the frequency of output writes (in the units specified by ${Output frequency_units})
 *    - frequency_units: the units of output frequency (nsteps, nmonths, nyears, nhours, ndays,...)
 *      snapshots have been written on a single nc file, the class will close the file, and open a new one
 *  - Checkpointing: parameters for checkpointing control
 *    - Frequency: the frequenct of checkpoints writes. This option is used/matters only if
 *      if Averaging Type is *not* Instant. A value of 0 is interpreted as 'no checkpointing'.
 *    - frequency_units: the units of restart history output.
 *  - Restart: parameters for history restart
 *    - filename_prefix: the history restart filename root.
 *    - Perform Restart: if this is a restarted run, and Averaging Type is not Instant, this flag
 *      determines whether we want to restart the output history or start from scrach. That is,
 *      you can set this to false to force a fresh new history, even in a restarted run.

 *  Notes:
 *   - you can specify lists with either of the two syntaxes:
 *
 *    LIST_NAME: [item_1,item_2,...,item_N]
 *    LIST_NAME:
 *      - item_1
 *      - item_2
 *        ...
 *      - item_N
 *
 *   - in case of single-grid tests, you can specify fields names by adding 'Field Names' directly
 *     in the top-level parameter list. In that case, you can also add 'IO Grid Name' in the top-level
 *     parameter list.
 *   - each instance of this class can only handle ONE grid, so if multiple grids are specified,
 *     you will need one instance per grid.
 *   - usage of this class is to create an output file, write data to the file and close the file.
 *   - this class keeps a temp array for all output fields to be used to perform averaging.
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 *  (2021-08-19) Luca Bertagna (SNL)
 *  (2021-10-14) Luca Bertagna (SNL)
 *  (2021-11-10) Luca Bertagna (SNL)
 */

namespace scream
{

class AtmosphereOutput
{
public:
  using fm_type       = FieldManager;
  using grid_type     = AbstractGrid;
  using gm_type       = GridsManager;
  using remapper_type = AbstractRemapper;
  using atm_diag_type = AtmosphereDiagnostic;

  using KT = KokkosTypes<DefaultDevice>;
  template<int N>
  using view_Nd_dev  = typename KT::template view_ND<Real,N>;
  template<int N>
  using view_Nd_host = typename KT::template view_ND<Real,N>::HostMirror;

  using view_1d_dev  = view_Nd_dev<1>;
  using view_1d_host = view_Nd_host<1>;

  virtual ~AtmosphereOutput () = default;

  // Constructor
  AtmosphereOutput(const ekat::Comm& comm, const ekat::ParameterList& params,
                   const std::shared_ptr<const fm_type>& field_mgr,
                   const std::shared_ptr<const gm_type>& grids_mgr);

  // Short version for outputing a list of fields (no remapping supported)
  AtmosphereOutput(const ekat::Comm& comm,
                   const std::vector<Field>& fields,
                   const std::shared_ptr<const grid_type>& grid);

  // Main Functions
  void restart (const std::string& filename);
  void init();
  void reset_dev_views();
  void update_avg_cnt_view(const Field&, view_1d_dev& dev_view);
  void setup_output_file (const std::string& filename, const std::string& fp_precision, const scorpio::FileMode mode);

  void init_timestep (const util::TimeStamp& start_of_step);
  void run (const std::string& filename,
            const bool output_step, const bool checkpoint_step,
            const int nsteps_since_last_output,
            const bool allow_invalid_fields = false);

  long long res_dep_memory_footprint () const;

  std::shared_ptr<const AbstractGrid> get_io_grid () const {
    return m_io_grid;
  }

  // Option to add a logger
  void set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& atm_logger) {
      m_atm_logger = atm_logger;
  }

protected:
  // Internal functions
  void set_grid (const std::shared_ptr<const AbstractGrid>& grid);
  void set_field_manager (const std::shared_ptr<const fm_type>& field_mgr, const std::string& mode);
  void set_field_manager (const std::shared_ptr<const fm_type>& field_mgr, const std::vector<std::string>& modes);

  std::shared_ptr<const fm_type> get_field_manager (const std::string& mode) const;

  void register_dimensions(const std::string& name);
  void register_variables(const std::string& filename, const std::string& fp_precision, const scorpio::FileMode mode);
  void set_decompositions(const std::string& filename);
  std::vector<scorpio::offset_t> get_var_dof_offsets (const FieldLayout& layout);
  void register_views();
  Field get_field(const std::string& name, const std::string& mode) const;
  void compute_diagnostic (const std::string& name, const bool allow_invalid_fields = false);
  void set_diagnostics();
  std::shared_ptr<AtmosphereDiagnostic>
  create_diagnostic (const std::string& diag_name);

  // Tracking the averaging of any filled values:
  void set_avg_cnt_tracking(const std::string& name, const FieldLayout& layout);

  // --- Internal variables --- //
  ekat::Comm                          m_comm;

  // We store two shared pointers for field managers:
  // io_field_manager stores the fields in the layout for output
  // sim_field_manager points to the simulation field manager
  // when remapping horizontally these two field managers may be different.
  std::map<std::string,std::shared_ptr<const fm_type>> m_field_mgrs;
  std::shared_ptr<const grid_type>            m_io_grid;
  std::shared_ptr<remapper_type>              m_horiz_remapper;
  std::shared_ptr<remapper_type>              m_vert_remapper;
  std::shared_ptr<const gm_type>              m_grids_manager;

  // How to combine multiple snapshots in the output: Instant, Max, Min, Average
  OutputAvgType     m_avg_type;
  Real              m_avg_coeff_threshold = 0.5; // % of unfilled values required to not just assign value as FillValue

  // Internal maps to the output fields, how the columns are distributed, the file dimensions and the global ids.
  std::vector<std::string>                              m_fields_names;
  std::vector<std::string>                              m_avg_cnt_names;
  std::map<std::string,std::string>                     m_field_to_avg_cnt_map;
  std::map<std::string,std::string>                     m_field_to_avg_cnt_suffix;
  std::map<std::string,FieldLayout>                     m_layouts;
  std::map<std::string,int>                             m_dims;
  std::map<std::string,std::shared_ptr<atm_diag_type>>  m_diagnostics;
  std::map<std::string,std::vector<std::string>>        m_diag_depends_on_diags;
  std::map<std::string,bool>                            m_diag_computed;
  DefaultMetadata                                       m_default_metadata;

  // Use float, so that if output fp_precision=float, this is a representable value.
  // Otherwise, you would get an error from Netcdf, like
  //   NetCDF: Numeric conversion not representable
  // Also, by default, don't pick max float, to avoid any overflow if the value
  // is used inside other calculation and/or remap.
  float m_fill_value = constants::DefaultFillValue<float>().value;

  // Local views of each field to be used for "averaging" output and writing to file.
  std::map<std::string,view_1d_host>    m_host_views_1d;
  std::map<std::string,view_1d_dev>     m_dev_views_1d;

  bool m_add_time_dim;
  bool m_track_avg_cnt = false;

  // The logger to be used throughout the ATM to log message
  std::shared_ptr<ekat::logger::LoggerBase> m_atm_logger;
};

} //namespace scream

#endif // SCREAM_SCORPIO_OUTPUT_HPP
