#ifndef SCREAM_SCORPIO_OUTPUT_HPP
#define SCREAM_SCORPIO_OUTPUT_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/manager/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"
#include "share/io/eamxx_io_utils.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/util/eamxx_utils.hpp"

#include <ekat_comm.hpp>
#include <ekat_parameter_list.hpp>

/*  The AtmosphereOutput class handles an output stream in SCREAM.
 *  Typical usage is to register an AtmosphereOutput object with the OutputManager (see
 eamxx_output_manager.hpp
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
 *  filename_prefix:                    STRING
 *  averaging_type:                     STRING
 *  max_snapshots_per_file:             INT                   (default: 1)
 *  fields:
 *     GRID_NAME_1:
 *        field_names:                  ARRAY OF STRINGS
 *     GRID_NAME_2:
 *        field_names:                  ARRAY OF STRINGS
 *        output_data_layout:           STRING                (default: 'default')
 *     ...
 *     GRID_NAME_N:
 *        field_names:                  ARRAY OF STRINGS
 *  output_control:
 *    frequency:                        INT
 *    frequency_units:                  STRING                (default: nsteps)
 *  restart:
 *    filename_prefix:                  STRING                (default: ${filename_prefix})
 *    skip_restart_if_rhist_not_found:  BOOL                  (default: false)
 *  -----
 *  The meaning of these parameters is the following:
 *  - filename_prefix: the output filename root.
 *  - averaging_type: a string that describes which type of output, current options are:
 *      instant - no averaging, output each snap as is.
 *      average - average of the field over some interval.
 *      min     - minimum value of the field over time interval.
 *      max     - maximum value of the field over time interval.
 *    Here, 'time interval' is described by ${Output frequency} and ${Output frequency_units}.
 *    E.g., with 'Output frequency'=10 and 'Output frequency_units'="Days", the time interval is 10
 days.
 *  - fields: parameters specifying fields to output
 *     - GRID_NAME: parameters specifyign fields to output from grid $GRID_NAME
 *        - field_names: names of fields defined on grid $grid_name that need to be outputed
 *        - output_data_layout: attempt to 'remap' fields to this data layout first.
 *          This option is mostly used to enable dyn->phys_gll remap (to save storage),
 *          but can be expanded to support other online layout changes (e.g., transpose
 *          (ncol,[...,]nlev) layouts). The default value 'default' will use the default
 *          for this grid. Usually, 'default' does the same thing as 'native' (which does
 *          no layout change, and outputs the fields as they appear on this grid), but for
 *          the dynamics grid the default behavior is to remap to phys_gll layout.
 *          Other values are available on a per grid basis. E.g., the dynamics grid
 *          supports 'dg' (same as 'native') and 'cg' (which implements the dyn->phys_gll
 *          remap, which is the same behavior as 'default').
 *  - max_snapshots_per_file: the maximum number of snapshots saved per file. After this many
 *    snapshots, the current files is closed and a new file created.
 *  - Output: parameters for output control
 *    - frequency: the frequency of output writes (in the units specified by ${Output
 frequency_units})
 *    - frequency_units: the units of output frequency (nsteps, nmonths, nyears, nhours, ndays,...)
 *      snapshots have been written on a single nc file, the class will close the file, and open a
 new one
 *  - Checkpointing: parameters for checkpointing control
 *    - frequency: the frequenct of checkpoints writes. This option is used/matters only if
 *      if averaging_type is *not* instant. A value of 0 is interpreted as 'no checkpointing'.
 *    - frequency_units: the units of restart history output.
 *  - restart: parameters for history restart
 *    - filename_prefix: the history restart filename root.
 *    - skip_restart_if_rhist_not_found: if this is a restarted run and this is true, skip the
 *      hist restart if the proper filename is not found in rpointer. Allows to add a new stream
 *      upon restart.

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
 *   - in case of single-grid tests, you can specify fields names by adding 'field_names' directly
 *     in the top-level parameter list.
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
  using diag_ptr_type = std::shared_ptr<AtmosphereDiagnostic>;

  ~AtmosphereOutput();

  // Constructor
  AtmosphereOutput(const ekat::Comm &comm, const ekat::ParameterList &params,
                   const std::shared_ptr<const fm_type> &field_mgr, const std::string &grid_name);

  // Short version for outputing a list of fields (no remapping supported)
  AtmosphereOutput(const ekat::Comm &comm, const std::vector<Field> &fields,
                   const std::shared_ptr<const grid_type> &grid);

  // Main Functions
  void restart(const std::string &filename);
  void init();
  void reset_scorpio_fields();
  void setup_output_file(const std::string &filename, const std::string &fp_precision,
                         const scorpio::FileMode mode);

  void init_timestep(const util::TimeStamp &start_of_step);
  void run(const std::string &filename, const bool output_step, const bool checkpoint_step,
           const int nsteps_since_last_output, const bool allow_invalid_fields = false);

  long long res_dep_memory_footprint() const;

  std::shared_ptr<const AbstractGrid>
  get_io_grid() const
  {
    return m_io_grid;
  }

  void set_logger(const std::shared_ptr<ekat::logger::LoggerBase> &atm_logger);

protected:
  template <typename T> using strmap_t = std::map<std::string, T>;

  using strvec_t = std::vector<std::string>;

  // Internal functions
  void register_variables(const std::string &filename, const std::string &fp_precision,
                          const scorpio::FileMode mode);
  void set_decompositions(const std::string &filename);
  void compute_diagnostics(const bool allow_invalid_fields);
  void process_requested_fields();
  strvec_t get_var_dimnames(const FieldLayout &layout) const;

  // Tracking the averaging of any filled values:
  void set_avg_cnt_tracking(const std::string &name, const FieldLayout &layout);

  // --- Internal variables --- //
  ekat::Comm m_comm;

  // We store separate shared pointers for field mgrs at different stages of IO:
  // More specifically, the order of operations is as follows:
  //  - compute diags (if any)
  //  - vert remap (if any)
  //  - horiz remap (if any)
  //  - call scorpio
  // and the field mgrs are connected by the following ops
  //         VERT_REMAP        HORIZ_REMAP       TALLY_UPDATE
  //  FromModel -> AfterVertRemap -> AfterHorizRemap -> Scorpio
  // The last 2 field mgrs contain DIFFERENT fields if 1+ of the following happens:
  //  - fields are padded: we have PackSize>1 during remaps, but Scorpio needs CONTIGUOUS memory
  //  - there's no remap (so the first 3 FM are the same), but the field is a subfield: again NOT
  //  CONTIGUOUS
  //  - the avg type is NOT instant: we need a separate Field to store the tallies
  // Also, FromModel is NOT the same field mgr as stored in the AD. In particular, it is a "clone"
  // of the AD field mgr but restricted to the grid that this object is handling, AND we stuff all
  // diags in this field mgr (so that we do not pollute the AD field mgr with output-only fields).
  // NOTE: if avg_type!=Instant, then ALL fields in the last two field mgrs are different, otherwise
  // SOME field
  //       MAY be the same. E.g., field that are NOT subfields and are NOT padded can be "soft
  //       copies", to reduce memory footprint and runtime costs.
  enum Phase {
    FromModel,       // Output fields as from the model (or diags computed from model fields)
    AfterVertRemap,  // Output fields after vertical remap
    AfterHorizRemap, // Output fields after horiz remap
    Scorpio // Output fields to pass to scorpio (may differ from the above in case of packing)
  };
  std::map<Phase, std::shared_ptr<fm_type>> m_field_mgrs;

  std::shared_ptr<const grid_type> m_io_grid;
  std::shared_ptr<remapper_type> m_horiz_remapper;
  std::shared_ptr<remapper_type> m_vert_remapper;

  // How to combine multiple snapshots in the output: instant, Max, Min, Average
  OutputAvgType m_avg_type;
  Real m_avg_coeff_threshold =
      0.5; // % of unfilled values required to not just assign value as FillValue

  // Internal maps to the output fields, how the columns are distributed, the file dimensions and
  // the global ids.
  strvec_t m_fields_names;
  strmap_t<Field> m_field_to_avg_count;
  std::vector<Field> m_avg_counts;
  strmap_t<std::string> m_field_to_avg_cnt_suffix;
  strmap_t<strvec_t> m_vars_dims;
  strmap_t<int> m_dims_len;
  std::list<diag_ptr_type> m_diagnostics;

  static strmap_t<diag_ptr_type> m_diag_repo;

  // Field aliasing support
  strmap_t<std::string> m_alias_to_orig; // Map from alias names to original names (used to set io attribute)

  DefaultMetadata m_default_metadata;

  bool m_add_time_dim;
  bool m_track_avg_cnt         = false;
  std::string m_decomp_dimname = "";

  std::shared_ptr<ekat::logger::LoggerBase> m_atm_logger =
      console_logger(ekat::logger::LogLevel::warn);

  std::string m_stream_name; // used in error msgs to help distinguish which stream this is
};

} // namespace scream

#endif // SCREAM_SCORPIO_OUTPUT_HPP
