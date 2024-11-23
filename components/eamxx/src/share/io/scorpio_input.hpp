#ifndef SCREAM_SCORPIO_INPUT_HPP
#define SCREAM_SCORPIO_INPUT_HPP

#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/logging/ekat_logger.hpp"

/*  The AtmosphereInput class handles all input streams to SCREAM.
 *  It is important to note that there does not exist an InputManager,
 *  like in the case of output.  So all input streams have to be managed
 *  directly by the process that requires it. The reason is that input
 *  operations are less convoluted than output ones (e.g., no averaging).
 *
 *  Note: the init, and finalize are separate routines that are outward facing in
 *  this class to facilitate cases where reading input over some number of simulation
 *  timesteps is possible.
 *
 *  The EKAT parameter list contains the following options to control input behavior:
 *  -----
 *  Input Parameters
 *    Filename: STRING
 *    Field Names:   ARRAY OF STRINGS
 *  -----
 *  The meaning of these parameters is the following:
 *   - Filename: the name of the input file to be read.
 *   - Field Names: list of names of fields to load from file. Should match the name in the file and the name in the field manager.
 *
 *  TODO: add a rename option if variable names differ in file and field manager.
 *
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 *  (2021-08-19) Luca Bertagna (SNL)
 */

namespace scream
{

class AtmosphereInput 
{
public:
  using fm_type       = FieldManager;
  using grid_type     = AbstractGrid;

  using KT = KokkosTypes<DefaultDevice>;
  template<int N>
  using view_Nd_host = typename KT::template view_ND<Real,N>::HostMirror;
  using view_1d_host = view_Nd_host<1>;

  // --- Constructor(s) & Destructor --- //
  // NOTE: non-trivial constructors simply call the corresponding init method
  AtmosphereInput () = default;
  AtmosphereInput (const ekat::ParameterList& params,
                   const std::shared_ptr<const fm_type>& field_mgr);
  AtmosphereInput (const ekat::ParameterList& params,
                   const std::shared_ptr<const grid_type>& grid,
                   const std::map<std::string,view_1d_host>& host_views_1d,
                   const std::map<std::string,FieldLayout>&  layouts);
  AtmosphereInput (const std::string& filename,
                   const std::shared_ptr<const grid_type>& grid,
                   const std::vector<Field>& fields,
                   const bool skip_grid_checks = false);
  // This constructor only sets the minimal info, deferring initialization
  // to when set_field_manager/reset_fields and reset_filename are called
  AtmosphereInput (const std::vector<std::string>& fields_names,
                   const std::shared_ptr<const grid_type>& grid);

  // Due to resource acquisition (in scorpio), avoid copies
  AtmosphereInput (const AtmosphereInput&) = delete;
  ~AtmosphereInput ();

  // Due to resource acquisition (in scorpio), avoid copies
  AtmosphereInput& operator= (const AtmosphereInput&) = delete;

  // --- Methods --- //
  // Initialize the class for reading into FieldManager-owned fields.
  //  - params: input parameters (must contain at least "Filename")
  //  - field_mgr: the FieldManager containing the Field's where the
  //               variables from the input filed will be read into.
  //               Fields can be padded/strided.
  void init (const ekat::ParameterList& params,
             const std::shared_ptr<const fm_type>& field_mgr);

  // Initialize the class for reading into user-provided flattened 1d host views.
  //  - params: input parameters (must contain at least "Filename")
  //  - grid: the grid where the variables live
  //  - host_views_1d: the 1d flattened views where data will be read into.
  //                   These views must be contiguous (no padding/striding).
  //  - layouts: the layout of the vars (used to reshape the views).
  void init (const ekat::ParameterList& params,
             const std::shared_ptr<const grid_type>& grid,
             const std::map<std::string,view_1d_host>& host_views_1d,
             const std::map<std::string,FieldLayout>&  layouts);

  // Read fields that were required via parameter list.
  void read_variables (const int time_index = -1);

  // Cleans up the class
  void finalize();

  // Getters
  std::string get_filename() { return m_filename; } // Simple getter to query the filename for this stream.

  // Expose the ability to set/reset fields/field_manager for cases like data interpolation,
  // where we swap pointers but all the scorpio data structures are unchanged.
  void set_field_manager (const std::shared_ptr<const fm_type>& field_mgr);
  void set_fields (const std::vector<Field>& fields);
  void reset_filename (const std::string& filename);

  // Option to add a logger
  void set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& atm_logger) {
      m_atm_logger = atm_logger;
  }
protected:

  void set_grid (const std::shared_ptr<const AbstractGrid>& grid);
  void set_views (const std::map<std::string,view_1d_host>& host_views_1d,
                  const std::map<std::string,FieldLayout>&  layouts);
  void init_scorpio_structures ();

  void set_decompositions();

  std::vector<std::string> get_vec_of_dims (const FieldLayout& layout);

  // Internal variables
  ekat::ParameterList   m_params;

  std::shared_ptr<const fm_type>        m_field_mgr;
  std::shared_ptr<const AbstractGrid>   m_io_grid;

  std::map<std::string, view_1d_host>   m_host_views_1d;
  std::map<std::string, FieldLayout>    m_layouts;
  
  std::string               m_filename;
  std::vector<std::string>  m_fields_names;

  bool m_inited_with_fields        = false;
  bool m_inited_with_views         = false;

  // The logger to be used throughout the ATM to log message
  std::shared_ptr<ekat::logger::LoggerBase> m_atm_logger;
}; // Class AtmosphereInput

} //namespace scream

#endif // SCREAM_SCORPIO_INPUT_HPP
