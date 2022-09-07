#ifndef SCREAM_SCORPIO_INPUT_HPP
#define SCREAM_SCORPIO_INPUT_HPP

#include "share/io/scream_scorpio_interface.hpp"
#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/grid/grids_manager.hpp"

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/mpi/ekat_comm.hpp"

/*  The AtmosphereInput class handles all input streams to SCREAM.
 *  It is important to note that there does not exist an InputManager,
 *  like in the case of output.  So all input streams have to be managed
 *  directly by the process that requires it.
 *
 *  The typical input call will be to the outward facing routine 'pull_input'.
 *
 *  Currently, input can only handle single timesnap input files.  In other words
 *  files that will be opened, read and closed within the same timestep and only
 *  store one timesnap of data.
 *
 *  Note: the init, and finalize are separate routines that are outward facing in
 *  this class to facilitate cases where reading input over some number of simulation
 *  timesteps is possible.
 *
 *  At construction time ALL output instances require at least an EKAT comm group and
 *  an EKAT parameter list. The following arguments depend on the input class use case.
 *  See constructors documentation for more info.
 *
 *  The EKAT parameter list contains the following options to control input behavior:
 *  -----
 *  Input Parameters
 *    Filename: STRING
 *    Fields:   ARRAY OF STRINGS
 *  -----
 *  The meaning of these parameters is the following:
 *   - Filename: the name of the input file to be read.
 *   - Fields: list of names of fields to load from file. Should match the name in the file and the name in the field manager.
 *  Note: you can specify lists (such as the 'Fields' list above) with either of the two syntaxes
 *    Fields: [field_name1, field_name2, ... , field_name_N]
 *    Fields:
 *      - field_name_1
 *      - field_name_2
 *        ...
 *      - field_name_N
 *  Note: an alternative way of specifying Fields names is to have
 *    Grid: STRING
 *    Fields:
 *      $GRID: [field_name1,...,field_name_N]
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
  using gm_type       = GridsManager;
  using remapper_type = AbstractRemapper;

  using KT = KokkosTypes<DefaultDevice>;
  template<int N>
  using view_Nd_host = typename KT::template view_ND<Real,N>::HostMirror;
  using view_1d_host = view_Nd_host<1>;

  // --- Constructor(s) & Destructor --- //
  // Creates bare input. Will require a call to one of the two 'init' methods.
  // Constructor inputs:
  //  - comm: the MPI comm used for I/O operations. Notice that the
  //          associated MPI group *must* be contained in the MPI
  //          group in the overall atm comm, but it *might* be smaller
  //          (and likely *will* be smaller, for large scale runs).
  //  - params: a parameter list containing info on the file/fields to load.
  AtmosphereInput (const ekat::Comm& comm,
                   const ekat::ParameterList& params);
  // Creates input to read into FieldManager-owned fields.
  // Constructor inputs:
  //  - comm: the MPI comm used for I/O operations. Notice that the
  //          associated MPI group *must* be contained in the MPI
  //          group in the overall atm comm, but it *might* be smaller
  //          (and likely *will* be smaller, for large scale runs).
  //  - params: a parameter list containing info on the file/fields to load.
  //  - field_mgr: the FieldManager containing the Field's where the
  //               variables from the input filed will be read into.
  //               Fields can be padded/strided.
  // It calls init(field_mgr) at the end.
  // TODO: is comm superfluous, considering we can get it from the grid in the field_mgr?
  AtmosphereInput (const ekat::ParameterList& params,
                   const std::shared_ptr<const fm_type>& field_mgr,
                   const std::shared_ptr<const gm_type>& grids_mgr = nullptr);

  // Creates input to read into user-provided flattened 1d host views.
  // Constructor inputs:
  //  - comm: the MPI comm used for I/O operations. Notice that the
  //          associated MPI group *must* be contained in the MPI
  //          group in the overall atm comm, but it *might* be smaller
  //          (and likely *will* be smaller, for large scale runs).
  //  - params: a parameter list containing info on the file/fields to load.
  //  - grid: the grid where the variables live
  //  - host_views_1d: the 1d flattened views where data will be read into.
  //                   These views must be contiguous (no padding/striding).
  //  - layouts: the layout of the vars (used to reshape the views).
  // It calls init(grid,host_views_1d,layouts) at the end.
  // TODO: do not require layouts, and read them from file.
  // TODO: is comm superfluous, considering we can get it from the grid?
  AtmosphereInput (const ekat::ParameterList& params,
                   const std::shared_ptr<const grid_type>& grid,
                   const std::map<std::string,view_1d_host>& host_views_1d,
                   const std::map<std::string,FieldLayout>&  layouts);

  virtual ~AtmosphereInput () = default;

  // --- Methods --- //
  // In case the class was constructed with the minimal ctor, these methods
  // allow to finalize initialization later.
  // NOTE: these two init methods are mutually exclusive
  void init (const std::shared_ptr<const fm_type>& field_mgr,
             const std::shared_ptr<const gm_type>& grids_mgr = nullptr);
  void init (const std::shared_ptr<const grid_type>& grid,
             const std::map<std::string,view_1d_host>& host_views_1d,
             const std::map<std::string,FieldLayout>&  layouts);

  // Read fields that were required via parameter list.
  void read_variables (const int time_index = -1);
  int read_int_scalar (const std::string& name);
  void finalize();

protected:

  void set_fields_and_grid_names (const std::vector<std::string>& grid_aliases);
  void build_remapper (const std::shared_ptr<const gm_type>& grids_mgr);
  void set_grid (const std::shared_ptr<const AbstractGrid>& grid);
  void set_field_manager (const std::shared_ptr<const fm_type>& field_mgr,
                          const std::shared_ptr<const gm_type>& grids_mgr);
  void set_views (const std::map<std::string,view_1d_host>& host_views_1d,
                  const std::map<std::string,FieldLayout>&  layouts);
  void init_scorpio_structures ();

  void register_fields_specs ();

  void register_variables();
  void set_degrees_of_freedom();

  std::vector<std::string> get_vec_of_dims (const FieldLayout& layout);
  std::string get_io_decomp (const FieldLayout& layout);
  std::vector<scorpio::offset_t> get_var_dof_offsets (const FieldLayout& layout);

  // Internal variables
  ekat::Comm            m_comm;
  ekat::ParameterList   m_params;

  std::shared_ptr<const fm_type>        m_field_mgr;
  std::shared_ptr<const AbstractGrid>   m_io_grid;
  std::shared_ptr<remapper_type>        m_remapper;

  std::map<std::string, view_1d_host>   m_host_views_1d;
  std::map<std::string, FieldLayout>    m_layouts;
  
  std::string               m_filename;
  std::string               m_io_grid_name;
  std::vector<std::string>  m_fields_names;

  bool m_inited_with_fields        = false;
  bool m_inited_with_views         = false;
}; // Class AtmosphereInput

} //namespace scream

#endif // SCREAM_SCORPIO_INPUT_HPP
