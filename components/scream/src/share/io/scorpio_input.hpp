#ifndef SCREAM_SCORPIO_INPUT_HPP
#define SCREAM_SCORPIO_INPUT_HPP

#include "share/io/scream_scorpio_interface.hpp"

#include "share/field/field_manager.hpp"
#include "share/grid/abstract_grid.hpp"
#include "scream_config.h"

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
 *  At construction time input instances require
 *  1. an EKAT comm group and
 *  2. a EKAT parameter list
 *  3. a shared pointer to the field manager
 *  The parameter list contains all of the essential data regarding
 *  the input file name and the set of variables to be read.  The parameter list can be
 *  created locally in the process needing input, see src/share/io/tests/ for examples of
 *  setting up the input parameter list.
 *
 *  A typical input parameter list looks like:
 *  -----
 *  Input Parameters
 *    FILENAME: STRING
 *    GRID: STRING
 *    FIELDS: [field_name_1, field_name_2, ..., field_name_N]
 *  -----
 *  where,
 *  FILENAME: is a string value of the name of the input file to be read.
 *  GRID: is a string of the grid the input file is written on, currently only "Physics" is supported.
 *  FIELDS: list of names of fields to load from file. Should match the name in the file and the name in the field manager.
 *  TODO: add a rename option if variable names differ in file and field manager.
 *
 * Usage:
 * 1. Construct an instance of the AtmosphereInput class:
 *    AtmosphereInput in_object(comm,params,field_mgr);
 * 2. Use read_variables to gather the registered variables
 *    in_object.read_variables();
 *
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 */

namespace scream
{

class AtmosphereInput 
{
public:
  using fm_type   = FieldManager<Real>;
  using grid_type = AbstractGrid;

  using KT = KokkosTypes<DefaultDevice>;
  template<int N>
  using view_Nd_host = typename KT::template view_ND<Real,N>::HostMirror;
  using view_1d_host = view_Nd_host<1>;

  // --- Constructor(s) & Destructor --- //
  AtmosphereInput (const ekat::Comm& comm,
                   const ekat::ParameterList& params);
  AtmosphereInput (const ekat::Comm& comm,
                   const ekat::ParameterList& params,
                   const std::shared_ptr<const fm_type>& field_mgr);
  AtmosphereInput (const ekat::Comm& comm,
                   const ekat::ParameterList& params,
                   const std::shared_ptr<const grid_type>& grid,
                   const std::map<std::string,view_1d_host>& host_views_1d,
                   const std::map<std::string,FieldLayout>&  layouts);

  virtual ~AtmosphereInput () = default;

  // --- Methods --- //
  // Sets up the scorpio metadata to preare for reading
  // The first version manually specifies grid, host 1d views, and layouts.
  // The second version grabs everything from the field manager.
  // NOTE: the first version cannot handle padded fields, while the latter
  //       can, and can also handle subfields.
  void init(const std::shared_ptr<const grid_type>& grid,
            const std::map<std::string,view_1d_host>& host_views_1d,
            const std::map<std::string,FieldLayout>&  layouts);
  void init(const std::shared_ptr<const fm_type>& field_mgr);

  // Read fields that were required via parameter list.
  void read_variables ();
  int read_int_scalar (const std::string& name);
  void finalize();

protected:
  // Internal functions
  void set_parameters (const ekat::ParameterList& params);
  void set_grid (const std::shared_ptr<const AbstractGrid>& grid);
  void set_field_manager (const std::shared_ptr<const fm_type>& field_mgr);

  void init_scorpio_structures ();
  void register_variables();
  void set_degrees_of_freedom();

  std::vector<std::string> get_vec_of_dims (const FieldLayout& layout);
  std::string get_io_decomp (const std::vector<std::string>& vec_of_dims);
  std::vector<int> get_var_dof_offsets (const FieldLayout& layout);

  // Internal variables
  ekat::ParameterList m_params;
  ekat::Comm          m_comm;

  std::shared_ptr<const fm_type>        m_field_mgr;
  std::shared_ptr<const AbstractGrid>   m_grid;

  std::map<std::string, view_1d_host>   m_host_views_1d;
  std::map<std::string, FieldLayout>    m_layouts;
  
  std::string               m_filename;
  std::vector<std::string>  m_fields_names;

  // Whether we are reading the history restart file of a model output (see output class for details)
  bool m_is_history_restart = false;
  bool m_is_inited = false;

}; // Class AtmosphereInput

} //namespace scream

#endif // SCREAM_SCORPIO_INPUT_HPP
