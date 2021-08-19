#ifndef SCREAM_OUTPUT_MANAGER_HPP
#define SCREAM_OUTPUT_MANAGER_HPP

#include "share/field/field_manager.hpp"
#include "share/grid/grids_manager.hpp"

#include "share/io/scorpio_output.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_parse_yaml_file.hpp"

namespace scream
{

/*
 * The OutputManager Class handles all individual output streams.
 * Rather than have the SCREAM-AD or unit tests call the scorpio_output
 * class objects directly, the OutputManager stores these in a vector,
 * m_output_streams.
 *
 * Similar to a typical atmospheric process object, the OutputManager
 * has an init, run and finalize routine which is called during the AD
 * during those respective steps.
 *
 * PROPER USAGE:
 * The output manager requires a communication group, a parameter list of control
 * variables, a grid manager and a field manager.
 * Each of these four things are set using the setter function 'set_X' where
 *   X = comm, for the EKAT comm group.
 *     = params, for the parameter list.  In typical SCREAM runs this parameter
 *       list is a sublist of the scream control yaml file.
 *     = grids, for the grids mananger
 *     = fm, for the field manager.
 * The setup of the output manager in the SCREAM-AD is one of the last steps to
 * ensure that all four of the above objects have already been constructed.
 * see /control/atmospheric_driver.cpp for an example.
 *
 * For UNIT TESTS:
 * The output manager does require a comm group, list of parameters, grids manager
 * and field manager to work.  If output is desired in a unit test then these
 * must be established.  There are examples in /src/share/io/tests of how to
 * establish a simple grids manager and field manager.  As well as how to
 * locally create a parameter list.
 *
 * Adding output streams mid-simulation:
 * It is possible to add an output stream after init has been called by calling
 * the internal function 'new_output' which takes an EKAT parameter list as input.
 * See comments in new_output below for more details.
 *
 * --------------------------------------------------------------------------------
 *  (2020-10-21) Aaron S. Donahue (LLNL)
 */
class OutputManager
{
public:
  using fm_type = FieldManager<Real>;

  // Constructor(s) & Destructor
  OutputManager () = default;
  virtual ~OutputManager () = default;

  void setup (const ekat::Comm& io_comm, const ekat::ParameterList& params,
              const std::shared_ptr<const fm_type>& field_mgr,
              const bool runtype_restart);
  void run(util::TimeStamp& current_ts);
  void finalize();

protected:
  // Add an output stream.
  void new_output(const ekat::ParameterList& params, const bool model_restart_output);

  // Craft the restart parameter list
  void make_restart_param_list(ekat::ParameterList& params);

  using output_type = AtmosphereOutput;
  using output_ptr_type = std::shared_ptr<output_type>;

  std::vector<output_ptr_type>        m_output_streams;
  ekat::Comm                          m_io_comm;
  ekat::ParameterList                 m_params;
  std::shared_ptr<const fm_type>      m_field_mgr;
  std::string                         m_ref_grid_name;

  bool m_runtype_restart  = false;

}; // class OutputManager

} // namespace scream

#endif // SCREAM_OUTPUT_MANAGER_HPP
