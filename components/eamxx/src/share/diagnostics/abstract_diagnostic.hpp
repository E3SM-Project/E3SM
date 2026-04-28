#ifndef SCREAM_ABSTRACT_DIAGNOSTIC_HPP
#define SCREAM_ABSTRACT_DIAGNOSTIC_HPP

#include "share/field/field.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_comm.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_factory.hpp>
#include <ekat_string_utils.hpp>

#include <list>
#include <map>
#include <memory>
#include <string>

namespace scream
{

/*
 * A class representing a diagnostic.
 *
 * A diagnostic takes zero or more fields as input and produces a single
 * output field. It is similar to AtmosphereProcess, but  does not participate
 * in the atmosphere time-stepping loop; it is only used for output/IO.
 */

class AbstractDiagnostic
{
public:
  AbstractDiagnostic (const ekat::Comm& comm,
                      const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid);

  virtual ~AbstractDiagnostic () = default;

  // The name of the diagnostic
  virtual std::string name () const = 0;

  // Initialize the diagnostic. Calls initialize_impl().
  void initialize ();

  // Returns true after initialize() has been called successfully.
  bool is_initialized () const { return m_is_initialized; }

  // Return the parameter list
  ekat::ParameterList& get_params () { return m_params; }

  // Return the list of input fields needed by this diagnostic
  const std::list<std::string>& get_input_fields_names () const { return m_field_in_names; }

  // Store input field in the map
  void set_input_field (const Field& f);

  // Allows the diagnostic to save some start-of-step quantity (e.g., in case
  // we need to compute tendencies, or accumulated stuff)
  virtual void init_timestep (const util::TimeStamp& /* start_of_step */) {}

  // Returns the diagnostic output field
  Field get_diagnostic () const;

  // Compute the diagnostic (skips if inputs have not changed since last call)
  void compute_diagnostic (const util::TimeStamp& ts);

protected:

  // Derived classes override this for any setup needed after fields are set.
  virtual void initialize_impl () { /* Nothing to do */ }

  // Derived classes implement this to compute the output from the inputs.
  virtual void compute_diagnostic_impl () = 0;

  // MPI communicator
  ekat::Comm              m_comm;

  // Parameter list
  ekat::ParameterList     m_params;

  // The grid on which this diagnostic is defined (set via set_grid())
  std::shared_ptr<const AbstractGrid> m_grid;

  // Diagnostics produce a single output field
  Field m_diagnostic_output;

  // Timestamp of the last diag evaluation
  util::TimeStamp m_last_eval_ts;

  // Input fields
  std::list<std::string>        m_field_in_names;
  std::map<std::string,Field>   m_fields_in;

  bool m_is_initialized = false;
};

// Factory for atmosphere diagnostics
using DiagnosticFactory =
    ekat::Factory<AbstractDiagnostic,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AbstractDiagnostic>,
                  const ekat::Comm&,const ekat::ParameterList&,
                  const std::shared_ptr<const AbstractGrid>&>;

// Convenience functions to create a diagnostic (with and without a grid)
template <typename AtmDiagType>
inline std::shared_ptr<AbstractDiagnostic>
create_diagnostic (const ekat::Comm& comm,
                   const ekat::ParameterList& p,
                   const std::shared_ptr<const AbstractGrid>& grid)
{
  return std::make_shared<AtmDiagType>(comm,p,grid);
}

} //namespace scream

#endif // SCREAM_ABSTRACT_DIAGNOSTIC_HPP
