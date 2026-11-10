#ifndef SCREAM_ABSTRACT_DIAGNOSTIC_HPP
#define SCREAM_ABSTRACT_DIAGNOSTIC_HPP

#include "share/field/field.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/util/eamxx_time_stamp.hpp"
#include "share/util/eamxx_bfbhash.hpp"

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
 * The base class for all diagnostic quantities.
 *
 * A diagnostic quantity is a Field that is not part of the model state,
 * nor it can impact it in any way. It is derived from the state, and usually
 * only computed for the sake of IO to reduce the footprint of what must
 * be written to file (as well as taking advantage of runtime compute capabilities,
 * which may make the calculation faster than in post-processing tools).
 *
 * The lifetime of a diagnostic is as follows:
 *  - constructor: create the diag following *exactly* the signature of the base class.
 *                 This is b/c diags are constructed at run time via a factory, which
 *                 requires a uniform signature for all constructors
 *  - set inputs: the diag provides the name of the input fields, which customers must
 *                retrieve and set in the diag before any attempt at computing it.
 *  - initialize: by the time this function is called, ALL input fields (if any) MUST
 *                have been set in the class (via the above method), while by the time
 *                it returns, the stored diagnostic field MUST be allocated. It calls the
 *                protected 'initialize_impl', which derived classes MUST implement.
 *  - init_timestep: this is optional (the base class has an empty impl), but provides
 *                   the diag the possibility of doing some start-of-step calculations.
 *                   E.g., the FieldPrev diag, which computes the field at the previous
 *                   timestep, can use this call to copy the field before the timestep begins.
 *  - compute: this is where the calculation is (usually) performed. It calls the private
 *             method 'compute_impl', which derived classes MUST implement.
 *  - get: this can be used to retrieve the diagnostic field (can be called once, as the
 *         stored field is allocated once at initialization and reused at every step).
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
  Field get () const;

  // Compute the diagnostic (skips if inputs have not changed since last call)
  void compute (const util::TimeStamp& ts);

protected:

  // Derived classes override this for any setup needed after fields are set.
  virtual void initialize_impl () { /* Nothing to do */ }

  // Derived classes implement this to compute the output from the inputs.
  virtual void compute_impl () = 0;

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

  bfbhash::HashType             m_last_eval_ts_hash = 0;

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
