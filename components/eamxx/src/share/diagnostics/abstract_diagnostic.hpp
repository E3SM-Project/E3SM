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
#include <vector>

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
 *
 * Most diagnostics compute exactly one field, which they store in m_diagnostic_output.
 * A diagnostic that computes several fields in a single pass (e.g. a satellite
 * simulator, which produces a whole suite of outputs from one expensive call)
 * declares the extra ones with add_output(), and they are then reachable by name
 * via get(name). Nothing changes for the single-output case: get() still returns
 * m_diagnostic_output, and a diag that never calls add_output() behaves exactly
 * as it did before.
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

  // Return the input fields themselves, keyed by name. Available once they have
  // all been set. Prefer this over looking the names up in a FieldManager: an
  // input may be an intermediate that no FieldManager knows about.
  const std::map<std::string,Field>& get_input_fields () const { return m_fields_in; }

  // Store input field in the map
  void set_input_field (const Field& f);

  // Allows the diagnostic to save some start-of-step quantity (e.g., in case
  // we need to compute tendencies, or accumulated stuff)
  virtual void init_timestep (const util::TimeStamp& /* start_of_step */) {}

  // Returns the diagnostic output field. For a diag computing several fields,
  // this is the primary one: m_diagnostic_output if it was set, and otherwise
  // the first output declared via add_output().
  Field get () const;

  // Returns the output field named 'name'. Throws if this diag does not compute
  // a field by that name.
  Field get (const std::string& name) const;

  // Whether this diag computes a field named 'name'.
  bool has_output (const std::string& name) const;

  // All the fields this diag computes, primary first.
  std::vector<Field> get_outputs () const;

  // Names of all the fields this diag computes, primary first.
  std::vector<std::string> get_output_names () const;

  // Compute the diagnostic (skips if inputs have not changed since last call)
  void compute (const util::TimeStamp& ts);

protected:

  // Derived classes override this for any setup needed after fields are set.
  virtual void initialize_impl () { /* Nothing to do */ }

  // Derived classes implement this to compute the output from the inputs.
  virtual void compute_impl () = 0;

  // Declare an output field beyond m_diagnostic_output. This is what a diag
  // computing several fields uses; single-output diags keep assigning
  // m_diagnostic_output and never call this. A diag with no natural primary
  // output may leave m_diagnostic_output unset and declare all of its outputs
  // here, in which case the first one declared plays the role of the primary.
  // Call this from initialize_impl(), and only with allocated fields.
  void add_output (const Field& f);

  // MPI communicator
  ekat::Comm              m_comm;

  // Parameter list
  ekat::ParameterList     m_params;

  // The grid on which this diagnostic is defined (set via set_grid())
  std::shared_ptr<const AbstractGrid> m_grid;

  // The primary output field. Most diagnostics produce only this one.
  Field m_diagnostic_output;

  // Any further output field, for the diagnostics that compute more than one.
  // Kept separate from m_diagnostic_output so that single-output diags, and
  // everything that consumes them, are untouched. See add_output().
  std::vector<Field> m_extra_outputs;

  // Timestamp of the last diag evaluation
  util::TimeStamp m_last_eval_ts;

  // Timestamp compute() was called with. Set before compute_impl() runs, so a
  // diag that needs to know when it is being evaluated can read it there.
  util::TimeStamp m_current_ts;

  // Input fields
  std::list<std::string>        m_field_in_names;
  std::map<std::string,Field>   m_fields_in;

  bfbhash::HashType             m_last_eval_ts_hash = 0;

  bool m_is_initialized = false;

private:
  // The stored output named 'name', or nullptr if there is none.
  const Field* find_output (const std::string& name) const;
};

// Factory for atmosphere diagnostics
using DiagnosticFactory =
    ekat::Factory<AbstractDiagnostic,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AbstractDiagnostic>,
                  const ekat::Comm&,const ekat::ParameterList&,
                  const std::shared_ptr<const AbstractGrid>&>;

/*
 * The fields computed by the diagnostics that compute more than one.
 *
 * A diagnostic is asked for by the name of the field it produces, so a diag
 * computing several fields is reachable under several names, none of which is
 * the key it is registered in the factory under. Which diag to build for a given
 * field name must therefore be known *before* anything is built, which rules out
 * asking an instance: it is declared here, next to the factory registration.
 *
 * Single-output diags never appear here. They are found in the factory under
 * their own name, exactly as they always have been.
 */
class DiagOutputsRegistry
{
public:
  static DiagOutputsRegistry& instance ();

  // Declare that the diag registered in the factory under 'factory_key'
  // computes the fields named in 'outputs'.
  void register_outputs (const std::string& factory_key,
                         const std::vector<std::string>& outputs);

  // The factory key of the diag computing 'output_name', or an empty string if
  // no registered diag computes it.
  std::string provider_of (const std::string& output_name) const;

private:
  DiagOutputsRegistry () = default;

  // output field name -> factory key
  std::map<std::string,std::string> m_provider;
};

// Convenience functions to create a diagnostic (with and without a grid)
template <typename AtmDiagType>
inline std::shared_ptr<AbstractDiagnostic>
create_diagnostic (const ekat::Comm& comm,
                   const ekat::ParameterList& p,
                   const std::shared_ptr<const AbstractGrid>& grid)
{
  return std::make_shared<AtmDiagType>(comm,p,grid);
}

// Register a diag that computes several fields: it goes in the factory like any
// other, and its outputs are declared so that a request for any of them finds
// it. Registering with the factory alone is enough for a single-output diag.
template <typename AtmDiagType>
inline void
register_multi_output_diagnostic (const std::string& factory_key,
                                  const std::vector<std::string>& outputs)
{
  DiagnosticFactory::instance().register_product(factory_key,&create_diagnostic<AtmDiagType>);
  DiagOutputsRegistry::instance().register_outputs(factory_key,outputs);
}

} //namespace scream

#endif // SCREAM_ABSTRACT_DIAGNOSTIC_HPP
