#ifndef SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP
#define SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_process.hpp"

#include <map>
#include <string>
#include <vector>

namespace scream
{

/*
 * A class representing an atmosphere diagnostic.
 * Diagnostics are similar to atmosphere processes in that
 * they take field(s) as input and produce some output.
 *
 * This subclass definition makes a few distinctions that are
 * unique to diagnostics.
 *
 * 1) Fields from the field manager can only be accessed as "Required".
 *    Meaning that an atmosphere diagnostic can not compute or update
 *    a field in the field manager.
 *    - note, if you want to update something in the field manager it
 *            should happen in an atmosphere process.
 * 2) Diagnostics produce one or more field outputs.
 *    Single-output diagnostics use m_diagnostic_output.
 *    Multi-output diagnostics use the m_diagnostic_outputs map.
 *    Typically the field output will be to buffered memory, meaning it
 *    could be overwritten later in the simulation.
 *    A diagnostic output is meant to be used for OUTPUT or a PROPERTY CHECK.
 */

// TODO: inheriting from AtmosphereProcess is conceptually wrong. It was done
//       out of convenience, but we should revisit that choice.

class AtmosphereDiagnostic : public AtmosphereProcess
{
public:

  // Constructor(s)
  AtmosphereDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  virtual ~AtmosphereDiagnostic () = default;

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // Most (all?) diagnostics will be defined on the physics grid, so we use that
  // by default. Derived classes can, of course, override this.
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert("physics");
    return s;
  }

  // Getting diagnostic output(s)
  // Single-output diagnostics use get_diagnostic() with no args.
  // Multi-output diagnostics use get_diagnostic(name) and get_diagnostic_names().
  Field get_diagnostic () const;
  Field get_diagnostic (const std::string& name) const;

  // Return the names of all diagnostic output fields.
  // Multi-output diagnostics should override this.
  virtual std::vector<std::string> get_diagnostic_names () const;

  // Whether this diagnostic produces multiple output fields
  bool is_multi_output () const {
    return m_diagnostic_outputs.size() > 1;
  }

  // Allows the diagnostic to save some start-of-step quantity (e.g., in case
  // we need to compute tendencies, or accumulated stuff)
  virtual void init_timestep (const util::TimeStamp& /* start_of_step */) {}

  void compute_diagnostic (const double dt = 0);
protected:

  void set_required_field_impl (const Field& f) final;
  void set_computed_field_impl (const Field& f) final;
  void set_computed_group_impl (const FieldGroup& group) final;

  virtual void compute_diagnostic_impl () = 0;

  // By default, diagnostic don't do any initialization/finalization stuff.
  // Derived classes can override, of course
  void initialize_impl (const RunType /*run_type*/) { /* Nothing to do */ }
  void run_impl (const double dt);
  void finalize_impl   () { /* Nothing to do */ }

  // Some diagnostics will need the timestep, store here.
  double m_dt;

  // Single-output diagnostics store their output here.
  // These two storages are mutually exclusive:
  // - Single-output diags populate m_diagnostic_output directly
  // - Multi-output diags populate m_diagnostic_outputs map in create_requests()
  // Public accessors check both.
  Field m_diagnostic_output;

  // Multi-output diagnostics store their outputs here (keyed by field name).
  std::map<std::string, Field> m_diagnostic_outputs;

  // Timestamp of the last diag evaluation
  util::TimeStamp m_last_eval_ts;
};

// A short name for the factory for atmosphere diagnostics
// WARNING: you do not need to write your own creator function to register your atmosphere diagnostic in the factory.
//          You could, but there's no need. You can simply register the common one right below, using your
//          atmosphere diagnostic subclass name as templated argument. If you roll your own creator function, you
//          *MUST* ensure that it correctly sets up the self pointer after creating the shared_ptr.
//          This is *necessary* until we can safely switch to std::enable_shared_from_this.
//          For more details, see the comments at the top of ekat_std_enable_shared_from_this.hpp.
using AtmosphereDiagnosticFactory =
    ekat::Factory<AtmosphereDiagnostic,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AtmosphereDiagnostic>,
                  const ekat::Comm&,const ekat::ParameterList&>;

// Create an atmosphere process, and correctly set up the (weak) pointer to self.
template <typename AtmDiagType>
inline std::shared_ptr<AtmosphereDiagnostic>
create_atmosphere_diagnostic (const ekat::Comm& comm, const ekat::ParameterList& p) {
  return std::make_shared<AtmDiagType>(comm,p);
}
} //namespace scream

#endif // SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP
