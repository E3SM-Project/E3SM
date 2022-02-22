#ifndef SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP
#define SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_process.hpp"

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
 * 2) Diagnostics produce a single field output.
 *    Typically the field output will be to buffered memory, meaning it
 *    could be overwritten later in the simulation.
 *    A diagnostic output is meant to be used for OUTPUT or a PROPERTY CHECK.
 */

class AtmosphereDiagnostic : public AtmosphereProcess
{
public:

  // Constructor(s)
  explicit AtmosphereDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  virtual ~AtmosphereDiagnostic () = default;

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // Getting the diagnostic output
  Field get_diagnostic (const Real dt);

  void set_computed_field (const Field& f) final;
  void set_computed_group (const FieldGroup& group) final;
protected:

  // Diagnostics are meant to return a field
  Field m_diagnostic_output;


};

} //namespace scream

#endif // SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP
