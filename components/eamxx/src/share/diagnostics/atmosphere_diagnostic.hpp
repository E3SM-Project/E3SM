#ifndef SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP
#define SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP

#include "share/field/field.hpp"
#include "share/field/field_identifier.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/data_managers/field_request.hpp"
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
 * A class representing an atmosphere diagnostic.
 *
 * A diagnostic takes one or more fields as input and produces a single
 * output field. Unlike AtmosphereProcess, this class does not participate
 * in the atmosphere time-stepping loop; it is only used for output/IO.
 *
 * This class does NOT inherit from AtmosphereProcess. It is a standalone
 * base class that depends only on eamxx_core, eamxx_utils, eamxx_field,
 * and eamxx_grid libraries.
 */

class AtmosphereDiagnostic
{
public:

  // Constructor(s)
  AtmosphereDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  virtual ~AtmosphereDiagnostic () = default;

  // The name of the diagnostic
  virtual std::string name () const = 0;

  // Set the grid used by this diagnostic. This also calls create_requests(),
  // which derived classes must implement to declare their required fields.
  void set_grid (const std::shared_ptr<const AbstractGrid>& grid);

  // Derived classes must override this to declare their required fields
  // (via add_field<Required>(...)) and set up m_diagnostic_output.
  virtual void create_requests () = 0;

  // Initialize the diagnostic. Calls initialize_impl().
  void initialize ();

  // Returns true after initialize() has been called successfully.
  bool is_initialized () const { return m_is_initialized; }

  // Return the parameter list
  ekat::ParameterList& get_params () { return m_params; }

  // Return the list of field requests declared by create_requests()
  const std::list<FieldRequest>& get_field_requests () const { return m_field_requests; }

  // Set a required input field (called after create_requests(), before initialize())
  void set_required_field (const Field& f);

  // Return all input fields
  const std::list<Field>& get_fields_in () const { return m_fields_in; }

  // Retrieve an input field by name (after set_required_field has been called)
  const Field& get_field_in (const std::string& field_name) const;
        Field& get_field_in (const std::string& field_name);

  // Allows the diagnostic to save some start-of-step quantity (e.g., in case
  // we need to compute tendencies, or accumulated stuff)
  virtual void init_timestep (const util::TimeStamp& /* start_of_step */) {}

  // Returns the diagnostic output field
  Field get_diagnostic () const;

  // Compute the diagnostic (skips if inputs have not changed since last call)
  void compute_diagnostic (const double dt = 0);

protected:

  // Derived classes override this for any setup needed after fields are set.
  // Called by initialize().
  virtual void initialize_impl () { /* Nothing to do */ }

  // Derived classes implement this to compute the output from the inputs.
  virtual void compute_diagnostic_impl () = 0;

  // Helper template methods to declare required fields in create_requests().
  // All overloads ultimately create a FieldRequest and append it to m_field_requests.

  // With full layout+units+grid (optional pack size)
  template<RequestType RT>
  void add_field (const std::string& name, const FieldLayout& layout,
                  const ekat::units::Units& u, const std::string& grid_name,
                  const int ps = 1)
  {
    static_assert(RT==Required, "Error! Diagnostics can only add 'Required' fields.\n");
    m_field_requests.emplace_back(FieldIdentifier(name,layout,u,grid_name),ps);
  }

  // With group name
  template<RequestType RT>
  void add_field (const std::string& name, const FieldLayout& layout,
                  const ekat::units::Units& u, const std::string& grid_name,
                  const std::string& group, const int ps = 1)
  {
    static_assert(RT==Required, "Error! Diagnostics can only add 'Required' fields.\n");
    m_field_requests.emplace_back(FieldIdentifier(name,layout,u,grid_name),group,ps);
  }

  // With a pre-built FieldIdentifier
  template<RequestType RT>
  void add_field (const FieldIdentifier& fid, const int ps = 1)
  {
    static_assert(RT==Required, "Error! Diagnostics can only add 'Required' fields.\n");
    m_field_requests.emplace_back(fid,ps);
  }

  // With just name + grid_name (incomplete request, resolved when the field is set)
  template<RequestType RT>
  void add_field (const std::string& name, const std::string& grid_name,
                  const int ps = 1)
  {
    static_assert(RT==Required, "Error! Diagnostics can only add 'Required' fields.\n");
    m_field_requests.emplace_back(name,grid_name,std::list<std::string>{},ps);
  }

  // MPI communicator
  ekat::Comm              m_comm;

  // Parameter list
  ekat::ParameterList     m_params;

  // The grid on which this diagnostic is defined (set via set_grid())
  std::shared_ptr<const AbstractGrid> m_grid;

  // Some diagnostics will need the timestep, store here
  double m_dt;

  // Diagnostics produce a single output field
  Field m_diagnostic_output;

  // Timestamp of the last diag evaluation
  util::TimeStamp m_last_eval_ts;

private:

  bool m_is_initialized = false;

  // Required input fields
  std::list<FieldRequest> m_field_requests;
  std::list<Field>        m_fields_in;

  // Map from field name to pointer into m_fields_in (for fast lookup)
  std::map<std::string, Field*> m_fields_in_pointers;
};

// Factory for atmosphere diagnostics
using AtmosphereDiagnosticFactory =
    ekat::Factory<AtmosphereDiagnostic,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AtmosphereDiagnostic>,
                  const ekat::Comm&,const ekat::ParameterList&>;

// Convenience function to create a diagnostic and register it in the factory
template <typename AtmDiagType>
inline std::shared_ptr<AtmosphereDiagnostic>
create_atmosphere_diagnostic (const ekat::Comm& comm, const ekat::ParameterList& p) {
  return std::make_shared<AtmDiagType>(comm,p);
}

} //namespace scream

#endif // SCREAM_ATMOSPHERE_DIAGNOSTIC_HPP
