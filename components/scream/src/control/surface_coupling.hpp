#ifndef SCREAM_SURFACE_COUPLING_HPP
#define SCREAM_SURFACE_COUPLING_HPP

#include <share/atmosphere_process.hpp>
#include <share/parameter_list.hpp>

#include <share/field/field_repository.hpp>

namespace scream {

// This class is responsible to import/export fields from/to the rest of
// E3SM, via the mct coupler.
class SurfaceCoupling : public AtmosphereProcess {
public:
  template<typename MS>
  using field_repo_type = FieldRepository<Real,MS>;

  SurfaceCoupling (const ParameterList& params);

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Coupling; }

  std::string name () const { return "surface_coupling"; }

  // The communicator associated with this atm process
  const Comm& get_comm () const { return m_comm; }

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize (const Comm& comm);

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run ( /* inputs? */ );

  // Clean up
  void finalize ( /* inputs */ );

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, ExecMemSpace>& field_repo) const;

  // Providing a list of required and computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_fields_to_import; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_fields_to_export; }

protected:

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real, ExecMemSpace, MemoryManaged>& f);
  void set_computed_field_impl (const Field<      Real, ExecMemSpace, MemoryManaged>& f);

  std::set<FieldIdentifier> m_fields_to_export;
  std::set<FieldIdentifier> m_fields_to_import;

  field_repo_type<HostMemSpace>   m_host_field_repo;
  field_repo_type<ExecMemSpace>   m_device_field_repo;

  Comm    m_comm;
};

} // namespace scream

#endif // SCREAM_SURFACE_COUPLING_HPP
