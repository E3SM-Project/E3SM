#ifndef SCREAM_SURFACE_COUPLING_HPP
#define SCREAM_SURFACE_COUPLING_HPP

#include "share/atmosphere_process.hpp"
#include "share/field_repository.hpp"

namespace scream {

// This class is responsible to import/export fields from/to the rest of
// E3SM, via the mct coupler.
class SurfaceCoupling : public AtmosphereProcess {
public:
  template<typename VT>
  using field_repo_type = FieldRepository<VT>;

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Coupling; }

  std::string name () const { return "surface_coupling"; }

  // The communicator associated with this atm process
  const Comm& get_comm () const { return m_surface_coupling_comm; }

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize ( /* inputs? */ );

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run ( /* inputs? */ );

  // Clean up
  void finalize ( /* inputs */ );

  // Providing a list of required and computed fields
  const std::list<std::shared_ptr<FieldHeader>>&  get_required_fields () const { return m_fields_to_import; }
  const std::list<std::shared_ptr<FieldHeader>>&  get_computed_fields () const { return m_fields_to_export; }

  // Setting the field in the atmosphere process
  void set_required_field (const Field<const Real*, ExecMemSpace, true>& /*f*/) { /* impl */ }
  void set_computed_field (const Field<      Real*, ExecMemSpace, true>& /*f*/) { /* impl */ }

protected:

  std::list<std::shared_ptr<FieldHeader>> m_fields_to_export;
  std::list<std::shared_ptr<FieldHeader>> m_fields_to_import;

  field_repo_type<HostMemSpace>   m_host_field_repo;

  Comm    m_surface_coupling_comm;
};

} // namespace scream

#endif // SCREAM_SURFACE_COUPLING_HPP
