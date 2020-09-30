#ifndef SCREAM_SURFACE_COUPLING_HPP
#define SCREAM_SURFACE_COUPLING_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "share/field/field_repository.hpp"

namespace scream {

// This class is responsible to import/export fields from/to the rest of
// E3SM, via the mct coupler.
class SurfaceCoupling : public AtmosphereProcess {
public:
  using host_device_type = HostDevice;

  template<typename MS>
  using field_repo_type = FieldRepository<Real,MS>;

  explicit SurfaceCoupling (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of the block (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Coupling; }

  // The type of grids required by the process
  std::set<std::string> get_required_grids () const {
    // TODO: define what grid the coupling runs on. Check with MOAB folks.
    static std::set<std::string> s;
    s.insert(e2str(GridType::Undefined));
    return s;
  }

  std::string name () const { return "surface_coupling"; }

  // The communicator associated with this atm process
  const ekat::Comm& get_comm () const { return m_comm; }

  // Get the grid from the grids manager
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize (const util::TimeStamp& t0);

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run (const Real dt);

  // Clean up
  void finalize ( /* inputs */ );

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const;

  // Providing a list of required and computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_fields_to_import; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_fields_to_export; }

protected:

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real, device_type>& f);
  void set_computed_field_impl (const Field<      Real, device_type>& f);

  std::set<FieldIdentifier> m_fields_to_export;
  std::set<FieldIdentifier> m_fields_to_import;

  field_repo_type<host_device_type> m_host_field_repo;
  field_repo_type<device_type>      m_device_field_repo;

  ekat::Comm    m_comm;
};

} // namespace scream

#endif // SCREAM_SURFACE_COUPLING_HPP
