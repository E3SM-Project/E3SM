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

  // The initialization method should prepare all stuff needed to import/export from/to
  // f90 structures.
  void initialize ( /* inputs? */ );

  // The run method is responsible for exporting atm states to the e3sm coupler, and
  // import surface states from the e3sm coupler.
  void run ( /* inputs? */ );

  // Clean up
  void finalize ( /* inputs */ );

protected:

  field_repo_type<ExecViewManaged<Real*>>   m_device_field_repo;
  field_repo_type<HostViewManaged<Real*>>   m_host_field_repo;
};

} // namespace scream

#endif // SCREAM_SURFACE_COUPLING_HPP
