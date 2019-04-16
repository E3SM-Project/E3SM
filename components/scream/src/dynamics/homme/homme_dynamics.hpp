#ifndef SCREAM_HOMME_DYNAMICS_HPP
#define SCREAM_HOMME_DYNAMICS_HPP

#include "share/atmosphere_process.hpp"
#include "share/parameter_list.hpp"

#include <string>

namespace scream
{

/*
 *  The class responsible to handle the atmosphere dynamics
 *
 *  The AD should store exactly ONE instance of this class stored
 *  in its list of subcomponents (the AD should make sure of this).
 *
 *  Note: for now, Scream is only going to accommodate HOMME as a dynamics
 *  dycore.
 */

class HommeDynamics : public AtmosphereProcess
{
public:

  // Constructor(s)
  explicit HommeDynamics (const ParameterList& params);

  // The type of the subcomponent (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Dynamics; }

  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(e2str(GridType::Dynamics));
    return s;
  }

  // The name of the subcomponent
  std::string name () const { return "dynamics"; }

  // The communicator used by the dynamics
  const Comm& get_comm () const { return m_dynamics_comm; }

  // These are the three main interfaces:
  void initialize (const Comm& comm, const std::shared_ptr<const GridsManager> grids_manager);
  void run        (/* what inputs? */);
  void finalize   (/* what inputs? */);

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const;

  // Get the set of required/computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_computed_fields; }

protected:

  // Setting the field in the atmosphere process
  void set_required_field_impl (const Field<const Real, device_type>& /*f*/) { /* impl */ }
  void set_computed_field_impl (const Field<      Real, device_type>& /*f*/) { /* impl */ }

  std::set<FieldIdentifier> m_required_fields;
  std::set<FieldIdentifier> m_computed_fields;

  Comm      m_dynamics_comm;

};

inline AtmosphereProcess* create_homme_dynamics(const ParameterList& p) {
  return new HommeDynamics(p);
}

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HPP
