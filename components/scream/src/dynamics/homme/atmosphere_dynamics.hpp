#ifndef SCREAM_HOMME_DYNAMICS_HPP
#define SCREAM_HOMME_DYNAMICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/scream_parameter_list.hpp"

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
  using field_type       = Field<      Real,device_type>;
  using const_field_type = Field<const Real,device_type>;

  // Constructor(s)
  HommeDynamics (const Comm& comm, const ParameterList& params);

  // The type of the subcomponent (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Dynamics; }

  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert("SE Dynamics");
    return s;
  }

  // The name of the subcomponent
  std::string name () const { return "Dynamics"; }

  // The communicator used by the dynamics
  const Comm& get_comm () const { return m_dynamics_comm; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // These are the three main interfaces:
  void initialize (const util::TimeStamp& t0);
  void run        (const Real dt);
  void finalize   (/* what inputs? */);

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const;

  // Get the set of required/computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_computed_fields; }

protected:

  // Setting the fields in the atmosphere process
  void set_required_field_impl (const Field<const Real, device_type>& f);
  void set_computed_field_impl (const Field<      Real, device_type>& f);

  std::set<FieldIdentifier> m_required_fields;
  std::set<FieldIdentifier> m_computed_fields;

  std::map<std::string,const_field_type>  m_dyn_fields_in;
  std::map<std::string,field_type>        m_dyn_fields_out;

  util::TimeStamp   m_current_ts;
  Comm              m_dynamics_comm;
};

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HPP
