#ifndef SCREAM_HOMME_DYNAMICS_HPP
#define SCREAM_HOMME_DYNAMICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"

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
  using field_type       = Field<      Real>;
  using const_field_type = Field<const Real>;

  // Constructor(s)
  HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of the subcomponent (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Dynamics; }

  std::set<std::string> get_required_grids () const {
    return std::set<std::string>{"Dynamics"};
  }

  // The name of the subcomponent
  std::string name () const { return "Dynamics"; }

  // The communicator used by the dynamics
  const ekat::Comm& get_comm () const { return m_dynamics_comm; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real>& field_repo) const;

  // Get the set of required/computed fields
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_computed_fields; }

protected:

  // These are the three main interfaces:
  void initialize_impl (const util::TimeStamp& t0);
  void run_impl        (const Real dt);
  void finalize_impl   (/* what inputs? */);

  // Setting the fields in the atmosphere process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

  std::set<FieldIdentifier> m_required_fields;
  std::set<FieldIdentifier> m_computed_fields;

  std::map<std::string,const_field_type>  m_dyn_fields_in;
  std::map<std::string,field_type>        m_dyn_fields_out;

  // For certain tests, dynamics needs to init states
  std::shared_ptr<FieldInitializer>       m_initializer;

  ekat::ParameterList     m_params;
  ekat::Comm              m_dynamics_comm;
};

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HPP
