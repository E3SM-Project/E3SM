#ifndef SCREAM_HOMME_DYNAMICS_HPP
#define SCREAM_HOMME_DYNAMICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/grid/remap/abstract_remapper.hpp"
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

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Register all fields in the proper field manager(s).
  // Note: field_mgrs[grid_name] is the FM on grid $grid_name
  void register_fields (const std::map<std::string,std::shared_ptr<FieldManager<Real>>>& field_mgrs) const;

  // Retrieves an internal field, given field name and grid name.
  const Field<Real>& get_internal_field (const std::string& name, const std::string& grid) const;

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif
  void homme_pre_process (const int dt);
  void homme_post_process ();

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif
  // Sets the scream-allocated views inside the Homme c++ structures
  void init_homme_views ();

  // Propagates initial conditions to homme
  void import_initial_conditions ();

  // Updates p_mid
  void update_pressure ();

  void initialize_impl ();

protected:
  void run_impl        (const int dt);
  void finalize_impl   ();

  // Dynamics updates the "tracers" group, and needs to do some extra checks on the group.
  void set_computed_group_impl (const FieldGroup<Real>& group);

  // Override the check computed fields impl so we can repair slightly negative tracer values.
  void check_computed_fields_impl ();

  // Computes total number of bytes needed for local variables
  int requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Creates an internal field, not to be shared with the AD's FieldManager
  void create_internal_field (const std::string& name,
                              const std::vector<FieldTag>& tags,
                              const std::vector<int>& dims,
                              const std::string& grid);

  // Some helper fields. WARNING: only one copy for each internal field!
  std::map<std::string,field_type>  m_internal_fields;

  // Remapper for inputs and outputs, plus a special one for initial conditions
  std::shared_ptr<AbstractRemapper<Real>>   m_p2d_remapper;
  std::shared_ptr<AbstractRemapper<Real>>   m_d2p_remapper;
  std::shared_ptr<AbstractRemapper<Real>>   m_ic_remapper_fwd;
  std::shared_ptr<AbstractRemapper<Real>>   m_ic_remapper_bwd;

  // The dynamics and reference grids
  std::shared_ptr<const AbstractGrid>  m_dyn_grid;
  std::shared_ptr<const AbstractGrid>  m_ref_grid;
};

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HPP
