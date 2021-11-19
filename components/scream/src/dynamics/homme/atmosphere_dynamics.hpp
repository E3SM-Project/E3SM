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
  template<int N>
  using view_Nd_type = typename field_type::view_ND_type<Real,N>;

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
  void initialize_homme_state ();
  // Restart homme
  void restart_homme_state ();

  // Updates p_mid
  void update_pressure ();

  // Copy initial states from n0 timelevel to other timelevels
  void copy_dyn_states_to_all_timelevels ();

  void initialize_impl (const RunType run_type);
protected:
  void run_impl        (const int dt);
  void finalize_impl   ();

  // For simplicity, it's best to store the size of the tracers as soon as it is available.
  // We can do it the first time that the 'tracers' group is set
  void set_computed_group_impl (const FieldGroup<Real>& group);

  // Override the check computed fields impl so we can repair slightly negative tracer values.
  void check_computed_fields_impl ();

  // Computes total number of bytes needed for local variables
  int requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Creates an helper field, not to be shared with the AD's FieldManager
  void create_helper_field (const std::string& name,
                            const std::vector<FieldTag>& tags,
                            const std::vector<int>& dims,
                            const std::string& grid);

  // Some helper fields.
  std::map<std::string,field_type>  m_helper_fields;

  // Remapper for inputs and outputs, plus a special one for initial conditions
  std::shared_ptr<AbstractRemapper<Real>>   m_p2d_remapper;
  std::shared_ptr<AbstractRemapper<Real>>   m_d2p_remapper;
  std::shared_ptr<AbstractRemapper<Real>>   m_ic_remapper;

  // The dynamics and reference grids
  std::shared_ptr<const AbstractGrid>  m_dyn_grid;
  std::shared_ptr<const AbstractGrid>  m_ref_grid;
};

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HPP
