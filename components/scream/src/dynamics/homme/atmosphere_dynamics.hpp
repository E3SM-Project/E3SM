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

  // The communicator used by the dynamics
  const ekat::Comm& get_comm () const { return m_dynamics_comm; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Register all fields in the proper field manager(s).
  // Note: field_mgrs[grid_name] is the FM on grid $grid_name
  void register_fields (const std::map<std::string,std::shared_ptr<FieldManager<Real>>>& field_mgrs) const;

  // Dynamics updates 'TRACERS'.
  void set_updated_group (const FieldGroup<Real>& group);

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif
  void homme_pre_process (const Real dt);
  void homme_post_process ();
  // These are the three main interfaces:

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif
  void initialize_impl (const util::TimeStamp& t0);
protected:
  void run_impl        (const Real dt);
  void finalize_impl   ();

  // Setting the fields in the atmosphere process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

  void create_dyn_field (const std::string& name,
                         const std::vector<FieldTag>& tags,
                         const std::vector<int>& dims);

  // Computes total number of bytes needed for local variables
  int requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Fields on reference and dynamics grid
  // NOTE: the dyn grid fields are *NOT* in the FieldManager. We still use
  //       scream Field's (rather than, e.g., raw views) cause we want to use
  //       the remapper infrastructure to remap from/to ref grid to/from dyn grid.
  std::map<std::string,field_type>  m_ref_grid_fields;
  std::map<std::string,field_type>  m_dyn_grid_fields;

  // Remapper for inputs and outputs, plus a special one for initial conditions
  std::shared_ptr<AbstractRemapper<Real>>   m_p2d_remapper;
  std::shared_ptr<AbstractRemapper<Real>>   m_d2p_remapper;
  std::shared_ptr<AbstractRemapper<Real>>   m_ic_remapper;

  // The dynamics and reference grids
  std::shared_ptr<const AbstractGrid>  m_dyn_grid;
  std::shared_ptr<const AbstractGrid>  m_ref_grid;

  // Homme dyn parameters
  ekat::ParameterList     m_params;

  // The MPI communicator associated witht his atm process
  ekat::Comm              m_dynamics_comm;
};

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HPP
