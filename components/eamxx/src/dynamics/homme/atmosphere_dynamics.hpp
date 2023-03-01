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

  // Constructor(s) and Destructor
  HommeDynamics (const ekat::Comm& comm, const ekat::ParameterList& params);
  ~HommeDynamics ();

  // The type of the subcomponent (dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Dynamics; }

  // The name of the subcomponent
  std::string name () const { return "Dynamics"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif
  void homme_pre_process (const double dt);
  void homme_post_process (const double dt);

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
  void update_pressure (const std::shared_ptr<const AbstractGrid>& grid);

  // Copy initial states from n0 timelevel to other timelevels
  void copy_dyn_states_to_all_timelevels ();

  void initialize_impl (const RunType run_type);

  // fv_phys refers to the horizontal finite volume (FV) grid for column
  // parameterizations nested inside the horizontal element grid. The grid names
  // are "Physics PGN", where N in practice is 2. The name of each routine is
  // fv_phys_X, where X is the name of an existing HommeDynamics routine. If
  // fv_phys is not being used, each of these routines does an immediate exit,
  // so it's OK to always call the routine.
  //   For the finite volume (FV) physics grid, sometimes referred to as
  // "physgrid", we use the remapper that Homme provides. Store N in pgN; if the
  // grid is not FV, then this variable is set to -1.
  int m_phys_grid_pgN;
  void fv_phys_set_grids();
  void fv_phys_requested_buffer_size_in_bytes() const;
  void fv_phys_initialize_impl();
  void fv_phys_dyn_to_fv_phys(const bool restart = false);
  void fv_phys_pre_process();
  void fv_phys_post_process();
  // See [rrtmgp active gases] in atmosphere_dynamics_fv_phys.cpp.
  void fv_phys_rrtmgp_active_gases_init(const std::shared_ptr<const GridsManager>& gm);
  void fv_phys_rrtmgp_active_gases_remap();

  // Rayleigh friction functions
  void rayleigh_friction_init ();
  void rayleigh_friction_apply (const Real dt) const;

public:
  // Fast boolean function returning whether Physics PGN is being used.
  bool fv_phys_active() const;
  struct GllFvRemapTmp;
  void remap_dyn_to_fv_phys(GllFvRemapTmp* t = nullptr) const;
  void remap_fv_phys_to_dyn() const;
  
protected:
  void run_impl        (const double dt);
  void finalize_impl   ();

  // We need to store the size of the tracers group as soon as it is available.
  // In particular, we absolutely need to store it *before* the call to
  // requested_buffer_size_in_bytes, where Homme::ForcingFunctor queries
  // Homme::Tracers for qsize.
  void set_computed_group_impl (const FieldGroup& group);

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Creates an helper field, not to be shared with the AD's FieldManager
  void create_helper_field (const std::string& name,
                            const std::vector<FieldTag>& tags,
                            const std::vector<int>& dims,
                            const std::string& grid);

  // Some helper fields.
  std::map<std::string,Field>  m_helper_fields;

  // Remapper for inputs and outputs, plus a special one for initial
  // conditions. These are used when the physics grid is the continuous GLL
  // point grid.
  std::shared_ptr<AbstractRemapper>   m_p2d_remapper;
  std::shared_ptr<AbstractRemapper>   m_d2p_remapper;
  std::shared_ptr<AbstractRemapper>   m_ic_remapper;

  // The dynamics and reference grids
  std::shared_ptr<const AbstractGrid> m_dyn_grid;  // Dynamics DGLL
  std::shared_ptr<const AbstractGrid> m_phys_grid; // Column parameterizations grid
  std::shared_ptr<const AbstractGrid> m_cgll_grid; // Unique CGLL

  template<int N>
  using RPack = ekat::Pack<Real,N>;

  using KT = KokkosTypes<DefaultDevice>;
  template<typename ScalarT>
  using view_1d = typename KT::template view_1d<ScalarT>;

  // Rayleigh friction decay rate profile
  view_1d<RPack<SCREAM_PACK_SIZE>> m_otau;

  // Rayleigh friction paramaters
  int m_rayk0;      // Vertical level at which rayleigh friction term is centered.
  Real m_raykrange; // Range of rayleigh friction profile.
  Real m_raytau0;   // Approximate value of decay time at model top (days)
                    // if set to 0, no rayleigh friction is applied
};

} // namespace scream

#endif // SCREAM_HOMME_DYNAMICS_HPP
