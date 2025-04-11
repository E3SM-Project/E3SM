#ifndef SCREAM_IOP_FORCING_HPP
#define SCREAM_IOP_FORCING_HPP

#include "ekat/ekat_parameter_list.hpp"
#include "ekat/ekat_workspace.hpp"

#include "share/atm_process/atmosphere_process.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
#include "share/util/eamxx_column_ops.hpp"

#include "physics/share/physics_constants.hpp"

#include <string>

namespace scream
{
/*
 * The class responsible for running EAMxx with an intensive
 * observation period (IOP).
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
 * Currently the only use case is the doubly
 * periodic model (DP-SCREAM).
 */

class IOPForcing : public scream::AtmosphereProcess
{
  // Typedefs for process
  using KT           = ekat::KokkosTypes<DefaultDevice>;
  using ESU          = ekat::ExeSpaceUtils<KT::ExeSpace>;
  using Pack         = ekat::Pack<Real, SCREAM_PACK_SIZE>;
  using IntPack      = ekat::Pack<int, Pack::n>;
  using Mask         = ekat::Mask<Pack::n>;
  using WorkspaceMgr = ekat::WorkspaceManager<Pack, KT::Device>;
  using Workspace    = WorkspaceMgr::Workspace;

  using MemberType = KT::MemberType;
  template <typename T>
  using view_1d = KT::view_1d<T>;
  template <typename T>
  using view_2d = KT::view_2d<T>;
  template <typename T>
  using uview_1d = ekat::Unmanaged<view_1d<T>>;
  template <typename T>
  using uview_2d = ekat::Unmanaged<view_2d<T>>;

  using ColOps = ColumnOps<DefaultDevice, Real>;
  using C      = physics::Constants<Real>;



public:

  // Constructors
  IOPForcing (const ekat::Comm& comm, const ekat::ParameterList& params)
    : AtmosphereProcess(comm, params) {}

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "iop"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif

  void initialize_impl (const RunType run_type);

  // Compute effects of large scale subsidence on T, q, u, and v.
  KOKKOS_FUNCTION
  static void advance_iop_subsidence(const KT::MemberType& team,
                                     const int nlevs,
                                     const Real dt,
                                     const Real ps,
                                     const view_1d<const Pack>& pmid,
                                     const view_1d<const Pack>& pint,
                                     const view_1d<const Pack>& pdel,
                                     const view_1d<const Pack>& omega,
                                     const Workspace& workspace,
                                     const view_1d<Pack>& u,
                                     const view_1d<Pack>& v,
                                     const view_1d<Pack>& T,
                                     const view_2d<Pack>& Q);

  // Apply large scale forcing for temperature and water vapor as provided by the IOP file
  KOKKOS_FUNCTION
  static void advance_iop_forcing(const KT::MemberType& team,
                                  const int nlevs,
                                  const Real dt,
                                  const view_1d<const Pack>& divT,
                                  const view_1d<const Pack>& divq,
                                  const view_1d<Pack>& T,
                                  const view_1d<Pack>& qv);

  // Provide coriolis forcing to u and v winds, using large scale winds specified in IOP forcing file.
  KOKKOS_FUNCTION
  static void iop_apply_coriolis(const KT::MemberType& team,
                                 const int nlevs,
                                 const Real dt,
                                 const Real lat,
                                 const view_1d<const Pack>& u_ls,
                                 const view_1d<const Pack>& v_ls,
                                 const view_1d<Pack>& u,
                                 const view_1d<Pack>& v);
  
  void run_impl        (const double dt);

protected:

  void finalize_impl   () {}

  // Creates an helper field, not to be shared with the AD's FieldManager
  void create_helper_field (const std::string& name,
                            const FieldLayout& layout,
                            const std::string& grid_name,
                            const int          ps = 1);

  void set_computed_group_impl (const FieldGroup& group);

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Keep track of field dimensions and other scalar values
  // needed in IOP
  Int m_num_cols;
  Int m_num_levs;
  Int m_num_tracers;

  struct Buffer {
    Pack* wsm_data;
  };

  // Some helper fields.
  std::map<std::string,Field> m_helper_fields;

  // Struct which contains local variables
  Buffer m_buffer;

  // WSM for internal local variables
  WorkspaceMgr m_workspace_mgr;

  std::shared_ptr<const AbstractGrid> m_grid;
}; // class IOPForcing

} // namespace scream

#endif // SCREAM_IOP_FORCING_HPP
