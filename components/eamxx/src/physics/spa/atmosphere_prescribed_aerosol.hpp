#ifndef SCREAM_PRESCRIBED_AEROSOL_HPP
#define SCREAM_PRESCRIBED_AEROSOL_HPP

#include "physics/spa/spa_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_lin_interp.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the calculation of the subgrid cloud fractions
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

class SPA : public AtmosphereProcess
{
public:
  using gid_type         = AbstractGrid::gid_type;

  using SPAFunc         = spa::SPAFunctions<Real, DefaultDevice>;
  using Spack           = SPAFunc::Spack;
  using Pack            = ekat::Pack<Real,Spack::n>;
  using KT              = ekat::KokkosTypes<DefaultDevice>;

  using view_1d         = typename SPAFunc::view_1d<Spack>;
  using view_2d         = typename SPAFunc::view_2d<Spack>;
  using view_3d         = typename SPAFunc::view_3d<Spack>;
  using view_1d_dof     = typename SPAFunc::view_1d<const gid_type>;

  template<typename ScalarT>
  using uview_1d = Unmanaged<typename KT::template view_1d<ScalarT>>;
  template<typename ScalarT>
  using uview_2d = Unmanaged<typename KT::template view_2d<ScalarT>>;

  // Constructors
  SPA (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Simple Prescribed Aerosols (SPA)"; }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    // Used to store temporary data during spa_main
    SPAFunc::SPAInput spa_temp;

    // Temporary to use
    uview_2d<Spack> p_mid_src;
  };
protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void initialize_spa_impl ();
  void run_impl        (const double dt);
  void finalize_impl   ();

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  // Keep track of field dimensions and the iteration count
  int m_num_cols; 
  int m_num_levs;
  int m_num_src_levs;
  int m_nk_pack;
  int m_nswbands = 14;
  int m_nlwbands = 16;

  // DOF information
  view_1d_dof m_dofs_gids;
  int         m_total_global_dofs; // Needed to make sure that remap data matches grid.
  gid_type    m_min_global_dof;

  // Struct which contains temporary variables used during spa_main
  Buffer m_buffer;

  // SPA specific files
  std::string m_spa_remap_file;
  std::string m_spa_data_file;

  // Structures to store the data used for interpolation
  SPAFunc::SPATimeState     SPATimeState;
  SPAFunc::SPAHorizInterp   SPAHorizInterp;
  SPAFunc::SPAInput         SPAData_start;
  SPAFunc::SPAInput         SPAData_end;
  SPAFunc::SPAOutput        SPAData_out;

  std::shared_ptr<const AbstractGrid>   m_grid;
}; // class SPA 

} // namespace scream

#endif // SCREAM_SPA_HPP
