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
  using field_type       = Field<      Real>;
  using const_field_type = Field<const Real>;

  using SPAFunc         = spa::SPAFunctions<Real, DefaultDevice>;
  using Spack           = SPAFunc::Spack;
  using Pack            = ekat::Pack<Real,Spack::n>;
  using KT              = ekat::KokkosTypes<DefaultDevice>;
  using WSM             = ekat::WorkspaceManager<Spack, KT::Device>;
  using LIV             = ekat::LinInterp<Real,Spack::n>;

  using view_1d         = typename SPAFunc::view_1d<Spack>;
  using view_2d         = typename SPAFunc::view_2d<Spack>;
  using uview_2d        = Unmanaged<view_2d>;

  // Constructors
  SPA (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Simple Prescribed Aerosols (SPA)"; }

  // The communicator used by subcomponent
  const ekat::Comm& get_comm () const { return m_spa_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_spa_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    // 1d view scalar, size (ncol)
    static constexpr int num_1d_scalar = 0;
    // 2d view packed, size (ncol, nlev_packs)
    static constexpr int num_2d_vector = 2;
    static constexpr int num_2dp1_vector = 0;

    uview_2d p_mid_src;
    uview_2d ccn_src;

    Spack* wsm_data;
  };
protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const util::TimeStamp& t0);
  void run_impl        (const Real dt);
  void finalize_impl   ();

  // Setting the fields in the atmospheric process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

  // Computes total number of bytes needed for local variables
  int requested_buffer_size_in_bytes() const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager &buffer_manager);

  std::map<std::string,const_field_type>  m_spa_fields_in;
  std::map<std::string,field_type>        m_spa_fields_out;

  ekat::Comm          m_spa_comm;
  ekat::ParameterList m_spa_params;

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_levs;
  Int m_nk_pack;

  // Struct which contains local variables
  Buffer m_buffer;

}; // class SPA 

} // namespace scream

#endif // SCREAM_SPA_HPP
