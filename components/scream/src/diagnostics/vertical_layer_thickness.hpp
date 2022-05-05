#ifndef EAMXX_VERTICAL_LAY_THICK_DIAGNOSTIC_HPP
#define EAMXX_VERTICAL_LAY_THICK_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class VerticalLayerThicknessDiagnostic : public AtmosphereDiagnostic
{
public:
  template <typename S>
  using SmallPack     = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;

  using Spack         = SmallPack<Real>;
  using Pack          = ekat::Pack<Real,Spack::n>;
  using PF            = scream::PhysicsFunctions<DefaultDevice>;
  using KT            = KokkosTypes<DefaultDevice>;
  using MemberType    = typename KT::MemberType;
  using view_1d       = typename KT::template view_1d<Spack>;
  using view_2d       = typename KT::template view_2d<Spack>;
  using view_2d_const = typename KT::template view_2d<const Spack>;

  template<typename ScalarT>
  using uview_1d = Unmanaged<typename KT::template view_1d<ScalarT>>;
  template<typename ScalarT>
  using uview_2d = Unmanaged<typename KT::template view_2d<ScalarT>>;

  // Constructors
  VerticalLayerThicknessDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Set type to diagnostic
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // The name of the diagnostic
  std::string name () const { return "Potential Temperature"; } 

  // Get the required grid for the diagnostic
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Actual diagnostic calculation 
  struct run_diagnostic_impl {
    run_diagnostic_impl() = default;
    // Functor for Kokkos loop
    KOKKOS_INLINE_FUNCTION
    void operator() (const MemberType& team) const {
      const int icol = team.league_rank();
    }
    int m_ncol, m_nlevs;
    view_2d_const        T_mid;
    view_2d_const        p_mid;
    view_2d_const        pseudo_density;
    view_2d_const        qv;
    view_2d              dz;
    // assign variables to this structure
    void set_variables(const int ncol, const int nlevs,
                       const view_2d_const& T_mid_, const view_2d_const& p_mid_, const view_2d_const& pseudo_density_,
                       const view_2d_const& qv_,
                       const view_2d& dz_
      )
    {
      m_ncol   = ncol;
      m_nlevs  = nlevs;
      // IN
      T_mid = T_mid_;
      p_mid = p_mid_;
      pseudo_density = pseudo_density_;
      qv = qv_;
      // OUT
      dz = dz_;
    }
  }; // struct run_diagnostic_impl

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const RunType run_type);
  void run_impl        (const int dt);
  void finalize_impl   ();

  // Keep track of field dimensions
  Int m_num_cols; 
  Int m_num_levs;

  // Structure to run the diagnostic
  run_diagnostic_impl  run_diagnostic;

}; // class VerticalLayerThicknessDiagnostic

} //namespace scream

#endif // EAMXX_VERTICAL_LAY_THICK_DIAGNOSTIC_HPP
