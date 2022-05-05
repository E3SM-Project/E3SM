#ifndef EAMXX_VERTICAL_LAY_INT_DIAGNOSTIC_HPP
#define EAMXX_VERTICAL_LAY_INT_DIAGNOSTIC_HPP

#include "share/atm_process/atmosphere_diagnostic.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

namespace scream
{

/*
 * This diagnostic will produce the potential temperature.
 */

class VerticalLayerInterfaceDiagnostic : public AtmosphereDiagnostic
{
public:
  template <typename S>
  using SmallPack     = ekat::Pack<S,SCREAM_SMALL_PACK_SIZE>;
  using IntSmallPack  = SmallPack<Int>;

  using Spack         = SmallPack<Real>;
  using Smask         = ekat::Mask<Spack::n>;
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
  VerticalLayerInterfaceDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params);

  // Set type to diagnostic
  AtmosphereProcessType type () const { return AtmosphereProcessType::Diagnostic; }

  // The name of the diagnostic
  std::string name () const { return "Vertical Layer Interface"; } 

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
      const int nlev_packs = ekat::npack<Spack>(m_nlevs);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
        const auto range = ekat::range<IntSmallPack>(k*Spack::n);
        const Smask in_nlev_range = (range < m_nlevs);

        dz(icol,k) = PF::calculate_dz(pseudo_density(icol,k), p_mid(icol,k), T_mid(icol,k), qv(icol,k));
      });
      team.team_barrier();
      const auto& dz_s    = ekat::subview(dz,    icol);
      const auto& z_int_s = ekat::subview(z_int, icol);
      PF::calculate_z_int(team,m_nlevs,dz_s,surf_geopotential,z_int_s);
    }
    Real surf_geopotential;
    int m_ncol, m_nlevs;
    view_2d_const        T_mid;
    view_2d_const        p_mid;
    view_2d_const        pseudo_density;
    view_2d_const        qv;
    view_2d              dz;
    view_2d              z_int;
    // assign variables to this structure
    void set_variables(const Real surf_geo, const int ncol, const int nlevs,
                       const view_2d_const& T_mid_, const view_2d_const& p_mid_, const view_2d_const& pseudo_density_,
                       const view_2d_const& qv_,
                       const view_2d& z_int_
      )
    {
      surf_geopotential = surf_geo;
      m_ncol   = ncol;
      m_nlevs  = nlevs;
      const int nlev_packs = ekat::npack<Spack>(nlevs);
      dz = view_2d("",ncol,nlev_packs);
      // IN
      T_mid = T_mid_;
      p_mid = p_mid_;
      pseudo_density = pseudo_density_;
      qv = qv_;
      // OUT
      z_int = z_int_;
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

}; // class VerticalLayerInterfaceDiagnostic

} //namespace scream

#endif // EAMXX_VERTICAL_LAY_INT_DIAGNOSTIC_HPP
