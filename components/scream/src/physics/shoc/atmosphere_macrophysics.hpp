#ifndef SCREAM_SHOC_MACROPHYSICS_HPP
#define SCREAM_SHOC_MACROPHYSICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "physics/shoc/shoc_main_impl.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "physics/share/physics_functions.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the atmosphere microphysics
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
 *  Note: for now, scream is only going to accommodate SHOC as macrophysics
*/

class SHOCMacrophysics : public scream::AtmosphereProcess
{
  using SHF          = shoc::Functions<Real, DefaultDevice>;
  using physics_fun  = scream::physics::Functions<Real, DefaultDevice>;
  using C            = physics::Constants<Real>;
  using KT           = ekat::KokkosTypes<DefaultDevice>;

  using Spack = typename SHF::Spack;
  using IntSmallPack = typename SHF::IntSmallPack;
  using Pack1d = typename SHF::Pack1d;
  using Smask = typename SHF::Smask;
  using view_1d  = typename SHF::view_1d<Pack1d>;
  using view_1d_const  = typename SHF::view_1d<const Pack1d>;
  using view_2d  = typename SHF::view_2d<SHF::Spack>;
  using view_2d_const  = typename SHF::view_2d<const Spack>;
  using sview_2d = typename KokkosTypes<DefaultDevice>::template view_2d<Real>;
  using view_3d = typename SHF::view_3d<Spack>;
  using view_3d_const = typename SHF::view_3d<const Spack>;

public:
  using field_type       = Field<      Real>;
  using const_field_type = Field<const Real>;

  // Constructors
  SHOCMacrophysics (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Macrophysics"; }

  // The communicator used by subcomponent
  const ekat::Comm& get_comm () const { return m_shoc_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_shoc_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Register all fields in the given repo
  void register_fields (const std::map<std::string,std::shared_ptr<FieldRepository<Real>>>& field_repos) const;

  // SHOC updates the 'TRACERS' group.
  void set_updated_group (const FieldGroup<Real>& group);

  /*--------------------------------------------------------------------------------------------*/
  // Most individual processes have a pre-processing step that constructs needed variables from
  // the set of fields stored in the field manager.  A structure like this defines those operations,
  // which can then be called during run_impl in the main .cpp code.
  // Structure to handle the local generation of data needed by shoc_main in run_impl
  struct SHOCPreprocess {
    SHOCPreprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank();

      const Real zvir = C::ZVIR;
      const Real latvap = C::LatVap;
      const Real cpair = C::Cpair;
      const Real ggr = C::gravit;

      // Transpose and repack tracer array
      // TODO: remove once shoc tracer traspose is implemented
      for (int k=0; k<nlev; ++k) {
        for (int q=0; q<num_tracers; ++q) {
          const int q_v = q/Spack::n;
          const int q_p = q%Spack::n;
          const int k_v = k/Spack::n;
          const int k_p = k%Spack::n;

          tracers(i,k,q_v)[q_p] = Q(i,q,k_v)[k_p];
        }
      }

      const int nlevi_v = nlev/Spack::n;
      const int nlevi_p = nlev%Spack::n;

      const auto sub_z_int = ekat::subview(z_int, i);
      const auto s_z_int = ekat::scalarize(sub_z_int);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
        // Exner
        const Spack p_mid_ik(p_mid(i,k));
        const Smask p_mid_mask(!isnan(p_mid_ik) and p_mid_ik>0.0);
        auto exner_ik = physics_fun::get_exner(p_mid_ik,p_mid_mask);
        exner(i,k) = exner_ik;

        tke(i,k) = ekat::max(sp(0.004), tke(i,k));

        // Tracers are updated as a group. These tracers act as seperate inputs to shoc_main
        // and should have different output as tracer group. We make a copy and pass these to
        // shoc_main, then copy back to tracer group in postprocessing.
        // TODO: remove *_copy views once SHOC can request a subset of tracers.
        tke_copy(i,k) = tke(i,k);
        qc_copy(i,k) = qc(i,k);
        qv_copy(i,k) = qv(i,k);

        // Temperature
        qw(i,k) = qv(i,k) + qc(i,k);

        const auto theta_zt = T_mid(i,k)/exner(i,k);
        thlm(i,k) = theta_zt - (latvap/cpair)*qc(i,k);
        thv(i,k)  = theta_zt*(1 + zvir*qv(i,k) - qc(i,k));

        // Dry static energy
        const Spack T_mid_ik(T_mid(i,k));
        const Spack z_mid_ik(z_mid(i,k));
        const Real  phis_i(phis(i)[0]);
        const Smask range_mask(!isnan(T_mid_ik) && T_mid_ik>0.0 &&
                               !isnan(z_mid_ik) && z_mid_ik>0.0 &&
                               !isnan(phis_i)   && phis_i>0.0);
        auto dse_ik = physics_fun::get_dse(T_mid_ik,z_mid_ik,phis_i,range_mask);
        shoc_s(i,k) = dse_ik;

        Spack zi_k, zi_kp1;
        auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
        auto range_pack2 = range_pack1;
        range_pack2.set(range_pack1 > nlev, 1);
        ekat::index_and_shift<1>(s_z_int, range_pack2, zi_k, zi_kp1);
        dz(i,k).set(range_pack1 < nlev, zi_k - zi_kp1);

        zt_grid(i,k) = z_mid(i,k) - z_int(i, nlevi_v)[nlevi_p];
        rrho(i,k) = (1/ggr)*(pseudo_density(i,k)/dz(i,k));

        wm_zt(i,k) = -1*omega(i,k)/(rrho(i,k)*ggr);

        zi_grid(i,k) = z_int(i,k) - z_int(i, nlevi_v)[nlevi_p];
      });
      zi_grid(i,nlevi_v)[nlevi_p] = 0;

      team.team_barrier();

      const auto zt_grid_s = ekat::subview(zt_grid, i);
      const auto zi_grid_s = ekat::subview(zi_grid, i);
      const auto rrho_s    = ekat::subview(rrho, i);
      const auto rrho_i_s  = ekat::subview(rrho_i, i);
      SHF::linear_interp(team,zt_grid_s,zi_grid_s,rrho_s,rrho_i_s,nlev,nlev+1,0);
      team.team_barrier();

      const int nlev_v = (nlev-1)/Spack::n;
      const int nlev_p = (nlev-1)%Spack::n;

      wpthlp_sfc(i)[0] = surf_sens_flux(i)[0]/(cpair*rrho_i(i,nlev_v)[nlev_p]);
      wprtp_sfc(i)[0]  = surf_latent_flux(i)[0]/rrho_i(i,nlev_v)[nlev_p];
      upwp_sfc(i)[0]   = surf_u_mom_flux(i)[0]/rrho_i(i,nlev_v)[nlev_p];
      vpwp_sfc(i)[0]   = surf_v_mom_flux(i)[0]/rrho_i(i,nlev_v)[nlev_p];

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_tracer_packs), [&] (const Int& q) {
        wtracer_sfc(i,q) = 0;
      });
    } // operator

    // Local variables
    int ncol, nlev, nlev_packs, num_tracers, num_tracer_packs;
    view_2d_const T_mid;
    view_2d_const z_int;
    view_2d_const z_mid;
    view_2d_const p_mid;
    view_2d_const pseudo_density;
    view_2d_const omega;
    view_1d_const phis;
    view_1d_const surf_sens_flux;
    view_1d_const surf_latent_flux;
    view_1d_const surf_u_mom_flux;
    view_1d_const surf_v_mom_flux;
    view_2d_const qv;
    view_2d       qv_copy;
    view_3d       Q;
    view_2d       qc;
    view_2d       qc_copy;
    view_2d       shoc_s;
    view_2d       tke;
    view_2d       tke_copy;
    view_2d       rrho;
    view_2d       rrho_i;
    view_2d       thv;
    view_2d       dz;
    view_2d       zt_grid;
    view_2d       zi_grid;
    view_1d       wpthlp_sfc;
    view_1d       wprtp_sfc;
    view_1d       upwp_sfc;
    view_1d       vpwp_sfc;
    view_2d       wtracer_sfc;
    view_2d       wm_zt;
    view_2d       exner;
    view_2d       thlm;
    view_2d       qw;
    view_2d       cloud_frac;
    view_3d       tracers;

    // Assigning local variables
    void set_variables(const int ncol_, const int nlev_, const int num_tracers_,
                       const int nlev_packs_, const int num_tracer_packs_,
                       view_2d_const T_mid_, view_2d_const z_int_,
                       view_2d_const z_mid_, view_2d_const p_mid_, view_2d_const pseudo_density_,
                       view_2d_const omega_,
                       view_1d_const phis_, view_1d_const surf_sens_flux_, view_1d_const surf_latent_flux_,
                       view_1d_const surf_u_mom_flux_, view_1d_const surf_v_mom_flux_,
                       view_2d_const qv_, view_2d qv_copy_, view_3d Q_, view_2d qc_, view_2d qc_copy_,
                       view_2d tke_, view_2d tke_copy_, view_2d s_, view_2d rrho_, view_2d rrho_i_,
                       view_2d thv_, view_2d dz_,view_2d zt_grid_,view_2d zi_grid_, view_1d wpthlp_sfc_,
                       view_1d wprtp_sfc_,view_1d upwp_sfc_,view_1d vpwp_sfc_, view_2d wtracer_sfc_,
                       view_2d wm_zt_,view_2d exner_,view_2d thlm_,view_2d qw_,view_3d tracers_)
    {
      ncol = ncol_;
      nlev = nlev_;
      nlev_packs = nlev_packs_;
      num_tracers = num_tracers_;
      num_tracer_packs = num_tracer_packs_;
      // IN
      T_mid = T_mid_;
      z_int = z_int_;
      z_mid = z_mid_;
      p_mid = p_mid_;
      pseudo_density = pseudo_density_;
      omega = omega_;
      phis = phis_;
      surf_sens_flux = surf_sens_flux_;
      surf_latent_flux = surf_latent_flux_;
      surf_u_mom_flux = surf_u_mom_flux_;
      surf_v_mom_flux = surf_v_mom_flux_;
      qv = qv_;
      qv_copy = qv_copy_;
      Q = Q_;
      // OUT
      qc = qc_;
      qc_copy = qc_copy_;
      shoc_s = s_;
      tke = tke_;
      tke_copy = tke_copy_;
      rrho = rrho_;
      rrho_i = rrho_i_;
      thv = thv_;
      dz = dz_;
      zt_grid = zt_grid_;
      zi_grid = zi_grid_;
      wpthlp_sfc = wpthlp_sfc_;
      wprtp_sfc = wprtp_sfc_;
      upwp_sfc = upwp_sfc_;
      vpwp_sfc = vpwp_sfc_;
      wtracer_sfc = wtracer_sfc_;
      wm_zt = wm_zt_;
      exner = exner_;
      thlm = thlm_;
      qw = qw_;
      tracers = tracers_;
    } // set_variables
  }; // SHOCPreprocess
  /* --------------------------------------------------------------------------------------------*/

  /*--------------------------------------------------------------------------------------------*/
  // Structure to handle the generation of data needed by the rest of the model based on output from
  // shoc_main.
  struct SHOCPostprocess {
    SHOCPostprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank();

      const Real cpair = C::Cpair;
      const Real inv_qc_relvar_max = 10;
      const Real inv_qc_relvar_min = 0.001;

      // Transpose and repack tracer array
      // TODO: remove once shoc tracer traspose is implemented
      for (int k=0; k<nlev; ++k) {
        for (int q=0; q<num_tracers; ++q) {
          const int q_v = q/Spack::n;
          const int q_p = q%Spack::n;
          const int k_v = k/Spack::n;
          const int k_p = k%Spack::n;

          Q(i,q,k_v)[k_p] = tracers(i,k,q_v)[q_p];
        }
      }

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
        // Tracers were updated as a group. Copy output from seperate inputs instead
        // of using tracer group output.
        tke(i,k) = tke_copy(i,k);
        qc(i,k) = qc_copy(i,k);
        qv(i,k) = qv_copy(i,k);

        cldfrac_liq(i,k) = ekat::min(cldfrac_liq(i,k), 1);
        sgs_buoy_flux(i,k) = sgs_buoy_flux(i,k)*rrho(i,k)*cpair;

        inv_qc_relvar(i,k) = 1;
        const auto condition = (qc(i,k) != 0 && qc2(i,k) != 0);
        inv_qc_relvar(i,k).set(condition,
                               ekat::min(inv_qc_relvar_max,
                                         ekat::max(inv_qc_relvar_min,
                                                   ekat::square(qc(i,k))/qc2(i,k))));
      });
    } // operator

    // Local variables
    int ncol, nlev, nlev_packs, num_tracers, num_tracer_packs;
    view_2d_const rrho;
    view_2d qv, qc, tke;
    view_2d_const qv_copy, qc_copy, tke_copy;
    view_2d_const qc2;
    view_3d Q;
    view_3d_const tracers;
    view_2d cldfrac_liq;
    view_2d sgs_buoy_flux;
    view_2d inv_qc_relvar;

    // Assigning local variables
    void set_variables(const int ncol_, const int nlev_, const int num_tracers_,
                       const int nlev_packs_, const int num_tracer_packs_,
                       view_2d_const rrho_,
                       view_2d qv_, view_2d_const qv_copy_, view_2d qc_, view_2d_const qc_copy_,
                       view_2d tke_, view_2d_const tke_copy_, view_2d_const qc2_,
                       view_3d Q_, view_3d_const tracers_,
                       view_2d cldfrac_liq_, view_2d sgs_buoy_flux_, view_2d inv_qc_relvar_)
    {
      ncol = ncol_;
      nlev = nlev_;
      nlev_packs = nlev_packs_;
      num_tracers = num_tracers_;
      num_tracer_packs = num_tracer_packs_;
      rrho = rrho_;
      qv = qv_;
      qv_copy = qv_copy_;
      qc = qc_;
      qc_copy = qc_copy_;
      tke = tke_;
      tke_copy = tke_copy_;
      qc2 = qc2_;
      Q = Q_;
      tracers = tracers_;
      cldfrac_liq = cldfrac_liq_;
      sgs_buoy_flux = sgs_buoy_flux_;
      inv_qc_relvar = inv_qc_relvar_;
    } // set_variables
  }; // SHOCPostprocess
  /* --------------------------------------------------------------------------------------------*/

protected:

  // The three main interfaces for the subcomponent
  void initialize_impl (const util::TimeStamp& t0);
  void run_impl        (const Real dt);
  void finalize_impl   ();

  // Setting the fields in the atmospheric process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

  std::map<std::string,const_field_type>  m_shoc_fields_in;
  std::map<std::string,field_type>        m_shoc_fields_out;

  util::TimeStamp     m_current_ts;
  ekat::Comm          m_shoc_comm;
  ekat::ParameterList m_shoc_params;

  // Keep track of field dimensions and other scalar values
  // needed in shoc_main
  Int m_num_cols;
  Int m_num_levs;
  Int m_npbl;
  Int m_nadv;
  Int m_num_tracers;
  Int hdtime;

  // Store the structures for each arguement to shoc_main;
  SHF::SHOCInput input;
  SHF::SHOCInputOutput input_output;
  SHF::SHOCOutput output;
  SHF::SHOCHistoryOutput history_output;

  // Structures which compute pre/post process
  SHOCPreprocess shoc_preprocess;
  SHOCPostprocess shoc_postprocess;
}; // class SHOCMacrophysics

} // namespace scream

#endif // SCREAM_SHOC_MACROPHYSICS_HPP
