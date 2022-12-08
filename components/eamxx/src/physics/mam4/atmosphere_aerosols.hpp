#ifndef SCREAM_MAM4_AEROSOLS_HPP
#define SCREAM_MAM4_AEROSOLS_HPP

#include <share/atm_process/atmosphere_process.hpp>
//#include <share/util/scream_common_physics_functions.hpp>
//#include <share/atm_process/ATMBufferManager.hpp>

#include <ekat/ekat_parameter_list.hpp>
#include <mam4xx/mam4.hpp>

#include <string>

#ifndef KOKKOS_ENABLE_CUDA
#define protected public
#define private public
#endif

namespace scream
{

/*
 * The class responsible for handling MAM4 aerosols
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
*/

class MAM4Aerosols final : public scream::AtmosphereProcess
{
  using PF           = scream::PhysicsFunctions<DefaultDevice>;
  using KT           = ekat::KokkosTypes<DefaultDevice>;

  using ColumnView   = mam4::ColumnView;
  using ThreadTeam   = mam4::ThreadTeam;

public:

  // Constructor
  MAM4Aerosols(const ekat::Comm& comm, const ekat::ParameterList& params);

protected:

  // --------------------------------------------------------------------------
  // AtmosphereProcess overrides (see share/atm_process/atmosphere_process.hpp)
  // --------------------------------------------------------------------------

  // process metadata
  AtmosphereProcessType type() const override;
  std::string name() const override;

  // grid
  void set_grids(const std::shared_ptr<const GridsManager> grids_manager) override;

  // management of common atm process memory
  size_t requested_buffer_size_in_bytes() const override;
  void init_buffers(const ATMBufferManager &buffer_manager) override;

  // process behavior
  void initialize_impl(const RunType run_type) override;
  void run_impl(const int dt) override;
  void finalize_impl() override;

  // MAM4xx updates the 'tracers' group.
  void set_computed_group_impl(const FieldGroup& group) override;

private:

  // Atmosphere processes often have a pre-processing step that constructs
  // required variables from the set of fields stored in the field manager.
  // This functor implements this step, which is called during run_impl.
  struct MAM4Preprocess {
    MAM4Preprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank();

      const Real zvir = C::ZVIR;
      const Real latvap = C::LatVap;
      const Real cpair = C::Cpair;
      const Real ggr = C::gravit;
      const Real inv_ggr = 1/ggr;

      const int nlev_packs = ekat::npack<Spack>(nlev);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
        //--------------------------
        // Wet to dry mixing ratios
        //--------------------------
        //
        // Since tracers from the host model (or AD) are wet mixing ratios, and
        // MAM4 expects these tracers in dry mixing ratios, we convert the wet
        // mixing ratios to dry mixing ratios for all the tracers.
        //
        // The function calculate_drymmr_from_wetmmr takes 2 arguments:
        // 1. wet mmr
        // 2. "wet" water vapor mixing ratio
        //
        // Units of all tracers become [kg/kg(dry-air)] for mass mixing ratios and
        // [#/kg(dry-air)] for number mixing ratios after the following
        // conversion. qv is converted to dry mmr in the next parallel for.
        for (int iq = 0; iq < num_qtracers-1; ++iq) {
          qtracers(i,convert_wet_dry_idx_d(iq),k) = PF::calculate_drymmr_from_wetmmr(qtracers(i,convert_wet_dry_idx_d(iq),k), qv(i,k));
        }

        // convert qv to dry mmr
        qv(i,k) = PF::calculate_drymmr_from_wetmmr(qv(i,k), qv(i,k));

        const Smask in_nlev_range = (range < nlev);

        // Tracers are updated as a group. The tracers tke and qc act as separate inputs to shoc_main()
        // and are therefore updated differently to the bundled tracers. Here, we make a copy if each
        // of these tracers and pass to shoc_main() so that changes to the tracer group does not alter
        // tke or qc  values. Then during post processing, we copy back correct values of tke and qc
        // to tracer group in postprocessing.
        // TODO: remove *_copy views once SHOC can request a subset of tracers.
        tke_copy(i,k) = tke(i,k);
        qc_copy(i,k)  = qc(i,k); //at this point, qc should be dry mmr [kg/kg(dry air)]


        qw(i,k) = qv(i,k) + qc(i,k);

        // atmospheric state
        const auto theta_zt = PF::calculate_theta_from_T(T_mid(i,k),p_mid(i,k));
        thlm(i,k) = theta_zt-(theta_zt/T_mid(i,k))*(latvap/cpair)*qc(i,k);
        thv(i,k)  = theta_zt*(1 + zvir*qv(i,k) - qc(i,k));

        // Vertical layer thickness
        dz(i,k) = PF::calculate_dz(pseudo_density(i,k), p_mid(i,k), T_mid(i,k), qv(i,k));

        rrho(i,k) = inv_ggr*(pseudo_density(i,k)/dz(i,k));
        wm_zt(i,k) = -1*omega(i,k)/(rrho(i,k)*ggr);
      });
      team.team_barrier();

      // Compute vertical layer heights
      const auto dz_s    = ekat::subview(dz,    i);
      const auto z_int_s = ekat::subview(z_int, i);
      const auto z_mid_s = ekat::subview(z_mid, i);
      PF::calculate_z_int(team,nlev,dz_s,z_surf,z_int_s);
      team.team_barrier();
      PF::calculate_z_mid(team,nlev,z_int_s,z_mid_s);
      team.team_barrier();

      const int nlevi_v = nlev/Spack::n;
      const int nlevi_p = nlev%Spack::n;
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
        zt_grid(i,k) = z_mid(i,k) - z_int(i, nlevi_v)[nlevi_p];
        zi_grid(i,k) = z_int(i,k) - z_int(i, nlevi_v)[nlevi_p];

        // Dry static energy
        shoc_s(i,k) = PF::calculate_dse(T_mid(i,k),z_mid(i,k),phis(i));
      });
      zi_grid(i,nlevi_v)[nlevi_p] = 0;
      team.team_barrier();

      const auto zt_grid_s = ekat::subview(zt_grid, i);
      const auto zi_grid_s = ekat::subview(zi_grid, i);
      const auto rrho_s    = ekat::subview(rrho, i);
      const auto rrho_i_s  = ekat::subview(rrho_i, i);
      SHF::linear_interp(team,zt_grid_s,zi_grid_s,rrho_s,rrho_i_s,nlev,nlev+1,0);
      team.team_barrier();

      // For now, we are considering dy=dx. Here, we
      // will need to compute dx/dy instead of cell_length
      // if we have dy!=dx.
      cell_length(i) = PF::calculate_dx_from_area(area(i),lat(i));

      const auto exner_int = PF::exner_function(p_int(i,nlevi_v)[nlevi_p]);
      const auto inv_exner_int_surf = 1/exner_int;

      wpthlp_sfc(i) = (surf_sens_flux(i)/(cpair*rrho_i(i,nlevi_v)[nlevi_p]))*inv_exner_int_surf;
      wprtp_sfc(i)  = surf_evap(i)/rrho_i(i,nlevi_v)[nlevi_p];
      upwp_sfc(i)   = surf_mom_flux(i,0)/rrho_i(i,nlevi_v)[nlevi_p];
      vpwp_sfc(i)   = surf_mom_flux(i,1)/rrho_i(i,nlevi_v)[nlevi_p];

      const int num_qtracer_packs = ekat::npack<Spack>(num_qtracers);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_qtracer_packs), [&] (const Int& q) {
        wtracer_sfc(i,q) = 0;
      });
    } // operator

    int ncol, nlev, num_qtracers;

    // used for converting between wet and dry mixing ratios
    view_1d_int convert_wet_dry_idx_d;

    // local atmospheric state column variables
    view_2d_const T_mid;  // temperature at grid midpoints [K]
    view_2d_const p_mid;  // total pressure at grid midpoints [Pa]
    view_2d       qv;     // water vapor mass mixing ratio, not const because it
                          // must be converted from wet to dry
                          // [kg vapor/kg dry air]
    view_2d       height; // height at grid interfaces [m]
    view_2d       pdel;   // hydrostatic "pressure thickness" at grid
                          // interfaces [Pa]
    view_1d_const pblh;   // planetary boundary layer height [m]

    // local aerosol-related gases
    view_2d       q_soag;  // secondary organic aerosol gas [kg gas/kg dry air]
    view_2d       q_h2so4; // H2SO3 gas [kg/kg dry air]
    view_2d       q_nh3;   // NH3 gas [kg/kg dry air]

    // local aerosols (more to appear as we improve this atm process)
    view_2d       q_aitken_so4; // SO4 aerosol in aitken mode [kg/kg dry air]

    // assigns local variables
    void set_variables(const int ncol_, const int nlev_, const int num_qtracers_,
                       const view_1d_int& convert_wet_dry_idx_d_,
                       const view_2d_const& T_mid_,
                       const view_2d_const& p_mid_,
                       const view_2d&       qv_,
                       const view_2d&       height_,
                       const view_2d&       pdel_,
                       const view_1d_const& pblh_,
                       const view_2d&       q_soag_,
                       const view_2d&       q_h2so4,
                       const view_2d&       q_nh3,
                       const view_2d&       q_aitken_so4) {
      ncol = ncol_;
      nlev = nlev_;
      num_qtracers = num_qtracers_;
      convert_wet_dry_idx_d = convert_wet_dry_idx_d_;
      T_mid = T_mid_;
      p_mid = p_mid_;
      qv = qv_;
      height = height_;
      pdel = pdel_;
      pblh = pblh_;
      q_soag = q_soag_;
      q_h2so4 = q_h2so4_;
      q_nh3 = q_nh3_;
      q_aitken_so4 = q_aitken_so4_;
    } // set_variables
  }; // MAM4Preprocess

  // Postprocessing functor
  struct MAM4Postprocess {
    MAM4Postprocess() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank();

      const Real inv_qc_relvar_max = 10;
      const Real inv_qc_relvar_min = 0.001;

      //In the following loop, tke, qc and qv are updated. All these variables are slices of "qtracers" array
      //After these updates, all tracers (except TKE) in the qtarcers array will be converted
      //from dry mmr to wet mmr
      const int nlev_packs = ekat::npack<Spack>(nlev);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
        // See comment in SHOCPreprocess::operator() about the necessity of *_copy views
        tke(i,k) = tke_copy(i,k);
        qc(i,k)  = qc_copy(i,k);

        qv(i,k) = qw(i,k) - qc(i,k);

        cldfrac_liq(i,k) = ekat::min(cldfrac_liq(i,k), 1);

        inv_qc_relvar(i,k) = 1;
        const auto condition = (qc(i,k) != 0 && qc2(i,k) != 0);
        if (condition.any()) {
          //inv_qc_relvar is used in P3 and is computed here using qc and qc2, which are in dry mmr
          inv_qc_relvar(i,k).set(condition,
                                 ekat::min(inv_qc_relvar_max,
                                           ekat::max(inv_qc_relvar_min,
                                                     ekat::square(qc(i,k))/qc2(i,k))));
        }

        // Temperature
        const Spack dse_ik(dse(i,k));
        const Spack z_mid_ik(z_mid(i,k));
        const Real  phis_i(phis(i));
        T_mid(i,k) = PF::calculate_temperature_from_dse(dse_ik,z_mid_ik,phis_i);


      /*--------------------------------------------------------------------------------
       *DRY-TO-WET MMRs:
       *-----------------
       *Since the host model (or AD) expects wet mixing ratios, we need to convert dry
       *mixing ratios from SHOC to wet mixing ratios except for qv[which will be converted
       *in the following "parallel for" after all the other tracers] and TKE [which is
       already in wet mmr].
       *---------------------------------------------------------------------------------
       */

          //NOTE:Function calculate_wetmmr_from_drymmr takes 2 arguments: ( dry mmr and "dry"
          //water vapor mixing ratio)
          //Units of all tracers (except TKE and qv) will become [kg/kg(wet-air)] for mass and
          //[#/kg(wet-air)] for number after the following conversion. qv will be converted
          //to wet mmr in the next parallel for
          for (Int iq = 0; iq < num_qtracers-2; ++iq)
            qtracers(i,convert_wet_dry_idx_d(iq),k) = PF::calculate_wetmmr_from_drymmr(qtracers(i,convert_wet_dry_idx_d(iq),k), qv(i,k));
          qv(i,k) = PF::calculate_wetmmr_from_drymmr(qv(i,k), qv(i,k));
      });

      // If necessary, set appropriate boundary fluxes for energy and mass conservation checks.
      // Any boundary fluxes not included in SHOC interface are set to 0.
      if (compute_mass_and_energy_fluxes) {
        vapor_flux(i) = surf_evap(i);
        water_flux(i) = 0.0;
        ice_flux(i)   = 0.0;
        heat_flux(i)  = surf_sens_flux(i);
      }
    } // operator

    // Local variables
    int ncol, nlev, num_qtracers;
    view_1d_int convert_wet_dry_idx_d;
    view_2d_const rrho;
    view_2d qv, qc, tke;
    view_2d_const tke_copy, qc_copy, qw;
    view_3d qtracers;
    view_2d_const qc2;
    view_2d cldfrac_liq;
    view_2d inv_qc_relvar;
    view_2d T_mid;
    view_2d_const dse,z_mid;
    view_1d_const phis;
    bool compute_mass_and_energy_fluxes = false;
    view_1d_const surf_evap;
    view_1d_const surf_sens_flux;
    view_1d vapor_flux;
    view_1d water_flux;
    view_1d ice_flux;
    view_1d heat_flux;

    // Assigning local variables
    void set_variables(const int ncol_, const int nlev_, const int num_qtracers_,
                       const view_1d_int& convert_wet_dry_idx_d_,
                       const view_2d_const& rrho_,
                       const view_2d& qv_, const view_2d_const& qw_, const view_2d& qc_, const view_2d_const& qc_copy_,
                       const view_2d& tke_, const view_2d_const& tke_copy_, const view_3d& qtracers_, const view_2d_const& qc2_,
                       const view_2d& cldfrac_liq_, const view_2d& inv_qc_relvar_,
                       const view_2d& T_mid_, const view_2d_const& dse_, const view_2d_const& z_mid_, const view_1d_const phis_)
    {
      ncol = ncol_;
      nlev = nlev_;
      num_qtracers = num_qtracers_;
      convert_wet_dry_idx_d = convert_wet_dry_idx_d_;
      rrho = rrho_;
      qv = qv_;
      qw = qw_;
      qc = qc_;
      qc_copy = qc_copy_;
      tke = tke_;
      tke_copy = tke_copy_;
      qtracers = qtracers_;
      qc2 = qc2_;
      cldfrac_liq = cldfrac_liq_;
      inv_qc_relvar = inv_qc_relvar_;
      T_mid = T_mid_;
      dse = dse_;
      z_mid = z_mid_;
      phis = phis_;
    } // set_variables

    void set_mass_and_energy_fluxes (const view_1d_const& surf_evap_, const view_1d_const surf_sens_flux_,
				     const view_1d& vapor_flux_, const view_1d& water_flux_,
                                     const view_1d& ice_flux_, const view_1d& heat_flux_)
    {
      compute_mass_and_energy_fluxes = true;
      surf_evap = surf_evap_;
      surf_sens_flux = surf_sens_flux_;
      vapor_flux = vapor_flux_;
      water_flux = water_flux_;
      ice_flux = ice_flux_;
      heat_flux = heat_flux_;
    }
  }; // MAM4Postprocess

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    static constexpr int num_1d_scalar_ncol = 18;
    static constexpr int num_1d_scalar_nlev = 1;
    static constexpr int num_2d_vector_mid  = 22;
    static constexpr int num_2d_vector_int  = 13;
    static constexpr int num_2d_vector_tr   = 1;

    uview_1d<Real> cell_length;
    uview_1d<Real> wpthlp_sfc;
    uview_1d<Real> wprtp_sfc;
    uview_1d<Real> upwp_sfc;
    uview_1d<Real> vpwp_sfc;

    uview_1d<Spack> pref_mid;

    uview_2d<Spack> z_mid;
    uview_2d<Spack> z_int;
    uview_2d<Spack> rrho;
    uview_2d<Spack> rrho_i;
    uview_2d<Spack> thv;
    uview_2d<Spack> dz;
    uview_2d<Spack> zt_grid;
    uview_2d<Spack> zi_grid;
    uview_2d<Spack> wtracer_sfc;
    uview_2d<Spack> wm_zt;
    uview_2d<Spack> inv_exner;
    uview_2d<Spack> thlm;
    uview_2d<Spack> qw;
    uview_2d<Spack> dse;
    uview_2d<Spack> tke_copy;
    uview_2d<Spack> qc_copy;
    uview_2d<Spack> shoc_ql2;
    uview_2d<Spack> shoc_mix;
    uview_2d<Spack> isotropy;
    uview_2d<Spack> w_sec;
    uview_2d<Spack> thl_sec;
    uview_2d<Spack> qw_sec;
    uview_2d<Spack> qwthl_sec;
    uview_2d<Spack> wthl_sec;
    uview_2d<Spack> wqw_sec;
    uview_2d<Spack> wtke_sec;
    uview_2d<Spack> uw_sec;
    uview_2d<Spack> vw_sec;
    uview_2d<Spack> w3;
    uview_2d<Spack> wqls_sec;
    uview_2d<Spack> brunt;

    Spack* wsm_data;
  };
  /* --------------------------------------------------------------------------------------------*/

  // MAM4 aerosol particle size description
  mam4::AeroConfig aero_config_;

  // Aerosol processes
  mam4::NucleationProcess nucleation_;

  // WSM for internal local variables
  ekat::WorkspaceManager<Spack, KT::Device> workspace_mgr_;

  // physics grid for column information
  std::shared_ptr<const AbstractGrid> grid_;
};

} // namespace scream

#endif // SCREAM_MAM4_AEROSOLS_HPP
