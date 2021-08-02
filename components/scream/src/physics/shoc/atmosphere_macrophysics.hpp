#ifndef SCREAM_SHOC_MACROPHYSICS_HPP
#define SCREAM_SHOC_MACROPHYSICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "physics/shoc/shoc_main_impl.hpp"
#include "physics/shoc/shoc_functions.hpp"
#include "share/util/scream_common_physics_functions.hpp"
#include "share/atm_process/ATMBufferManager.hpp"

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
  using PF           = scream::PhysicsFunctions<DefaultDevice>;
  using C            = physics::Constants<Real>;
  using KT           = ekat::KokkosTypes<DefaultDevice>;

  using Spack                = typename SHF::Spack;
  using IntSmallPack         = typename SHF::IntSmallPack;
  using Smask                = typename SHF::Smask;
  using view_1d              = typename SHF::view_1d<Real>;
  using view_1d_const        = typename SHF::view_1d<const Real>;
  using view_2d              = typename SHF::view_2d<SHF::Spack>;
  using view_2d_const        = typename SHF::view_2d<const Spack>;
  using sview_2d             = typename KokkosTypes<DefaultDevice>::template view_2d<Real>;
  using view_3d              = typename SHF::view_3d<Spack>;
  using view_3d_const        = typename SHF::view_3d<const Spack>;

  using WSM = ekat::WorkspaceManager<Spack, KT::Device>;

  template<typename ScalarT>
  using uview_1d = Unmanaged<typename KT::template view_1d<ScalarT>>;
  template<typename ScalarT>
  using uview_2d = Unmanaged<typename KT::template view_2d<ScalarT>>;

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

      const int nlev_packs = ekat::npack<Spack>(nlev);
      const int num_qtracer_packs = ekat::npack<Spack>(num_qtracers);

      const int nlevi_v = nlev/Spack::n;
      const int nlevi_p = nlev%Spack::n;

      const auto sub_z_int = ekat::subview(z_int, i);
      const auto s_z_int = ekat::scalarize(sub_z_int);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
        const auto range = ekat::range<IntSmallPack>(k*Spack::n);
        const Smask in_nlev_range = (range < nlev);

        // Inverse of Exner. Assert that exner != 0 when in range before computing.
        const Spack p_mid_ik(p_mid(i,k));
        const Spack exner_ik = PF::exner_function(p_mid_ik);
        const Smask nonzero = (exner_ik != 0);
        EKAT_KERNEL_ASSERT((nonzero || !in_nlev_range).any());
        inv_exner(i,k).set(nonzero, 1/exner_ik);

        tke(i,k) = ekat::max(sp(0.004), tke(i,k));

        // Tracers are updated as a group. These tracers act as seperate inputs to shoc_main
        // and should have different output as tracer group. We make a copy and pass these to
        // shoc_main, then copy back to tracer group in postprocessing.
        // TODO: remove *_copy views once SHOC can request a subset of tracers.
        tke_copy(i,k) = tke(i,k);
        qc_copy(i,k) = qc(i,k);
        qv_copy(i,k) = qv(i,k);

        // Temperature
        const Spack T_mid_ik(T_mid(i,k));
        qw(i,k) = qv(i,k) + qc(i,k);

        const auto& theta_zt = PF::calculate_theta_from_T(T_mid_ik,p_mid_ik);
        thlm(i,k) = theta_zt - (latvap/cpair)*qc(i,k);
        thv(i,k)  = theta_zt*(1 + zvir*qv(i,k) - qc(i,k));

        // Dry static energy
        const Spack z_mid_ik(z_mid(i,k));
        const Real  phis_i(phis(i));
        shoc_s(i,k) = PF::calculate_dse(T_mid_ik,z_mid_ik,phis_i);

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

      const auto& zt_grid_s = ekat::subview(zt_grid, i);
      const auto& zi_grid_s = ekat::subview(zi_grid, i);
      const auto& rrho_s    = ekat::subview(rrho, i);
      const auto& rrho_i_s  = ekat::subview(rrho_i, i);
      SHF::linear_interp(team,zt_grid_s,zi_grid_s,rrho_s,rrho_i_s,nlev,nlev+1,0);
      team.team_barrier();

      const int nlev_v = (nlev-1)/Spack::n;
      const int nlev_p = (nlev-1)%Spack::n;

      // For now, we are considering dy=dx. Here, we
      // will need to compute dx/dy instead of cell_length
      // if we have dy!=dx.
      cell_length(i) = sqrt(area(i));

      wpthlp_sfc(i) = surf_sens_flux(i)/(cpair*rrho_i(i,nlev_v)[nlev_p]);
      wprtp_sfc(i)  = surf_latent_flux(i)/rrho_i(i,nlev_v)[nlev_p];
      upwp_sfc(i)   = surf_u_mom_flux(i)/rrho_i(i,nlev_v)[nlev_p];
      vpwp_sfc(i)   = surf_v_mom_flux(i)/rrho_i(i,nlev_v)[nlev_p];

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_qtracer_packs), [&] (const Int& q) {
        wtracer_sfc(i,q) = 0;
      });
    } // operator

    // Local variables
    int ncol, nlev, num_qtracers;
    view_1d_const        area;
    view_2d_const        T_mid;
    view_2d_const        z_int;
    view_2d_const        z_mid;
    view_2d_const        p_mid;
    view_2d_const        pseudo_density;
    view_2d_const        omega;
    view_1d_const        phis;
    view_1d_const        surf_sens_flux;
    view_1d_const        surf_latent_flux;
    view_1d_const        surf_u_mom_flux;
    view_1d_const        surf_v_mom_flux;
    view_2d_const        qv;
    view_2d              qv_copy;
    view_1d              cell_length;
    view_2d              qc;
    view_2d              qc_copy;
    view_2d              shoc_s;
    view_2d              tke;
    view_2d              tke_copy;
    view_2d              rrho;
    view_2d              rrho_i;
    view_2d              thv;
    view_2d              dz;
    view_2d              zt_grid;
    view_2d              zi_grid;
    view_1d              wpthlp_sfc;
    view_1d              wprtp_sfc;
    view_1d              upwp_sfc;
    view_1d              vpwp_sfc;
    view_2d              wtracer_sfc;
    view_2d              wm_zt;
    view_2d              inv_exner;
    view_2d              thlm;
    view_2d              qw;
    view_2d              cloud_frac;

    // Assigning local variables
    void set_variables(const int ncol_, const int nlev_, const int num_qtracers_,
                       const view_1d_const& area_,
                       const view_2d_const& T_mid_, const view_2d_const& z_int_,
                       const view_2d_const& z_mid_, const view_2d_const& p_mid_, const view_2d_const& pseudo_density_,
                       const view_2d_const& omega_,
                       const view_1d_const& phis_, const view_1d_const& surf_sens_flux_, const view_1d_const& surf_latent_flux_,
                       const view_1d_const& surf_u_mom_flux_, const view_1d_const& surf_v_mom_flux_,
                       const view_2d_const& qv_, const view_2d& qv_copy_, const view_2d& qc_, const view_2d& qc_copy_,
                       const view_2d& tke_, const view_2d& tke_copy_, const view_1d& cell_length_,
                       const view_2d& s_, const view_2d& rrho_, const view_2d& rrho_i_,
                       const view_2d& thv_, const view_2d& dz_,const view_2d& zt_grid_,const view_2d& zi_grid_, const view_1d& wpthlp_sfc_,
                       const view_1d& wprtp_sfc_,const view_1d& upwp_sfc_,const view_1d& vpwp_sfc_, const view_2d& wtracer_sfc_,
                       const view_2d& wm_zt_,const view_2d& inv_exner_,const view_2d& thlm_,const view_2d& qw_)
    {
      ncol = ncol_;
      nlev = nlev_;
      num_qtracers = num_qtracers_;
      // IN
      area = area_;
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
      // OUT
      qc = qc_;
      qc_copy = qc_copy_;
      shoc_s = s_;
      tke = tke_;
      tke_copy = tke_copy_;
      cell_length = cell_length_;
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
      inv_exner = inv_exner_;
      thlm = thlm_;
      qw = qw_;
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

      const int nlev_packs = ekat::npack<Spack>(nlev);
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
        if (condition.any()) {
          inv_qc_relvar(i,k).set(condition,
                                 ekat::min(inv_qc_relvar_max,
                                           ekat::max(inv_qc_relvar_min,
                                                     ekat::square(qc(i,k))/qc2(i,k))));
        }
      });
    } // operator

    // Local variables
    int ncol, nlev;
    view_2d_const rrho;
    view_2d qv, qc, tke;
    view_2d_const qv_copy, qc_copy, tke_copy;
    view_2d_const qc2;
    view_2d cldfrac_liq;
    view_2d sgs_buoy_flux;
    view_2d inv_qc_relvar;

    // Assigning local variables
    void set_variables(const int ncol_, const int nlev_,
                       const view_2d_const& rrho_,
                       const view_2d& qv_, const view_2d_const& qv_copy_, const view_2d& qc_, const view_2d_const& qc_copy_,
                       const view_2d& tke_, const view_2d_const& tke_copy_, const view_2d_const& qc2_,
                       const view_2d& cldfrac_liq_, const view_2d& sgs_buoy_flux_, const view_2d& inv_qc_relvar_)
    {
      ncol = ncol_;
      nlev = nlev_;
      rrho = rrho_;
      qv = qv_;
      qv_copy = qv_copy_;
      qc = qc_;
      qc_copy = qc_copy_;
      tke = tke_;
      tke_copy = tke_copy_;
      qc2 = qc2_;
      cldfrac_liq = cldfrac_liq_;
      sgs_buoy_flux = sgs_buoy_flux_;
      inv_qc_relvar = inv_qc_relvar_;
    } // set_variables
  }; // SHOCPostprocess
  /* --------------------------------------------------------------------------------------------*/

  // Structure for storing local variables initialized using the ATMBufferManager
  struct Buffer {
    static constexpr int num_1d_scalar     = 5;
    static constexpr int num_2d_vector_mid = 18;
    static constexpr int num_2d_vector_int = 11;
    static constexpr int num_2d_vector_tr  = 1;

    uview_1d<Real> cell_length;
    uview_1d<Real> wpthlp_sfc;
    uview_1d<Real> wprtp_sfc;
    uview_1d<Real> upwp_sfc;
    uview_1d<Real> vpwp_sfc;

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
    uview_2d<Spack> s;
    uview_2d<Spack> qv_copy;
    uview_2d<Spack> qc_copy;
    uview_2d<Spack> tke_copy;
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


protected:

  // The three main interfaces for the subcomponent
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

  KokkosTypes<DefaultDevice>::view_1d<Real> m_cell_area;

  // Struct which contains local variables
  Buffer m_buffer;

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
