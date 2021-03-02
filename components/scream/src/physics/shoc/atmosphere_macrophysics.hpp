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
  void register_fields (FieldRepository<Real>& field_repo) const;

  // SHOC updates the 'TRACERS' group.
  void set_updated_group (const FieldGroup<Real>& group);

  // Get the set of required/computed fields and groups
  const std::set<FieldIdentifier>& get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>& get_computed_fields () const { return m_computed_fields; }
  std::set<GroupRequest> get_updated_groups () const { return m_inout_groups_req; }

  /*--------------------------------------------------------------------------------------------*/
  // Most individual processes have a pre-processing step that constructs needed variables from
  // the set of fields stored in the field manager.  A structure like this defines those operations,
  // which can then be called during run_impl in the main .cpp code.
  // Structure to handle the local generation of data needed by shoc_main in run_impl
  struct SHOCPreamble {
    SHOCPreamble() = default;

    KOKKOS_INLINE_FUNCTION
    void operator()(const Kokkos::TeamPolicy<KT::ExeSpace>::member_type& team) const {
      const int i = team.league_rank();

      const Real zvir = C::ZVIR;
      const Real latvap = C::LatVap;
      const Real cpair = C::Cpair;
      const Real ggr = C::gravit;

      // Transpose and repack tracer array
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

      const auto sub_zi = ekat::subview(zi, i);
      const auto s_zi = ekat::scalarize(sub_zi);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
        // Exner
        const Spack opmid(pmid(i,k));
        const Smask opmid_mask(!isnan(opmid) and opmid>0.0);
        auto oexner = physics_fun::get_exner(opmid,opmid_mask);
        exner(i,k) = oexner;

        // Temperature
        qw(i,k) = shoc_qv(i,k) + shoc_ql(i,k);

        const auto theta_zt = t(i,k)/exner(i,k);
        thlm(i,k) = theta_zt - (latvap/cpair)*shoc_ql(i,k);
        thv(i,k)  = theta_zt*(1 + zvir*shoc_qv(i,k) - shoc_ql(i,k));

        // TKE
        tke_zt(i,k) = ekat::max(sp(0.004), tke_zt(i,k));

        // Cloudfrac
        cloud_frac(i,k) = alst(i,k);

        Spack zi_k, zi_kp1;
        auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
        auto range_pack2 = range_pack1;
        range_pack2.set(range_pack1 > nlev, 1);
        ekat::index_and_shift<1>(s_zi, range_pack2, zi_k, zi_kp1);
        dz(i,k).set(range_pack1 < nlev, zi_k - zi_kp1);

        zt_grid(i,k) = zm(i,k) - zi(i, nlevi_v)[nlevi_p];
        rrho(i,k) = (1/ggr)*(pdel(i,k)/dz(i,k));
        wm_zt(i,k) = -1*omega(i,k)/(rrho(i,k)*ggr);

        zi_grid(i,k) = zi(i,k) - zi(i, nlevi_v)[nlevi_p];
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

      wpthlp_sfc(i)[0] = shf(i)[0]/(cpair*rrho_i(i,nlev_v)[nlev_p]);
      wprtp_sfc(i)[0]  = cflx(i)[0]/rrho_i(i,nlev_v)[nlev_p];
      upwp_sfc(i)[0]   = wsx(i)[0]/rrho_i(i,nlev_v)[nlev_p];
      vpwp_sfc(i)[0]   = wsy(i)[0]/rrho_i(i,nlev_v)[nlev_p];

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, num_tracer_packs), [&] (const Int& q) {
        wtracer_sfc(i,q) = 0;
      });
    } // operator

    // Local variables
    int ncol, nlev, nlev_packs, num_tracers, num_tracer_packs;
    view_2d_const t;
    view_2d_const alst;
    view_2d_const zi;
    view_2d_const zm;
    view_2d_const pmid;
    view_2d_const pdel;
    view_2d_const omega;
    view_1d_const shf;
    view_1d_const cflx;
    view_1d_const wsx;
    view_1d_const wsy;
    view_2d_const shoc_qv;
    view_3d       Q;
    view_2d       shoc_ql;
    view_2d       shoc_s;
    view_2d       tke_zt;
    view_2d       um;
    view_2d       vm;
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
                       view_2d_const t_, view_2d_const alst_, view_2d_const zi_,
                       view_2d_const zm_, view_2d_const pmid_, view_2d_const pdel_,
                       view_2d_const omega_,
                       view_1d_const shf_, view_1d_const cflx_, view_1d_const wsx_, view_1d_const wsy_,
                       view_2d_const shoc_qv_, view_3d Q_, view_2d shoc_ql_, view_2d tke_,
                       view_2d s_, view_2d u_, view_2d v_, view_2d rrho_, view_2d rrho_i_,view_2d thv_,
                       view_2d dz_,view_2d zt_grid_,view_2d zi_grid_, view_1d wpthlp_sfc_,
                       view_1d wprtp_sfc_,view_1d upwp_sfc_,view_1d vpwp_sfc_, view_2d wtracer_sfc_,
                       view_2d wm_zt_,view_2d exner_,view_2d thlm_,view_2d qw_,view_2d cloud_frac_,view_3d tracers_)
    {
      ncol = ncol_;
      nlev = nlev_;
      nlev_packs = nlev_packs_;
      num_tracers = num_tracers_;
      num_tracer_packs = num_tracer_packs_;
      // IN
      t = t_;
      alst = alst_;
      zi = zi_;
      zm = zm_;
      pmid = pmid_;
      pdel = pdel_;
      omega = omega_;
      shf = shf_;
      cflx = cflx_;
      wsx = wsx_;
      wsy = wsy_;
      shoc_qv = shoc_qv_;
      Q = Q_;
      // OUT
      shoc_ql = shoc_ql_;
      shoc_s = s_;
      tke_zt = tke_;
      um = u_;
      vm = v_;
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
      cloud_frac = cloud_frac_;
      tracers = tracers_;
    } // set_variables
  }; // SHOCPreamble
  /* --------------------------------------------------------------------------------------------*/

protected:

  // The three main interfaces for the subcomponent
  void initialize_impl (const util::TimeStamp& t0);
  void run_impl        (const Real dt);
  void finalize_impl   ();

  // Setting the fields in the atmospheric process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

  std::set<FieldIdentifier> m_required_fields;
  std::set<FieldIdentifier> m_computed_fields;
  std::set<GroupRequest>    m_inout_groups_req;

  std::map<std::string,const_field_type>  m_shoc_fields_in;
  std::map<std::string,field_type>        m_shoc_fields_out;

  // Used to init some fields. For now, only needed for stand-alone shoc runs
  std::shared_ptr<FieldInitializer>  m_initializer;

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
  SHOCPreamble shoc_preamble;
}; // class SHOCMacrophysics

} // namespace scream

#endif // SCREAM_SHOC_MACROPHYSICS_HPP
