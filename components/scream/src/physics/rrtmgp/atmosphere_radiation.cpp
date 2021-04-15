#include "ekat/ekat_assert.hpp"
#include "physics/rrtmgp/atmosphere_radiation.hpp"
#include "physics/rrtmgp/rrtmgp_heating_rate.hpp"
#include "cpp/rrtmgp/mo_gas_concentrations.h"
#include "YAKL/YAKL.h"

namespace scream {

  using KT = KokkosTypes<DefaultDevice>;
  using ExeSpace = KT::ExeSpace;
  using MemberType = KT::MemberType;

RRTMGPRadiation::RRTMGPRadiation (const ekat::Comm& comm, const ekat::ParameterList& params) 
  : AtmosphereProcess::AtmosphereProcess(), m_rrtmgp_comm (comm), m_rrtmgp_params (params) {
}  // RRTMGPRadiation::RRTMGPRadiation

void RRTMGPRadiation::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {

  using namespace ekat::units;

  auto kgkg = kg/kg;
  kgkg.set_string("kg/kg");
  auto m3 = m * m * m;
  auto Wm2 = W / m / m;
  Wm2.set_string("Wm2");
  auto nondim = m/m;  // dummy unit for non-dimensional fields
  auto micron = m / 1000000;

  using namespace ShortFieldTagsNames;

  auto grid = grids_manager->get_grid("Physics");
  m_ncol = grid->get_num_local_dofs();
  m_nlay = grid->get_num_vertical_levels();

  // Set up dimension layouts
  FieldLayout scalar2d_layout     { {COL   }, {m_ncol    } };
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_ncol,m_nlay} };
  FieldLayout scalar3d_layout_int { {COL,LEV}, {m_ncol,m_nlay+1} };
  // Use VAR field tag for gases for now; consider adding a tag?
  FieldLayout gas_layout          { {COL,LEV,NGAS}, {m_ncol,m_nlay,m_ngas} };
  FieldLayout scalar2d_swband_layout { {COL,SWBND}, {m_ncol,m_nswbands} };

  // Set required (input) fields here
  add_field<Required>("p_mid" , scalar3d_layout_mid, Pa, grid->name());
  add_field<Required>("p_int", scalar3d_layout_int, Pa, grid->name());
  add_field<Required>("pseudo_density", scalar3d_layout_mid, Pa, grid->name());
  add_field<Required>("t_int" , scalar3d_layout_int, K , grid->name());
  add_field<Required>("gas_vmr", gas_layout, kgkg, grid->name());
  add_field<Required>("surf_alb_direct", scalar2d_swband_layout, nondim, grid->name());
  add_field<Required>("surf_alb_diffuse", scalar2d_swband_layout, nondim, grid->name());
  add_field<Required>("cos_zenith", scalar2d_layout, nondim, grid->name());
  add_field<Required>("qc", scalar3d_layout_mid, kg/kg, grid->name());
  add_field<Required>("qi", scalar3d_layout_mid, kg/kg, grid->name());
  add_field<Required>("cldfrac_tot", scalar3d_layout_mid, nondim, grid->name());
  add_field<Required>("eff_radius_qc", scalar3d_layout_mid, micron, grid->name());
  add_field<Required>("eff_radius_qi", scalar3d_layout_mid, micron, grid->name());

  // Set computed (output) fields
  add_field<Updated >("T_mid"     , scalar3d_layout_mid, K  , grid->name());
  add_field<Computed>("SW_flux_dn", scalar3d_layout_int, Wm2, grid->name());
  add_field<Computed>("SW_flux_up", scalar3d_layout_int, Wm2, grid->name());
  add_field<Computed>("SW_flux_dn_dir", scalar3d_layout_int, Wm2, grid->name());
  add_field<Computed>("LW_flux_up", scalar3d_layout_int, Wm2, grid->name());
  add_field<Computed>("LW_flux_dn", scalar3d_layout_int, Wm2, grid->name());

}  // RRTMGPRadiation::set_grids

int RRTMGPRadiation::requested_buffer_size_in_bytes() const
{
  const int interface_request = Buffer::num_1d_ncol*m_ncol*sizeof(Real) +
                                Buffer::num_2d_nlay*m_ncol*m_nlay*sizeof(Real) +
                                Buffer::num_2d_nlay_p1*m_ncol*(m_nlay+1)*sizeof(Real) +
                                Buffer::num_2d_nswbands*m_ncol*m_nswbands*sizeof(Real) +
                                Buffer::num_3d_ngas*m_ncol*m_nlay*m_ngas*sizeof(Real);

  return interface_request;
} // RRTMGPRadiation::requested_buffer_size

void RRTMGPRadiation::init_buffers(const ATMBufferManager &buffer_manager)
{
  EKAT_REQUIRE_MSG(buffer_manager.allocated_bytes() >= requested_buffer_size_in_bytes(), "Error! Buffers size not sufficient.\n");

  Real* mem = reinterpret_cast<Real*>(buffer_manager.get_memory());

  // 1d array
  m_buffer.mu0 = decltype(m_buffer.mu0)("mu0", mem, m_ncol);
  mem += m_buffer.mu0.totElems();

  // 2d arrays
  m_buffer.p_lay = decltype(m_buffer.p_lay)("p_lay", mem, m_ncol, m_nlay);
  mem += m_buffer.p_lay.totElems();
  m_buffer.t_lay = decltype(m_buffer.t_lay)("t_lay", mem, m_ncol, m_nlay);
  mem += m_buffer.t_lay.totElems();
  m_buffer.p_del = decltype(m_buffer.p_del)("p_del", mem, m_ncol, m_nlay);
  mem += m_buffer.p_del.totElems();
  m_buffer.qc = decltype(m_buffer.qc)("qc", mem, m_ncol, m_nlay);
  mem += m_buffer.qc.totElems();
  m_buffer.qi = decltype(m_buffer.qi)("qi", mem, m_ncol, m_nlay);
  mem += m_buffer.qi.totElems();
  m_buffer.cldfrac_tot = decltype(m_buffer.cldfrac_tot)("cldfrac_tot", mem, m_ncol, m_nlay);
  mem += m_buffer.cldfrac_tot.totElems();
  m_buffer.eff_radius_qc = decltype(m_buffer.eff_radius_qc)("eff_radius_qc", mem, m_ncol, m_nlay);
  mem += m_buffer.eff_radius_qc.totElems();
  m_buffer.eff_radius_qi = decltype(m_buffer.eff_radius_qi)("eff_radius_qi", mem, m_ncol, m_nlay);
  mem += m_buffer.eff_radius_qi.totElems();
  m_buffer.tmp2d = decltype(m_buffer.tmp2d)("tmp2d", mem, m_ncol, m_nlay);
  mem += m_buffer.tmp2d.totElems();
  m_buffer.lwp = decltype(m_buffer.lwp)("lwp", mem, m_ncol, m_nlay);
  mem += m_buffer.lwp.totElems();
  m_buffer.iwp = decltype(m_buffer.iwp)("iwp", mem, m_ncol, m_nlay);
  mem += m_buffer.iwp.totElems();
  m_buffer.sw_heating = decltype(m_buffer.sw_heating)("sw_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.sw_heating.totElems();
  m_buffer.lw_heating = decltype(m_buffer.lw_heating)("lw_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.lw_heating.totElems();
  m_buffer.rad_heating = decltype(m_buffer.rad_heating)("rad_heating", mem, m_ncol, m_nlay);
  mem += m_buffer.rad_heating.totElems();

  m_buffer.p_lev = decltype(m_buffer.p_lev)("p_lev", mem, m_ncol, m_nlay+1);
  mem += m_buffer.p_lev.totElems();
  m_buffer.t_lev = decltype(m_buffer.t_lev)("t_lev", mem, m_ncol, m_nlay+1);
  mem += m_buffer.t_lev.totElems();
  m_buffer.sw_flux_up = decltype(m_buffer.sw_flux_up)("sw_flux_up", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_up.totElems();
  m_buffer.sw_flux_dn = decltype(m_buffer.sw_flux_dn)("sw_flux_dn", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_dn.totElems();
  m_buffer.sw_flux_dn_dir = decltype(m_buffer.sw_flux_dn_dir)("sw_flux_dn_dir", mem, m_ncol, m_nlay+1);
  mem += m_buffer.sw_flux_dn_dir.totElems();
  m_buffer.lw_flux_up = decltype(m_buffer.lw_flux_up)("lw_flux_up", mem, m_ncol, m_nlay+1);
  mem += m_buffer.lw_flux_up.totElems();
  m_buffer.lw_flux_dn = decltype(m_buffer.lw_flux_dn)("lw_flux_dn", mem, m_ncol, m_nlay+1);
  mem += m_buffer.lw_flux_dn.totElems();

  m_buffer.sfc_alb_dir = decltype(m_buffer.sfc_alb_dir)("surf_alb_direct", mem, m_ncol, m_nswbands);
  mem += m_buffer.sfc_alb_dir.totElems();
  m_buffer.sfc_alb_dif = decltype(m_buffer.sfc_alb_dif)("surf_alb_diffuse", mem, m_ncol, m_nswbands);
  mem += m_buffer.sfc_alb_dif.totElems();

  // 3d array
  m_buffer.gas_vmr = decltype(m_buffer.gas_vmr)("gas_vmr", mem, m_ncol, m_nlay, m_ngas);
  mem += m_buffer.gas_vmr.totElems();

  int used_mem = (reinterpret_cast<Real*>(mem) - buffer_manager.get_memory())*sizeof(Real);
  EKAT_REQUIRE_MSG(used_mem==requested_buffer_size_in_bytes(), "Error! Used memory != requested memory for RRTMGPRadiation.");
} // RRTMGPRadiation::init_buffers

void RRTMGPRadiation::initialize_impl(const util::TimeStamp& /* t0 */) {
  // Names of active gases
  // TODO: this needs to be not hard-coded...I wanted to get these from
  // input files, but need to get around the rrtmgp_initializer logic.
  // Maybe can store the gas names somewhere else?
  string1d gas_names_1d("gas_names",m_ngas);
  for (int igas = 1; igas <= m_ngas; igas++) {
    gas_names_1d(igas) = m_gas_names[igas-1];
  }

  // Initialize GasConcs object to pass to RRTMGP initializer;
  // This is just to provide gas names
  // Make GasConcs from gas_vmr and gas_names_1d
  GasConcs gas_concs;
  gas_concs.init(gas_names_1d,1,1);
  rrtmgp::rrtmgp_initialize(gas_concs);
}

void RRTMGPRadiation::run_impl (const Real dt) {
  // Get data from the FieldManager
  auto d_pmid = m_rrtmgp_fields_in.at("p_mid").get_reshaped_view<const Real**>();
  auto d_pint = m_rrtmgp_fields_in.at("p_int").get_reshaped_view<const Real**>();
  auto d_pdel = m_rrtmgp_fields_in.at("pseudo_density").get_reshaped_view<const Real**>();
  auto d_tint = m_rrtmgp_fields_in.at("t_int").get_reshaped_view<const Real**>();
  auto d_gas_vmr = m_rrtmgp_fields_in.at("gas_vmr").get_reshaped_view<const Real***>();
  auto d_sfc_alb_dir = m_rrtmgp_fields_in.at("surf_alb_direct").get_reshaped_view<const Real**>();
  auto d_sfc_alb_dif = m_rrtmgp_fields_in.at("surf_alb_diffuse").get_reshaped_view<const Real**>();
  auto d_mu0 = m_rrtmgp_fields_in.at("cos_zenith").get_reshaped_view<const Real*>();
  auto d_qc = m_rrtmgp_fields_in.at("qc").get_reshaped_view<const Real**>();
  auto d_qi = m_rrtmgp_fields_in.at("qi").get_reshaped_view<const Real**>();
  auto d_cldfrac_tot = m_rrtmgp_fields_in.at("cldfrac_tot").get_reshaped_view<const Real**>();
  auto d_rel = m_rrtmgp_fields_in.at("eff_radius_qc").get_reshaped_view<const Real**>();
  auto d_rei = m_rrtmgp_fields_in.at("eff_radius_qi").get_reshaped_view<const Real**>();
  auto d_tmid = m_rrtmgp_fields_out.at("T_mid").get_reshaped_view<Real**>();
  auto d_sw_flux_up = m_rrtmgp_fields_out.at("SW_flux_up").get_reshaped_view<Real**>();
  auto d_sw_flux_dn = m_rrtmgp_fields_out.at("SW_flux_dn").get_reshaped_view<Real**>();
  auto d_sw_flux_dn_dir = m_rrtmgp_fields_out.at("SW_flux_dn_dir").get_reshaped_view<Real**>();
  auto d_lw_flux_up = m_rrtmgp_fields_out.at("LW_flux_up").get_reshaped_view<Real**>();
  auto d_lw_flux_dn = m_rrtmgp_fields_out.at("LW_flux_dn").get_reshaped_view<Real**>();

  // Create YAKL arrays. RRTMGP expects YAKL arrays with styleFortran, i.e., data has ncol
  // as the fastest index. For this reason we must copy the data.
  auto p_lay          = m_buffer.p_lay;
  auto t_lay          = m_buffer.t_lay;
  auto p_lev          = m_buffer.p_lev;
  auto p_del          = m_buffer.p_del;
  auto t_lev          = m_buffer.t_lev;
  auto gas_vmr        = m_buffer.gas_vmr;
  auto sfc_alb_dir    = m_buffer.sfc_alb_dir;
  auto sfc_alb_dif    = m_buffer.sfc_alb_dif;
  auto mu0            = m_buffer.mu0;
  auto qc             = m_buffer.qc;
  auto qi             = m_buffer.qi;
  auto cldfrac_tot    = m_buffer.cldfrac_tot;
  auto rel            = m_buffer.eff_radius_qc;
  auto rei            = m_buffer.eff_radius_qi;
  auto sw_flux_up     = m_buffer.sw_flux_up;
  auto sw_flux_dn     = m_buffer.sw_flux_dn;
  auto sw_flux_dn_dir = m_buffer.sw_flux_dn_dir;
  auto lw_flux_up     = m_buffer.lw_flux_up;
  auto lw_flux_dn     = m_buffer.lw_flux_dn;

  // Copy data from the FieldManager to the YAKL arrays
  {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();

      mu0(i+1) = d_mu0(i);
      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, m_nlay), [&] (const int& k) {
        p_lay(i+1,k+1)       = d_pmid(i,k);
        t_lay(i+1,k+1)       = d_tmid(i,k);
        p_del(i+1,k+1)       = d_pdel(i,k);
        qc(i+1,k+1)          = d_qc(i,k);
        qi(i+1,k+1)          = d_qi(i,k);
        cldfrac_tot(i+1,k+1) = d_cldfrac_tot(i,k);
        rel(i+1,k+1)         = d_rel(i,k);
        rei(i+1,k+1)         = d_rei(i,k);
        p_lev(i+1,k+1)       = d_pint(i,k);
        t_lev(i+1,k+1)       = d_tint(i,k);

        Kokkos::parallel_for(Kokkos::ThreadVectorRange(team, m_ngas), [&] (const int& g) {
          gas_vmr(i+1,k+1,g+1) = d_gas_vmr(i,k,g);
        });
      });

      p_lev(i+1,m_nlay+1) = d_pint(i,m_nlay);
      t_lev(i+1,m_nlay+1) = d_tint(i,m_nlay);

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, m_nswbands), [&] (const int& k) {
        sfc_alb_dir(i+1,k+1) = d_sfc_alb_dir(i,k);
        sfc_alb_dif(i+1,k+1) = d_sfc_alb_dif(i,k);
      });
    });
  }
  Kokkos::fence();
  
  // Make GasConcs from gas_vmr and gas_names
  // TODO: only allocate at initialization and
  // just update values here
  string1d gas_names_1d("gas_names",m_ngas);
  for (int igas = 1; igas <= m_ngas; igas++) {
    gas_names_1d(igas) = m_gas_names[igas-1];
  }

  // Create and populate a GasConcs object to pass to RRTMGP driver
  GasConcs gas_concs;
  gas_concs.init(gas_names_1d,m_ncol,m_nlay);
  auto tmp2d = m_buffer.tmp2d;
  for (int igas = 1; igas <= m_ngas; igas++) {
    parallel_for(Bounds<2>(m_ncol,m_nlay), YAKL_LAMBDA(int icol, int ilay) {
        tmp2d(icol,ilay) = gas_vmr(icol,ilay,igas);
    });
    gas_concs.set_vmr(gas_names_1d(igas), tmp2d);
  }

  // Compute layer cloud mass (per unit area)
  auto lwp = m_buffer.lwp;
  auto iwp = m_buffer.iwp;
  scream::rrtmgp::mixing_ratio_to_cloud_mass(qc, cldfrac_tot, p_del, lwp);
  scream::rrtmgp::mixing_ratio_to_cloud_mass(qi, cldfrac_tot, p_del, iwp);
  // Convert to g/m2 (needed by RRTMGP)
  parallel_for(Bounds<2>(m_nlay,m_ncol), YAKL_LAMBDA(int ilay, int icol) {
      lwp(icol,ilay) = 1e3 * lwp(icol,ilay);
      iwp(icol,ilay) = 1e3 * iwp(icol,ilay);
  });

  // Run RRTMGP driver
  rrtmgp::rrtmgp_main(
    m_ncol, m_nlay,
    p_lay, t_lay, p_lev, t_lev,
    gas_concs,
    sfc_alb_dir, sfc_alb_dif, mu0,
    lwp, iwp, rel, rei,
    sw_flux_up, sw_flux_dn, sw_flux_dn_dir,
    lw_flux_up, lw_flux_dn
  );

  // Compute and apply heating rates
  auto sw_heating  = m_buffer.sw_heating;
  auto lw_heating  = m_buffer.lw_heating;
  auto rad_heating = m_buffer.rad_heating;
  rrtmgp::compute_heating_rate(
    sw_flux_up, sw_flux_dn, p_del, sw_heating
  );
  rrtmgp::compute_heating_rate(
    lw_flux_up, lw_flux_dn, p_del, lw_heating
  );
  parallel_for(Bounds<2>(m_nlay,m_ncol), YAKL_LAMBDA(int ilay, int icol) {
    rad_heating(icol,ilay) = sw_heating(icol,ilay) + lw_heating(icol,ilay);
    t_lay(icol,ilay) = t_lay(icol,ilay) + rad_heating(icol,ilay) * dt;
  });

  // Copy ouput data back to FieldManager
  {
    const auto policy = ekat::ExeSpaceUtils<ExeSpace>::get_default_team_policy(m_ncol, m_nlay);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
      const int i = team.league_rank();

      Kokkos::parallel_for(Kokkos::TeamThreadRange(team, m_nlay+1), [&] (const int& k) {
        if (k < m_nlay) d_tmid(i,k) = t_lay(i+1,k+1);

        d_sw_flux_up(i,k)     = sw_flux_up(i+1,k+1);
        d_sw_flux_dn(i,k)     = sw_flux_dn(i+1,k+1);
        d_sw_flux_dn_dir(i,k) = sw_flux_dn_dir(i+1,k+1);
        d_lw_flux_up(i,k)     = lw_flux_up(i+1,k+1);
        d_lw_flux_dn(i,k)     = lw_flux_dn(i+1,k+1);
      });
    });
  }
}

void RRTMGPRadiation::finalize_impl  () {
  rrtmgp::rrtmgp_finalize();
}

void RRTMGPRadiation::set_required_field_impl(const Field<const Real>& f) {
  const auto& name = f.get_header().get_identifier().name();
  m_rrtmgp_fields_in.emplace(name,f);
  m_rrtmgp_host_views_in[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_in[name] = m_rrtmgp_host_views_in[name].data();

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void RRTMGPRadiation::set_computed_field_impl(const Field<      Real>& f) {
  const auto& name = f.get_header().get_identifier().name();
  m_rrtmgp_fields_out.emplace(name,f);
  m_rrtmgp_host_views_out[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_out[name] = m_rrtmgp_host_views_out[name].data();

  // Add myself as provider for the field
  add_me_as_provider(f);
}

}  // namespace scream
