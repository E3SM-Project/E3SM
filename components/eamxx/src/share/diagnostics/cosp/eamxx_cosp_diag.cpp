#include "eamxx_cosp_diag.hpp"
#include "cosp_functions.hpp"

#include "share/physics/physics_constants.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"
#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_assert.hpp>
#include <ekat_units.hpp>

namespace scream
{

// =========================================================================================
std::vector<std::string> CospDiag::output_names ()
{
  return {"isccp_cldtot", "isccp_ctptau", "modis_ctptau", "misr_cthtau"};
}

// =========================================================================================
CospDiag::
CospDiag (const ekat::Comm& comm, const ekat::ParameterList& params,
          const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  // How many subcolumns to use for COSP. The process took this from its yaml
  // section; a diag is built from an expression, which carries no such knob, so
  // in practice this is the default. It is still read from the params, so that
  // it works if a caller ever does supply them.
  m_num_subcols = m_params.get<Int>("cosp_subcolumns", 10);

  m_num_cols = m_grid->get_num_local_dofs();
  m_num_levs = m_grid->get_num_vertical_levels();

  // Same inputs the process declared as Required. They are resolved against the
  // FieldManager (or built as diags in their own right) by whoever wires us up.
  for (const auto& n : {"surf_radiative_T", "sunlit_mask", "p_mid", "p_int",
                        "T_mid", "phis", "pseudo_density", "cldfrac_rad",
                        "qv", "qc", "qi", "dtau067", "dtau105",
                        "eff_radius_qc", "eff_radius_qi"}) {
    m_field_in_names.push_back(n);
  }
}

// =========================================================================================
CospDiag::~CospDiag ()
{
  // Mirrors the process' finalize_impl. Only if we got as far as initializing,
  // since finalizing the Fortran side twice (or without an init) is not safe.
  if (m_is_initialized) {
    CospFunc::finalize();
  }
}

// =========================================================================================
void CospDiag::initialize_impl ()
{
  using namespace ekat::units;
  using namespace ekat::prefixes;
  using namespace ShortFieldTagsNames;

  CospFunc::initialize(m_num_cols, m_num_subcols, m_num_levs);

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto percent = none.rename("%");

  const auto& grid_name = m_grid->name();

  FieldLayout scalar2d     = m_grid->get_2d_scalar_layout();
  FieldLayout scalar3d_mid = m_grid->get_3d_scalar_layout(LEV);
  FieldLayout scalar3d_int = m_grid->get_3d_scalar_layout(ILEV);
  FieldLayout scalar4d_ctptau ( {COL,CMP,CMP},
                                {m_num_cols,m_num_tau,m_num_ctp},
                                {e2str(COL), "cosp_tau", "cosp_prs"});
  FieldLayout scalar4d_cthtau ( {COL,CMP,CMP},
                                {m_num_cols,m_num_tau,m_num_cth},
                                {e2str(COL), "cosp_tau", "cosp_cth"});

  // The outputs. The first one declared is the primary, which is what a caller
  // asking for the diag rather than for a specific field gets.
  const std::map<std::string,FieldLayout> out_layouts = {
    {"isccp_cldtot", scalar2d},
    {"isccp_ctptau", scalar4d_ctptau},
    {"modis_ctptau", scalar4d_ctptau},
    {"misr_cthtau" , scalar4d_cthtau},
  };
  for (const auto& n : output_names()) {
    Field f (FieldIdentifier(n,out_layouts.at(n),percent,grid_name));
    f.allocate_view();
    add_output(f);
  }

  // Masks for the 4d fields: the sunlit mask, broadcast over tau/prs and
  // tau/cth. They are filled in compute_impl, once sunlit_mask has a value.
  FieldIdentifier mctp_fid ("sunlit_mask_ctptau", scalar4d_ctptau, none, grid_name, DataType::IntType);
  FieldIdentifier mcth_fid ("sunlit_mask_cthtau", scalar4d_cthtau, none, grid_name, DataType::IntType);
  Field mctp(mctp_fid,true);
  Field mcth(mcth_fid,true);
  const std::map<std::string,Field> masks = {
    {"isccp_cldtot", m_fields_in.at("sunlit_mask")},
    {"isccp_ctptau", mctp},
    {"modis_ctptau", mctp},
    {"misr_cthtau" , mcth},
  };
  for (const auto& n : output_names()) {
    auto f = get(n);
    f.set_valid_mask(masks.at(n));
    // Night columns are fill_value, so averaging must track a per-cell count
    // rather than divide by nsteps.
    f.get_header().set_may_be_filled(true);
  }

  // Temporaries
  m_z_mid = Field(FieldIdentifier("z_mid",scalar3d_mid,m,grid_name),true);
  m_z_int = Field(FieldIdentifier("z_int",scalar3d_int,m,grid_name),true);
  m_sunlit_real = Field(FieldIdentifier("sunlit_mask_real",scalar2d,none,grid_name),true);

  create_cosp_geometry_data();
}

// =========================================================================================
void CospDiag::create_cosp_geometry_data () const
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The grid may already carry these: the diag can be rebuilt, and (during the
  // transition) the Cosp process may have defined them first. Defining them
  // twice throws, so ask.
  if (m_grid->has_geometry_data("cosp_tau")) {
    return;
  }

  // Layouts for 1D coordinate variables
  FieldLayout cosp_tau_layout ({CMP}, {m_num_tau}, {"cosp_tau"});
  FieldLayout cosp_prs_layout ({CMP}, {m_num_ctp}, {"cosp_prs"});
  FieldLayout cosp_cth_layout ({CMP}, {m_num_cth}, {"cosp_cth"});
  // Layouts for 2D bounds variables (dim x 2)
  FieldLayout cosp_tau_bnds_layout ({CMP,CMP}, {m_num_tau,2}, {"cosp_tau","nbnd"});
  FieldLayout cosp_prs_bnds_layout ({CMP,CMP}, {m_num_ctp,2}, {"cosp_prs","nbnd"});
  FieldLayout cosp_cth_bnds_layout ({CMP,CMP}, {m_num_cth,2}, {"cosp_cth","nbnd"});

  auto cosp_tau_f      = m_grid->create_geometry_data("cosp_tau",      cosp_tau_layout,     none);
  auto cosp_prs_f      = m_grid->create_geometry_data("cosp_prs",      cosp_prs_layout,     Pa);
  auto cosp_cth_f      = m_grid->create_geometry_data("cosp_cth",      cosp_cth_layout,     m);
  auto cosp_tau_bnds_f = m_grid->create_geometry_data("cosp_tau_bnds", cosp_tau_bnds_layout, none);
  auto cosp_prs_bnds_f = m_grid->create_geometry_data("cosp_prs_bnds", cosp_prs_bnds_layout, Pa);
  auto cosp_cth_bnds_f = m_grid->create_geometry_data("cosp_cth_bnds", cosp_cth_bnds_layout, m);

  // In IO, make the output of this geo data conditional to vars using it being present in the stream
  cosp_tau_f.get_header().set_extra_data("io_output_if_dim_exists", std::string("cosp_tau"));
  cosp_prs_f.get_header().set_extra_data("io_output_if_dim_exists", std::string("cosp_prs"));
  cosp_cth_f.get_header().set_extra_data("io_output_if_dim_exists", std::string("cosp_cth"));
  cosp_tau_bnds_f.get_header().set_extra_data("io_output_if_dim_exists", std::string("cosp_tau"));
  cosp_prs_bnds_f.get_header().set_extra_data("io_output_if_dim_exists", std::string("cosp_prs"));
  cosp_cth_bnds_f.get_header().set_extra_data("io_output_if_dim_exists", std::string("cosp_cth"));

  // Retrieve bin centers and edges from the Fortran COSP interface (mod_cosp_config).
  // The F90 arrays for edges have shape (2, nbins) in Fortran column-major order, so the
  // flat memory layout is [lower_0, upper_0, lower_1, upper_1, ...], which maps directly
  // onto the (nbins, 2) host views below.
  auto tau_h      = cosp_tau_f.get_view<Real*,Host>();
  auto tau_bnds_h = cosp_tau_bnds_f.get_view<Real**,Host>();
  auto prs_h      = cosp_prs_f.get_view<Real*,Host>();
  auto prs_bnds_h = cosp_prs_bnds_f.get_view<Real**,Host>();
  auto cth_h      = cosp_cth_f.get_view<Real*,Host>();
  auto cth_bnds_h = cosp_cth_bnds_f.get_view<Real**,Host>();

  CospFunc::get_bins(m_num_tau, m_num_ctp, m_num_cth,
                         tau_h.data(), tau_bnds_h.data(),
                         prs_h.data(), prs_bnds_h.data(),
                         cth_h.data(), cth_bnds_h.data());

  cosp_tau_f.sync_to_dev();
  cosp_tau_bnds_f.sync_to_dev();
  cosp_prs_f.sync_to_dev();
  cosp_prs_bnds_f.sync_to_dev();
  cosp_cth_f.sync_to_dev();
  cosp_cth_bnds_f.sync_to_dev();

  // Set "bounds" CF attribute on each coordinate variable to link to its bounds variable
  using stratts_t = std::map<std::string,std::string>;
  cosp_tau_f.get_header().get_extra_data<stratts_t>("io: string attributes")["bounds"] = "cosp_tau_bnds";
  cosp_prs_f.get_header().get_extra_data<stratts_t>("io: string attributes")["bounds"] = "cosp_prs_bnds";
  cosp_cth_f.get_header().get_extra_data<stratts_t>("io: string attributes")["bounds"] = "cosp_cth_bnds";
}

// =========================================================================================
void CospDiag::compute_impl ()
{
  // Get fields from the input map; note that we get host views because this
  // interface serves primarily as a wrapper to a c++ to f90 bridge for the COSP
  // all then need to be copied to layoutLeft views to permute the indices for
  // F90.

  // Ensure host data of input fields is up to date
  for (const auto& n : {"qv", "qc", "qi", "surf_radiative_T", "T_mid", "p_mid",
                        "p_int", "cldfrac_rad", "eff_radius_qc",
                        "eff_radius_qi", "dtau067", "dtau105"}) {
    m_fields_in.at(n).sync_to_host();
  }

  m_sunlit_real.deep_copy(m_fields_in.at("sunlit_mask"));
  m_sunlit_real.sync_to_host();

  // Compute z_mid
  const auto T_mid_d = m_fields_in.at("T_mid").get_view<const Real**>();
  const auto qv_d    = m_fields_in.at("qv").get_view<const Real**>();
  const auto p_mid_d = m_fields_in.at("p_mid").get_view<const Real**>();
  const auto phis_d  = m_fields_in.at("phis").get_view<const Real*>();
  const auto pseudo_density_d = m_fields_in.at("pseudo_density").get_view<const Real**>();
  const auto z_mid_d = m_z_mid.get_view<Real**>();
  const auto z_int_d = m_z_int.get_view<Real**>();
  const auto ncol = m_num_cols;
  const auto nlev = m_num_levs;

  using KT       = KokkosTypes<DefaultDevice>;
  using ExeSpace = typename KT::ExeSpace;
  using TPF      = ekat::TeamPolicyFactory<ExeSpace>;
  using PF       = scream::PhysicsFunctions<DefaultDevice>;

  const auto scan_policy = TPF::get_thread_range_parallel_scan_team_policy(ncol, nlev);
  const Real g = physics::Constants<Real>::gravit.value;
  Kokkos::parallel_for(scan_policy, KOKKOS_LAMBDA (const KT::MemberType& team) {
      const int i = team.league_rank();
      const auto p_mid_s = ekat::subview(p_mid_d, i);
      const auto T_mid_s = ekat::subview(T_mid_d, i);
      const auto qv_s = ekat::subview(qv_d, i);
      const auto z_int_s = ekat::subview(z_int_d, i);
      const auto z_mid_s = ekat::subview(z_mid_d, i);
      const Real z_surf  = phis_d(i) / g;
      const auto pseudo_density_s = ekat::subview(pseudo_density_d, i);

      // 1. Compute dz (recycle z_mid_s as a temporary)
      const auto dz_s = z_mid_s; //
      PF::calculate_dz(team, pseudo_density_s, p_mid_s, T_mid_s, qv_s, dz_s);
      team.team_barrier();

      // 2. Compute z_int (vertical scan)
      PF::calculate_z_int(team,nlev,dz_s,z_surf,z_int_s);
      team.team_barrier();

      // 3. Compute z_mid (int->mid interpolation)
      PF::calculate_z_mid(team,nlev,z_int_s,z_mid_s);
      team.team_barrier();
  });
  Kokkos::fence();

  m_z_mid.sync_to_host();
  const auto z_mid_h   = m_z_mid.get_view<const Real**,Host>();
  const auto T_mid_h   = m_fields_in.at("T_mid").get_view<const Real**, Host>();
  const auto qv_h      = m_fields_in.at("qv").get_view<const Real**, Host>();
  const auto p_mid_h   = m_fields_in.at("p_mid").get_view<const Real**,Host>();
  const auto qc_h      = m_fields_in.at("qc").get_view<const Real**, Host>();
  const auto qi_h      = m_fields_in.at("qi").get_view<const Real**, Host>();
  const auto sunlit_h  = m_sunlit_real.get_view<const Real*, Host>();
  const auto skt_h     = m_fields_in.at("surf_radiative_T").get_view<const Real*, Host>();
  const auto p_int_h   = m_fields_in.at("p_int").get_view<const Real**, Host>();
  const auto cldfrac_h = m_fields_in.at("cldfrac_rad").get_view<const Real**, Host>();
  const auto reff_qc_h = m_fields_in.at("eff_radius_qc").get_view<const Real**, Host>();
  const auto reff_qi_h = m_fields_in.at("eff_radius_qi").get_view<const Real**, Host>();
  const auto dtau067_h = m_fields_in.at("dtau067").get_view<const Real**, Host>();
  const auto dtau105_h = m_fields_in.at("dtau105").get_view<const Real**, Host>();

  auto isccp_cldtot = get("isccp_cldtot");
  auto isccp_ctptau = get("isccp_ctptau");
  auto modis_ctptau = get("modis_ctptau");
  auto misr_cthtau  = get("misr_cthtau");

  auto isccp_cldtot_h = isccp_cldtot.get_view<Real*, Host>();
  auto isccp_ctptau_h = isccp_ctptau.get_view<Real***, Host>();
  auto modis_ctptau_h = modis_ctptau.get_view<Real***, Host>();
  auto misr_cthtau_h  = misr_cthtau. get_view<Real***, Host>();

  Real emsfc_lw = 0.99;
  CospFunc::main(
          m_num_cols, m_num_subcols, m_num_levs, m_num_tau, m_num_ctp, m_num_cth, emsfc_lw,
          sunlit_h, skt_h, T_mid_h, p_mid_h, p_int_h, z_mid_h, qv_h, qc_h, qi_h,
          cldfrac_h, reff_qc_h, reff_qi_h, dtau067_h, dtau105_h,
          isccp_cldtot_h, isccp_ctptau_h, modis_ctptau_h, misr_cthtau_h
  );

  // Mask night values
  constexpr auto fill_value = constants::fill_value<Real>;
  for (int i = 0; i < m_num_cols; i++) {
    if (sunlit_h(i) == 0) {
      // if night, set to fill val
      isccp_cldtot_h(i) = fill_value;
      for (int j = 0; j < m_num_tau; j++) {
        for (int k = 0; k < m_num_ctp; k++) {
          isccp_ctptau_h(i,j,k) = fill_value;
          modis_ctptau_h(i,j,k) = fill_value;
        }
        for (int k = 0; k < m_num_cth; k++) {
          misr_cthtau_h (i,j,k) = fill_value;
        }
      }
    }
  }

  // Make sure dev data is up to date
  isccp_cldtot.sync_to_dev();
  isccp_ctptau.sync_to_dev();
  modis_ctptau.sync_to_dev();
  misr_cthtau.sync_to_dev();

  // Update the ctptau and cthtau masks by broadcasting sunlit mask
  const auto& sunlit = m_fields_in.at("sunlit_mask");
  auto ctptau = isccp_ctptau.get_valid_mask();
  auto cthtau = misr_cthtau.get_valid_mask();

  auto sunlit_v = sunlit.get_view<const int*>();
  auto ctptau_v = ctptau.get_view<int***>();
  auto cthtau_v = cthtau.get_view<int***>();
  auto do_ctp = KOKKOS_LAMBDA (int icol, int itau, int ictp) {
    ctptau_v(icol,itau,ictp) = sunlit_v(icol);
  };
  auto do_cth = KOKKOS_LAMBDA (int icol, int itau, int icth) {
    cthtau_v(icol,itau,icth) = sunlit_v(icol);
  };
  using exec_space = typename DefaultDevice::execution_space;
  using policy_t = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<3>>;
  policy_t policy_ctp({0,0,0},{m_num_cols,m_num_tau,m_num_ctp});
  policy_t policy_cth({0,0,0},{m_num_cols,m_num_tau,m_num_cth});
  Kokkos::parallel_for(policy_ctp,do_ctp);
  Kokkos::parallel_for(policy_cth,do_cth);
}

} // namespace scream
