#include "physics/shoc/shoc_inputs_initializer.hpp"
#include "physics/shoc/shoc_main_impl.hpp"
#include "ekat/util/ekat_file_utils.hpp"
#include "physics/share/physics_functions.hpp"

#include <array>
#include <fstream>

namespace scream
{

void SHOCInputsInitializer::add_field (const field_type &f)
{
  const auto& id = f.get_header().get_identifier();
  
  m_fields.emplace(id.name(),f);
  m_fields_id.insert(id);
}

void SHOCInputsInitializer::
add_field (const field_type &f, const field_type& f_ref,
           const remapper_ptr_type& remapper)
{
  if (m_remapper) {
    // Sanity check
    EKAT_REQUIRE_MSG (m_remapper->get_src_grid()->name()==remapper->get_src_grid()->name(),
      "Error! A remapper was already set in SHOCInputsInitializer, but its src grid differs from"
      "       the grid of the input remapper of this call.\n");
  } else {
    m_remapper = remapper;
    m_remapper->registration_begins();
  }

  const auto& id = f.get_header().get_identifier();
  const auto& id_ref = f_ref.get_header().get_identifier();

  // To the AD, we only expose the fact that we init f_ref...
  m_fields_id.insert(id_ref);

  // ...but SHOC only knows how to init f...
  m_fields.emplace(id.name(),f);

  // ...hence, we remap to f_ref.
  m_remapper->register_field(f, f_ref);
}

// =========================================================================================
void SHOCInputsInitializer::initialize_fields ()
{
  using namespace shoc;
  using SHF          = Functions<Real, DefaultDevice>;
  using C            = scream::physics::Constants<Real>;
  using PhysicsFun   = scream::physics::Functions<Real, DefaultDevice>;
  using Scalar       = typename SHF::Scalar;
  using Spack        = typename SHF::Spack;
  using Smask        = typename SHF::Smask;
  using Pack1d       = typename SHF::Pack1d;
  using IntSmallPack = typename SHF::IntSmallPack;
  using view_2d      = typename SHF::view_2d<Spack>;
  using view_1d_scalar      = typename SHF::view_1d<Real>;
  using MemberType  = typename SHF::MemberType;

  const Real cpair = C::Cpair;
  const Real ggr = C::gravit;
  const Real zvir = C::ZVIR;
  const Real latvap = C::LatVap;
  const Real p0 = C::P0;
  const Real rair = C::Rair;

  // Safety check: if we're asked to init anything at all,
  // To simplify the initializer we first define all the fields we expect to have to initialize.
  std::vector<std::string> fields_to_init;
  fields_to_init.push_back("pref_mid");
  fields_to_init.push_back("t");
  fields_to_init.push_back("alst");
  fields_to_init.push_back("zi");
  fields_to_init.push_back("zm");
  fields_to_init.push_back("omega");
  fields_to_init.push_back("shf");
  fields_to_init.push_back("cflx_k0");
  fields_to_init.push_back("wsx");
  fields_to_init.push_back("wsy");
  fields_to_init.push_back("shoc_qv");
  fields_to_init.push_back("host_dx");
  fields_to_init.push_back("host_dy");
  fields_to_init.push_back("pmid");
  fields_to_init.push_back("pint");
  fields_to_init.push_back("pdel");
  fields_to_init.push_back("phis");
  fields_to_init.push_back("s");
  fields_to_init.push_back("tke");
  fields_to_init.push_back("u");
  fields_to_init.push_back("v");
  fields_to_init.push_back("Q");
  fields_to_init.push_back("wthv_sec");
  fields_to_init.push_back("tkh");
  fields_to_init.push_back("tk");
  fields_to_init.push_back("shoc_ql");
  int count = 0;
  std::string list_of_fields = "";
  for (auto name : fields_to_init)
  {
    list_of_fields += name;
    list_of_fields += ", ";
    count += m_fields.count(name);
  }
 
  EKAT_REQUIRE_MSG(count!=0,"Error in shoc_inputs_initializer: no fields have declared this initializer.  Check shoc interface.");

  EKAT_REQUIRE_MSG (count==(int)fields_to_init.size(),
    "Error! SHOCInputsInitializer is expected to init " + std::to_string(fields_to_init.size()) + " fields:\n"
    "       " + list_of_fields + "\n"
    "       Instead found " + std::to_string(count) + " fields.\n"
    "       Please, check the atmosphere processes you are using,\n"
    "       and make sure they agree on who's initializing each field.\n");

  // Initialize the fields that we expect.
  // Get device views
  auto d_pref_mid = m_fields.at("pref_mid").get_reshaped_view<Spack*>();
  auto d_t        = m_fields.at("t").get_reshaped_view<Spack**>();
  auto d_alst     = m_fields.at("alst").get_reshaped_view<Spack**>();
  auto d_zi       = m_fields.at("zi").get_reshaped_view<Spack**>();
  auto d_zm       = m_fields.at("zm").get_reshaped_view<Spack**>();
  auto d_omega    = m_fields.at("omega").get_reshaped_view<Spack**>();
  auto d_shf      = m_fields.at("shf").get_reshaped_view<Pack1d*>();
  auto d_cflx_k0  = m_fields.at("cflx_k0").get_reshaped_view<Pack1d*>();
  auto d_wsx      = m_fields.at("wsx").get_reshaped_view<Pack1d*>();
  auto d_wsy      = m_fields.at("wsy").get_reshaped_view<Pack1d*>();
  auto d_shoc_qv  = m_fields.at("shoc_qv").get_reshaped_view<Spack**>();
  auto d_host_dx  = m_fields.at("host_dx").get_reshaped_view<Pack1d*>();
  auto d_host_dy  = m_fields.at("host_dy").get_reshaped_view<Pack1d*>();
  auto d_pmid     = m_fields.at("pmid").get_reshaped_view<Spack**>();
  auto d_pint     = m_fields.at("pint").get_reshaped_view<Spack**>();
  auto d_pdel     = m_fields.at("pdel").get_reshaped_view<Spack**>();
  auto d_phis     = m_fields.at("phis").get_reshaped_view<Pack1d*>();
  auto d_s        = m_fields.at("s").get_reshaped_view<Spack**>();
  auto d_tke      = m_fields.at("tke").get_reshaped_view<Spack**>();
  auto d_u        = m_fields.at("u").get_reshaped_view<Spack**>();
  auto d_v        = m_fields.at("v").get_reshaped_view<Spack**>();
  auto d_Q        = m_fields.at("Q").get_reshaped_view<Spack***>();
  auto d_wthv_sec = m_fields.at("wthv_sec").get_reshaped_view<Spack**>();
  auto d_tkh      = m_fields.at("tkh").get_reshaped_view<Spack**>();
  auto d_tk       = m_fields.at("tk").get_reshaped_view<Spack**>();
  auto d_shoc_ql  = m_fields.at("shoc_ql").get_reshaped_view<Spack**>();

  // Reference elevations for data interpolation (bottom-to-top).
  const std::vector<Real> z_ref{0.0, 520.0, 1480.0, 2000.0, 3000.0};
  const std::vector<Real> qw_ref{sp(1e-3*17), sp(1e-3*16.3), sp(1e-3*10.7), sp(1e-3*4.2), sp(1e-3*3)};
  const std::vector<Real> ql_ref{0.0, sp(1e-3*5), sp(1e-3*7), sp(1e-3*6), 0.0};
  const std::vector<Real> theta_ref{sp(299.7), sp(298.7), sp(302.4), sp(308.2), sp(312.85)};

  // Wind speed interpolation data.
  const std::vector<Real> wind_z_ref{0.0, 700.0, 3000.0};
  const std::vector<Real> u_ref{-7.75, -8.75, -4.61};
  const std::vector<Real> v_ref{0.0, 0.1, 0.0};
  const std::vector<Real> w_ref{0.1, 0.1, 0.0};

  // Copy interpolation data to views
  view_1d_scalar
    d_z_ref("z_ref",5),
    d_qw_ref("qw_ref",5),
    d_ql_ref("ql_ref",5),
    d_theta_ref("theta_ref",5),
    d_wind_z_ref("wind_z_ref",3),
    d_u_ref("u_ref",3),
    d_v_ref("v_ref",3),
    d_w_ref("w_ref",3);

  auto h_z_ref = Kokkos::create_mirror_view(d_z_ref);
  auto h_qw_ref = Kokkos::create_mirror_view(d_qw_ref);
  auto h_ql_ref = Kokkos::create_mirror_view(d_ql_ref);
  auto h_theta_ref = Kokkos::create_mirror_view(d_theta_ref);
  auto h_wind_z_ref = Kokkos::create_mirror_view(d_wind_z_ref);
  auto h_u_ref = Kokkos::create_mirror_view(d_u_ref);
  auto h_v_ref = Kokkos::create_mirror_view(d_v_ref);
  auto h_w_ref = Kokkos::create_mirror_view(d_w_ref);

  for (int i=0; i<5; ++i) {
    if (i < 3) {
      h_wind_z_ref(i) = wind_z_ref[i];

      h_u_ref(i) = u_ref[i];
      h_v_ref(i) = v_ref[i];
      h_w_ref(i) = w_ref[i];
    }
    h_z_ref(i) = z_ref[i];
    h_qw_ref(i) = qw_ref[i];
    h_ql_ref(i) = ql_ref[i];
    h_theta_ref(i) = theta_ref[i];
  }

  Kokkos::deep_copy(d_z_ref,h_z_ref);
  Kokkos::deep_copy(d_qw_ref,h_qw_ref);
  Kokkos::deep_copy(d_ql_ref,h_ql_ref);
  Kokkos::deep_copy(d_theta_ref,h_theta_ref);
  Kokkos::deep_copy(d_wind_z_ref,h_wind_z_ref);
  Kokkos::deep_copy(d_u_ref,h_u_ref);
  Kokkos::deep_copy(d_v_ref,h_v_ref);
  Kokkos::deep_copy(d_w_ref,h_w_ref);

  auto mdims = m_fields.at("t").get_header().get_identifier().get_layout();
  int ncol = mdims.dim(0);
  int nlev   = mdims.dim(1);
  int nlevi = m_fields.at("zi").get_header().get_identifier().get_layout().dim(1);

  const Int nlev_packs = ekat::npack<Spack>(nlev);
  const Int nlevi_packs = ekat::npack<Spack>(nlevi);

  view_2d rrho("rrho", ncol, nlev_packs);
  view_2d rrhoi("rrhoi", ncol, nlevi_packs);
  view_2d qw("qw", ncol, nlev_packs);


  const int nlev_v = (nlev-1)/Spack::n;
  const int nlev_p = (nlev-1)%Spack::n;
  const int nlev_p1_v = nlev/Spack::n;
  const int nlev_p1_p = nlev%Spack::n;

  // Initalize inputs. Values are taken from shoc_ic_cases.cpp. In some cases,
  // input values here are reversed engineered from shoc_ic_cases.cpp values
  const Scalar ztop = 2400.0;
  const Scalar dz = ztop/nlev;

  const auto policy = ekat::ExeSpaceUtils<>::get_default_team_policy(ncol, nlev_packs);
  Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const MemberType& team) {
    const Int i = team.league_rank();

    // host_dx/y
    d_host_dx(i)[0] = 5300;
    d_host_dy(i)[0] = 5300;

    // phis
    d_phis(i)[0] = 0;

    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {

      const auto indices = ekat::range<IntSmallPack>(k*Spack::n);
      const auto inv_indices = (nlev-1) - indices;

      // Some inputs initialize to 0
      d_alst(i,k).set(indices < nlev, 0);
      d_wthv_sec(i,k).set(indices < nlev, 0);
      d_tkh(i,k).set(indices < nlev, 0);
      d_tk(i,k).set(indices < nlev, 0);
      d_tke(i,k).set(indices < nlev, 0);

      // Set zi
      d_zi(i,k).set(inv_indices >= 0 && indices < nlev, (Spack)(inv_indices+1)*dz);
    });
    team.team_barrier();
    d_zi(i,nlev_p1_v)[nlev_p1_p] = 0;

    const auto sub_zi= ekat::subview(d_zi, i);
    const auto s_zi = ekat::scalarize(sub_zi);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {
      // Set zm (zt in run_and_cmp test
      Spack zi_k, zi_kp1;
      auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
      auto range_pack2 = range_pack1;

      // Original code was: range_pack2.set(range_pack1 > nlev, 1); but that caused mysterious test
      // failures on blake.
      vector_simd
      for (int s = 0; s < Spack::n; ++s) {
        if (range_pack1[s] > nlev) {
          range_pack2[s] = 1;
        }
      }
      ekat::index_and_shift<1>(s_zi, range_pack2, zi_k, zi_kp1);

      d_zm(i,k).set(range_pack1 < nlev, 0.5*(zi_k + zi_kp1));

      // Set cell-centered wind speeds and compute shoc_ql and shoc_qv
      for (int p=0; p<Spack::n && k*Spack::n+p<nlev; ++p) {
        d_u(i, k)[p] = interpolate_data(3, d_wind_z_ref, d_u_ref, d_zm(i,k)[p]);
        d_v(i, k)[p] = interpolate_data(3, d_wind_z_ref, d_v_ref, d_zm(i,k)[p]);

        d_shoc_ql(i,k)[p] = interpolate_data(5, d_z_ref, d_ql_ref, d_zm(i,k)[p]);
        qw(i,k)[p] = interpolate_data(5, d_z_ref, d_qw_ref, d_zm(i,k)[p]);
      }
      d_shoc_qv(i,k).set(range_pack1 < nlev, qw(i,k) - d_shoc_ql(i,k));

      // Set tracer array. We only consider the 3 tracers
      // required elsewhere in the code.
      for (int p=0; p<Spack::n && k*Spack::n+p < nlev; ++p) {
        d_Q(i,0,k)[p] = d_shoc_ql(i,k)[p];
        d_Q(i,1,k)[p] = d_shoc_qv(i,k)[p];
        d_Q(i,2,k)[p] = d_tke(i,k)[p];
      }
    });
    team.team_barrier();

    // Inputs related to pressure and temperature
    compute_column_pressure(i,nlev,rair,cpair,ggr,p0,d_z_ref,d_theta_ref,d_zm,d_pmid);
    compute_column_pressure(i,nlevi,rair,cpair,ggr,p0,d_z_ref,d_theta_ref,d_zi,d_pint);

    team.team_barrier();
    const auto sub_pint = ekat::subview(d_pint, i);
    const auto s_pint = ekat::scalarize(sub_pint);
    Kokkos::parallel_for(Kokkos::TeamThreadRange(team, nlev_packs), [&] (const Int& k) {

      if (i==0) {
        d_pref_mid(k) = d_pmid(i,k);
      }

      Spack pint_k, pint_kp1;
      auto range_pack1 = ekat::range<IntSmallPack>(k*Spack::n);
      auto range_pack2 = range_pack1;
      range_pack2.set(range_pack1 > nlev, 1);
      ekat::index_and_shift<1>(s_pint, range_pack2, pint_k, pint_kp1);
      d_pdel(i,k).set(range_pack1<nlev, abs(pint_k - pint_kp1));

      rrho(i,k).set(range_pack1<nlev, (1/ggr)*(d_pdel(i,k)/dz));

      Spack theta_zt;
      Spack w_field;
      for (int p=0; p<Spack::n && k*Spack::n+p<nlev; ++p) {
        theta_zt[p] = interpolate_data(5, d_z_ref, d_theta_ref, d_zm(i,k)[p]);
        if (i > 0) {
          theta_zt[p] += ((i%3)-sp(0.5))/nlev*(nlev-1 - (k*Spack::n+p));
        }
        w_field[p] = interpolate_data(3, d_wind_z_ref, d_w_ref, d_zm(i,k)[p]);
      }

      const Spack pmid_k(d_pmid(i,k));
      const Smask pmid_mask(!isnan(pmid_k) and pmid_k>0.0);
      const Spack exner = PhysicsFun::get_exner(pmid_k,pmid_mask);

      Spack th_atm = theta_zt + (latvap/cpair)*d_shoc_ql(i,k);
      d_t(i,k).set(range_pack1<nlev, th_atm*exner);

      d_omega(i,k).set(range_pack1<nlev, -w_field*(ggr*rrho(i,k)));

      d_s(i,k).set(range_pack1<nlev, cpair*exner*(theta_zt*(1+zvir*qw(i,k))) +
                                      ggr*d_zm(i,k));
    });
    team.team_barrier();

    const auto zt_grid_s = ekat::subview(d_zm, i);
    const auto zi_grid_s = ekat::subview(d_zi, i);
    const auto rrho_s    = ekat::subview(rrho, i);
    const auto rrho_i_s  = ekat::subview(rrhoi, i);
    SHF::linear_interp(team,zt_grid_s,zi_grid_s,rrho_s,rrho_i_s,nlev,nlev+1,0);

    d_shf(i)[0] = 1e-4*(cpair*rrhoi(i,nlev_v)[nlev_p]);
    d_cflx_k0(i)[0] = 1e-6*rrhoi(i,nlev_v)[nlev_p];
    d_wsx(i)[0] = 1e-2*rrhoi(i,nlev_v)[nlev_p];
    d_wsy(i)[0] = 1e-4*rrhoi(i,nlev_v)[nlev_p];

  }); // col loop

  if (m_remapper) {
    m_remapper->registration_ends();

    m_remapper->remap(true);

    // Now we can destroy the remapper
    m_remapper = nullptr;
  }
}

} // namespace scream
