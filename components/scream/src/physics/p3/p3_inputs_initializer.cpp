#include "physics/p3/p3_inputs_initializer.hpp"
#include "physics/p3/p3_main_impl.hpp"
#include "ekat/util/ekat_file_utils.hpp"

#include <array>
#include <fstream>

namespace scream
{

void P3InputsInitializer::add_field (const field_type &f)
{
  const auto& id = f.get_header().get_identifier();
  
  m_fields.emplace(id.name(),f);
  m_fields_id.insert(id);
}

void P3InputsInitializer::
add_field (const field_type &f, const field_type& f_ref,
           const remapper_ptr_type& remapper)
{
  if (m_remapper) {
    // Sanity check
    EKAT_REQUIRE_MSG (m_remapper->get_src_grid()->name()==remapper->get_src_grid()->name(),
      "Error! A remapper was already set in P3InputsInitializer, but its src grid differs from"
      "       the grid of the input remapper of this call.\n");
  } else {
    m_remapper = remapper;
    m_remapper->registration_begins();
  }

  const auto& id = f.get_header().get_identifier();
  const auto& id_ref = f_ref.get_header().get_identifier();

  // To the AD, we only expose the fact that we init f_ref...
  m_fields_id.insert(id_ref);

  // ...but P3 only knows how to init f...
  m_fields.emplace(id.name(),f);

  // ...hence, we remap to f_ref.
  m_remapper->register_field(f, f_ref);
}


// =========================================================================================
void P3InputsInitializer::initialize_fields ()
{
  using namespace p3;
  using P3F             = Functions<Real, DefaultDevice>;
  using Spack           = typename P3F::Spack;
  using Pack            = ekat::Pack<Real,Spack::n>;

  // Safety check: if we're asked to init anything at all,
  // To simplify the initializer we first define all the fields we expect to have to initialize.
  std::vector<std::string> fields_to_init;
  fields_to_init.push_back("T_atm");
  fields_to_init.push_back("ast");
  fields_to_init.push_back("ni_activated");
  fields_to_init.push_back("nc_nuceat_tend");
  fields_to_init.push_back("pmid");
  fields_to_init.push_back("dp");
  fields_to_init.push_back("zi");
  fields_to_init.push_back("qv_prev");
  fields_to_init.push_back("T_prev");
  fields_to_init.push_back("qv");
  fields_to_init.push_back("qc");
  fields_to_init.push_back("qr");
  fields_to_init.push_back("qi");
  fields_to_init.push_back("qm");
  fields_to_init.push_back("nc");
  fields_to_init.push_back("nr");
  fields_to_init.push_back("ni");
  fields_to_init.push_back("bm");
  fields_to_init.push_back("nccn_prescribed");
  fields_to_init.push_back("inv_qc_relvar");
  // TODO: Delete eventually, should instead be set in run interface: 
  fields_to_init.push_back("th_atm");  
  fields_to_init.push_back("dz");  
  fields_to_init.push_back("exner");  
  fields_to_init.push_back("cld_frac_l");  
  fields_to_init.push_back("cld_frac_i");  
  fields_to_init.push_back("cld_frac_r");  
  // then we should have been asked to init 20 + 6 (to be deleted) fields.
  int count = 0;
  std::string list_of_fields = "";
  for (auto name : fields_to_init)
  {
    list_of_fields += name;
    list_of_fields += ", ";
    count += m_fields.count(name);
  }
  
  if (count==0) {
    return;
  }

  EKAT_REQUIRE_MSG (count==fields_to_init.size(),
    "Error! P3InputsInitializer is expected to init " + std::to_string(fields_to_init.size()) + " fields:\n"
    "       " + list_of_fields + "\n"
    "       Instead found " + std::to_string(count) + " fields.\n"
    "       Please, check the atmosphere processes you are using,\n"
    "       and make sure they agree on who's initializing each field.\n");

  // Initialize the fields that we expect.
  // Get device views
  auto d_T_atm           = m_fields.at("T_atm").get_reshaped_view<Pack**>();
  auto d_ast             = m_fields.at("ast").get_reshaped_view<Pack**>();
  auto d_ni_activated    = m_fields.at("ni_activated").get_reshaped_view<Pack**>();
  auto d_nc_nuceat_tend  = m_fields.at("nc_nuceat_tend").get_reshaped_view<Pack**>();
  auto d_pmid            = m_fields.at("pmid").get_reshaped_view<Pack**>();
  auto d_dp              = m_fields.at("dp").get_reshaped_view<Pack**>();
  auto d_zi              = m_fields.at("zi").get_reshaped_view<Pack**>();
  auto d_qv_prev         = m_fields.at("qv_prev").get_reshaped_view<Pack**>();
  auto d_T_prev          = m_fields.at("T_prev").get_reshaped_view<Pack**>();
  auto d_qv              = m_fields.at("qv").get_reshaped_view<Pack**>();
  auto d_qc              = m_fields.at("qc").get_reshaped_view<Pack**>();
  auto d_qr              = m_fields.at("qr").get_reshaped_view<Pack**>();
  auto d_qi              = m_fields.at("qi").get_reshaped_view<Pack**>();
  auto d_qm              = m_fields.at("qm").get_reshaped_view<Pack**>();
  auto d_nc              = m_fields.at("nc").get_reshaped_view<Pack**>();
  auto d_nr              = m_fields.at("nr").get_reshaped_view<Pack**>();
  auto d_ni              = m_fields.at("ni").get_reshaped_view<Pack**>();
  auto d_bm              = m_fields.at("bm").get_reshaped_view<Pack**>();
  auto d_nccn_prescribed = m_fields.at("nccn_prescribed").get_reshaped_view<Pack**>();
  auto d_inv_qc_relvar   = m_fields.at("inv_qc_relvar").get_reshaped_view<Pack**>();
  // TODO: Delete eventually, should instead be set in run interface: 
  auto d_th_atm          = m_fields.at("th_atm").get_reshaped_view<Pack**>();  
  auto d_dz              = m_fields.at("dz").get_reshaped_view<Pack**>();  
  auto d_exner           = m_fields.at("exner").get_reshaped_view<Pack**>();  
  auto d_cld_frac_l      = m_fields.at("cld_frac_l").get_reshaped_view<Pack**>();  
  auto d_cld_frac_i      = m_fields.at("cld_frac_i").get_reshaped_view<Pack**>();  
  auto d_cld_frac_r      = m_fields.at("cld_frac_r").get_reshaped_view<Pack**>();

  
  // Create host mirrors
  auto h_T_atm            = Kokkos::create_mirror_view(d_T_atm          ); 
  auto h_ast              = Kokkos::create_mirror_view(d_ast            ); 
  auto h_ni_activated     = Kokkos::create_mirror_view(d_ni_activated   ); 
  auto h_nc_nuceat_tend   = Kokkos::create_mirror_view(d_nc_nuceat_tend ); 
  auto h_pmid             = Kokkos::create_mirror_view(d_pmid           ); 
  auto h_dp               = Kokkos::create_mirror_view(d_dp             ); 
  auto h_zi               = Kokkos::create_mirror_view(d_zi             ); 
  auto h_qv_prev          = Kokkos::create_mirror_view(d_qv_prev        ); 
  auto h_T_prev           = Kokkos::create_mirror_view(d_T_prev         ); 
  auto h_qv               = Kokkos::create_mirror_view(d_qv             ); 
  auto h_qc               = Kokkos::create_mirror_view(d_qc             ); 
  auto h_qr               = Kokkos::create_mirror_view(d_qr             ); 
  auto h_qi               = Kokkos::create_mirror_view(d_qi             ); 
  auto h_qm               = Kokkos::create_mirror_view(d_qm             ); 
  auto h_nc               = Kokkos::create_mirror_view(d_nc             ); 
  auto h_nr               = Kokkos::create_mirror_view(d_nr             ); 
  auto h_ni               = Kokkos::create_mirror_view(d_ni             ); 
  auto h_bm               = Kokkos::create_mirror_view(d_bm             ); 
  auto h_nccn_prescribed  = Kokkos::create_mirror_view(d_nccn_prescribed); 
  auto h_inv_qc_relvar    = Kokkos::create_mirror_view(d_inv_qc_relvar  ); 
  // TODO: Delete eventually, should instead be set in run interface: 
  auto h_th_atm           = Kokkos::create_mirror_view(d_th_atm         ); 
  auto h_dz               = Kokkos::create_mirror_view(d_dz             ); 
  auto h_exner            = Kokkos::create_mirror_view(d_exner          ); 
  auto h_cld_frac_l       = Kokkos::create_mirror_view(d_cld_frac_l     ); 
  auto h_cld_frac_i       = Kokkos::create_mirror_view(d_cld_frac_i     ); 
  auto h_cld_frac_r       = Kokkos::create_mirror_view(d_cld_frac_r     ); 

 // Initalize from text file 
  std::ifstream fid("p3_init_vals.txt", std::ifstream::in);
  std::string tmp_line;
  int icol_in_max = 0;
  while(getline(fid,tmp_line))
  {
    std::stringstream s(tmp_line);
    std::string field;
    std::vector<Real> field_vals;
    while (getline(s,field,' '))
    {
      field_vals.push_back(std::stod(field));
    }
    int icol  = (int)field_vals[0];
    int ipack = (int)field_vals[1] / Spack::n;
    int ivec  = (int)field_vals[1] % Spack::n;
    icol_in_max = std::max(icol,icol_in_max);
    int cnt = 2;
    h_qv(icol,ipack)[ivec]=field_vals[cnt++];
    h_th_atm(icol,ipack)[ivec]=field_vals[cnt++];
    h_pmid(icol,ipack)[ivec]=field_vals[cnt++];
    h_dz(icol,ipack)[ivec]=field_vals[cnt++];
    h_nc_nuceat_tend(icol,ipack)[ivec]=field_vals[cnt++];
    h_nccn_prescribed(icol,ipack)[ivec]=field_vals[cnt++];
    h_ni_activated(icol,ipack)[ivec]=field_vals[cnt++];
    h_inv_qc_relvar(icol,ipack)[ivec]=field_vals[cnt++];
    h_qc(icol,ipack)[ivec]=field_vals[cnt++];
    h_nc(icol,ipack)[ivec]=field_vals[cnt++];
    h_qr(icol,ipack)[ivec]=field_vals[cnt++];
    h_nr(icol,ipack)[ivec]=field_vals[cnt++];
    h_qi(icol,ipack)[ivec]=field_vals[cnt++];
    h_ni(icol,ipack)[ivec]=field_vals[cnt++];
    h_qm(icol,ipack)[ivec]=field_vals[cnt++];
    h_bm(icol,ipack)[ivec]=field_vals[cnt++];
    cnt+=5;
    //h_precip_liq_surf(icol,ipack)[ivec]=field_vals[cnt];    //
    //h_precip_ice_surf(icol,ipack)[ivec]=field_vals[cnt];    //
    //h_diag_eff_radius_qc(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_diag_eff_radius_qi(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_rho_qi(icol,ipack)[ivec]=field_vals[cnt]; //
    h_dp(icol,ipack)[ivec]=field_vals[cnt++];
    h_exner(icol,ipack)[ivec]=field_vals[cnt++];
    cnt+=6;
    //h_qv2qi_depos_tend(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_precip_total_tend(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_nevapr(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_qr_evap_tend(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_precip_liq_flux(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_precip_ice_flux(icol,ipack)[ivec]=field_vals[cnt]; //
    h_cld_frac_r(icol,ipack)[ivec]=field_vals[cnt++];
    h_cld_frac_l(icol,ipack)[ivec]=field_vals[cnt++];
    h_cld_frac_i(icol,ipack)[ivec]=field_vals[cnt++];
    cnt+=5;
    //h_mu_c(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_lamc(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_liq_ice_exchange(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_vap_liq_exchange(icol,ipack)[ivec]=field_vals[cnt]; //
    //h_vap_ice_exchange(icol,ipack)[ivec]=field_vals[cnt]; //
    h_qv_prev(icol,ipack)[ivec]=field_vals[cnt++]; 
    h_T_prev(icol,ipack)[ivec]=field_vals[cnt];

  } // while getline(fid,tmp_line)
  // For now use dummy values copied from `p3_ic_cases.cpp`, which is loaded from a file.
  // That file only has data for 3 columns, need to expand to >3 columns.
  auto mdims = m_fields.at("qc").get_header().get_identifier().get_layout();
  Int ncol = mdims.dim(0); 
  Int nk   = mdims.dim(1);
  for (int icol_i = icol_in_max;icol_i<ncol;icol_i++)
  {
    for (int k = 0;k<nk;k++)
    {
    int icol  = icol_i % icol_in_max;
    int ipack = k / Spack::n;
    int ivec  = k % Spack::n;
    h_qv(icol_i,ipack)[ivec]              = h_qv(icol,ipack)[ivec]             ;
    h_th_atm(icol_i,ipack)[ivec]          = h_th_atm(icol,ipack)[ivec]         ;
    h_pmid(icol_i,ipack)[ivec]            = h_pmid(icol,ipack)[ivec]           ;
    h_dz(icol_i,ipack)[ivec]              = h_dz(icol,ipack)[ivec]             ;
    h_nc_nuceat_tend(icol_i,ipack)[ivec]  = h_nc_nuceat_tend(icol,ipack)[ivec] ;
    h_nccn_prescribed(icol_i,ipack)[ivec] = h_nccn_prescribed(icol,ipack)[ivec];
    h_ni_activated(icol_i,ipack)[ivec]    = h_ni_activated(icol,ipack)[ivec]   ;
    h_inv_qc_relvar(icol_i,ipack)[ivec]   = h_inv_qc_relvar(icol,ipack)[ivec]  ;
    h_qc(icol_i,ipack)[ivec]              = h_qc(icol,ipack)[ivec]             ;
    h_nc(icol_i,ipack)[ivec]              = h_nc(icol,ipack)[ivec]             ;
    h_qr(icol_i,ipack)[ivec]              = h_qr(icol,ipack)[ivec]             ;
    h_nr(icol_i,ipack)[ivec]              = h_nr(icol,ipack)[ivec]             ;
    h_qi(icol_i,ipack)[ivec]              = h_qi(icol,ipack)[ivec]             ;
    h_ni(icol_i,ipack)[ivec]              = h_ni(icol,ipack)[ivec]             ;
    h_qm(icol_i,ipack)[ivec]              = h_qm(icol,ipack)[ivec]             ;
    h_bm(icol_i,ipack)[ivec]              = h_bm(icol,ipack)[ivec]             ;
    h_dp(icol_i,ipack)[ivec]              = h_dp(icol,ipack)[ivec]             ;
    h_exner(icol_i,ipack)[ivec]           = h_exner(icol,ipack)[ivec]          ;
    h_cld_frac_r(icol_i,ipack)[ivec]      = h_cld_frac_r(icol,ipack)[ivec]     ;
    h_cld_frac_l(icol_i,ipack)[ivec]      = h_cld_frac_l(icol,ipack)[ivec]     ;
    h_cld_frac_i(icol_i,ipack)[ivec]      = h_cld_frac_i(icol,ipack)[ivec]     ;
    h_qv_prev(icol_i,ipack)[ivec]         = h_qv_prev(icol,ipack)[ivec]        ;
    h_T_prev(icol_i,ipack)[ivec]          = h_T_prev(icol,ipack)[ivec]         ;
    }
  }
  printf("ASD - new version p3, %d %d, %d\n",ncol,nk,icol_in_max);

  // Deep copy back to device
  Kokkos::deep_copy(d_T_atm          , h_T_atm          ); 
  Kokkos::deep_copy(d_ast            , h_ast            ); 
  Kokkos::deep_copy(d_ni_activated   , h_ni_activated   ); 
  Kokkos::deep_copy(d_nc_nuceat_tend , h_nc_nuceat_tend ); 
  Kokkos::deep_copy(d_pmid           , h_pmid           ); 
  Kokkos::deep_copy(d_dp             , h_dp             ); 
  Kokkos::deep_copy(d_zi             , h_zi             ); 
  Kokkos::deep_copy(d_qv_prev        , h_qv_prev        ); 
  Kokkos::deep_copy(d_T_prev         , h_T_prev         ); 
  Kokkos::deep_copy(d_qv             , h_qv             ); 
  Kokkos::deep_copy(d_qc             , h_qc             ); 
  Kokkos::deep_copy(d_qr             , h_qr             ); 
  Kokkos::deep_copy(d_qi             , h_qi             ); 
  Kokkos::deep_copy(d_qm             , h_qm             ); 
  Kokkos::deep_copy(d_nc             , h_nc             ); 
  Kokkos::deep_copy(d_nr             , h_nr             ); 
  Kokkos::deep_copy(d_ni             , h_ni             ); 
  Kokkos::deep_copy(d_bm             , h_bm             ); 
  Kokkos::deep_copy(d_nccn_prescribed, h_nccn_prescribed); 
  Kokkos::deep_copy(d_inv_qc_relvar  , h_inv_qc_relvar  ); 
  // TODO: Delete eventually, should instead be set in run interface: 
  Kokkos::deep_copy(d_th_atm         , h_th_atm    ); 
  Kokkos::deep_copy(d_dz             , h_dz        ); 
  Kokkos::deep_copy(d_exner          , h_exner     ); 
  Kokkos::deep_copy(d_cld_frac_l     , h_cld_frac_l); 
  Kokkos::deep_copy(d_cld_frac_i     , h_cld_frac_i); 
  Kokkos::deep_copy(d_cld_frac_r     , h_cld_frac_r); 

  if (m_remapper) {
    m_remapper->registration_ends();

    m_remapper->remap(true);

    // Now we can destroy the remapper
    m_remapper = nullptr;
  }
}

} // namespace scream
