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
  using view_2d  = typename P3F::view_2d<Spack>;

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
  size_t count = 0;
  std::string list_of_fields = "";
  for (const auto& name : fields_to_init)
  {
    list_of_fields += name;
    list_of_fields += ", ";
    count += m_fields.count(name);
  }
  list_of_fields.erase(list_of_fields.size()-2); // Erase last ', '
  std::string inited_fields ="";
  for (const auto& it : m_fields) {
    inited_fields += it.first;
    inited_fields += ", ";
  }
  inited_fields.erase(inited_fields.size()-2); // Erase last ', '
 
  EKAT_REQUIRE_MSG(count!=0,"Error in p3_inputs_initializer: no fields have declared this initializer.  Check p3 interface."); 

  EKAT_REQUIRE_MSG (count==fields_to_init.size() || (count+1)==fields_to_init.size(),
    "Error! P3InputsInitializer is expected to init the following fields:\n"
    "       " + list_of_fields + "\n"
    "       (possibly with the exception of dp). Instead, it is initializing these fields:\n"
    "       " + inited_fields + "\n"
    "       Please, check the atmosphere processes you are using,\n"
    "       and make sure they agree on who's initializing each field.\n");

  // Initialize the fields that we expect.
  // Ask directly for host mirrors.
  auto h_T_atm           = m_fields.at("T_atm").get_reshaped_view<Pack**,Host>();
  auto h_ast             = m_fields.at("ast").get_reshaped_view<Pack**,Host>();
  auto h_ni_activated    = m_fields.at("ni_activated").get_reshaped_view<Pack**,Host>();
  auto h_nc_nuceat_tend  = m_fields.at("nc_nuceat_tend").get_reshaped_view<Pack**,Host>();
  auto h_pmid            = m_fields.at("pmid").get_reshaped_view<Pack**,Host>();
  auto h_zi              = m_fields.at("zi").get_reshaped_view<Pack**,Host>();
  auto h_qv_prev         = m_fields.at("qv_prev").get_reshaped_view<Pack**,Host>();
  auto h_T_prev          = m_fields.at("T_prev").get_reshaped_view<Pack**,Host>();
  auto h_qv              = m_fields.at("qv").get_reshaped_view<Pack**,Host>();
  auto h_qc              = m_fields.at("qc").get_reshaped_view<Pack**,Host>();
  auto h_qr              = m_fields.at("qr").get_reshaped_view<Pack**,Host>();
  auto h_qi              = m_fields.at("qi").get_reshaped_view<Pack**,Host>();
  auto h_qm              = m_fields.at("qm").get_reshaped_view<Pack**,Host>();
  auto h_nc              = m_fields.at("nc").get_reshaped_view<Pack**,Host>();
  auto h_nr              = m_fields.at("nr").get_reshaped_view<Pack**,Host>();
  auto h_ni              = m_fields.at("ni").get_reshaped_view<Pack**,Host>();
  auto h_bm              = m_fields.at("bm").get_reshaped_view<Pack**,Host>();
  auto h_nccn_prescribed = m_fields.at("nccn_prescribed").get_reshaped_view<Pack**,Host>();
  auto h_inv_qc_relvar   = m_fields.at("inv_qc_relvar").get_reshaped_view<Pack**,Host>();

  // dp might be inited by Homme.
  decltype(h_pmid) h_dp;
  if (m_fields.count("dp")==1) {
    h_dp = m_fields.at("dp").get_reshaped_view<Pack**,Host>();
  }

  // Set local views which are used in some of the initialization
  auto mdims = m_fields.at("qc").get_header().get_identifier().get_layout();
  Int ncol = mdims.dim(0); 
  Int nk   = mdims.dim(1);
  const Int nk_pack = ekat::npack<Spack>(nk);
  view_2d::HostMirror h_th_atm("th_atm",ncol,nk_pack);        
  view_2d::HostMirror h_dz("dz",ncol,nk_pack);                
  view_2d::HostMirror h_exner("exner",ncol,nk_pack);          
  
  // Initalize from text file 
  std::ifstream fid("p3_init_vals.txt", std::ifstream::in);
  EKAT_REQUIRE_MSG(!fid.fail(),"Error in p3_inputs_initializer.cpp loading p3_init_vals.txt file");
  std::string tmp_line;
  getline(fid,tmp_line);  // Read header and discard.
  int icol_in_max = 0;
  while(getline(fid,tmp_line))
  {
    std::stringstream s(tmp_line);
    std::string field;
    std::vector<Real> field_vals;
    std::string skip;
    int icol, k, ipack, ivec;
    s >> icol >> k;
    ipack = k / Spack::n;
    ivec  = k % Spack::n;
    icol_in_max = std::max(icol,icol_in_max);
    s >> h_qv(icol,ipack)[ivec]             ;
    s >> h_th_atm(icol,ipack)[ivec]         ;
    s >> h_pmid(icol,ipack)[ivec]           ;
    s >> h_dz(icol,ipack)[ivec]             ;
    s >> h_nc_nuceat_tend(icol,ipack)[ivec] ;
    s >> h_nccn_prescribed(icol,ipack)[ivec];
    s >> h_ni_activated(icol,ipack)[ivec]   ;
    s >> h_inv_qc_relvar(icol,ipack)[ivec]  ;
    s >> h_qc(icol,ipack)[ivec]             ;
    s >> h_nc(icol,ipack)[ivec]             ;
    s >> h_qr(icol,ipack)[ivec]             ;
    s >> h_nr(icol,ipack)[ivec]             ;
    s >> h_qi(icol,ipack)[ivec]             ;
    s >> h_ni(icol,ipack)[ivec]             ;
    s >> h_qm(icol,ipack)[ivec]             ;
    s >> h_bm(icol,ipack)[ivec]             ;
    for (int skp=0;skp<5;skp++) { s >> skip; }
    if (h_dp.data()!=nullptr) {
      s >> h_dp(icol,ipack)[ivec]   ;
    } else {
      s >> skip;
    }
    s >> h_exner(icol,ipack)[ivec];
    for (int skp=0;skp<7;skp++) { s >> skip; }
    s >> h_ast(icol,ipack)[ivec];
    for (int skp=0;skp<6;skp++) { s >> skip; }
    s >> h_qv_prev(icol,ipack)[ivec];
    s >> h_T_prev(icol,ipack)[ivec] ;

  } // while getline(fid,tmp_line)
  // For now use dummy values copied from `p3_ic_cases.cpp`, which is loaded from a file.
  // That file only has data for 3 columns, need to expand to >3 columns.
  for (int icol_i = 0;icol_i<ncol;icol_i++)
  {
    for (int k = nk-1;k>=0;k--)
    {
      int icol  = icol_i % (icol_in_max+1);
      int ipack = k / Spack::n;
      int ivec  = k % Spack::n;
      int ipack_p1 = (k+1) / Spack::n;
      int ivec_p1  = (k+1) % Spack::n;
      h_qv(icol_i,ipack)[ivec]              = h_qv(icol,ipack)[ivec]                              ;
      h_th_atm(icol_i,ipack)[ivec]          = h_th_atm(icol,ipack)[ivec]                          ;
      h_T_atm(icol_i,ipack)[ivec]           = h_th_atm(icol,ipack)[ivec]*h_exner(icol,ipack)[ivec];
      h_pmid(icol_i,ipack)[ivec]            = h_pmid(icol,ipack)[ivec]                            ;
      h_dz(icol_i,ipack)[ivec]              = h_dz(icol,ipack)[ivec]                              ;
      h_nc_nuceat_tend(icol_i,ipack)[ivec]  = h_nc_nuceat_tend(icol,ipack)[ivec]                  ;
      h_nccn_prescribed(icol_i,ipack)[ivec] = h_nccn_prescribed(icol,ipack)[ivec]                 ;
      h_ni_activated(icol_i,ipack)[ivec]    = h_ni_activated(icol,ipack)[ivec]                    ;
      h_inv_qc_relvar(icol_i,ipack)[ivec]   = h_inv_qc_relvar(icol,ipack)[ivec]                   ;
      h_qc(icol_i,ipack)[ivec]              = h_qc(icol,ipack)[ivec]                              ;
      h_nc(icol_i,ipack)[ivec]              = h_nc(icol,ipack)[ivec]                              ;
      h_qr(icol_i,ipack)[ivec]              = h_qr(icol,ipack)[ivec]                              ;
      h_nr(icol_i,ipack)[ivec]              = h_nr(icol,ipack)[ivec]                              ;
      h_qi(icol_i,ipack)[ivec]              = h_qi(icol,ipack)[ivec]                              ;
      h_ni(icol_i,ipack)[ivec]              = h_ni(icol,ipack)[ivec]                              ;
      h_qm(icol_i,ipack)[ivec]              = h_qm(icol,ipack)[ivec]                              ;
      h_bm(icol_i,ipack)[ivec]              = h_bm(icol,ipack)[ivec]                              ;
      if (h_dp.data()!=nullptr) {
        h_dp(icol_i,ipack)[ivec]              = h_dp(icol,ipack)[ivec]                              ;
      }
      h_exner(icol_i,ipack)[ivec]           = h_exner(icol,ipack)[ivec]                           ;
      h_ast(icol_i,ipack)[ivec]             = h_ast(icol,ipack)[ivec]                             ;
      h_qv_prev(icol_i,ipack)[ivec]         = h_qv_prev(icol,ipack)[ivec]                         ;
      h_T_prev(icol_i,ipack)[ivec]          = h_T_prev(icol,ipack)[ivec]                          ;
      // Initialize arrays not provided in input file
      if (k==nk-1) { h_zi(icol_i,ipack_p1)[ivec_p1] = 0.0; } // If this is the bottom layer then initialize zi[nk]
      h_zi(icol_i,ipack)[ivec] = h_zi(icol_i,ipack_p1)[ivec_p1] + h_dz(icol_i,ipack)[ivec];
    } // for k
  } // for icol_i

  // Use Field interface to deep copy back to device
  m_fields.at("T_atm").sync_to_dev();
  m_fields.at("ast").sync_to_dev();
  m_fields.at("ni_activated").sync_to_dev();
  m_fields.at("nc_nuceat_tend").sync_to_dev();
  m_fields.at("pmid").sync_to_dev();
  m_fields.at("zi").sync_to_dev();
  m_fields.at("qv_prev").sync_to_dev();
  m_fields.at("T_prev").sync_to_dev();
  m_fields.at("qv").sync_to_dev();
  m_fields.at("qc").sync_to_dev();
  m_fields.at("qr").sync_to_dev();
  m_fields.at("qi").sync_to_dev();
  m_fields.at("qm").sync_to_dev();
  m_fields.at("nc").sync_to_dev();
  m_fields.at("nr").sync_to_dev();
  m_fields.at("ni").sync_to_dev();
  m_fields.at("bm").sync_to_dev();
  m_fields.at("nccn_prescribed").sync_to_dev();
  m_fields.at("inv_qc_relvar").sync_to_dev();
  if (h_dp.data()!=nullptr) {
    m_fields.at("dp").sync_to_dev();
  }

  if (m_remapper) {
    m_remapper->registration_ends();

    m_remapper->remap(true);

    // Now we can destroy the remapper
    m_remapper = nullptr;
  }
}

} // namespace scream
