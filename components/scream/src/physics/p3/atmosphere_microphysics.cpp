#include "share/scream_assert.hpp"
#include "physics/p3/scream_p3_interface.hpp"
#include "physics/p3/atmosphere_microphysics.hpp"

namespace scream
{

// =========================================================================================
P3Microphysics::P3Microphysics (const Comm& comm,const ParameterList& params)
 : m_p3_comm (comm)
 , m_p3_params (params)
{
/* Anything that can be initialized without grid information can be initialized here.
 * Like universal constants, table lookups, p3 options.
*/
  if (!m_p3_params.isParameter("Grid")) {
    m_p3_params.set("Grid",std::string("SE Physics"));
  }
}

// =========================================================================================
void P3Microphysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{

  using namespace units;
  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  auto Qdp = Q * Pa;
  Q.set_string("kg/kg");
  Qdp.set_string("kg/kg Pa");

  constexpr int NVL = 72;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */
  constexpr int QSZ =  35;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */

  const auto& grid_name = m_p3_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);
  const int num_dofs = grid->get_num_local_dofs();
  const int nc = num_dofs;

  using namespace ShortFieldTagsNames;

  FieldLayout scalar3d_layout_mid { {COL,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout scalar3d_layout_int { {COL,VL}, {nc,NVL+1} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout_mid{ {COL,CMP,VL}, {nc,QSZ,NVL} };
  FieldLayout tracers_layout { {COL,VAR,VL}, {nc,QSZ,NVL} };

  // set requirements
//  m_required_fields.emplace("P3_req_test",  scalar3d_layout,   units::one, grid_name);
//  m_required_fields.emplace("ASD_req_test",  scalar3d_layout,   units::one, grid_name);
  // Inputs from Dynamics
  m_required_fields.emplace("dp", scalar3d_layout_mid,    Pa,             grid_name);
  m_required_fields.emplace("qdp", tracers_layout,        Qdp,          grid_name);  
  // Inputs from Physics  NOTE: Making them "computed" fields now because we will initialize them in P3 for testing purposes.
  m_computed_fields.emplace("ast",    scalar3d_layout_mid,            A,             grid_name);
  m_computed_fields.emplace("naai",   scalar3d_layout_mid,   pow(kg,-1),           grid_name);
  m_computed_fields.emplace("npccn",  scalar3d_layout_mid,   pow(kg,-1)*pow(s,-1), grid_name);
  // Inputs that are STATE variables
  m_computed_fields.emplace("pmid",  scalar3d_layout_mid,    Pa,                   grid_name);
  m_computed_fields.emplace("pdel",  scalar3d_layout_mid,    Pa,                   grid_name);
  m_computed_fields.emplace("zi",    scalar3d_layout_int, m,                   grid_name);

  // set computed
//  m_computed_fields.emplace("P3_comq_test", scalar3d_layout, units::one, grid_name);
  m_computed_fields.emplace("FQ", tracers_layout,      Q, grid_name);
  // INPUT/OUTPUT from STATE
  m_computed_fields.emplace("T",  scalar3d_layout_mid, K, grid_name);
  m_computed_fields.emplace("q",  vector3d_layout_mid, Q, grid_name);

}
// =========================================================================================
void P3Microphysics::initialize (const util::TimeStamp& t0)
{
  m_current_ts = t0;
  auto q_ptr = m_p3_fields_out.at("q").get_view().data();
  auto ast_ptr = m_p3_fields_out.at("ast").get_view().data();
  auto naai_ptr = m_p3_fields_out.at("naai").get_view().data();
  auto npccn_ptr = m_p3_fields_out.at("npccn").get_view().data();
  auto pmid_ptr = m_p3_fields_out.at("pmid").get_view().data();
  auto pdel_ptr = m_p3_fields_out.at("pdel").get_view().data();
  auto zi_ptr = m_p3_fields_out.at("zi").get_view().data();
  auto T_ptr = m_p3_fields_out.at("T").get_view().data();
  const std::string dimnames[2] = {"column", "level"};
  const std::string filename = "p3_output";
  // const int  dimrng[2]   = {218,72}; // TODO make this based on actual col,levels data
  // const int  dimnum = 2;

  p3_init_f90 (q_ptr,T_ptr,zi_ptr,pmid_ptr,pdel_ptr,
     ast_ptr, naai_ptr, npccn_ptr);
}

// =========================================================================================
void P3Microphysics::run (const Real dt)
{
  auto q_ptr = m_p3_fields_out.at("q").get_view().data();
  auto FQ_ptr = m_p3_fields_out.at("FQ").get_view().data();
  // auto dp_ptr = m_p3_fields_in.at("dp").get_view().data();
  auto qdp_ptr = m_p3_fields_in.at("qdp").get_view().data();
  auto ast_ptr = m_p3_fields_out.at("ast").get_view().data();
  auto naai_ptr = m_p3_fields_out.at("naai").get_view().data();
  auto npccn_ptr = m_p3_fields_out.at("npccn").get_view().data();
  auto pmid_ptr = m_p3_fields_out.at("pmid").get_view().data();
  auto pdel_ptr = m_p3_fields_out.at("pdel").get_view().data();
  auto zi_ptr = m_p3_fields_out.at("zi").get_view().data();
  auto T_ptr = m_p3_fields_out.at("T").get_view().data();

  p3_main_f90 (dt,q_ptr,FQ_ptr,qdp_ptr,T_ptr,zi_ptr,pmid_ptr,pdel_ptr,
     ast_ptr, naai_ptr, npccn_ptr);

  m_current_ts += dt;
  m_p3_fields_out.at("q").get_header().get_tracking().update_time_stamp(m_current_ts);
  m_p3_fields_out.at("FQ").get_header().get_tracking().update_time_stamp(m_current_ts);
  m_p3_fields_out.at("T").get_header().get_tracking().update_time_stamp(m_current_ts);
}
// =========================================================================================
void P3Microphysics::finalize()
{
  p3_finalize_f90 ();
}
// =========================================================================================

void P3Microphysics::register_fields (FieldRepository<Real, device_type>& field_repo) const {
  for (auto& fid : m_required_fields) {
    field_repo.register_field(fid);
  }
  for (auto& fid : m_computed_fields) {
    field_repo.register_field(fid);
  }
}

void P3Microphysics::set_required_field_impl (const Field<const Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_p3_fields_in.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as customer to the field
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void P3Microphysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_p3_fields_out.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as provider for the field
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

} // namespace scream
