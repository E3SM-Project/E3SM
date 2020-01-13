#include "share/scream_assert.hpp"
#include "physics/p3/scream_p3_interface.hpp"
#include "physics/p3/atmosphere_microphysics.hpp"

#include <array>

namespace scream
{
/* 
 * P3 Microphysics routines
*/

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

  // Inputs
  auto nondim = m/m;
  m_required_fields.emplace("qdp",    tracers_layout,           Qdp, grid_name);  
  m_required_fields.emplace("ast",    scalar3d_layout_mid,   nondim, grid_name);
  m_required_fields.emplace("naai",   scalar3d_layout_mid,     1/kg, grid_name);
  m_required_fields.emplace("npccn",  scalar3d_layout_mid, 1/(kg*s), grid_name);
  m_required_fields.emplace("pmid",   scalar3d_layout_mid,       Pa, grid_name);
  m_required_fields.emplace("pdel",   scalar3d_layout_mid,       Pa, grid_name);
  m_required_fields.emplace("zi",     scalar3d_layout_int,        m, grid_name);

  // Outputs
  m_computed_fields.emplace("FQ", tracers_layout,      Q, grid_name);
  m_computed_fields.emplace("T",  scalar3d_layout_mid, K, grid_name);
  m_computed_fields.emplace("q",  vector3d_layout_mid, Q, grid_name);
}
// =========================================================================================
void P3Microphysics::initialize (const util::TimeStamp& t0)
{
  m_current_ts = t0;

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_p3_fields_in) {
    Kokkos::deep_copy(m_p3_host_views_in.at(it.first),it.second.get_view());
  }
  for (auto& it : m_p3_fields_out) {
    Kokkos::deep_copy(m_p3_host_views_out.at(it.first),it.second.get_view());
  }

  // Note: in order to init some fields for tests, we need to write them,
  //       even if they are marked as inputs.
  std::map<std::string,Real*> ptrs_in;
  for (auto it : m_raw_ptrs_in) {
    ptrs_in[it.first] = const_cast<Real*>(it.second);
  }

  // Call f90 routine
  p3_init_f90 (m_raw_ptrs_out["q"], m_raw_ptrs_out["T"], ptrs_in["zi"], ptrs_in["pmid"], ptrs_in["pdel"], ptrs_in["ast"], ptrs_in["naai"], ptrs_in["npccn"]);

  // Copy outputs back to device
  for (auto& it : m_p3_fields_out) {
    Kokkos::deep_copy(it.second.get_view(),m_p3_host_views_out.at(it.first));
  }
}

// =========================================================================================
void P3Microphysics::run (const Real dt)
{
  // std::array<const char*, num_views> view_names = {"q", "FQ", "qdp", "T", "zi", "pmid", "pdel", "ast", "naai", "npccn"};

  std::vector<const Real*> in;
  std::vector<Real*> out;

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_p3_fields_in) {
    Kokkos::deep_copy(m_p3_host_views_in.at(it.first),it.second.get_view());
  }
  for (auto& it : m_p3_fields_out) {
    Kokkos::deep_copy(m_p3_host_views_out.at(it.first),it.second.get_view());
  }

  // Call f90 routine
  p3_main_f90 (dt, m_raw_ptrs_in["qdp"], m_raw_ptrs_in["zi"], m_raw_ptrs_in["pmid"], m_raw_ptrs_in["pdel"], m_raw_ptrs_in["ast"], m_raw_ptrs_in["naai"], m_raw_ptrs_in["npccn"], m_raw_ptrs_out["q"], m_raw_ptrs_out["FQ"], m_raw_ptrs_out["T"]);

  // Copy outputs back to device
  for (auto& it : m_p3_fields_out) {
    Kokkos::deep_copy(it.second.get_view(),m_p3_host_views_out.at(it.first));
  }
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
  const auto& name = f.get_header().get_identifier().name();
  m_p3_fields_in.emplace(name,f);
  m_p3_host_views_in[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_in[name] = m_p3_host_views_in[name].data();

  // Add myself as customer to the field
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void P3Microphysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_p3_fields_out.emplace(name,f);
  m_p3_host_views_out[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_out[name] = m_p3_host_views_out[name].data();

  // Add myself as provider for the field
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

/* 
 * Stand-alone Microphysics stub routines
*/

// =========================================================================================
SAMicrophysics::SAMicrophysics (const Comm& comm,const ParameterList& params)
 : m_sa_comm (comm)
 , m_sa_params (params)
{
/* Anything that can be initialized without grid information can be initialized here.
 * Like universal constants, table lookups, sa options.
*/
  if (!m_sa_params.isParameter("Grid")) {
    m_sa_params.set("Grid",std::string("SE Physics"));
  }
}

// =========================================================================================
void SAMicrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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

  const auto& grid_name = m_sa_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);
  const int num_dofs = grid->get_num_local_dofs();
  const int nc = num_dofs;

  using namespace ShortFieldTagsNames;

  FieldLayout scalar3d_layout_mid { {COL,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout scalar3d_layout_int { {COL,VL}, {nc,NVL+1} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout_mid{ {COL,CMP,VL}, {nc,QSZ,NVL} };
  FieldLayout tracers_layout { {COL,VAR,VL}, {nc,QSZ,NVL} };

  // Since this is a stand-alone stub for initialization, all fields are computed fields
  auto nondim = m/m;
  m_computed_fields.emplace("qdp",    tracers_layout,           Qdp, grid_name);  
  m_computed_fields.emplace("ast",    scalar3d_layout_mid,   nondim, grid_name);
  m_computed_fields.emplace("naai",   scalar3d_layout_mid,     1/kg, grid_name);
  m_computed_fields.emplace("npccn",  scalar3d_layout_mid, 1/(kg*s), grid_name);
  m_computed_fields.emplace("pmid",   scalar3d_layout_mid,       Pa, grid_name);
  m_computed_fields.emplace("pdel",   scalar3d_layout_mid,       Pa, grid_name);
  m_computed_fields.emplace("zi",     scalar3d_layout_int,        m, grid_name);

  // Outputs
  m_computed_fields.emplace("FQ", tracers_layout,      Q, grid_name);
  m_computed_fields.emplace("T",  scalar3d_layout_mid, K, grid_name);
  m_computed_fields.emplace("q",  vector3d_layout_mid, Q, grid_name);
}
// =========================================================================================
void SAMicrophysics::initialize (const util::TimeStamp& t0)
{
  m_current_ts = t0;

  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_sa_fields_in) {
    Kokkos::deep_copy(m_sa_host_views_in.at(it.first),it.second.get_view());
  }
  for (auto& it : m_sa_fields_out) {
    Kokkos::deep_copy(m_sa_host_views_out.at(it.first),it.second.get_view());
  }

  // Call f90 routine
  sa_init_f90 (m_raw_ptrs_out["q"], m_raw_ptrs_out["T"], m_raw_ptrs_out["zi"], m_raw_ptrs_out["pmid"], m_raw_ptrs_out["pdel"], m_raw_ptrs_out["ast"], m_raw_ptrs_out["naai"], m_raw_ptrs_out["npccn"]);

  // Copy outputs back to device
  for (auto& it : m_sa_fields_out) {
    Kokkos::deep_copy(it.second.get_view(),m_sa_host_views_out.at(it.first));
  }
}

// =========================================================================================
void SAMicrophysics::run (const Real dt)
{
  // Since it's just for initialization ,do nothing here 
}
// =========================================================================================
void SAMicrophysics::finalize()
{
  // Since it's just for initialization ,do nothing here 
}
// =========================================================================================

void SAMicrophysics::register_fields (FieldRepository<Real, device_type>& field_repo) const {
  for (auto& fid : m_required_fields) {
    field_repo.register_field(fid);
  }
  for (auto& fid : m_computed_fields) {
    field_repo.register_field(fid);
  }
}

void SAMicrophysics::set_required_field_impl (const Field<const Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_sa_fields_in.emplace(name,f);
  m_sa_host_views_in[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_in[name] = m_sa_host_views_in[name].data();

  // Add myself as customer to the field
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void SAMicrophysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_sa_fields_out.emplace(name,f);
  m_sa_host_views_out[name] = Kokkos::create_mirror_view(f.get_view());
  m_raw_ptrs_out[name] = m_sa_host_views_out[name].data();

  // Add myself as provider for the field
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}
} // namespace scream
