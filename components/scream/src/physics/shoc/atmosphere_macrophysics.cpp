#include "ekat/ekat_assert.hpp"
#include "physics/shoc/scream_shoc_interface.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"

namespace scream
{

// =========================================================================================
SHOCMacrophysics::SHOCMacrophysics (const ekat::Comm& comm,const ekat::ParameterList& /* params */)
  : m_shoc_comm (comm)
{
/* Anything that can be initialized without grid information can be initialized here.
 * Like universal constants, shoc options.
*/
}

// =========================================================================================
void SHOCMacrophysics::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");

  constexpr int NVL = SCREAM_NUM_VERTICAL_LEV;
  constexpr int QSZ =  35;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */

  auto grid = grids_manager->get_grid("Physics");
  const int num_dofs = grid->get_num_local_dofs();
  const int nc = num_dofs;

  auto VL = FieldTag::VerticalLevel;
  auto CO = FieldTag::Column;
  auto VR = FieldTag::Variable;

  FieldLayout scalar3d_layout { {CO,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout { {CO,VR,VL}, {nc,QSZ,NVL} };
  FieldLayout q_forcing_layout  { {CO,VR,VL}, {nc,QSZ,NVL} };

  // Input
  m_required_fields.emplace("dp", scalar3d_layout,  Pa, grid->name());
  // Input-Output
  m_required_fields.emplace("q",  vector3d_layout,  Q,  grid->name());
  m_required_fields.emplace("FQ", q_forcing_layout, Q,  grid->name());

  m_computed_fields.emplace("q" , vector3d_layout,  Q,  grid->name());
  m_computed_fields.emplace("FQ", q_forcing_layout, Q,  grid->name());

}
// =========================================================================================
void SHOCMacrophysics::initialize_impl (const util::TimeStamp& t0)
{
  // TODO: create the host mirrors once.
  auto q_dev = m_shoc_fields_out.at("q").get_view();
  auto q_host = Kokkos::create_mirror_view(q_dev);
  Kokkos::deep_copy(q_host,q_dev);
  auto q_ptr = q_host.data();

  shoc_init_f90 (q_ptr);

  Kokkos::deep_copy(q_dev,q_host);
}

// =========================================================================================
void SHOCMacrophysics::run_impl (const Real dt)
{
  // TODO: create the host mirrors once.
  auto q_dev   = m_shoc_fields_out.at("q").get_view();
  auto fq_dev  = m_shoc_fields_out.at("FQ").get_view();

  auto q_host   = Kokkos::create_mirror_view(q_dev);
  auto fq_host  = Kokkos::create_mirror_view(fq_dev);

  auto q_ptr   = q_host.data();
  auto fq_ptr  = fq_host.data();

  shoc_main_f90 (dt,q_ptr,fq_ptr);

  Kokkos::deep_copy(q_dev,q_host);
  Kokkos::deep_copy(fq_dev,fq_host);

  // Get the beginning-of-step timestamp and advance it to update our fields.
  auto ts = timestamp();
  ts += dt;
  m_shoc_fields_out.at("q").get_header().get_tracking().update_time_stamp(ts);
  m_shoc_fields_out.at("FQ").get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
void SHOCMacrophysics::finalize_impl()
{
  shoc_finalize_f90 ();
}
// =========================================================================================

void SHOCMacrophysics::register_fields (FieldRepository<Real, device_type>& field_repo) const {
  for (auto& fid : m_required_fields) {
    field_repo.register_field(fid);
  }
  for (auto& fid : m_computed_fields) {
    field_repo.register_field(fid);
  }
}

void SHOCMacrophysics::set_required_field_impl (const Field<const Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_shoc_fields_in.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as customer to the field
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void SHOCMacrophysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_shoc_fields_out.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as provider for the field
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

} // namespace scream
