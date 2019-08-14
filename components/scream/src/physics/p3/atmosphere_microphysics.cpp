#include "share/scream_assert.hpp"
#include "physics/p3/scream_p3_interface.hpp"
#include "physics/p3/atmosphere_microphysics.hpp"

namespace scream
{

// =========================================================================================
P3Microphysics::P3Microphysics (const Comm& comm,const ParameterList& /* params */)
  : m_p3_comm (comm)
{
/* Anything that can be initialized without grid information can be initialized here.
 * Like universal constants, table lookups, p3 options.
*/
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
  constexpr int QSZ =  9;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */

  auto grid = grids_manager->get_grid("Physics");
  const int num_dofs = grid->num_dofs();
  const int nc = num_dofs;

  auto VL = FieldTag::VerticalLevel;
  auto CO = FieldTag::Column;
  auto VR = FieldTag::Variable;
  auto TL = FieldTag::TimeLevel;

  FieldLayout scalar3d_layout { {CO,VL}, {nc,NVL} }; // Note that C++ and Fortran read array dimensions in reverse
  FieldLayout vector3d_layout { {CO,VR,VL}, {nc,QSZ,NVL} };
  FieldLayout tracers_state_layout { {CO,TL,VR,VL}, {nc,2,4,NVL} };
  FieldLayout scalar_state_3d_mid_layout { {CO,TL,VL} , {nc,2,NVL}};

  // set requirements
//  m_required_fields.emplace("P3_req_test",  scalar3d_layout,   units::one, "Physics");
//  m_required_fields.emplace("ASD_req_test",  scalar3d_layout,   units::one, "Physics");
  m_required_fields.emplace("dp"         , scalar_state_3d_mid_layout,      Pa, "Physics");
  m_required_fields.emplace("qdp"        , tracers_state_layout,      Qdp, "Physics");
  // set computed
  m_computed_fields.emplace("P3_comq_test", scalar3d_layout, units::one, "Physics");
  m_computed_fields.emplace("q"           , vector3d_layout,          Q, "Physics");

}
// =========================================================================================
void P3Microphysics::initialize (const util::TimeStamp& t0)
{
  m_time = t0;
  auto q_ptr = m_p3_fields_out.at("q").get_view().data();

  p3_init_f90 (q_ptr);

}

// =========================================================================================
void P3Microphysics::run (const double dt)
{
  auto q_ptr = m_p3_fields_out.at("q").get_view().data();
  auto qdp_ptr = m_p3_fields_in.at("qdp").get_view().data();
  auto dp_ptr = m_p3_fields_in.at("dp").get_view().data();

  p3_main_f90 (dt,q_ptr,qdp_ptr);

  m_time += dt;
  m_p3_fields_out.at("q").get_header().get_tracking().update_time_stamp(m_time);
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
}

void P3Microphysics::set_computed_field_impl (const Field<      Real, device_type>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_p3_fields_out.emplace(f.get_header().get_identifier().name(),f);
}

} // namespace scream
