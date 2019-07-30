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
void P3Microphysics::set_grid(const std::shared_ptr<const GridsManager> grids_manager) 
{

  constexpr int NVL = 72;  /* TODO THIS NEEDS TO BE CHANGED TO A CONFIGURABLE */

  auto grid = grids_manager->get_grid("Physics");
  const int num_dofs = grid->num_dofs();
  const int nc = num_dofs;

  auto VL = FieldTag::VerticalLevel;
  auto CO = FieldTag::Column;

  FieldLayout scalar3d_layout { {CO,VL}, {nc,NVL} };

  // set requirements
  m_required_fields.emplace("P3_req_test",  scalar3d_layout, "Physics");
  // set computed
  m_computed_fields.emplace("P3_comq_test", scalar3d_layout, "Physics");

}
// =========================================================================================
void P3Microphysics::initialize ()
{
  p3_init_f90 ();
}

// =========================================================================================
void P3Microphysics::run ()
{
  p3_main_f90 ();
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
