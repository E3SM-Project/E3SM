#include "ekat/ekat_assert.hpp"
#include "physics/shoc/scream_shoc_interface.hpp"
#include "physics/shoc/atmosphere_macrophysics.hpp"

namespace scream
{

// =========================================================================================
SHOCMacrophysics::SHOCMacrophysics (const ekat::Comm& comm,const ekat::ParameterList& params)
  : m_shoc_comm   (comm)
  , m_shoc_params (params)
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

  const auto& grid_name = m_shoc_params.get<std::string>("Grid");
  auto grid = grids_manager->get_grid(grid_name);

  const int ncols = grid->get_num_local_dofs();
  const int nlevs = grid->get_num_vertical_levels();

  auto VL = FieldTag::VerticalLevel;
  auto CO = FieldTag::Column;

  FieldLayout scalar3d_layout { {CO,VL}, {ncols,nlevs} }; // Note that C++ and Fortran read array dimensions in reverse

  // Input
  m_required_fields.emplace("dp", scalar3d_layout,  Pa, grid->name());

  // Input-output
  m_inout_groups_req.emplace("TRACERS",grid->name());
}
// =========================================================================================
void SHOCMacrophysics::
set_updated_group (const FieldGroup<Real>& group)
{
  EKAT_REQUIRE_MSG(group.m_info->size()>0,
    "Error! We were not expecting an empty field group.\n");

  const auto& name = group.m_info->m_group_name;

  EKAT_REQUIRE_MSG(name=="TRACERS",
    "Error! We were not expecting a field group called '" << name << "\n");

  EKAT_REQUIRE_MSG(group.m_info->m_bundled,
      "Error! Shoc expects bundled fields for tracers.\n");

  m_shoc_groups_inout.emplace(name,group);
}

// =========================================================================================
void SHOCMacrophysics::initialize_impl (const util::TimeStamp& /* t0 */)
{
  // TODO: create the host mirrors once.
  auto Q = m_shoc_groups_inout.at("TRACERS").m_bundle;
  Q->sync_to_host();

  auto q_ptr = Q->get_view<Host>().data();

  shoc_init_f90 (q_ptr);

  Q->sync_to_dev();
}

// =========================================================================================
void SHOCMacrophysics::run_impl (const Real dt)
{
  // Get bundled tracers
  auto Q  = m_shoc_groups_inout.at("TRACERS").m_bundle;

  Q->sync_to_host();

  auto q_ptr  = Q->get_view<Host>().data();

  shoc_main_f90 (dt,q_ptr);

  Q->sync_to_dev();

  // Get the beginning-of-step timestamp and advance it to update our fields.
  auto ts = timestamp();
  ts += dt;
  Q->get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
void SHOCMacrophysics::finalize_impl()
{
  shoc_finalize_f90 ();
}
// =========================================================================================

void SHOCMacrophysics::register_fields (FieldRepository<Real>& field_repo) const {
  for (auto& fid : m_required_fields) {
    field_repo.register_field(fid);
  }
  for (auto& fid : m_computed_fields) {
    field_repo.register_field(fid);
  }
}

void SHOCMacrophysics::set_required_field_impl (const Field<const Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_shoc_fields_in.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as customer to the field
  f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
}

void SHOCMacrophysics::set_computed_field_impl (const Field<Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  m_shoc_fields_out.emplace(f.get_header().get_identifier().name(),f);

  // Add myself as provider for the field
  f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
}

} // namespace scream
