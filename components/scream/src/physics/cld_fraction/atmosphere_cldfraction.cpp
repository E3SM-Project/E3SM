#include "atmosphere_cldfraction.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
  using namespace cldfrac;
// =========================================================================================
CldFraction::CldFraction (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_cldfraction_comm (comm)
 , m_cldfraction_params (params)
{
  // Nothing to do here
}

// =========================================================================================
void CldFraction::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg/kg;
  Q.set_string("kg/kg");
  Units nondim(0,0,0,0,0,0,0);

  const auto& grid_name = m_cldfraction_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column
  m_nk_pack  = ekat::npack<Spack>(m_num_levs);

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // Set of fields used strictly as input
  m_required_fields.emplace("qi",   scalar3d_layout_mid, Q,      grid_name);
  m_required_fields.emplace("alst", scalar3d_layout_mid, nondim, grid_name);

  // Set of fields used strictly as output
  m_computed_fields.emplace("ast",   scalar3d_layout_mid, nondim, grid_name);
  m_computed_fields.emplace("aist",  scalar3d_layout_mid, nondim, grid_name);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.

}

// =========================================================================================
void CldFraction::initialize_impl (const util::TimeStamp& t0)
{
  m_current_ts = t0;
}

// =========================================================================================
void CldFraction::run_impl (const Real dt)
{
/* Luca, do we need any of the following code?
  // Copy inputs to host. Copy also outputs, cause we might "update" them, rather than overwrite them.
  for (auto& it : m_cldfraction_fields_in) {
    it.second.sync_to_host();
  }
  for (auto& it : m_cldfraction_fields_out) {
    it.second.sync_to_host();
  }

  // Copy outputs back to device
  // LB: why?!?
  for (auto& it : m_cldfraction_fields_out) {
    it.second.sync_to_dev();
  }
*/

  
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio. 
  auto qi   = m_cldfraction_fields_in["qi"].get_reshaped_view<const Pack**>();
  auto alst = m_cldfraction_fields_in["alst"].get_reshaped_view<const Pack**>();
  auto aist = m_cldfraction_fields_out["aist"].get_reshaped_view<Pack**>();
  auto ast  = m_cldfraction_fields_out["ast"].get_reshaped_view<Pack**>();

  cldfracF::cldfraction_main(m_num_cols,m_num_levs,qi,alst,aist,ast);

  // Get a copy of the current timestamp (at the beginning of the step) and
  // advance it,
  auto ts = timestamp();
  ts += dt;
  for (auto& f : m_cldfraction_fields_out) {
    f.second.get_header().get_tracking().update_time_stamp(ts);
  }

}

// =========================================================================================
void CldFraction::finalize_impl()
{
  // Do nothing
}

// =========================================================================================
void CldFraction::register_fields (FieldRepository<Real>& field_repo) const {
  std::set<ci_string> q_names =
    { "qv","qc","qr","qi","qm","nc","nr","ni","bm"};

  for (const auto& fid : m_required_fields) {
    const auto& name = fid.name();
    if (q_names.count(name)>0) {
      field_repo.register_field<Pack>(fid,"TRACERS");
    } else {
      field_repo.register_field<Pack>(fid);
    }
  }
  for (const auto& fid : m_computed_fields) {
    const auto& name = fid.name();
    if (q_names.count(name)>0) {
      field_repo.register_field<Pack>(fid,"TRACERS");
    } else {
      field_repo.register_field<Pack>(fid);
    }
  }
}

void CldFraction::set_required_field_impl (const Field<const Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking checks
  // at the beginning of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_cldfraction_fields_in.emplace(name,f);
  m_cldfraction_host_views_in[name] = f.get_view<Host>();
  m_raw_ptrs_in[name] = m_cldfraction_host_views_in[name].data();

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void CldFraction::set_computed_field_impl (const Field<      Real>& f) {
  // Store a copy of the field. We need this in order to do some tracking updates
  // at the end of the run call. Other than that, there would be really
  // no need to store a scream field here; we could simply set the view ptr
  // in the Homme's view, and be done with it.
  const auto& name = f.get_header().get_identifier().name();
  m_cldfraction_fields_out.emplace(name,f);
  m_cldfraction_host_views_out[name] = f.get_view<Host>();
  m_raw_ptrs_out[name] = m_cldfraction_host_views_out[name].data();

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream
