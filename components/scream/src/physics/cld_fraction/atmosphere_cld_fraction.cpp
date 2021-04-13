#include "atmosphere_cld_fraction.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
  using namespace cld_fraction;
// =========================================================================================
CldFraction::CldFraction (const ekat::Comm& comm, const ekat::ParameterList& params)
 : m_cldfraction_comm (comm)
 , m_cld_fraction_params (params)
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

  const auto& grid_name = m_cld_fraction_params.get<std::string>("Grid");
  auto grid  = grids_manager->get_grid(grid_name);
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process
  using namespace ShortFieldTagsNames;

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // Set of fields used strictly as input
  add_field<Required>("qi",   scalar3d_layout_mid, Q,      grid_name);
  add_field<Required>("cldfrac_liq", scalar3d_layout_mid, nondim, grid_name);

  // Set of fields used strictly as output
  add_field<Computed>("cldfrac_tot",   scalar3d_layout_mid, nondim, grid_name);
  add_field<Computed>("cldfrac_ice",  scalar3d_layout_mid, nondim, grid_name);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.
}

// =========================================================================================
void CldFraction::initialize_impl (const util::TimeStamp& /* t0 */)
{
}

// =========================================================================================
void CldFraction::run_impl (const Real dt)
{
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio. 
  auto qi   = m_cld_fraction_fields_in["qi"].get_reshaped_view<const Pack**>();
  auto liq_cld_frac = m_cld_fraction_fields_in["cldfrac_liq"].get_reshaped_view<const Pack**>();
  auto ice_cld_frac = m_cld_fraction_fields_out["cldfrac_ice"].get_reshaped_view<Pack**>();
  auto tot_cld_frac = m_cld_fraction_fields_out["cldfrac_tot"].get_reshaped_view<Pack**>();

  CldFractionFunc::main(m_num_cols,m_num_levs,qi,liq_cld_frac,ice_cld_frac,tot_cld_frac);

  // Get a copy of the current timestamp (at the beginning of the step) and
  // advance it,
  auto ts = timestamp();
  ts += dt;
  for (auto& f : m_cld_fraction_fields_out) {
    f.second.get_header().get_tracking().update_time_stamp(ts);
  }
}

// =========================================================================================
void CldFraction::finalize_impl()
{
  // Do nothing
}

// =========================================================================================
void CldFraction::
register_fields (const std::map<std::string,std::shared_ptr<FieldManager<Real>>>& field_mgrs) const {
  const auto& grid_name = m_cld_fraction_params.get<std::string>("Grid");
  auto& field_mgr = *field_mgrs.at(grid_name);
  for (const auto& req : get_required_fields()) {
    const auto& fid = req.fid;
    const auto& name = fid.name();
    if (name == "qi") {
      field_mgr.register_field<Pack>(fid,"TRACERS");
    } else {
      field_mgr.register_field<Pack>(fid);
    }
  }
  for (const auto& req : get_computed_fields()) {
    const auto& fid = req.fid;
    const auto& name = fid.name();
    if (name == "qi") {
      field_mgr.register_field<Pack>(fid,"TRACERS");
    } else {
      field_mgr.register_field<Pack>(fid);
    }
  }
}

void CldFraction::set_required_field_impl (const Field<const Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_cld_fraction_fields_in.emplace(name,f);

  // Add myself as customer to the field
  add_me_as_customer(f);
}

void CldFraction::set_computed_field_impl (const Field<      Real>& f) {

  const auto& name = f.get_header().get_identifier().name();
  m_cld_fraction_fields_out.emplace(name,f);

  // Add myself as provider for the field
  add_me_as_provider(f);
}

} // namespace scream
