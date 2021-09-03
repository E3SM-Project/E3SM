#include "atmosphere_cld_fraction.hpp"
#include "share/field/field_property_checks/field_within_interval_check.hpp"

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
  using namespace ShortFieldTagsNames;

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

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and interfaces 
  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // Set of fields used strictly as input
  constexpr int ps = Pack::n;
  add_field<Required>("qi",          scalar3d_layout_mid, Q,      grid_name,"tracers",ps);
  add_field<Required>("cldfrac_liq", scalar3d_layout_mid, nondim, grid_name,ps);

  // Set of fields used strictly as output
  add_field<Computed>("cldfrac_tot",  scalar3d_layout_mid, nondim, grid_name,ps);
  add_field<Computed>("cldfrac_ice",  scalar3d_layout_mid, nondim, grid_name,ps);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.
}

// =========================================================================================
void CldFraction::initialize_impl (const util::TimeStamp& /* t0 */)
{
  // Set property checks for fields in this process
  auto frac_interval_check = std::make_shared<FieldWithinIntervalCheck<Real> >(0,1);
  m_fields_out["cldfrac_ice"].add_property_check(frac_interval_check);
  m_fields_out["cldfrac_tot"].add_property_check(frac_interval_check);
}

// =========================================================================================
void CldFraction::run_impl (const Real dt)
{
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio. 
  auto qi   = m_fields_in["qi"].get_view<const Pack**>();
  auto liq_cld_frac = m_fields_in["cldfrac_liq"].get_view<const Pack**>();
  auto ice_cld_frac = m_fields_out["cldfrac_ice"].get_view<Pack**>();
  auto tot_cld_frac = m_fields_out["cldfrac_tot"].get_view<Pack**>();

  CldFractionFunc::main(m_num_cols,m_num_levs,qi,liq_cld_frac,ice_cld_frac,tot_cld_frac);

  // Get a copy of the current timestamp (at the beginning of the step) and
  // advance it,
  auto ts = timestamp();
  ts += dt;
  for (auto& f : m_fields_out) {
    f.second.get_header().get_tracking().update_time_stamp(ts);
  }
}

// =========================================================================================
void CldFraction::finalize_impl()
{
  // Do nothing
}
// =========================================================================================

} // namespace scream
