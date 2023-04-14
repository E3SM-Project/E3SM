#include "atmosphere_cld_fraction.hpp"
#include "share/property_checks/field_within_interval_check.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
  using namespace cld_fraction;
// =========================================================================================
CldFraction::CldFraction (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
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
  auto nondim = Units::nondimensional();

  m_grid = grids_manager->get_grid("Physics");
  const auto& grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = m_grid->get_num_vertical_levels();  // Number of levels per column

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
void CldFraction::initialize_impl (const RunType /* run_type */)
{
  // Set property checks for fields in this process
  using Interval = FieldWithinIntervalCheck;
  add_postcondition_check<Interval>(get_field_out("cldfrac_ice"),m_grid,0.0,1.0,false);
  add_postcondition_check<Interval>(get_field_out("cldfrac_tot"),m_grid,0.0,1.0,false);
}

// =========================================================================================
void CldFraction::run_impl (const double /* dt */)
{
  // Calculate ice cloud fraction and total cloud fraction given the liquid cloud fraction
  // and the ice mass mixing ratio. 
  auto qi   = get_field_in("qi").get_view<const Pack**>();
  auto liq_cld_frac = get_field_in("cldfrac_liq").get_view<const Pack**>();
  auto ice_cld_frac = get_field_out("cldfrac_ice").get_view<Pack**>();
  auto tot_cld_frac = get_field_out("cldfrac_tot").get_view<Pack**>();

  CldFractionFunc::main(m_num_cols,m_num_levs,qi,liq_cld_frac,ice_cld_frac,tot_cld_frac);
}

// =========================================================================================
void CldFraction::finalize_impl()
{
  // Do nothing
}
// =========================================================================================

} // namespace scream
