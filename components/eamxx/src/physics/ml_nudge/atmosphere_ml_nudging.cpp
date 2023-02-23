#include "atmosphere_ml_nudging.hpp"

#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

#include <array>

namespace scream
{
  using namespace cld_fraction;
// =========================================================================================
MLNudging::MLNudging (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

// =========================================================================================
void MLNudging::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
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
  add_field<Required>("qv",          scalar3d_layout_mid, Q,      grid_name,"tracers",ps);

  // Set of fields used strictly as output
  add_field<Computed>("qv_nudging_tend",  scalar3d_layout_mid, Q, grid_name,ps);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.
}

// =========================================================================================
void MLNudging::initialize_impl (const RunType /* run_type */)
{
  // Set property checks for fields in this process
  using Interval = FieldWithinIntervalCheck;
  add_postcondition_check<Interval>(get_field_out("qv_nudging"),m_grid,0.0,1e-3,false);
}

// =========================================================================================
void MLNudging::run_impl (const double /* dt */)
{
  auto qv          = get_field_in("qv").get_view<const Pack**>();
  auto qv_nudging = get_field_out("qv_nudging").get_view<Pack**>();

// STRATEGY:: 
// 1. Try calling any type of python code here, ignoring the SCREAM variables and data
//    structures.
// 2. Then try passing something from SCREAM, in this case the view to `qv_nudging`, a
//    made-up variable for this test.  See if we can change the value of qv_nudging in
//    python.
// 3. Once we've accomplished 1 & 2 we can start experimenting with real SCREAM values,
//    updating the state in python, grabbing all the variables we want to nudge.

}

// =========================================================================================
void MLNudging::finalize_impl()
{
  // Do nothing
}
// =========================================================================================

} // namespace scream
