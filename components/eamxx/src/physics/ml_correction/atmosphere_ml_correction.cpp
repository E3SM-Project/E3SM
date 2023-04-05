#include <array>

#include "atmosphere_ml_correction.hpp"
#include "ekat/ekat_assert.hpp"
#include "ekat/util/ekat_units.hpp"

namespace scream {
// =========================================================================================
MLCorrection::MLCorrection(const ekat::Comm &comm,
                           const ekat::ParameterList &params)
    : AtmosphereProcess(comm, params) {
  // Nothing to do here
}

// =========================================================================================
void MLCorrection::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  // The units of mixing ratio Q are technically non-dimensional.
  // Nevertheless, for output reasons, we like to see 'kg/kg'.
  auto Q = kg / kg;
  Q.set_string("kg/kg");

  m_grid                = grids_manager->get_grid("Physics");
  const auto &grid_name = m_grid->name();
  m_num_cols = m_grid->get_num_local_dofs();  // Number of columns on this rank
  m_num_levs =
      m_grid->get_num_vertical_levels();  // Number of levels per column

  // Define the different field layouts that will be used for this process

  // Layout for 3D (2d horiz X 1d vertical) variable defined at mid-level and
  // interfaces
  FieldLayout scalar3d_layout_mid{{COL, LEV}, {m_num_cols, m_num_levs}};

  // Set of fields used strictly as input
  add_field<Required>("qv", scalar3d_layout_mid, Q, grid_name, "tracers");

  // Set of fields used strictly as output
  add_field<Computed>("qv_nudging_tend", scalar3d_layout_mid, Q, grid_name);

  // Set of fields used as input and output
  // - There are no fields used as both input and output.
}

// =========================================================================================
void MLCorrection::initialize_impl(const RunType /* run_type */) {
  // Nothing to do
}

// =========================================================================================
void MLCorrection::run_impl(const double /* dt */) {
  auto qv         = get_field_in("qv").get_view<const Real **>();
  auto qv_nudging = get_field_out("qv_nudging").get_view<Real **>();

  // ML correction proceess is not yet implemented
}

// =========================================================================================
void MLCorrection::finalize_impl() {
  // Do nothing
}
// =========================================================================================

}  // namespace scream
