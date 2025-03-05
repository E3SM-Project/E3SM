#include "diagnostics/exner.hpp"
#include "share/util/eamxx_common_physics_functions.hpp"

namespace scream
{

// =========================================================================================
ExnerDiagnostic::
ExnerDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
 : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void ExnerDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  auto nondim = Units::nondimensional();

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  auto scalar3d = grid->get_3d_scalar_layout(true);

  // The fields required for this diagnostic to be computed
  add_field<Required>("p_mid", scalar3d, Pa, grid_name);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d, nondim, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void ExnerDiagnostic::compute_diagnostic_impl()
{
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto& exner = m_diagnostic_output.get_view<Real**>();
  const auto& p_mid = get_field_in("p_mid").get_view<const Real**>();

  int nlevs = m_num_levs;
  Kokkos::parallel_for("ExnerDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*nlevs),
                       KOKKOS_LAMBDA(const int& idx) {
      const int icol = idx / nlevs;
      const int ilev = idx % nlevs;
      exner(icol,ilev) = PF::exner_function(p_mid(icol,ilev));
  });
  Kokkos::fence();
}
// =========================================================================================
} //namespace scream
