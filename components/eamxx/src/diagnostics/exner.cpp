#include "diagnostics/exner.hpp"
#include "share/util/scream_common_physics_functions.hpp"

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
  using namespace ShortFieldTagsNames;

  auto nondim = Units::nondimensional();

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };

  // The fields required for this diagnostic to be computed
  add_field<Required>("p_mid", scalar3d_layout_mid, Pa, grid_name, SCREAM_PACK_SIZE);

  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d_layout_mid, nondim, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(SCREAM_PACK_SIZE);
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void ExnerDiagnostic::compute_diagnostic_impl()
{
  using Pack = ekat::Pack<Real,SCREAM_PACK_SIZE>;
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto npacks  = ekat::npack<Pack>(m_num_levs);
  const auto& exner = m_diagnostic_output.get_view<Pack**>();
  const auto& p_mid = get_field_in("p_mid").get_view<const Pack**>();

  Kokkos::parallel_for("ExnerDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*npacks),
                       KOKKOS_LAMBDA(const int& idx) {
      const int icol  = idx / npacks;
      const int jpack = idx % npacks;
      exner(icol,jpack) = PF::exner_function(p_mid(icol,jpack));
  });
  Kokkos::fence();
}
// =========================================================================================
} //namespace scream
