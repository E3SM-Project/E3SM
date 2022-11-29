#include "diagnostics/ice_cloud_mask.hpp"
#include "physics/share/physics_functions.hpp" // also for ETI not on GPUs

namespace scream
{


// =========================================================================================
IceCloudMaskDiagnostic::IceCloudMaskDiagnostic (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

// =========================================================================================
void IceCloudMaskDiagnostic::set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;


  auto Q = kg/kg;
  Q.set_string("kg/kg");

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();
  m_num_cols = grid->get_num_local_dofs(); // Number of columns on this rank
  m_num_levs = grid->get_num_vertical_levels();  // Number of levels per column

  FieldLayout scalar3d_layout_mid { {COL,LEV}, {m_num_cols,m_num_levs} };
  constexpr int ps = Pack::n;

  // The fields required for this diagnostic to be computed
  add_field<Required>("qi",          scalar3d_layout_mid, Q,  grid_name, "tracers", ps);


  // Construct and allocate the diagnostic field
  FieldIdentifier fid (name(), scalar3d_layout_mid, K, grid_name);
  m_diagnostic_output = Field(fid);
  auto& C_ap = m_diagnostic_output.get_header().get_alloc_properties();
  C_ap.request_allocation(ps);
  m_diagnostic_output.allocate_view();
}
// =========================================================================================
void IceCloudMaskDiagnostic::initialize_impl(const RunType /* run_type */)
{
  auto ts = timestamp(); 
  m_diagnostic_output.get_header().get_tracking().update_time_stamp(ts);
}
// =========================================================================================
void IceCloudMaskDiagnostic::compute_diagnostic_impl()
{
  const auto npacks  = ekat::npack<Pack>(m_num_levs);
  auto qi_mid    = get_field_in("qi").get_view<const Pack**>();
  const auto& ice_cld_mask = m_diagnostic_output.get_view<Pack**>();

  Kokkos::parallel_for("IceCloudMaskDiagnostic",
                       Kokkos::RangePolicy<>(0,m_num_cols*npacks),
                       KOKKOS_LAMBDA (const int& idx) {
      const int icol  = idx / npacks;
      const int jpack = idx % npacks;
      const Real ice_frac_threshold = 1e-5;
      auto icecld = qi_mid(icol,jpack) > ice_frac_threshold;
      ice_cld_mask(icol,jpack) = 0.0;
      ice_cld_mask(icol,jpack).set(icecld, 1.0);
  });
  Kokkos::fence();
}
// =========================================================================================

} //namespace scream
