#include "shortwave_cloud_forcing.hpp"

namespace scream
{

ShortwaveCloudForcing::
ShortwaveCloudForcing (const ekat::Comm& comm, const ekat::ParameterList& params,
                                  const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ekat::units;

  m_field_in_names.push_back("SW_flux_dn");
  m_field_in_names.push_back("SW_flux_up");
  m_field_in_names.push_back("SW_clrsky_flux_dn");
  m_field_in_names.push_back("SW_clrsky_flux_up");

  auto m2 = (m*m).rename("m2");

  auto diag_layout = m_grid->get_2d_scalar_layout();
  FieldIdentifier diag_fid (name(),diag_layout,W/m2,m_grid->name());

  m_diagnostic_output = Field(diag_fid);
  m_diagnostic_output.allocate_view();
}

void ShortwaveCloudForcing::compute_diagnostic_impl()
{
  using KT = KokkosTypes<DefaultDevice>;

  const int ncols = m_grid->get_num_local_dofs();

  const auto& SWCF              = m_diagnostic_output.get_view<Real*>();
  const auto& SW_flux_dn        = m_fields_in.at("SW_flux_dn").get_view<const Real**>();
  const auto& SW_flux_up        = m_fields_in.at("SW_flux_up").get_view<const Real**>();
  const auto& SW_clrsky_flux_dn = m_fields_in.at("SW_clrsky_flux_dn").get_view<const Real**>();
  const auto& SW_clrsky_flux_up = m_fields_in.at("SW_clrsky_flux_up").get_view<const Real**>();

  typename KT::RangePolicy policy (0,ncols);
  auto lambda = KOKKOS_LAMBDA(const int icol) {
    SWCF(icol) = (SW_flux_dn(icol,0) - SW_flux_up(icol,0)) - (SW_clrsky_flux_dn(icol,0) - SW_clrsky_flux_up(icol,0));
  };
  Kokkos::parallel_for("ShortwaveCloudForcing",policy,lambda);
}

} //namespace scream
