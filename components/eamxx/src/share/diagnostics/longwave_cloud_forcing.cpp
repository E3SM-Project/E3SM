#include "longwave_cloud_forcing.hpp"

namespace scream
{

LongwaveCloudForcing::
LongwaveCloudForcing (const ekat::Comm& comm, const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid)
  : AbstractDiagnostic(comm,params,grid)
{
  using namespace ekat::units;

  // The fields required for this diagnostic to be computed
  m_field_in_names.push_back("LW_flux_up");
  m_field_in_names.push_back("LW_clrsky_flux_up");

  auto m2 = (m*m).rename("m2");
  auto diag_layout = m_grid->get_2d_scalar_layout();
  FieldIdentifier diag_fid(name(),diag_layout,W/m2,m_grid->name());

  m_diagnostic_output = Field(diag_fid);
  m_diagnostic_output.allocate_view();
}

void LongwaveCloudForcing::compute_diagnostic_impl()
{
  using KT = KokkosTypes<DefaultDevice>;

  const int ncols = m_grid->get_num_local_dofs();

  const auto& LWCF              = m_diagnostic_output.get_view<Real*>();
  const auto& LW_flux_up        = m_fields_in.at("LW_flux_up").get_view<const Real**>();
  const auto& LW_clrsky_flux_up = m_fields_in.at("LW_clrsky_flux_up").get_view<const Real**>();

  typename KT::RangePolicy policy (0,ncols);
  auto lambda = KOKKOS_LAMBDA(const int icol) {
    LWCF(icol) =  LW_clrsky_flux_up(icol,0) -  LW_flux_up(icol,0) ;
  };
  Kokkos::parallel_for("LongwaveCloudForcing",policy,lambda);
}

} //namespace scream
