#include "wind_speed.hpp"

namespace scream
{

WindSpeed::
WindSpeed (const ekat::Comm& comm, const ekat::ParameterList& params,
           const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ekat::units;
  using namespace ShortFieldTagsNames;

  const auto& grid_name = m_grid->name();

  auto scalar3d = m_grid->get_3d_scalar_layout(LEV);

  m_field_in_names.push_back("horiz_winds");

  FieldIdentifier fid (name(), scalar3d, m/s, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void WindSpeed::compute_diagnostic_impl()
{
  using KT      = KokkosTypes<DefaultDevice>;
  using MDRange = Kokkos::MDRangePolicy<typename KT::ExeSpace,Kokkos::Rank<2>>;

  const auto uv = m_fields_in.at("horiz_winds").get_view<const Real***>();
  const auto ws = m_diagnostic_output.get_view<Real**>();

  const int ncols = m_grid->get_num_local_dofs();
  const int nlevs = m_grid->get_num_vertical_levels();
  auto lambda = KOKKOS_LAMBDA(const int icol, const int ilev) {
    const auto& u = uv(icol,0,ilev);
    const auto& v = uv(icol,1,ilev);
    ws (icol,ilev) = sqrt(u*u + v*v);
  };
  MDRange policy ({0,0},{ncols,nlevs});
  Kokkos::parallel_for("Compute " + name(), policy, lambda);
}

} //namespace scream
