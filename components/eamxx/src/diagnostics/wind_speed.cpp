#include "diagnostics/wind_speed.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>

namespace scream
{

WindSpeed::
WindSpeed (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  // Nothing to do here
}

void WindSpeed::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  using namespace ekat::units;

  auto grid  = grids_manager->get_grid("Physics");
  const auto& grid_name = grid->name();

  m_ncols = grid->get_num_local_dofs();
  m_nlevs = grid->get_num_vertical_levels();

  auto scalar3d = grid->get_3d_scalar_layout(true);
  auto vector3d = grid->get_3d_vector_layout(true,2);

  // The fields required for this diagnostic to be computed
  add_field<Required>("horiz_winds", vector3d, Pa, grid_name);

  // Construct and allocate the 3d wind_speed field
  FieldIdentifier fid (name(), scalar3d, m/s, grid_name);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void WindSpeed::compute_diagnostic_impl()
{
  using KT = KokkosTypes<DefaultDevice>;
  using RP = typename KT::RangePolicy;

  const auto uv = get_field_in("horiz_winds").get_view<const Real***>();
  const auto ws = m_diagnostic_output.get_view<Real**>();

  const int nlevs = m_nlevs;
  Kokkos::parallel_for("Compute " + name(), RP(0,m_nlevs*m_ncols),
                       KOKKOS_LAMBDA(const int& idx) {
    const int icol = idx / nlevs;
    const int ilev = idx % nlevs;
    const auto& u = uv(icol,0,ilev);
    const auto& v = uv(icol,1,ilev);
    ws (icol,ilev) = sqrt(u*u + v*v);
  });
}

} //namespace scream
