#include "atm_density.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

namespace scream
{

AtmDensity::
AtmDensity(const ekat::Comm& comm, const ekat::ParameterList& params,
           const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  using namespace ShortFieldTagsNames;

  // The fields required for this diagnostic to be computed
  m_field_in_names.push_back("T_mid");
  m_field_in_names.push_back("pseudo_density");
  m_field_in_names.push_back("p_mid");
  m_field_in_names.push_back("qv");
}

void AtmDensity::initialize_impl()
{
  using namespace ekat::units;

  // Construct and allocate the diagnostic field
  const auto& T = m_fields_in.at("T_mid");
  auto fid = T.get_header().get_identifier().clone(name()).reset_units(kg/(m*m*m));
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void AtmDensity::compute_diagnostic_impl()
{
  using PF = scream::PhysicsFunctions<DefaultDevice>;

  const auto atm_dens = m_diagnostic_output.get_view<Real**>();
  const auto T    = m_fields_in.at("T_mid").get_view<const Real**>();
  const auto p    = m_fields_in.at("p_mid").get_view<const Real**>();
  const auto qv   = m_fields_in.at("qv").get_view<const Real**>();
  const auto pseudo_density = m_fields_in.at("pseudo_density").get_view<const Real**>();

  int nlevs = m_grid->get_num_vertical_levels();
  int ncols = m_grid->get_num_local_dofs();
  Kokkos::parallel_for("AtmosphereDensityDiagnostic",
                       Kokkos::RangePolicy<>(0,ncols*nlevs),
                       KOKKOS_LAMBDA(const int& idx) {
      const int icol  = idx / nlevs;
      const int ilev = idx % nlevs;
      auto dz = PF::calculate_dz(pseudo_density(icol,ilev),p(icol,ilev),T(icol,ilev),qv(icol,ilev));
      atm_dens(icol,ilev) = PF::calculate_density(pseudo_density(icol,ilev),dz);
  });
}

} //namespace scream
