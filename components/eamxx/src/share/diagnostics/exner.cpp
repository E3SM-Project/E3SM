#include "exner.hpp"
#include "share/physics/eamxx_common_physics_functions.hpp"

namespace scream
{

Exner::Exner (const ekat::Comm& comm, const ekat::ParameterList& params,
              const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm,params,grid)
{
  // The fields required for this diagnostic to be computed
  m_field_in_names.push_back("p_mid");
}

void Exner::initialize_impl()
{
  const auto&p = m_fields_in.at("p_mid");

  // Construct and allocate the diagnostic field
  auto fid = p.get_header().get_identifier().clone(name()).reset_units(ekat::units::none);
  m_diagnostic_output = Field(fid);
  m_diagnostic_output.allocate_view();
}

void Exner::compute_impl()
{
  using PF = PhysicsFunctions<DefaultDevice>;

  const auto& exner = m_diagnostic_output.get_view<Real**>();
  const auto& p_mid = m_fields_in.at("p_mid").get_view<const Real**>();

  int nlevs = m_grid->get_num_vertical_levels();
  int ncols = m_grid->get_num_local_dofs();
  Kokkos::parallel_for("Exner",
                       Kokkos::RangePolicy<>(0,ncols*nlevs),
                       KOKKOS_LAMBDA(const int& idx) {
      const int icol = idx / nlevs;
      const int ilev = idx % nlevs;
      exner(icol,ilev) = PF::exner_function(p_mid(icol,ilev));
  });
}

} //namespace scream
