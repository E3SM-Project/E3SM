#include "diagnostics/vert_contract.hpp"

#include "physics/share/physics_constants.hpp"
#include "share/field/field_utils.hpp"

namespace scream {

VertContractDiag::VertContractDiag(const ekat::Comm &comm,
                                   const ekat::ParameterList &params)
    : AtmosphereDiagnostic(comm, params) {}

void VertContractDiag::set_grids(
    const std::shared_ptr<const GridsManager> grids_manager) {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &fn = m_params.get<std::string>("field_name");
  const auto &gn = m_params.get<std::string>("grid_name");
  const auto g   = grids_manager->get_grid("Physics");

  add_field<Required>(fn, gn);

  // We support either sum or avg 
  m_contract_method = m_params.get<std::string>("contract_method");

  m_diag_name       = fn + m_contract_method;

  auto scalar3d = g->get_3d_scalar_layout(true);
  add_field<Required>("pseudo_density", scalar3d, Pa, gn);
}

void VertContractDiag::initialize_impl(const RunType /*run_type*/) {
  using namespace ShortFieldTagsNames;
  using PC           = scream::physics::Constants<Real>;
  const auto &f      = get_fields_in().front();
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  constexpr Real g = PC::gravit;

  EKAT_REQUIRE_MSG(layout.rank() >= 1 && layout.rank() <= 3,
                   "Error! Field rank not supported by VertContractDiag.\n"
                   " - field name: " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags().back() == LEV,
                   "Error! VertContractDiag diagnostic expects a layout ending "
                   "with the 'LEV' tag.\n"
                   " - field name  : " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");

  FieldIdentifier d_fid(m_diag_name, layout.clone().strip_dim(LEV),
                        fid.get_units(), fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // scale the weighting field (pseudo density) by g
  m_weighting = get_field_in("pseudo_density");
  m_weighting.scale(sp(1.0) / g);
  // if "avg" is in the method name, we need to scale the weighting by its sum
  if(m_contract_method.find("avg") != std::string::npos) {
    auto sum = field_sum<Real>(m_weighting, &m_comm);
    m_weighting.scale(sp(1.0) / sum);
  }
}

void VertContractDiag::compute_diagnostic_impl() {
  const auto &f = get_fields_in().front();
  const auto &d = m_diagnostic_output;
  // Call the vert_contraction impl that will take care of everything
  vert_contraction<Real>(d, f, m_weighting, &m_comm);
}

}  // namespace scream
