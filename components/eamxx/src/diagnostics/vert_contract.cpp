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
  const auto g   = grids_manager->get_grid(gn);

  add_field<Required>(fn, gn);

  // we support either sum or avg
  m_contract_method = m_params.get<std::string>("contract_method");
  // we support either dp or dz weighting
  m_weighting_method = m_params.get<std::string>("weighting_method");

  m_diag_name = fn + m_contract_method + "_" + m_weighting_method;

  auto scalar3d = g->get_3d_scalar_layout(true);
  if(m_weighting_method == "dp") {
    add_field<Required>("pseudo_density", scalar3d, Pa, gn);
  } else if(m_weighting_method == "dz") {
    add_field<Required>("dz", scalar3d, m, gn);
  } else {
    EKAT_ERROR_MSG("Unsupported weighting method: " + m_weighting_method + "\n");
  }
}

void VertContractDiag::initialize_impl(const RunType /*run_type*/) {
  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  const auto &f      = get_fields_in().front();
  const auto &fid    = f.get_header().get_identifier();
  const auto &layout = fid.get_layout();

  EKAT_REQUIRE_MSG(layout.rank() >= 1 && layout.rank() <= 3,
                   "Error! Field rank not supported by VertContractDiag.\n"
                   " - field name: " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(layout.tags().back() == LEV,
                   "Error! VertContractDiag diagnostic expects a layout ending "
                   "with the 'LEV' tag.\n"
                   " - field name  : " + fid.name() + "\n"
                   " - field layout: " + layout.to_string() + "\n");

  ekat::units::Units diag_units = ekat::units::Units::nondimensional();
  if (m_weighting_method == "dp") {
    m_weighting = get_field_in("pseudo_density").clone("weighting");
    if (m_contract_method == "avg") {
      // we scale by the weighting/sum(weighting), so we return the unit of the field
      diag_units = fid.get_units();
    } else { // m_contract_method = sum
      // we scale by the weighting, so we use fid units * weighting (but we scale by 1/g for dp)
      diag_units = fid.get_units() * m_weighting.get_header().get_identifier().get_units() / (m / (s*s));
    }
  } else { // m_weighting_method = "dz"
    m_weighting = get_field_in("dz").clone("weighting");
    if (m_contract_method == "avg") {
      // we scale by the weighting/sum(weighting), so we return the unit of the field
      diag_units = fid.get_units();
    } else { // m_contract_method = sum
      // we scale by the weighting, so we use fid units * weighting units
      diag_units = fid.get_units() * m_weighting.get_header().get_identifier().get_units();
    }
  }

  FieldIdentifier d_fid(m_diag_name, layout.clone().strip_dim(LEV), diag_units,
                        fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();
}

void VertContractDiag::compute_diagnostic_impl() {
  const auto &f = get_fields_in().front();
  const auto &d = m_diagnostic_output;

  // update the weights; if weighting by dp, we need to scale by 1/g
  if (m_weighting_method == "dp") {
    auto g = scream::physics::Constants<Real>::gravit;
    m_weighting.update(get_field_in("pseudo_density"), 1 / g, sp(0));
  } else if(m_weighting_method == "dz") {
    m_weighting.update(get_field_in("dz"), 1, 0);
  }

  // if "avg" is in the method name, we need to scale the weighting by its 1/sum
  if (m_contract_method.find("avg") != std::string::npos) {
    auto sum = field_sum<Real>(m_weighting, &m_comm);
    m_weighting.scale(sp(1.0) / sum);
  }

  // call the vert_contraction impl that will take care of everything
  vert_contraction<Real>(d, f, m_weighting, &m_comm);
}

}  // namespace scream
