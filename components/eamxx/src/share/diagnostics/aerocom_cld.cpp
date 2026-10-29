#include "aerocom_cld.hpp"
#include "aerocom_weights.hpp"

namespace scream {

AeroComCld::
AeroComCld(const ekat::Comm &comm,
           const ekat::ParameterList &params,
           const std::shared_ptr<const AbstractGrid>& grid)
 : AbstractDiagnostic(comm, params, grid)
{
  m_field_name = m_params.get<std::string>("field_name");
  const auto position = m_params.get<std::string>("position");
  const auto phase    = m_params.get<std::string>("phase","none");

  // Depend on the (shared) weights, so that all fields sampled at the same
  // position and phase use one computation of them
  m_weights_name = AeroComWeights::expression(position,phase);
  m_diag_name = m_field_name + ".aerocom('" + position + "'" +
                (phase=="none" ? std::string("") : ",phase='" + phase + "'") + ")";

  if (m_field_name!="mask") {
    m_field_in_names.push_back(m_field_name);
  }
  m_field_in_names.push_back(m_weights_name);
}

void AeroComCld::initialize_impl()
{
  using namespace ShortFieldTagsNames;

  auto units = ekat::units::none;
  if (m_field_name!="mask") {
    const auto& f   = m_fields_in.at(m_field_name);
    const auto& fid = f.get_header().get_identifier();
    EKAT_REQUIRE_MSG (fid.get_layout().congruent(m_grid->get_3d_scalar_layout(LEV)) and
                      f.data_type()==DataType::RealType,
        "Error! AeroComCld needs a Real (COL,LEV) field.\n"
        " - field name: " + m_field_name + "\n"
        " - field layout: " + fid.get_layout().to_string() + "\n"
        " - field data type: " + e2str(f.data_type()) + "\n");
    units = fid.get_units();
  }

  FieldIdentifier fid(m_diag_name, m_grid->get_2d_scalar_layout(), units, m_grid->name());
  m_diagnostic_output = Field(fid,true);
}

void AeroComCld::compute_impl()
{
  const auto w = m_fields_in.at(m_weights_name).get_view<const Real**>();
  const auto x = m_diagnostic_output.get_view<Real*>();

  const int ncols = m_grid->get_num_local_dofs();
  const int nlevs = m_grid->get_num_vertical_levels();

  // For 'mask', X is 1
  const bool is_mask = m_field_name=="mask";
  Field::view_dev_t<const Real**> X;
  Field::view_dev_t<const int**> m;
  bool masked = false;
  if (not is_mask) {
    const auto& f = m_fields_in.at(m_field_name);
    X = f.get_view<const Real**>();
    masked = f.has_valid_mask();
    if (masked) {
      m = f.get_valid_mask().get_view<const int**>();
    }
  }

  using RangePolicy = Kokkos::RangePolicy<typename KokkosTypes<DefaultDevice>::ExeSpace>;
  Kokkos::parallel_for("AeroComCld", RangePolicy(0,ncols),
                       KOKKOS_LAMBDA(const int icol) {
    Real sum = 0;
    for (int k=0; k<nlevs; ++k) {
      if (w(icol,k)!=0 and (not masked or m(icol,k)!=0)) {
        sum += (is_mask ? 1 : X(icol,k))*w(icol,k);
      }
    }
    x(icol) = sum;
  });
}

}  // namespace scream
