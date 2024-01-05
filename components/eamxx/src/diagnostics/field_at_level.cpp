#include "diagnostics/field_at_level.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

namespace scream
{

// =========================================================================================
FieldAtLevel::FieldAtLevel (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
{
  const auto& f = params.get<Field>("Field");
  const auto& fid = f.get_header().get_identifier();

  // Sanity checks
  using namespace ShortFieldTagsNames;
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (layout.rank()>1 && layout.rank()<=6,
      "Error! Field rank not supported by FieldAtLevel.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + to_string(layout) + "\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + to_string(layout) + "\n");

  // Note: you may ask why we can't just store f and be done, rather than go through the
  // add_field infrastructure. Unfortunately, there are some checks in the base classes
  // that require the diagnostic to have 1+ required fields. So we have to do this.
  // TODO: one day we may make atm diags *not* inherit from atm process...
  add_field<Required>(fid);

  // We can also create the diagnostic already!
  const auto& location = params.get<std::string>("Field Level Location");
  const auto diag_field_name = f.name() + "_at_" + location;

  FieldIdentifier d_fid (diag_field_name,layout.strip_dim(tag),fid.get_units(),fid.get_grid_name());
  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();

  // Figure out the level
  if (ekat::starts_with(location,"lev")) {
    const auto& lev = location.substr(4);
    EKAT_REQUIRE_MSG (lev.find_first_not_of("0123456789")==std::string::npos,
        "Error! Invalid level specification for FieldAtLevel diagnostic.\n"
        "  - input value: '" + location + "'\n"
        "  - expected: 'levN', with N an integer.\n");
    m_field_level = std::stoi(lev);
    EKAT_REQUIRE_MSG (m_field_level<layout.dims().back(),
        "Error! Invalid level specification. Level index out of bounds.\n"
        "  - input level: " + std::to_string(m_field_level) + "\n"
        "  - number of levels: " + std::to_string(layout.dims().back()) + "\n");
  } else if (location=="model_top") {
    m_field_level = 0;
  } else if (location=="model_bot") {
    m_field_level = layout.dims().back()-1;
  } else {
    EKAT_ERROR_MSG (
        "Error! Invalid level specification for FieldAtLevel diagnostic.\n"
        "  - input value: '" + location + "'\n"
        "  - expected: 'model_top','model_bot', or 'levN', with N an integer.\n");
  }
}

// =========================================================================================
void FieldAtLevel::compute_diagnostic_impl()
{
  const auto& f = get_fields_in().front();
  const auto& diag_layout = m_diagnostic_output.get_header().get_identifier().get_layout();
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  RangePolicy policy(0,diag_layout.size());
  auto level  = m_field_level;
  switch (diag_layout.rank()) {
    case 1:
      {
        auto f_view = f.get_view<const Real**>();
        auto d_view = m_diagnostic_output.get_view<      Real*>();
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            d_view(idx) = f_view(idx,level);
        });
      }
      break;
    case 2:
      {
        auto f_view = f.get_view<const Real***>();
        auto d_view = m_diagnostic_output.get_view<      Real**>();
        const int dim1 = diag_layout.dims()[1];
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            const int i = idx / dim1;
            const int j = idx % dim1;
            d_view(i,j) = f_view(i,j,level);
        });
      }
      break;
    case 3:
      {
        auto f_view = f.get_view<const Real****>();
        auto d_view = m_diagnostic_output.get_view<      Real***>();
        const int dim1 = diag_layout.dims()[1];
        const int dim2 = diag_layout.dims()[2];
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            const int i = (idx / dim2) / dim1;
            const int j = (idx / dim2) % dim1;
            const int k =  idx % dim2;
            d_view(i,j,k) = f_view(i,j,k,level);
        });
      }
      break;
    case 4:
      {
        auto f_view = f.get_view<const Real*****>();
        auto d_view = m_diagnostic_output.get_view<      Real****>();
        const int dim1 = diag_layout.dims()[1];
        const int dim2 = diag_layout.dims()[2];
        const int dim3 = diag_layout.dims()[3];
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            const int i = ((idx / dim3) / dim2) / dim1;
            const int j = ((idx / dim3) / dim2) % dim1;
            const int k =  (idx / dim3) % dim2;
            const int l =   idx % dim3;
            d_view(i,j,k,l) = f_view(i,j,k,l,level);
        });
      }
      break;
    case 5:
      {
        auto f_view = f.get_view<const Real******>();
        auto d_view = m_diagnostic_output.get_view<      Real*****>();
        const int dim1 = diag_layout.dims()[1];
        const int dim2 = diag_layout.dims()[2];
        const int dim3 = diag_layout.dims()[3];
        const int dim4 = diag_layout.dims()[4];
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            const int i = (((idx / dim4) / dim3) / dim2) / dim1;
            const int j = (((idx / dim4) / dim3) / dim2) % dim1;
            const int k =  ((idx / dim4) / dim3) % dim2;
            const int l =   (idx / dim4) % dim3;
            const int m =    idx / dim4;
            d_view(i,j,k,l,m) = f_view(i,j,k,l,m,level);
        });
      }
      break;
  }
  Kokkos::fence();
}

} //namespace scream
