#include "diagnostics/field_at_level.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

namespace scream
{

// =========================================================================================
FieldAtLevel::FieldAtLevel (const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereDiagnostic(comm,params)
  , m_field_layout(m_params.get<FieldLayout>("Field Layout"))
  , m_field_units(m_params.get<ekat::units::Units>("Field Units"))
{
  m_field_name  = m_params.get<std::string>("Field Name");

  EKAT_REQUIRE_MSG (m_field_layout.rank()>1 && m_field_layout.rank()<=6,
      "Error! Field rank not supported by FieldAtLevel.\n"
      " - field name: " + m_field_name + "\n"
      " - field layout: " + to_string(m_field_layout) + "\n");
  using namespace ShortFieldTagsNames;
  EKAT_REQUIRE_MSG (ekat::contains(std::vector<FieldTag>{LEV,ILEV},m_field_layout.tags().back()),
      "Error! FieldAtLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + m_field_name + "\n"
      " - field layout: " + to_string(m_field_layout) + "\n");
}

// =========================================================================================
void FieldAtLevel::set_grids(const std::shared_ptr<const GridsManager> /* grids_manager */)
{
  const auto& gname  = m_params.get<std::string>("Grid Name");

  add_field<Required>(m_field_name, m_field_layout, m_field_units, gname);

  // Calculate the actual level to slice field
  const auto& lev_str = m_params.get<std::string>("Field Level");
  if (lev_str=="tom") {
    m_field_level = 0;
  } else if (lev_str=="bot") {
    m_field_level = m_field_layout.dims().back()-1;
  } else {
    auto is_int = [](const std::string& s) -> bool {
      return s.find_first_not_of("0123456789")==std::string::npos;
    };
    EKAT_REQUIRE_MSG (is_int(lev_str),
        "Error! Entry 'Field Level' must be 'tom', 'bot', or a string representation of an integer.\n");

    m_field_level  = std::stoi(lev_str);
  }
  EKAT_REQUIRE_MSG (m_field_level>=0 && m_field_level<m_field_layout.dims().back(),
      "Error! Invalid value for 'Field Level' in FieldAtLevel diagnostic.\n"
      "  - field name  : " + m_field_name + "\n"
      "  - field layout: " + to_string(m_field_layout) + "\n"
      "  - field level : " + std::to_string(m_field_level) + "\n");

}
// =========================================================================================
void FieldAtLevel::compute_diagnostic_impl()
{
  const auto& f = get_field_in(m_field_name);
  const int dsize = m_field_layout.size () / m_field_layout.dims().back();
  using RangePolicy = Kokkos::RangePolicy<Field::device_t::execution_space>;
  RangePolicy policy(0,dsize);
  auto level  = m_field_level;
  switch (m_field_layout.rank()) {
    case 2:
      {
        auto f_view = f.get_view<const Real**>();
        auto d_view = m_diagnostic_output.get_view<      Real*>();
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            d_view(idx) = f_view(idx,level);
        });
      }
      break;
    case 3:
      {
        auto f_view = f.get_view<const Real***>();
        auto d_view = m_diagnostic_output.get_view<      Real**>();
        const int dim1 = m_field_layout.dims()[1];
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            const int i = idx / dim1;
            const int j = idx % dim1;
            d_view(i,j) = f_view(i,j,level);
        });
      }
      break;
    case 4:
      {
        auto f_view = f.get_view<const Real****>();
        auto d_view = m_diagnostic_output.get_view<      Real***>();
        const int dim1 = m_field_layout.dims()[1];
        const int dim2 = m_field_layout.dims()[2];
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            const int i = (idx / dim2) / dim1;
            const int j = (idx / dim2) % dim1;
            const int k =  idx % dim2;
            d_view(i,j,k) = f_view(i,j,k,level);
        });
      }
      break;
    case 5:
      {
        auto f_view = f.get_view<const Real*****>();
        auto d_view = m_diagnostic_output.get_view<      Real****>();
        const int dim1 = m_field_layout.dims()[1];
        const int dim2 = m_field_layout.dims()[2];
        const int dim3 = m_field_layout.dims()[3];
        Kokkos::parallel_for(m_diagnostic_output.name(),policy,KOKKOS_LAMBDA(const int idx) {
            const int i = ((idx / dim3) / dim2) / dim1;
            const int j = ((idx / dim3) / dim2) % dim1;
            const int k =  (idx / dim3) % dim2;
            const int l =   idx % dim3;
            d_view(i,j,k,l) = f_view(i,j,k,l,level);
        });
      }
      break;
    case 6:
      {
        auto f_view = f.get_view<const Real******>();
        auto d_view = m_diagnostic_output.get_view<      Real*****>();
        const int dim1 = m_field_layout.dims()[1];
        const int dim2 = m_field_layout.dims()[2];
        const int dim3 = m_field_layout.dims()[3];
        const int dim4 = m_field_layout.dims()[4];
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

void FieldAtLevel::set_required_field_impl (const Field& f)
{
  const auto& f_fid = f.get_header().get_identifier();

  // Now that we have the exact units of f, we can build the diagnostic field
  FieldIdentifier d_fid (f_fid.name() + "@lev_" + std::to_string(m_field_level),
                         m_field_layout.strip_dim(m_field_layout.rank()-1),
                         f_fid.get_units(),
                         f_fid.get_grid_name());

  m_diagnostic_output = Field(d_fid);
  m_diagnostic_output.allocate_view();
}

} //namespace scream
