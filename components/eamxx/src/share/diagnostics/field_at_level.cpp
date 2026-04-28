#include "field_at_level.hpp"

#include "share/util/eamxx_universal_constants.hpp"

#include <ekat_std_utils.hpp>

namespace scream
{

// =========================================================================================
FieldAtLevel::FieldAtLevel (const ekat::Comm& comm, const ekat::ParameterList& params,
                            const std::shared_ptr<const AbstractGrid>& grid)
  : AbstractDiagnostic(comm,params,grid)
{
  m_field_name = m_params.get<std::string>("field_name");
  m_field_in_names.push_back(m_field_name);
}

void FieldAtLevel::
initialize_impl ()
{
  const auto& f = m_fields_in.at(m_field_in_names.front());

  // Sanity checks
  using namespace ShortFieldTagsNames;
  EKAT_REQUIRE_MSG (f.data_type()==DataType::IntType or f.data_type()==DataType::RealType,
      "[FieldAtLevel] Error! Unsupported field data type.\n"
      " - field name: " + f.name() + "\n"
      " - data type : " + e2str(f.data_type()) + "\n");

  const auto& fid    = f.get_header().get_identifier();
  const auto& layout = fid.get_layout();
  EKAT_REQUIRE_MSG (layout.rank()>=2 && layout.rank()<=6,
      "Error! Field rank not supported by FieldAtLevel.\n"
      " - field name: " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n"
      "NOTE: if you requested something like 'field_horiz_avg_at_Y',\n"
      "      you can avoid this error by requesting 'fieldX_at_Y_horiz_avg' instead.\n");
  const auto tag = layout.tags().back();
  EKAT_REQUIRE_MSG (tag==LEV || tag==ILEV,
      "Error! FieldAtLevel diagnostic expects a layout ending with 'LEV'/'ILEV' tag.\n"
      " - field name  : " + fid.name() + "\n"
      " - field layout: " + layout.to_string() + "\n");

  // Figure out the level
  const auto& location = m_params.get<std::string>("vertical_location");
  if (ekat::starts_with(location,"lev_")) {
    const auto& lev = location.substr(4);
    EKAT_REQUIRE_MSG (lev.find_first_not_of("0123456789")==std::string::npos,
        "Error! Invalid level specification for FieldAtLevel diagnostic.\n"
        "  - input value: '" + location + "'\n"
        "  - expected: 'lev_N', with N an integer.\n");
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

  // All good, create the diag output
  auto diag_name = m_field_name + "_at_" + location;
  auto d_fid = fid.clone(diag_name).reset_layout(layout.clone().strip_dim(tag));
  m_diagnostic_output = Field(d_fid,true);
  if (f.has_valid_mask()) {
    m_diagnostic_output.create_valid_mask();
    m_diagnostic_output.get_header().set_may_be_filled(true);
  }

  using stratts_t = std::map<std::string,std::string>;

  // Propagate any io string attribute from input field to diag field
  const auto& src_atts = f.get_header().get_extra_data<stratts_t>("io: string attributes");
        auto& dst_atts = m_diagnostic_output.get_header().get_extra_data<stratts_t>("io: string attributes");
  for (const auto& [name, val] : src_atts) {
    dst_atts[name] = val;
  }
}

// =========================================================================================
void FieldAtLevel::compute_diagnostic_impl()
{
  const auto& f = m_fields_in.at(m_field_name);

  auto ALL = Kokkos::ALL;
  bool masked = m_diagnostic_output.has_valid_mask();
  switch (f.rank()) {
    case 1:
      {
        if (f.data_type()==DataType::IntType) {
          auto fv = f.get_view<const int*>();
          auto dv = m_diagnostic_output.get_view<int>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,m_field_level));
        } else if (f.data_type()==DataType::RealType) {
          auto fv = f.get_view<const Real*>();
          auto dv = m_diagnostic_output.get_view<Real>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,m_field_level));
        }

        if (masked) {
          auto fmv = f.get_valid_mask().get_view<const int*>();
          auto dmv = m_diagnostic_output.get_valid_mask().get_view<int>();
          Kokkos::deep_copy(dmv,Kokkos::subview(fmv,m_field_level));
        }
      } break;
    case 2:
      {
        if (f.data_type()==DataType::IntType) {
          auto fv = f.get_view<const int**>();
          auto dv = m_diagnostic_output.get_view<int*>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,m_field_level));
        } else if (f.data_type()==DataType::RealType) {
          auto fv = f.get_view<const Real**>();
          auto dv = m_diagnostic_output.get_view<Real*>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,m_field_level));
        }

        if (masked) {
          auto fmv = f.get_valid_mask().get_view<const int**>();
          auto dmv = m_diagnostic_output.get_valid_mask().get_view<int*>();
          Kokkos::deep_copy(dmv,Kokkos::subview(fmv,ALL,m_field_level));
        }
      } break;
    case 3:
      {
        if (f.data_type()==DataType::IntType) {
          auto fv = f.get_view<const int***>();
          auto dv = m_diagnostic_output.get_view<int**>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,ALL,m_field_level));
        } else if (f.data_type()==DataType::RealType) {
          auto fv = f.get_view<const Real***>();
          auto dv = m_diagnostic_output.get_view<Real**>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,ALL,m_field_level));
        }

        if (masked) {
          auto fmv = f.get_valid_mask().get_view<const int***>();
          auto dmv = m_diagnostic_output.get_valid_mask().get_view<int**>();
          Kokkos::deep_copy(dmv,Kokkos::subview(fmv,ALL,ALL,m_field_level));
        }
      } break;
    case 4:
      {
        if (f.data_type()==DataType::IntType) {
          auto fv = f.get_view<const int****>();
          auto dv = m_diagnostic_output.get_view<int***>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,ALL,ALL,m_field_level));
        } else if (f.data_type()==DataType::RealType) {
          auto fv = f.get_view<const Real****>();
          auto dv = m_diagnostic_output.get_view<Real***>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,ALL,ALL,m_field_level));
        }

        if (masked) {
          auto fmv = f.get_valid_mask().get_view<const int****>();
          auto dmv = m_diagnostic_output.get_valid_mask().get_view<int***>();
          Kokkos::deep_copy(dmv,Kokkos::subview(fmv,ALL,ALL,ALL,m_field_level));
        }
      } break;
    case 5:
      {
        if (f.data_type()==DataType::IntType) {
          auto fv = f.get_view<const int*****>();
          auto dv = m_diagnostic_output.get_view<int****>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,ALL,ALL,ALL,m_field_level));
        } else if (f.data_type()==DataType::RealType) {
          auto fv = f.get_view<const Real*****>();
          auto dv = m_diagnostic_output.get_view<Real****>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,ALL,ALL,ALL,m_field_level));
        }

        if (masked) {
          auto fmv = f.get_valid_mask().get_view<const int*****>();
          auto dmv = m_diagnostic_output.get_valid_mask().get_view<int****>();
          Kokkos::deep_copy(dmv,Kokkos::subview(fmv,ALL,ALL,ALL,ALL,m_field_level));
        }
      } break;
    case 6:
      {
        if (f.data_type()==DataType::IntType) {
          auto fv = f.get_view<const int******>();
          auto dv = m_diagnostic_output.get_view<int*****>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,ALL,ALL,ALL,ALL,m_field_level));
        } else if (f.data_type()==DataType::RealType) {
          auto fv = f.get_view<const Real******>();
          auto dv = m_diagnostic_output.get_view<Real*****>();
          Kokkos::deep_copy(dv,Kokkos::subview(fv,ALL,ALL,ALL,ALL,ALL,m_field_level));
        }

        if (masked) {
          auto fmv = f.get_valid_mask().get_view<const int******>();
          auto dmv = m_diagnostic_output.get_valid_mask().get_view<int*****>();
          Kokkos::deep_copy(dmv,Kokkos::subview(fmv,ALL,ALL,ALL,ALL,ALL,m_field_level));
        }
      } break;
    default:
      EKAT_ERROR_MSG (
          "[FieldAtLevel] Unexpected field rank. You should have gotten an error before though.\n"
          " - field name: " + f.name() + "\n"
          " - field rank: " + std::to_string(f.rank()) + "\n");
  }

  // TODO: remove when IO stops relying on mask=0 entries being already set to FillValue
  if (masked) {
    const auto fv = f.data_type()==DataType::RealType ? constants::fill_value<Real>
                                                      : constants::fill_value<int>;
    m_diagnostic_output.deep_copy(fv,m_diagnostic_output.get_valid_mask(),true);
  }
}

} //namespace scream
