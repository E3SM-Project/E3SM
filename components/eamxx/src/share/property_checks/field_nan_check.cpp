#include "share/property_checks/field_nan_check.hpp"
#include "share/field/field_utils.hpp"
#include "share/util//eamxx_array_utils.hpp"

#include "ekat/util/ekat_math_utils.hpp"

namespace scream
{

FieldNaNCheck::
FieldNaNCheck (const Field& f,
               const std::shared_ptr<const AbstractGrid>& grid)
 : m_grid (grid)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (f.rank()<=6,
      "Error in FieldNaNCheck constructor: unsupported field rank.\n"
      "  - Field name: " + f.name() << "\n"
      "  - Field rank: " + std::to_string(f.rank()) + "\n");
  EKAT_REQUIRE_MSG (field_valid_data_types().has_v(f.data_type()),
      "Error in FieldNaNCheck constructor: field data type not supported.\n"
      "  - Field name: " + f.name() << "\n"
      "  - Field rank: " + std::to_string(f.rank()) + "\n");

  EKAT_REQUIRE_MSG (grid==nullptr || f.get_header().get_identifier().get_grid_name()==grid->name(),
      "Error! The name of the input grid does not match the grid name stored in the field identifier.\n"
      "  - Field name: " + f.name() + "\n"
      "  - Field grid name: " + f.get_header().get_identifier().get_grid_name() + "\n"
      "  - Input grid name: " + grid->name() + "\n");

  // We can't repair NaN's.
  set_fields ({f},{false});
}

template<typename ST>
PropertyCheck::ResultAndMsg FieldNaNCheck::check_impl() const {
  using const_ST    = typename std::add_const<ST>::type;

  const auto& f = fields().front();

  const auto& layout = f.get_header().get_identifier().get_layout();
  const auto extents = layout.extents();
  const auto size = layout.size();

  int invalid_idx = -1;
  using max_t = Kokkos::Max<int>;
  // below, we can't be sure the field we consider has a continuous allocation,
  // so we use get_strided_view()
  switch (layout.rank()) {
    case 1:
      {
        auto v = f.template get_strided_view<const_ST*>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int i, int& result) {
          if (ekat::is_invalid(v(i))) {
            result = i;
          }
        }, max_t(invalid_idx));
      }
      break;
    case 2:
      {
        auto v = f.template get_strided_view<const_ST**>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j;
          unflatten_idx(idx,extents,i,j);
          if (ekat::is_invalid(v(i,j))) {
            result = idx;
          }
        }, max_t(invalid_idx));
      }
      break;
    case 3:
      {
        auto v = f.template get_strided_view<const_ST***>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j,k;
          unflatten_idx(idx,extents,i,j,k);
          if (ekat::is_invalid(v(i,j,k))) {
            result = idx;
          }
        }, max_t(invalid_idx));
      }
      break;
    case 4:
      {
        auto v = f.template get_strided_view<const_ST****>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j,k,l;
          unflatten_idx(idx,extents,i,j,k,l);
          if (ekat::is_invalid(v(i,j,k,l))) {
            result = idx;
          }
        }, max_t(invalid_idx));
      }
      break;
    case 5:
      {
        auto v = f.template get_strided_view<const_ST*****>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j,k,l,m;
          unflatten_idx(idx,extents,i,j,k,l,m);
          if (ekat::is_invalid(v(i,j,k,l,m))) {
            result = idx;
          }
        }, max_t(invalid_idx));
      }
      break;
    case 6:
      {
        auto v = f.template get_strided_view<const_ST******>();
        Kokkos::parallel_reduce(size, KOKKOS_LAMBDA(int idx, int& result) {
          int i,j,k,l,m,n;
          unflatten_idx(idx,extents,i,j,k,l,m,n);
          if (ekat::is_invalid(v(i,j,k,l,m,n))) {
            result = idx;
          }
        }, max_t(invalid_idx));
      }
      break;
    default:
      EKAT_ERROR_MSG (
          "Internal error in FieldNaNCheck: unsupported field rank.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }

  PropertyCheck::ResultAndMsg res_and_msg;
  res_and_msg.result = invalid_idx<0 ? CheckResult::Pass : CheckResult::Fail;
  res_and_msg.msg = "";
  if (res_and_msg.result==CheckResult::Fail) {
    auto& indices = res_and_msg.fail_loc_indices = unflatten_idx(layout.dims(),invalid_idx);
    res_and_msg.fail_loc_tags = layout.tags();
    res_and_msg.msg  = "FieldNaNCheck failed.\n";
    res_and_msg.msg += "  - field id: " + f.get_header().get_identifier().get_id_string() + "\n";
    using namespace ShortFieldTagsNames;

    int col_lid;

    if (m_grid) {
      // We are storing grid info, and the field is over columns. Get col id and coords.
      col_lid = indices[0];
      auto gids = m_grid->get_dofs_gids().get_strided_view<const AbstractGrid::gid_type*,Host>();

      res_and_msg.msg += "  - indices (w/ global column index): (" + std::to_string(gids(col_lid));
      for (size_t i=1; i<indices.size(); ++i) {
        res_and_msg.msg += "," + std::to_string(indices[i]);
      }
      res_and_msg.msg += ")\n";
      const bool has_latlon = m_grid->has_geometry_data("lat") && m_grid->has_geometry_data("lon");
      if (has_latlon) {
        auto lat = m_grid->get_geometry_data("lat").get_internal_view_data<const Real,Host>();
        auto lon = m_grid->get_geometry_data("lon").get_internal_view_data<const Real,Host>();
        res_and_msg.msg += "  - lat/lon: (" + std::to_string(lat[col_lid]) + ", " + std::to_string(lon[col_lid]) + ")\n";
      }
      bool has_additional_col_info = not additional_data_fields().empty();
      if (has_additional_col_info) {
        std::stringstream msg;
        msg << "  - additional data (w/ local column index):\n";
        for (auto& f : additional_data_fields()) {
          f.sync_to_host();
          msg << "\n";
          print_field_hyperslab(f, {COL}, {col_lid}, msg);
        }
        msg << "\n  END OF ADDITIONAL DATA\n";
        res_and_msg.msg += msg.str();
      }
    }

  } else {
    res_and_msg.msg = "FieldNaNCheck passed.\n";
  }
  return res_and_msg;
}

PropertyCheck::ResultAndMsg FieldNaNCheck::check() const {
  const auto& f = fields().front();

  switch (f.data_type()) {
    case DataType::IntType:
      return check_impl<int>();
    case DataType::FloatType:
      return check_impl<float>();
    case DataType::DoubleType:
      return check_impl<double>();
    default:
      EKAT_ERROR_MSG (
          "Internal error in FieldNaNCheck: unsupported field data type.\n"
          "You should not have reached this line. Please, contact developers.\n");
  }
  return ResultAndMsg{};
}

} // namespace scream
