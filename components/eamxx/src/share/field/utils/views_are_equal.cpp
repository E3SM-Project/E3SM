#include "share/field/field_utils.hpp"

namespace scream {

namespace impl {

// NOTE: we implement on host, as this is just used in unit tests anyways
template<typename ST>
bool views_are_equal(const Field& f1, const Field& f2, const ekat::Comm* comm)
{
  // Get physical layout (shoudl be the same for both fields)
  const auto& l1 = f1.get_header().get_identifier().get_layout();
  const auto& l2 = f2.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (l1==l2,
      "Error! Input fields have different layouts.\n");

  // For simplicity, we perform the check on Host only. This is not a big
  // limitation, since this code is likely used only in testing.
  f1.sync_to_host();
  f2.sync_to_host();

  // Reshape based on the rank, then loop over all entries.
  // NOTE: because views_are_equal() is only used for testing, we generalize by
  // always calling get_strided_view(), even if it could be contiguous.
  // This allows us to call this function on any type of field,
  // including multi-slice subfields
  bool same_locally = true;
  const auto& dims = l1.dims();
  switch (l1.rank()) {
    case 0:
      {
        auto v1 = f1.template get_strided_view<const ST,Host>();
        auto v2 = f2.template get_strided_view<const ST,Host>();
        if (v1() != v2()) {
          same_locally = false;
          break;
        }
        break;
      }
    case 1:
      {
        auto v1 = f1.template get_strided_view<const ST*,Host>();
        auto v2 = f2.template get_strided_view<const ST*,Host>();
        for (int i=0; i<dims[0]; ++i) {
          if (v1(i) != v2(i)) {
            same_locally = false;
            break;
          }
        }
      }
      break;
    case 2:
      {
        auto v1 = f1.template get_strided_view<const ST**,Host>();
        auto v2 = f2.template get_strided_view<const ST**,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; j<dims[1]; ++j) {
            if (v1(i,j) != v2(i,j)) {
              same_locally = false;
              break;
            }
        }}
      }
      break;
    case 3:
      {
        auto v1 = f1.template get_strided_view<const ST***,Host>();
        auto v2 = f2.template get_strided_view<const ST***,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; same_locally && j<dims[1]; ++j) {
            for (int k=0; k<dims[2]; ++k) {
              if (v1(i,j,k) != v2(i,j,k)) {
                same_locally = false;
                break;
              }
        }}}
      }
      break;
    case 4:
      {
        auto v1 = f1.template get_strided_view<const ST****,Host>();
        auto v2 = f2.template get_strided_view<const ST****,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; same_locally && j<dims[1]; ++j) {
            for (int k=0; same_locally && k<dims[2]; ++k) {
              for (int l=0; l<dims[3]; ++l) {
                if (v1(i,j,k,l) != v2(i,j,k,l)) {
                  same_locally = false;
                  break;
                }
        }}}}
      }
      break;
    case 5:
      {
        auto v1 = f1.template get_strided_view<const ST*****,Host>();
        auto v2 = f2.template get_strided_view<const ST*****,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; same_locally && j<dims[1]; ++j) {
            for (int k=0; same_locally && k<dims[2]; ++k) {
              for (int l=0; same_locally && l<dims[3]; ++l) {
                for (int m=0; m<dims[4]; ++m) {
                  if (v1(i,j,k,l,m) != v2(i,j,k,l,m)) {
                    same_locally = false;
                    break;
                  }
        }}}}}
      }
      break;
    case 6:
      {
        auto v1 = f1.template get_strided_view<const ST******,Host>();
        auto v2 = f2.template get_strided_view<const ST******,Host>();
        for (int i=0; same_locally && i<dims[0]; ++i) {
          for (int j=0; same_locally && j<dims[1]; ++j) {
            for (int k=0; same_locally && k<dims[2]; ++k) {
              for (int l=0; same_locally && l<dims[3]; ++l) {
                for (int m=0; same_locally && m<dims[4]; ++m) {
                  for (int n=0; n<dims[5]; ++n) {
                    if (v1(i,j,k,l,m,n) != v2(i,j,k,l,m,n)) {
                      same_locally = false;
                      break;
                    }
        }}}}}}
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank.\n");
  }

  if (comm) {
    bool same_globally;
    comm->all_reduce(&same_locally,&same_globally,1,MPI_LAND);
    return same_globally;
  } else {
    return same_locally;
  }
}

} // namespace impl

bool views_are_equal(const Field& f1, const Field& f2, const ekat::Comm* comm) {
  EKAT_REQUIRE_MSG (f1.data_type()==f2.data_type(),
      "Error! Views have different data types.\n");

  bool ret = false;
  switch (f1.data_type()) {
    case DataType::IntType:
      ret = impl::views_are_equal<int>(f1,f2,comm); break;
    case DataType::FloatType:
      ret = impl::views_are_equal<float>(f1,f2,comm); break;
    case DataType::DoubleType:
      ret = impl::views_are_equal<double>(f1,f2,comm); break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized field data type.\n");
  }
  return ret;
}

} // namespace scream
