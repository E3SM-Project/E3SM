#include "share/field/field_utils.hpp"

#include <ekat_kokkos_types.hpp>

namespace scream {

namespace impl {

template<typename ST>
void transpose (const Field& src, Field& tgt) {
  using exec_space = Field::device_t::execution_space;
  constexpr auto Right = Kokkos::Iterate::Right;

  const auto& fl = src.get_header().get_identifier().get_layout();
  const auto& d = fl.dims();
  switch(fl.rank()) {
    case 0: [[ fallthrough ]];
    case 1:
      tgt.deep_copy(src);
      break;
    case 2:
    {
      using mdpolicy_t = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<2,Right,Right>>;
      auto src_v = src.get_strided_view<const ST**>();
      auto tgt_v = tgt.get_strided_view<ST**>();
      auto lambda = KOKKOS_LAMBDA(int i, int j) {
        tgt_v(j,i) = src_v(i,j);
      };
      Kokkos::parallel_for(mdpolicy_t({0,0},{d[0],d[1]}),lambda);
      Kokkos::fence();
      break;
    }
    case 3:
    {
      using mdpolicy_t = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<3,Right,Right>>;
      auto src_v = src.get_strided_view<const ST***>();
      auto tgt_v = tgt.get_strided_view<ST***>();
      auto lambda = KOKKOS_LAMBDA(int i, int j, int k) {
        tgt_v(k,j,i) = src_v(i,j,k);
      };
      Kokkos::parallel_for(mdpolicy_t({0,0,0},{d[0],d[1],d[2]}),lambda);
      Kokkos::fence();
      break;
    }
    case 4:
    {
      using mdpolicy_t = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<4,Right,Right>>;
      auto src_v = src.get_strided_view<const ST****>();
      auto tgt_v = tgt.get_strided_view<ST****>();
      auto lambda = KOKKOS_LAMBDA(int i, int j, int k, int l) {
        tgt_v(l,k,j,i) = src_v(i,j,k,l);
      };
      Kokkos::parallel_for(mdpolicy_t({0,0,0,0},{d[0],d[1],d[2],d[3]}),lambda);
      Kokkos::fence();
      break;
    }
    default:
      EKAT_ERROR_MSG("Unsupported rank (" + std::to_string(fl.rank()) + ") for field transposition.\n");
  }
}

} // namespace impl

void transpose (const Field& src, Field& tgt)
{
  // Check tgt field has the right layout
  const auto& src_id = src.get_header().get_identifier();
  const auto& tgt_id = tgt.get_header().get_identifier();
  EKAT_REQUIRE_MSG (src_id.get_layout().congruent(tgt_id.get_layout().transpose()),
      "Error! Input transpose field layout is incompatible with src field.\n"
      " - src field name: " + src.name() + "\n"
      " - tgt field name: " + tgt.name() + "\n"
      " - src field layout: " + src_id.get_layout().to_string() + "\n"
      " - tgt field layout: " + tgt_id.get_layout().to_string() + "\n");

  EKAT_REQUIRE_MSG (src_id.get_units()==tgt_id.get_units(),
      "Error! Input transpose field units are incompatible with src field.\n"
      " - src field name: " + src.name() + "\n"
      " - tgt field name: " + tgt.name() + "\n"
      " - src field units: " + src_id.get_units().get_si_string() + "\n"
      " - tgt field units: " + tgt_id.get_units().get_si_string() + "\n");

  EKAT_REQUIRE_MSG (src_id.data_type()==tgt_id.data_type(),
      "Error! Input transpose field data type is incompatible with src field.\n"
      " - src field name: " + src.name() + "\n"
      " - tgt field name: " + tgt.name() + "\n"
      " - src field type: " + e2str(src_id.data_type()) + "\n"
      " - tgt field type: " + e2str(tgt_id.data_type()) + "\n");

  EKAT_REQUIRE_MSG (src.is_allocated(),
      "Error! Input src field is not allocated yet.\n"
      " - src field name: " + src.name() + "\n");
  EKAT_REQUIRE_MSG (tgt.is_allocated(),
      "Error! Input tgt field is not allocated yet.\n"
      " - tgt field name: " + tgt.name() + "\n");

  const auto dt = src.data_type();
  if (dt==DataType::DoubleType) {
    impl::transpose<double>(src,tgt);
  } else if (dt==DataType::FloatType) {
    impl::transpose<float>(src,tgt);
  } else if (dt==DataType::IntType) {
    impl::transpose<int>(src,tgt);
  } else {
    EKAT_ERROR_MSG (
        "Error! Unsupported data type for field transposition.\n"
        " - src field name: " + src.name() + "\n"
        " - tgt field name: " + tgt.name() + "\n"
        " - data type: " + e2str(src_id.data_type()) + "\n");
  }
}

Field transpose (const Field& src)
{
  const auto& src_id = src.get_header().get_identifier();
  FieldIdentifier id(src_id.name()+"_transpose", src_id.get_layout().transpose(),
                     src_id.get_units(), src_id.get_grid_name(),
                     src_id.data_type());
  Field ft (id,true);
  transpose(src,ft);
  return ft;
}

} // namespace scream
