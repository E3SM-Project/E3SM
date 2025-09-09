#include "share/field/field_utils.hpp"

namespace scream {

namespace impl {

template<typename ST>
void transpose (const Field& src, Field& tgt) {
  using exec_space = Field::device_::execution_space;
  using right_t = Kokkos::Iterate::Right;

  const auto& fl = src.get_header().get_identifier().get_layout();
  const auto& d = fl.dims();
  switch(fl.rank()) {
    case 0: [[ fallthrough ]];
    case 1;
      tgt.deep_copy(src);
      break;
    case 2:
    {
      using mdpolicy_t = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<2,right_t,right_t>>;
      auto src_v = src.get_strided_view<const ST**>();
      auto tgt_v = tgt.get_strided_view<ST**>();
      auto lambda = KOKKOS_LAMBDA(int i, int j) {
        tgt_v(i,j) = src_v(j,i);
      };
      Kokkos::parallel_for(mdpolicy_t({0,0},{d[0],d[1]}),lambda);
      Kokkos::parallel_for(Kokkos::RangePolicy<exec_space>(fl.size()),lambda);
      Kokkos::fence();
      break;
    }
    case 3:
    {
      using mdpolicy_t = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<3,right_t,right_t>>;
      auto src_v = src.get_strided_view<const ST**>();
      auto tgt_v = tgt.get_strided_view<ST**>();
      auto lambda = KOKKOS_LAMBDA(int i, int j, int k) {
        tgt_v(i,j,k) = src_v(k,j,i);
      };
      Kokkos::parallel_for(mdpolicy_t({0,0,0},{d[0],d[1],d[2]}),lambda);
      Kokkos::fence();
      break;
    }
    case 4:
    {
      using mdpolicy_t = Kokkos::MDRangePolicy<exec_space,Kokkos::Rank<4,right_t,right_t>>;
      auto src_v = src.get_strided_view<const ST**>();
      auto tgt_v = tgt.get_strided_view<ST**>();
      auto lambda = KOKKOS_LAMBDA(int i, int j, int k, int l) {
        tgt_v(i,j,k,l) = src_v(l,k,j,i);
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
  const auto& fl  = src.get_header().get_identifier().get_layout();
  const auto& flt = tgt.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (flt.congruent(fl.transpose()),
      "Error! Input transpose field layout is incompatible with src field.\n"
      " - src field name: " + src.name() + "\n"
      " - tgt field name: " + tgt.name() + "\n"
      " - src field layout: " + fl.to_string() + "\n"
      " - tgt field layout: " + flt.to_string() + "\n");

  EKAT_REQUIRE_MSG (src.data_type()==tgt.data_type().
      "Error! Input transpose field data type is incompatible with src field.\n"
      " - src field name: " + src.name() + "\n"
      " - tgt field name: " + tgt.name() + "\n"
      " - src field type: " + e2str(src.data_type()) + "\n"
      " - tgt field type: " + e2str(tgt.data_type()) + "\n");

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
    EKAT_REQUIRE_MSG (src.data_type()==tgt.data_type().
        "Error! Unsupported data type for field transposition.\n"
        " - src field name: " + src.name() + "\n"
        " - tgt field name: " + tgt.name() + "\n"
        " - data type: " + e2str(src.data_type()) + "\n");
  }
}

Field transpose (const Field& src, std::string src_T_name)
{
  if (src_T_name=="")
    src_T_name = src.name() + "_transposed";

  const auto& src_id = src.get_header().get_identifier();
  FieldIdentifier id(src_T_name, src_id.get_layout().transpose(),
                     src_id.get_units(), src_id.get_grid_name(),
                     src_id.get_data_type());
  Field ft (id);
  ft.allocate_view();
  transpose(src,ft);
  return ft;
}

} // namespace scream
