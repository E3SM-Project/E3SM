#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include "field.hpp"

namespace scream {

template<typename RT1, typename RT2>
bool views_are_equal(const Field<RT1>& f1, const Field<RT2>& f2) {
  static_assert(std::is_same<typename std::remove_cv<RT1>::type,
                             typename std::remove_cv<RT2>::type>::value,
                "Error! Real types must be the same (except possibly for cv qualifiers).\n");
  // Get physical layout (same for both fields)
  const auto& layout = f1.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (f2.get_header().get_identifier().get_layout()==layout,
      "Error! You should only call 'views_are_equal' with two fields with the same physical layout.\n");

  // The alloc props might be different, so we need to get both, and handle views separately
  const auto& p1 = f1.get_header().get_alloc_properties();
  const auto& p2 = f2.get_header().get_alloc_properties();

  // Get physical and allocation sizes
  const int size = layout.size();
  const int last_phys_dim = layout.dims().back();
  const int ext1 = p1.get_last_extent();
  const int ext2 = p2.get_last_extent();

  // Get views
  const auto v1 = f1.get_view();
  const auto v2 = f2.get_view();

  // Simple range policy for a 1d view
  using exec_space = typename decltype(v1)::execution_space;
  Kokkos::RangePolicy<exec_space> policy(0,size);

  int num_diffs = 0;
  Kokkos::parallel_reduce(policy,KOKKOS_LAMBDA(const int i, int& sum){
      // Separate last (physical) dimension from the others (if any)
      const int outer_dims = i / last_phys_dim;
      const int inner_dim  = i % last_phys_dim;
      const int idx1 = outer_dims*ext1 + inner_dim;
      const int idx2 = outer_dims*ext2 + inner_dim;

      if (v1(idx1)!=v2(idx2)) {
        printf("iter,i,j,idx1,idx2:%d,%d,%d,%d,%d\n",i,outer_dims,inner_dim,idx1,idx2);
        printf("v1(idx) = %3.15f\nv2(idx) = %3.15f\n",v1(idx1),v2(idx2));
        ++sum;
      }
    },num_diffs
  );
  return (num_diffs == 0);
}

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
