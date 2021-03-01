#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include "field_tag.hpp"
#include "field.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

namespace scream {

// The type of the layout, that is, the kind of field it represent.
enum class LayoutType {
  Invalid,
  Scalar2D,
  Vector2D,
  Tensor2D,
  Scalar3D,
  Vector3D,
  Tensor3D
};

inline LayoutType get_layout_type (const std::vector<FieldTag>& field_tags) {
  using ekat::erase;
  using ekat::count;
  using namespace ShortFieldTagsNames;

  auto tags = field_tags;

  const int n_element = count(tags,EL);
  const int n_column  = count(tags,COL);
  const int ngp       = count(tags,GP);

  // Start from undefined/invalid
  LayoutType result = LayoutType::Invalid;

  if (n_element>0 && ngp==2 && n_column==0) {
    // A Dynamics layout

    // Remove the Element and the two GaussPoint tags
    erase(tags,EL);
    erase(tags,GP);
    erase(tags,GP);
  } else if (n_element==0 && ngp==0 && n_column>0) {
    // A Physics layout

    // Remove the column tag
    erase(tags,COL);
  } else {
    // Not a supported layout.
    return result;
  }

  // Get the size of what's left
  const auto size = tags.size();
  switch (size) {
    case 0:
      result = LayoutType::Scalar2D;
      break;
    case 1:
      // The only tag left should be 'CMP', 'TL', 'VAR', or 'LEV'/'ILEV'
      if (tags[0]==CMP || tags[0]==TL || tags[0]==VAR) {
        result = LayoutType::Vector2D;
      } else if (tags[0]==LEV || tags[0]==ILEV) {
        result = LayoutType::Scalar3D;
      }
      break;
    case 2:
      // Possible scenarios:
      //  1) <CMP|VAR|TL,LEV|ILEV>
      //  2) <VAR|TL,CMP>
      //  3) <TL,VAR>
      //  4) <CMP1,CMP2>
      if ( (tags[1]==LEV || tags[1]==ILEV) && (tags[0]==CMP || tags[0]==VAR || tags[0]==TL)) {
        result = LayoutType::Vector3D;
      } else if ((tags[0]==CMP1 && tags[1]==CMP2) ||
                 (tags[1]==CMP && (tags[0]==VAR || tags[0]==TL)) ||
                 (tags[0]==TL && tags[1]==VAR)) {
        result = LayoutType::Tensor2D;
      }
      break;
    case 3:
      // The only scenarios are:
      //  1) <CMP1, CMP2, LEV|ILEV>
      //  2) <TL,  VAR, LEV|ILEV>
      //  3) <TL|VAR,  CMP, LEV|ILEV>
      if (tags[2]==LEV || tags[2]==ILEV) {
        if ((tags[0]==CMP1 && tags[1]==CMP2) ||
            (tags[0]==TL && (tags[1]==VAR || tags[1]==CMP)) ||
            (tags[0]==VAR && tags[1]==CMP)) {
          result = LayoutType::Tensor3D;
        }
      }
  }
  
  return result;
}

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
