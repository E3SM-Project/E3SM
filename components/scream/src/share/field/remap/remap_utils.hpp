#ifndef SCREAM_REMAP_UTILS_HPP
#define SCREAM_REMAP_UTILS_HPP

#include "share/field/field_layout.hpp"

namespace scream
{

class Comm;

// The type of the layout, that is, the kind of field it represent.
// Note: the following example of layouts are only indicative. The actual tags
//       can in principle be shuffled. Usually, Dynamics will have Element first,
//       and VerticalLevel last (if present), though this is not really necessary
//       from the point of view of detecting the layout type.
enum class LayoutType {
  Invalid,    // Not supported
              //        Dynamics                                            Physics
  Scalar2D,   // <Element, GaussPoint, GaussPoint>                    <Column>
  Vector2D,   // <Element, Component, GaussPoint, GaussPoint>         <Column,Component>
  Scalar3D,   // <Element, GaussPoint, GaussPoint, VerticalLevel>                  <Column,VerticalLevel>
  Tensor2D,   // <Element, ComponentX, ComponentY, GaussPoint, GaussPoint>         <Column,ComponentX,ComponentY>
  Vector3D,   // [num_elems, num_comp, np, np, num_lev].              [num_comp, num_cols, num_lev] 
  Tensor3D,   // [num_elems, num_comp1, num_comp2, np, np, num_lev].  [num_comp1, num_comp2, num_cols, num_lev] 
};

// An enum to mark fields' layouts as physics or dynamics
// Note: to mark a layout as Dyn, it must contain the tag 'Element' (exactly once),
//       while not containing the tag 'Column'.
//       On the other hand, for a layout to be marked as 'Phys' it must contain
//       the tag 'Column' (exactly once), while not containing the tag 'Element' or
//       the tag 'GaussPoint'
enum class PhysDyn {
  Undefined,
  Phys,
  Dyn
};

std::pair<PhysDyn,LayoutType> get_layout_specs (const FieldLayout& layout);

bool is_valid_layout (const FieldLayout& layout);
bool compatible_layout_types (const LayoutType& src, const LayoutType& tgt);
bool is_tags_permutation (const FieldLayout& src, const FieldLayout& tgt);
bool is_phys_dyn_remap (const FieldLayout& src, const FieldLayout& tgt);
bool is_mpi_remap (const FieldLayout& src, const FieldLayout& tgt, const Comm& comm);

} // namespace scream

#endif // SCREAM_REMAP_UTILS_HPP
