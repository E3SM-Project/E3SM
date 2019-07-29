#ifndef SCREAM_REMAP_UTILS_HPP
#define SCREAM_REMAP_UTILS_HPP

#include "share/field/field_layout.hpp"
#include "share/grid/grid_utils.hpp"

namespace scream
{

// The type of the layout, that is, the kind of field it represent.
// Note: the following example of layouts are only indicative. The actual tags
//       can in principle be shuffled. Usually, Dynamics will have Element first,
//       and VerticalLevel last (if present), though this is not really necessary
//       from the point of view of detecting the layout type.
enum class LayoutType {
  Invalid,    // Not supported
              //        Dynamics                                                                        Physics
  Scalar2D,   // <Element,GaussPoint,GaussPoint>                                          <Column>
  Vector2D,   // <Element,Component,GaussPoint,GaussPoint>                                <Column,Component>
  Tensor2D,   // <Element,ComponentX,ComponentY,GaussPoint,GaussPoint>                    <Column,ComponentX,ComponentY>
  Scalar3D,   // <Element,GaussPoint,GaussPoint,VerticalLevel>                            <Column,VerticalLevel>
  Vector3D,   // <Element,Component,GaussPoint,GaussPoint,VerticalLevel>                  <Column,Component,VerticalLevel>
  Tensor3D,   // <Element,ComponentX,ComponentY,GaussPoint,GaussPoint,VerticalLevel>      <Column,ComponentX,ComponentY,VerticalLevel>
};

LayoutType get_layout_type (const FieldLayout& layout);

GridType get_grid_type (const FieldLayout& layout);

} // namespace scream

#endif // SCREAM_REMAP_UTILS_HPP
