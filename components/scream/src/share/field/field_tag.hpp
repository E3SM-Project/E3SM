#ifndef SCREAM_FIELD_TAG_HPP
#define SCREAM_FIELD_TAG_HPP

#include "share/scream_assert.hpp"

#include <string>

namespace scream
{

/*
 *   An enum to tag fields dimensions
 *
 *   Field tags are meant to be used to determine what each dimension
 *   in a field refers to, and help distinguish fields.
 *   For instance, say there are two classes A and B, storing a field
 *   called 'tracers', but they expect the layout to be different.
 *   Namely, A expects to index the field as (element, dim, gauss point, gauss point),
 *   while B expects to index it as (element, gauss point, gauss point, dim).
 *   A check on the field name and rank is not enough to establish
 *   that the two fields are not the same. If the number of dimensions
 *   is equal to the number of points, even a check on the extents
 *   would not distinguish the two. At this point, using a tag for each
 *   dimension is the only way to distiguish the two.
 */

enum class FieldTag {
  Invalid,
  Element,
  VerticalLevel,
  Column,
  GaussPoint,
  Component,
  ComponentX,
  ComponentY,
  TimeLevel,
  Variable
};

inline std::string tag2string (const FieldTag ft) {
  std::string name = "";
  switch(ft) {
    case FieldTag::Invalid:
      name = "Invalid";
      break;
    case FieldTag::Element:
      name = "Element";
      break;
    case FieldTag::VerticalLevel:
      name = "VerticalLevel";
      break;
    case FieldTag::TimeLevel:
      name = "TimeLevel";
      break;
    case FieldTag::Column:
      name = "Column";
      break;
    case FieldTag::GaussPoint:
      name = "GaussPoint";
      break;
    case FieldTag::Component:
      name = "Component";
      break;
    case FieldTag::ComponentX:
      name = "ComponentX";
      break;
    case FieldTag::ComponentY:
      name = "ComponentY";
      break;
    case FieldTag::Variable:
      name = "Variable";
      break;
    default:
      scream_error_msg("Error! Unrecognized field tag.");
  }

  return name;
}

// If using tags a lot, consider adding 'using namespace ShortFieldTagsNames' locally to your function or cpp file.
namespace ShortFieldTagsNames {

  constexpr auto EL   = FieldTag::Element;
  constexpr auto COL  = FieldTag::Column;
  constexpr auto GP   = FieldTag::GaussPoint;
  constexpr auto TL   = FieldTag::TimeLevel;
  constexpr auto VAR  = FieldTag::Variable;
  constexpr auto VL   = FieldTag::VerticalLevel;
  constexpr auto CMP  = FieldTag::Component;
  constexpr auto CMPX = FieldTag::ComponentX;
  constexpr auto CMPY = FieldTag::ComponentY;
}

} // namespace scream

#endif // SCREAM_FIELD_TAG_HPP
