#ifndef SCREAM_FIELD_TAG_HPP
#define SCREAM_FIELD_TAG_HPP

#include "ekat/ekat_assert.hpp"

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
      name = "EL";
      break;
    case FieldTag::VerticalLevel:
      name = "VL";
      break;
    case FieldTag::TimeLevel:
      name = "TL";
      break;
    case FieldTag::Column:
      name = "COL";
      break;
    case FieldTag::GaussPoint:
      name = "GP";
      break;
    case FieldTag::Component:
      name = "CMP";
      break;
    case FieldTag::ComponentX:
      name = "CMPX";
      break;
    case FieldTag::ComponentY:
      name = "CMPY";
      break;
    case FieldTag::Variable:
      name = "VAR";
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized field tag.");
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
