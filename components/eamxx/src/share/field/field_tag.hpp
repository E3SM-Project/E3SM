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
  LevelMidPoint,
  LevelInterface,
  Column,
  GaussPoint,
  Component,
  TimeLevel,
};

// If using tags a lot, consider adding 'using namespace ShortFieldTagsNames'
// locally to your function or cpp file.
// TODO: if/when we require std=c++20, this can be removed, and user can do
//   using enum FieldTag;
namespace ShortFieldTagsNames {

  constexpr auto INV  = FieldTag::Invalid;
  constexpr auto EL   = FieldTag::Element;
  constexpr auto COL  = FieldTag::Column;
  constexpr auto GP   = FieldTag::GaussPoint;
  constexpr auto TL   = FieldTag::TimeLevel;
  constexpr auto LEV  = FieldTag::LevelMidPoint;
  constexpr auto ILEV = FieldTag::LevelInterface;
  constexpr auto CMP  = FieldTag::Component;
}

inline std::string e2str (const FieldTag ft) {
  using namespace ShortFieldTagsNames;
  std::string name = "";
  switch(ft) {
    case FieldTag::Invalid:
      name = "Invalid";
      break;
    case EL:
      name = "elem";
      break;
    case LEV:
      name = "lev";
      break;
    case ILEV:
      name = "ilev";
      break;
    case FieldTag::TimeLevel:
      name = "tl";
      break;
    case FieldTag::Column:
      name = "ncol";
      break;
    case FieldTag::GaussPoint:
      name = "gp";
      break;
    case FieldTag::Component:
      name = "dim";
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized field tag.");
  }
  return name;
}

inline std::vector<std::string> tags2str (const std::vector<FieldTag>& tags) {
  std::vector<std::string> names;
  names.reserve(tags.size());
  for (auto t : tags) {
    names.push_back(e2str(t));
  }
  return names;
}

// Allow to stream FieldTag values as strings.
inline std::ostream& operator<< (std::ostream& out, const FieldTag t) {
  out << e2str(t);
  return out;
}

} // namespace scream

#endif // SCREAM_FIELD_TAG_HPP
