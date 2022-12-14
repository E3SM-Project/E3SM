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
  Component1,
  Component2,
  Component3,
  TimeLevel,
  // Added for RRTMGP, TODO: Revisit this approach, is there a better way than adding more field tags?
  Gases,
  ShortWaveBand,
  ShortWaveGpoint,
  LongWaveBand,
  LongWaveGpoint
};

inline std::string e2str (const FieldTag ft) {
  std::string name = "";
  switch(ft) {
    case FieldTag::Invalid:
      name = "Invalid";
      break;
    case FieldTag::Element:
      name = "EL";
      break;
    case FieldTag::LevelMidPoint:
      name = "LEV";
      break;
    case FieldTag::LevelInterface:
      name = "ILEV";
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
    case FieldTag::Component1:
      name = "CMP1";
      break;
    case FieldTag::Component2:
      name = "CMP2";
      break;
    case FieldTag::Component3:
      name = "CMP3";
      break;
    // Added for rrtmgp - see TODO item above
    case FieldTag::Gases:
      name = "NGAS";
      break;
    case FieldTag::ShortWaveBand:
      name = "SWBND";
      break;
    case FieldTag::ShortWaveGpoint:
      name = "SWGPT";
      break;
    case FieldTag::LongWaveBand:
      name = "LWBND";
      break;
    case FieldTag::LongWaveGpoint:
      name = "LWGPT";
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized field tag.");
  }

  return name;
}

// If using tags a lot, consider adding 'using namespace ShortFieldTagsNames'
// locally to your function or cpp file.
// TODO: if/when we require std=c++20, this can be removed, and user can do
//   using enum FieldTag;
namespace ShortFieldTagsNames {

  constexpr auto EL   = FieldTag::Element;
  constexpr auto COL  = FieldTag::Column;
  constexpr auto GP   = FieldTag::GaussPoint;
  constexpr auto TL   = FieldTag::TimeLevel;
  constexpr auto LEV  = FieldTag::LevelMidPoint;
  constexpr auto ILEV = FieldTag::LevelInterface;
  constexpr auto CMP  = FieldTag::Component;
  constexpr auto CMP1 = FieldTag::Component1;
  constexpr auto CMP2 = FieldTag::Component2;
  constexpr auto CMP3 = FieldTag::Component3;
  // Added for rrtmgp - see TODO item above
  constexpr auto NGAS = FieldTag::Gases;
  constexpr auto SWBND = FieldTag::ShortWaveBand;
  constexpr auto LWBND = FieldTag::LongWaveBand;
  constexpr auto SWGPT = FieldTag::ShortWaveGpoint;
  constexpr auto LWGPT = FieldTag::LongWaveGpoint;
}

// Allow to stream FieldTag values as strings.
inline std::ostream& operator<< (std::ostream& out, const FieldTag t) {
  out << e2str(t);
  return out;
}

} // namespace scream

#endif // SCREAM_FIELD_TAG_HPP
