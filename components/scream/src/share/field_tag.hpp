#ifndef SCREAM_FIELD_TAG_HPP
#define SCREAM_FIELD_TAG_HPP

#include "error_defs.hpp"
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
 *   
 *   This can also be used (by the driver, or by the atm processes) to
 *   determine whether an input/output field has the same layout in the
 *   field repository as in the atm process, and therefore whether it
 *   needs 'remapping' or not (here 'remapping' is in the programming
 *   sense, not in the mathematical (e.g. eulerian-lagrangian) sense).
 */

enum class FieldTag {
  Element,
  Level,
  GaussPoint,
  Component,
  Tracer,
  State
};

inline std::string tag2string (const FieldTag ft) {
  std::string name = "";
  switch(ft) {
    case FieldTag::Element:
      name = "Element";
      break;
    case FieldTag::Level:
      name = "Level";
      break;
    case FieldTag::GaussPoint:
      name = "GaussPoint";
      break;
    case FieldTag::Component:
      name = "Component";
      break;
    case FieldTag::Tracer:
      name = "Tracer";
      break;
    case FieldTag::State:
      name = "State";
      break;
    default:
      error::runtime_abort("Error! Unrecognized fielt tag.",-1);
  }

  return name;
}

} // namespace scream

#endif // SCREAM_FIELD_TAG_HPP
