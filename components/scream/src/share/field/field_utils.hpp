#ifndef SCREAM_FIELD_UTILS_HPP
#define SCREAM_FIELD_UTILS_HPP

#include "field_tag.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"

namespace scream {

// How a field has to be init-ed
// Note: internally, Value still creates a FieldInitializer object, but can be done
//       behind the scenes by the infrastructure
enum class InitType {
  // NotNeeded,    // No initialization is needed for this field
  Value,        // Should be inited to a specific value
  Initializer,  // A FieldInitializer object should take care of this
  None          // No initialization is needed/expected
};

inline std::string e2str (const InitType e) {
  std::string s;
  switch (e) {
    // case InitType::NotNeeded:   s = "NotNeeded";    break;
    case InitType::None:        s = "None";         break;
    case InitType::Value:       s = "Value";        break;
    case InitType::Initializer: s = "Initializer";  break;
    default: s = "INVALID";
  }

  return s;
}

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

  auto tags = field_tags;

  constexpr auto EL   = FieldTag::Element;
  constexpr auto COL  = FieldTag::Column;
  constexpr auto GP   = FieldTag::GaussPoint;
  constexpr auto TL   = FieldTag::TimeLevel;
  constexpr auto CMP  = FieldTag::Component;
  constexpr auto CMPX = FieldTag::ComponentX;
  constexpr auto CMPY = FieldTag::ComponentY;
  constexpr auto VAR  = FieldTag::Variable;
  constexpr auto VL   = FieldTag::VerticalLevel;

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
      // The only tag left should be 'Component', 'TimeLevel', 'Variable', or 'VerticalLevel
      if (tags[0]==CMP || tags[0]==TL || tags[0]==VAR) {
        result = LayoutType::Vector2D;
      } else if (tags[0]==VL) {
        result = LayoutType::Scalar3D;
      }
      break;
    case 2:
      // Two possible scenarios:
      //  1) <Component|Variable|TimeLevel,VerticalLevel>
      //  2) <Variable|TimeLevel,Component>
      //  3) <TimeLevel,Variable>
      //  4) <ComponentX,ComponentY>
      if (tags[1]==VL && (tags[0]==CMP || tags[0]==VAR || tags[0]==TL)) {
        result = LayoutType::Vector3D;
      } else if ((tags[0]==CMPX && tags[1]==CMPY) ||
                 (tags[1]==CMP && (tags[0]==VAR || tags[0]==TL)) ||
                 (tags[0]==TL && tags[1]==VAR)) {
        result = LayoutType::Tensor2D;
      }
      break;
    case 3:
      // The only scenarios are:
      //  1) <ComponentX, ComponentY, VerticalLevel>
      //  2) <TimeLevel,  Variable, VerticalLevel>
      //  3) <TimeLevel|Variable,  Component, VerticalLevel>
      if (tags[2]==VL) {
        if ((tags[0]==CMPX && tags[1]==CMPY) ||
            (tags[0]==TL && (tags[1]==VAR || tags[1]==CMP)) ||
            (tags[0]==VAR && tags[1]==CMP)) {
          result = LayoutType::Tensor3D;
        }
      }
  }
  
  return result;
}

} // namespace scream

#endif // SCREAM_FIELD_UTILS_HPP
