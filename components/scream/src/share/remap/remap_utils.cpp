#include "share/util/scream_std_utils.hpp"
#include "share/remap/remap_utils.hpp"
#include "share/mpi/scream_comm.hpp"

#include <algorithm>

namespace scream
{

// Shortcuts to avoid using long iterators syntax (like v.begin() and v.end())
template<typename T>
int count (const std::vector<T>& vector, const T& value) {
  return std::count(vector.begin(), vector.end(), value);
}

// =================== Remap utilities ================== //

LayoutType get_layout_type (const FieldLayout& layout) {
  using util::erase;

  auto tags = layout.tags();

  const int n_element = count(tags,FieldTag::Element);
  const int n_column  = count(tags,FieldTag::Column);
  const int ngp = count(tags,FieldTag::GaussPoint);

  // Start from undefined/invalid
  LayoutType result = LayoutType::Invalid;

  if (n_element>0 && ngp==2 && n_column==0) {
    // A Dynamics layout

    // Remove the Element and the two GaussPoint tags
    erase(tags,FieldTag::Element);
    erase(tags,FieldTag::GaussPoint);
    erase(tags,FieldTag::GaussPoint);
  } else if (n_element==0 && ngp==0 && n_column>0) {
    // A Physics layout

    // Remove the column tag
    erase(tags,FieldTag::Column);
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
      if (tags[0]==FieldTag::Component || tags[0]==FieldTag::TimeLevel || tags[0]==FieldTag::Variable) {
        result = LayoutType::Vector2D;
      } else if (tags[0]==FieldTag::VerticalLevel) {
        result = LayoutType::Scalar3D;
      }
      break;
    case 2:
      // Possible scenarios:
      //  1) <Component|TimeLevel|Variable,VerticalLevel>
      //  2) <ComponentX,ComponentY>
      //  3) <Component,TimeLevel|Variable>
      //  4) <TimeLevel|Variable,TimeLevel|Variable>
      if (erase(tags,FieldTag::VerticalLevel)) {
        if (tags[0]==FieldTag::Component || tags[0]==FieldTag::TimeLevel || tags[0]==FieldTag::Variable) {
          result = LayoutType::Vector3D;
        }
      } else if (erase(tags,FieldTag::ComponentX)) {
        if (tags[0]==FieldTag::ComponentY) {
          result = LayoutType::Tensor2D;
        }
      } else if (erase(tags,FieldTag::Component)) {
        if (tags[0]==FieldTag::Variable || tags[0]==FieldTag::TimeLevel) {
          result = LayoutType::Tensor2D;
        }
      } else if (erase(tags,FieldTag::Variable)) {
        if (tags[0]==FieldTag::Variable || tags[0]==FieldTag::TimeLevel) {
          result = LayoutType::Tensor2D;
        }
      } else if (erase(tags,FieldTag::TimeLevel)) {
        if (tags[0]==FieldTag::TimeLevel) {
          result = LayoutType::Tensor2D;
        }
      }
      break;
    case 3:
      if (erase(tags,FieldTag::VerticalLevel)) {
        if (erase(tags,FieldTag::ComponentX)) {
          if (tags[0]==FieldTag::ComponentY) {
            result = LayoutType::Tensor3D;
          }
        } else if (erase(tags,FieldTag::Component)) {
          if (tags[0]==FieldTag::Variable || tags[0]==FieldTag::TimeLevel) {
            result = LayoutType::Tensor3D;
          }
        } else if (erase(tags,FieldTag::Variable)) {
          if (tags[0]==FieldTag::Variable || tags[0]==FieldTag::TimeLevel) {
            result = LayoutType::Tensor3D;
          }
        } else if (erase(tags,FieldTag::TimeLevel)) {
          if (tags[0]==FieldTag::TimeLevel) {
            result = LayoutType::Tensor3D;
          }
        }
      }
  }
  
  return result;
}

} // namespace scream
