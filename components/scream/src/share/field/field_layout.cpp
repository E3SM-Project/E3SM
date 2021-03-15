#include "field_layout.hpp"

namespace scream
{

FieldLayout::FieldLayout (const std::initializer_list<FieldTag>& tags)
 : m_rank(tags.size())
 , m_tags(tags)
{
  m_dims.resize(m_rank,-1);
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags)
 : m_rank(tags.size())
 , m_tags(tags)
{
  m_dims.resize(m_rank,-1);
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags,
                          const std::vector<int>& dims)
 : m_rank(tags.size())
 , m_tags(tags)
{
  m_dims.resize(m_rank,-1);
  set_dimensions(dims);
}

void FieldLayout::set_dimension (const int idim, const int dimension) {
  EKAT_REQUIRE_MSG(idim>=0 && idim<m_rank, "Error! Index out of bounds.");
  EKAT_REQUIRE_MSG(dimension>0, "Error! Dimensions must be positive.");
  EKAT_REQUIRE_MSG(m_dims[idim] == -1, "Error! You cannot reset field dimensions once set.\n");
  m_dims[idim] = dimension;
}

void FieldLayout::set_dimensions (const std::vector<int>& dims) {
  // Check, then set dims
  EKAT_REQUIRE_MSG(dims.size()==static_cast<size_t>(m_rank),
                     "Error! Input dimensions vector not properly sized.");
  for (int idim=0; idim<m_rank; ++idim) {
    set_dimension(idim,dims[idim]);
  }
}

LayoutType get_layout_type (const std::vector<FieldTag>& field_tags) {
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


} // namespace scream
