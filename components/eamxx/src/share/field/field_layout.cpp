#include "field_layout.hpp"

#include <ekat/util/ekat_string_utils.hpp>

namespace scream
{

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags,
                          const std::vector<int>& dims)
 : m_rank(tags.size())
 , m_tags(tags)
{
  m_dims.resize(m_rank,-1);
  m_extents = decltype(m_extents)("",m_rank);
  for (int idim=0; idim<m_rank; ++idim) {
    set_dimension(idim,dims[idim]);
  }
}

bool FieldLayout::is_vector_layout () const {
  const auto lt = get_layout_type (m_tags);
  return lt==LayoutType::Vector2D || lt==LayoutType::Vector3D;
}

bool FieldLayout::is_tensor_layout () const {
  const auto lt = get_layout_type (m_tags);
  return lt==LayoutType::Tensor2D || lt==LayoutType::Tensor3D;
}

// get the index of the CMP (Components) tag in the FieldLayout
// e.g., for FieldLayout f({COL, CMP, LEV}, {...});
// we have get_component_idx(f) = 1
int FieldLayout::get_vector_component_idx () const {
  EKAT_REQUIRE_MSG (is_vector_layout(),
      "Error! 'get_vector_dim' available only for vector layouts.\n"
      "       Current layout: " + e2str(get_layout_type(m_tags)) + "\n");

  using namespace ShortFieldTagsNames;
  std::vector<FieldTag> vec_tags = {CMP,NGAS,SWBND,LWBND,SWGPT,ISCCPTAU,ISCCPPRS};
  auto it = std::find_first_of (m_tags.cbegin(),m_tags.cend(),vec_tags.cbegin(),vec_tags.cend());

  EKAT_REQUIRE_MSG (it!=m_tags.cend(),
    "Error! Could not find a vector tag in the layout.\n"
    " - layout: " + to_string(*this) + "\n");

  return std::distance(m_tags.cbegin(),it);
}

// get the extent of the CMP (Components) tag in the FieldLayout
// e.g., for FieldLayout f({COL, CMP, LEV}, {ncol, ncmp, nlev});
// we have get_component_idx(f) = ncmp
int FieldLayout::get_vector_dim () const {
  // since we immediately call get_vector_component_idx(), the error checking
  // there should be sufficient
  return dim(get_vector_component_idx());
}

FieldTag FieldLayout::get_vector_tag () const {
  return m_tags[get_vector_component_idx()];
}

std::vector<int> FieldLayout::get_tensor_dims () const {
  EKAT_REQUIRE_MSG (is_tensor_layout(),
      "Error! 'get_tensor_dims' available only for tensor layouts.\n"
      "       Current layout: " + to_string(*this) + "\n"
      "       Layout type   : " + e2str(get_layout_type(m_tags)) + "\n");

  using namespace ShortFieldTagsNames;
  std::vector<FieldTag> cmp_tags = {CMP,NGAS,SWBND,LWBND,SWGPT,ISCCPTAU,ISCCPPRS};

  std::vector<int> idx;
  auto it = m_tags.begin();
  do {
    it = std::find_first_of (it,m_tags.cend(),cmp_tags.cbegin(),cmp_tags.cend());
    if (it!=m_tags.end()) {
      idx.push_back(std::distance(m_tags.begin(),it));
      ++it;
    }
  } while (it!=m_tags.end());

  EKAT_REQUIRE_MSG (idx.size()==2,
    "Error! Could not find a two tensor tags in the layout.\n"
    " - layout: " + to_string(*this) + "\n"
    " - detected tags indices: " + ekat::join(idx,",") + "\n");

  return idx;
}

std::vector<FieldTag> FieldLayout::get_tensor_tags () const {
  auto idx = get_tensor_dims();
  return {m_tags[idx[0]], m_tags[idx[1]]};
}

FieldLayout FieldLayout::strip_dim (const FieldTag tag) const {
  auto it = ekat::find(m_tags,tag);

  // Check if found
  EKAT_REQUIRE_MSG(it!=m_tags.end(), "Error! Tag '" + e2str(tag) + "' not found.\n");

  // Check only one tag (no ambiguity)
  EKAT_REQUIRE_MSG(ekat::count(m_tags,tag)==1,
                     "Error! Tag '" + e2str(tag) + "' appears multiple times.\n"
                     "       You must inspect tags() and dims() manually.\n");

  return strip_dim (std::distance(m_tags.begin(),it));
}

FieldLayout FieldLayout::strip_dim (const int idim) const {
  EKAT_REQUIRE_MSG (idim>=0 and idim<m_rank,
      "Error! Cannot strip dimension, because it is out of bounds.\n"
      "  - input dim index: " + std::to_string(idim) + "\n"
      "  - layout rank    : " + std::to_string(m_rank) + "\n");
  std::vector<FieldTag> t = tags();
  std::vector<int>      d = dims();
  t.erase(t.begin()+idim);
  d.erase(d.begin()+idim);
  return FieldLayout (t,d);
}

FieldLayout FieldLayout::clone_with_different_extent (const int idim, const int extent) const
{
  FieldLayout copy(m_tags,m_dims);
  copy.set_dimension(idim,extent);

  return copy;
}

void FieldLayout::set_dimension (const int idim, const int dimension) {
  EKAT_REQUIRE_MSG(idim>=0 && idim<m_rank, "Error! Index out of bounds.");
  EKAT_REQUIRE_MSG(dimension>=0, "Error! Dimensions must be non-negative.");
  m_dims[idim] = dimension;

  // Recompute device extents
  auto extents_h = Kokkos::create_mirror_view(m_extents);
  std::copy_n(m_dims.begin(),m_rank,extents_h.data());
  Kokkos::deep_copy(m_extents,extents_h);
}

LayoutType get_layout_type (const std::vector<FieldTag>& field_tags) {
  using ekat::erase;
  using ekat::count;
  using namespace ShortFieldTagsNames;

  auto tags = field_tags;

  const int n_element = count(tags,EL);
  const int n_column  = count(tags,COL);
  const int ngp       = count(tags,GP);
  const int nvlevs    = count(tags,LEV) + count(tags,ILEV);
  const int ncomps    = count(tags,CMP);

  // Start from undefined/invalid
  LayoutType result = LayoutType::Invalid;

  // We don't care about TimeLevel
  erase (tags,TL);

  if (n_element==1 && ngp==2 && n_column==0) {
    // A Dynamics layout

    // Remove the Element and the two GaussPoint tags
    erase(tags,EL);
    erase(tags,GP);
    erase(tags,GP);
  } else if (n_element==0 && ngp==0 && n_column==1) {
    // A Physics layout

    // Remove the column tag
    erase(tags,COL);
  } else if (tags.size()==0) {
    return LayoutType::Scalar0D;
  } else if (tags.size()==1 and tags[0]==CMP) {
    return LayoutType::Vector0D;
  } else if (tags.size()==1 and nvlevs==1) {
    return LayoutType::Scalar1D;
  } else if (tags.size()==2 and ncomps==1 and nvlevs==1) {
    return LayoutType::Vector1D;
  } else {
    // Not a supported layout.
    return result;
  }

  // Get the size of what's left
  const auto size = tags.size();
  auto is_lev_tag = [](const FieldTag t) {
    std::vector<FieldTag> lev_tags = {LEV,ILEV};
    return ekat::contains(lev_tags,t);
  };
  auto is_cmp_tag = [](const FieldTag t) {
    std::vector<FieldTag> cmp_tags = {CMP,NGAS,SWBND,LWBND,SWGPT,ISCCPTAU,ISCCPPRS};
    return ekat::contains(cmp_tags,t);
  };
  switch (size) {
    case 0:
      result = LayoutType::Scalar2D;
      break;
    case 1:
      // The only tag left should be a cmp tag or a lev tag
      if (is_cmp_tag(tags[0])) {
        result = LayoutType::Vector2D;
      } else if (is_lev_tag(tags[0])) {
        result = LayoutType::Scalar3D;
      }
      break;
    case 2:
      // Possible supported scenarios:
      //  1) <CMP,LEV|ILEV>
      //  3) <CMP1,CMP2>
      // where CMP,CMP1,CMP2 are any tag in cmp_tags
      if ( is_cmp_tag(tags[0]) and is_lev_tag(tags[1]) ) {
        result = LayoutType::Vector3D;
      } else if (is_cmp_tag(tags[0]) and is_cmp_tag(tags[1])) {
        result = LayoutType::Tensor2D;
      }
      break;
    case 3:
      // The only supported scenario is:
      //  1) <CMP1, CMP2, LEV|ILEV>
      // where CMP1,CMP2 are any tag in cmp_tags
      if (is_cmp_tag(tags[0]) and is_cmp_tag(tags[1]) and is_lev_tag(tags[2])) {
        result = LayoutType::Tensor3D;
      }
  }

  return result;
}

std::string to_string (const FieldLayout& layout)
{
  if (layout.rank()==0) {
    return "<>()";
  }

  std::string s;
  s += "<" + e2str(layout.tags()[0]);
  for (int dim=1; dim<layout.rank(); ++dim) {
    s += "," + e2str(layout.tags()[dim]);
  }
  s += ">(" + ekat::join(layout.dims(),",") + ")";

  return s;
}

} // namespace scream
