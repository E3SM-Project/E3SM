#include "field_layout.hpp"

#include <ekat/util/ekat_string_utils.hpp>

namespace scream
{

FieldLayout::
FieldLayout ()
 : FieldLayout({},{})
{
  // Nothing to do here
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags,
                          const std::vector<int>& dims)
 : FieldLayout (tags,dims,tags2str(tags))
{
  // Nothing to do here
}

FieldLayout::FieldLayout (const std::vector<FieldTag>& tags,
                          const std::vector<int>& dims,
                          const std::vector<std::string>& names)
 : m_rank (tags.size())
 , m_tags (tags)
 , m_names(names)
 , m_dims (dims)
{
  EKAT_REQUIRE_MSG (dims.size()==tags.size(),
      "Error! Tags and dims vectors dimensions mismatch.\n"
      "  tags size: " + std::to_string(tags.size()) + "\n"
      "  dims size: " + std::to_string(dims.size()) + "\n");

  EKAT_REQUIRE_MSG (names.size()==tags.size(),
      "Error! Tags and names vectors dimensions mismatch.\n"
      "  tags size : " + std::to_string(tags.size()) + "\n"
      "  names size: " + std::to_string(names.size()) + "\n");

  set_extents ();
  compute_type ();
}

bool FieldLayout::is_vector_layout () const {
  using namespace ShortFieldTagsNames;
  return ekat::count(m_tags,CMP)==1;
}

bool FieldLayout::is_tensor_layout () const {
  return m_type==LayoutType::Tensor2D || m_type==LayoutType::Tensor3D;
}

// get the index of the CMP (Components) tag in the FieldLayout
// e.g., for FieldLayout f({COL, CMP, LEV}, {...});
// we have get_vector_component_idx(f) = 1
int FieldLayout::get_vector_component_idx () const {
  EKAT_REQUIRE_MSG (is_vector_layout(),
      "Error! 'get_vector_dim' available only for vector layouts.\n"
      "       Current layout: " + e2str(type()) + "\n");

  using namespace ShortFieldTagsNames;
  return std::distance(m_tags.begin(),ekat::find(m_tags,CMP));
}

// get the extent of the CMP (Components) tag in the FieldLayout
// e.g., for FieldLayout f({COL, CMP, LEV}, {ncol, ncmp, nlev});
// we have get_vector_dim(f) = ncmp
int FieldLayout::get_vector_dim () const {
  // since we immediately call get_vector_component_idx(), the error checking
  // there should be sufficient
  return dim(get_vector_component_idx());
}

FieldTag FieldLayout::get_vector_tag () const {
  return m_tags[get_vector_component_idx()];
}

std::vector<int> FieldLayout::get_tensor_components_ids () const {
  EKAT_REQUIRE_MSG (is_tensor_layout(),
      "Error! 'get_tensor_dims' available only for tensor layouts.\n"
      "       Current layout: " + to_string() + "\n"
      "       Layout type   : " + e2str(m_type) + "\n");

  using namespace ShortFieldTagsNames;

  std::vector<int> idx;
  auto it = m_tags.begin();
  do {
    it = std::find(it,m_tags.cend(),CMP);
    if (it!=m_tags.end()) {
      idx.push_back(std::distance(m_tags.begin(),it));
      ++it;
    }
  } while (it!=m_tags.end());

  EKAT_REQUIRE_MSG (idx.size()==2,
    "Error! Could not find a two tensor tags in the layout.\n"
    " - layout: " + to_string() + "\n"
    " - detected tags indices: " + ekat::join(idx,",") + "\n");

  return idx;
}

std::vector<int> FieldLayout::get_tensor_dims () const {
  auto idx = get_tensor_components_ids();
  for (auto& i : idx) {
    i = m_dims[i];
  }
  return idx;
}

std::vector<FieldTag> FieldLayout::get_tensor_tags () const {
  auto idx = get_tensor_components_ids();
  return {m_tags[idx[0]], m_tags[idx[1]]};
}

FieldLayout& FieldLayout::strip_dim (const FieldTag t, const bool throw_if_not_found) {
  bool any_match = false;
  for (int i=0; i<rank(); ++i) {
    if (m_tags[i] == t) {
      strip_dim(i);
      --i;
      any_match = true;
    }
  }

  EKAT_REQUIRE_MSG (any_match or not throw_if_not_found,
      "Error! No match found when attempting to strip dim.\n"
      "  - input tag: " + e2str(t) + "\n"
      "  - layout   : " + to_string() + "\n");

  return *this;
}

FieldLayout& FieldLayout::strip_dim (const int idim) {
  EKAT_REQUIRE_MSG (idim>=0 and idim<m_rank,
      "Error! Cannot strip dimension, because it is out of bounds.\n"
      "  - input dim index: " + std::to_string(idim) + "\n"
      "  - layout rank    : " + std::to_string(m_rank) + "\n");

  m_tags.erase(m_tags.begin()+idim);
  m_names.erase(m_names.begin()+idim);
  m_dims.erase(m_dims.begin()+idim);
  --m_rank;

  set_extents ();
  compute_type ();
  return *this;
}

FieldLayout&
FieldLayout::append_dim (const FieldTag t, const int extent)
{
  return append_dim(t,extent,e2str(t));
}

FieldLayout&
FieldLayout::append_dim (const FieldTag t, const int extent, const std::string& name)
{
  m_tags.push_back(t);
  m_names.push_back(name);
  m_dims.push_back(extent);

  ++m_rank;
  set_extents();
  compute_type();
  return *this;
}

FieldLayout FieldLayout::clone() const
{
  return *this;
}

FieldLayout& FieldLayout::rename_dim (const int idim, const std::string& n)
{
  EKAT_REQUIRE_MSG(idim>=0 && idim<m_rank, "Error! Index out of bounds.");

  m_names[idim] = n;
  return *this;
}

FieldLayout& FieldLayout::rename_dim (const FieldTag t, const std::string& n, const bool throw_if_not_found)
{
  bool any_match = false;
  for (int i=0; i<rank(); ++i) {
    if (m_tags[i] == t) {
      m_names[i] = n;
      any_match = true;
    }
  }

  EKAT_REQUIRE_MSG (any_match or not throw_if_not_found,
      "Error! No match found when attempting to rename dim.\n"
      "  - input tag: " + e2str(t) + "\n"
      "  - layout   : " + to_string() + "\n");

  return *this;
}

FieldLayout& FieldLayout::strip_dims (const std::vector<FieldTag>& tags)
{
  for (auto t : tags) {
    if (has_tag(t)) {
      strip_dim(t,false);
    }
  }
  return *this;
}

FieldLayout& FieldLayout::rename_dims (const std::map<FieldTag,std::string>& new_names)
{
  for (const auto& it : new_names) {
    rename_dim(it.first,it.second,false);
  }

  return *this;
}

FieldLayout& FieldLayout::reset_dim (const int idim, const int extent)
{
  EKAT_REQUIRE_MSG(idim>=0 && idim<m_rank, "Error! Index out of bounds.");

  m_dims[idim] = extent;
  set_extents();
  return *this;
}

FieldLayout& FieldLayout::reset_dim (const FieldTag t, const int extent, const bool throw_if_not_found)
{
  bool any_match = false;
  for (int i=0; i<rank(); ++i) {
    if (m_tags[i] == t) {
      m_dims[i] = extent;
      any_match = true;
    }
  }

  if (any_match) {
    set_extents();
  } else {
    EKAT_REQUIRE_MSG (not throw_if_not_found,
        "Error! No match found when attempting to reset dim.\n"
        "  - input tag: " + e2str(t) + "\n"
        "  - layout   : " + to_string() + "\n");
  }

  return *this;
}

void FieldLayout::set_extents () {
  m_extents = decltype(m_extents)("",m_rank);
  m_extents_h = Kokkos::create_mirror_view(m_extents);
  std::copy_n(m_dims.begin(),m_rank,m_extents_h.data());
  Kokkos::deep_copy(m_extents,m_extents_h);
}

void FieldLayout::compute_type () {
  if (m_rank==0) {
    m_type = LayoutType::Scalar0D;
    return;
  }

  using namespace ShortFieldTagsNames;

  using ekat::erase;
  using ekat::count;

  auto tags = this->tags();

  const int n_element = count(tags,EL);
  const int n_column  = count(tags,COL);
  const int ngp       = count(tags,GP);
  const int nvlevs    = count(tags,LEV) + count(tags,ILEV);
  const int ncomps    = count(tags,CMP);

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
    m_type = LayoutType::Scalar0D; return;
  } else if (tags.size()==1 and tags[0]==CMP) {
    m_type = LayoutType::Vector0D; return;
  } else if (tags.size()==1 and nvlevs==1) {
    m_type = LayoutType::Scalar1D; return;
  } else if (tags[0]==CMP and tags[1]==CMP) {
    m_type = LayoutType::Tensor0D; return;
  } else if (tags.size()==2 and ncomps==1 and nvlevs==1) {
    m_type = LayoutType::Vector1D; return;
  } else {
    // Not a supported layout.
    m_type = LayoutType::Invalid; return;
  }

  // Get the size of what's left
  const auto size = tags.size();
  auto is_lev_tag = [](const FieldTag t) {
    return t==LEV or t==ILEV;
  };
  switch (size) {
    case 0:
      m_type = LayoutType::Scalar2D;
      break;
    case 1:
      // The only tag left should be 'CMP', 'TL', or 'LEV'/'ILEV'
      if (tags[0]==CMP || tags[0]==TL) {
        m_type = LayoutType::Vector2D;
      } else if (is_lev_tag(tags[0])) {
        m_type = LayoutType::Scalar3D;
      }
      break;
    case 2:
      // Possible supported scenarios:
      //  1) <CMP|TL,LEV|ILEV>
      //  2) <TL,CMP>
      //  3) <CMP,CMP>
      if ( is_lev_tag(tags[1]) && (tags[0]==CMP || tags[0]==TL)) {
        m_type = LayoutType::Vector3D;
      } else if (tags[0]==TL && tags[1]==CMP ) {
        m_type = LayoutType::Tensor2D;
      } else if (tags[0]==CMP and tags[1]==CMP) {
        m_type = LayoutType::Tensor2D;
      }
      break;
    case 3:
      // The only supported scenarios are:
      //  1) <TL,  CMP, LEV|ILEV>
      //  2) <CMP, CMP, LEV|ILEV>
      if ( tags[0]==TL && tags[1]==CMP && is_lev_tag(tags[2])) {
        m_type = LayoutType::Tensor3D;
      } else if ( tags[0]==CMP && tags[1]==CMP && is_lev_tag(tags[2])) {
        m_type = LayoutType::Tensor3D;
      }

      break;
    default:
      // If nothing worked, this type is not recognized
      m_type = LayoutType::Invalid;
  }
}

std::string FieldLayout::to_string () const
{
  std::string s;
  // Note: I don't know how to pass e2str directly to ekat::join,
  //       so wrap it in a lambda
  s += "<" + ekat::join(m_names,",") + ">";
  s += "(" + ekat::join(m_dims,",") + ")";

  return s;
}

} // namespace scream
