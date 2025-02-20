#ifndef SCREAM_FIELD_LAYOUT_HPP
#define SCREAM_FIELD_LAYOUT_HPP

#include "share/field/field_tag.hpp"
#include "share/eamxx_types.hpp"

#include <ekat/std_meta/ekat_std_utils.hpp>
#include <ekat/util/ekat_string_utils.hpp>
#include <ekat/ekat_assert.hpp>

#include <string>
#include <vector>

namespace scream
{

// The type of the layout, that is, the kind of field it represents.
enum class LayoutType {
  Invalid,
  Scalar0D,
  Vector0D,
  Tensor0D,
  Scalar1D,
  Vector1D,
  Scalar2D,
  Vector2D,
  Tensor2D,
  Scalar3D,
  Vector3D,
  Tensor3D
};

inline std::string e2str (const LayoutType lt) {
  std::string name;
  switch (lt) {
    case LayoutType::Scalar0D: name = "Scalar0D"; break;
    case LayoutType::Vector0D: name = "Vector0D"; break;
    case LayoutType::Scalar1D: name = "Scalar1D"; break;
    case LayoutType::Vector1D: name = "Vector1D"; break;
    case LayoutType::Scalar2D: name = "Scalar2D"; break;
    case LayoutType::Vector2D: name = "Vector2D"; break;
    case LayoutType::Tensor2D: name = "Tensor2D"; break;
    case LayoutType::Scalar3D: name = "Scalar3D"; break;
    case LayoutType::Vector3D: name = "Vector3D"; break;
    case LayoutType::Tensor3D: name = "Tensor3D"; break;
    case LayoutType::Invalid:  name = "INVALID" ; break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized LayoutType.\n");
  }
  return name;
}

/*
 *  A small class to hold basic info about a field layout
 *
 *  Note: the content of extents() is the same as the content of dims().
 *        The difference is that the former returns a device view, while
 *        the latter returns a std::vector.
 */

class FieldLayout {
public:
  using extents_type = typename KokkosTypes<DefaultDevice>::view_1d<int>;

  // Constructor(s)
  FieldLayout ();
  FieldLayout (const FieldLayout&) = default;
  FieldLayout (const std::vector<FieldTag>& tags,
               const std::vector<int>& dims);
  FieldLayout (const std::vector<FieldTag>& tags,
               const std::vector<int>& dims,
               const std::vector<std::string>& names);

  // Assignment (defaulted)
  FieldLayout& operator= (const FieldLayout&) = default;

  // Create invalid layout
  static FieldLayout invalid () { return FieldLayout({FieldTag::Invalid},{0}); }

  // ----- Getters ----- //

  LayoutType type () const { return m_type; }

  // Name and layout informations
  const std::vector<std::string>& names () const { return m_names; }
  const std::vector<FieldTag>& tags () const { return m_tags; }
  FieldTag tag  (const int idim) const;
  const std::string& name (const int idim) const;
  bool has_tag (const FieldTag t) const { return ekat::contains(m_tags,t); }
  bool has_tags (const std::vector<FieldTag>& tags) const;

  // The rank is the number of tags associated to this field.
  int rank () const  { return m_rank; }

  int dim_idx (const FieldTag t) const;

  int dim (const std::string& name) const;
  int dim (const FieldTag tag) const;
  int dim (const int idim) const;
  const std::vector<int>& dims () const { return m_dims; }
  const extents_type& extents () const { return m_extents; }
  const extents_type::HostMirror& extents_h () const { return m_extents_h; }

  long long  size () const;

  bool are_dimensions_set () const;

  // Check if this layout is that of a vector/tensor field
  bool is_vector_layout () const;
  bool is_tensor_layout () const;

  // If this is the layout of a vector field, get the idx of the
  // vector (CMP, Component) dimension
  // Note: throws if is_vector_layout()==false.
  int get_vector_component_idx () const;
  // get the dimension (extent) of the vector (CMP, Component) dimension
  // calls get_vector_component_idx()
  int get_vector_dim () const;
  FieldTag get_vector_tag () const;

  // If this is the layout of a tensor field, get the idx of the tensor dimensions
  // Note: throws if is_tensor_layout()==false.
  std::vector<int> get_tensor_components_ids () const;
  // Get the dimension (extent) of the tensor components. Calls get_tensor_components_ids
  std::vector<int> get_tensor_dims () const;
  std::vector<FieldTag> get_tensor_tags () const;

  // Change this layout by adding/removing a dimension or changing its extent/name
  // NOTE: the strip_dim/rename_dim/reset_dim overloads with FieldTag will alter *all*
  //       dimension matching the input tag
  FieldLayout& strip_dim (const FieldTag tag, const bool throw_if_not_found = true);
  FieldLayout& strip_dim (const int idim);
  FieldLayout& append_dim (const FieldTag t, const int extent);
  FieldLayout& append_dim (const FieldTag t, const int extent, const std::string& name);
  FieldLayout& rename_dim (const int idim, const std::string& n);
  FieldLayout& rename_dim (const FieldTag tag, const std::string& n, const bool throw_if_not_found = true);
  FieldLayout& reset_dim (const int idim, const int extent);
  FieldLayout& reset_dim (const FieldTag t, const int extent, const bool throw_if_not_found = true);

  // These overload allow to remove/rename dims *if found*. They won't throw if layout does not have them
  FieldLayout& strip_dims (const std::vector<FieldTag>& tags); // Does not throw if not found
  FieldLayout& rename_dims (const std::map<FieldTag,std::string>& new_names); // Does not throw if not found

  FieldLayout clone() const;

  // NOTE: congruent does not check the tags names. It only checks
  //       rank, m_tags, and m_dims. Use operator== if names are important
  bool congruent (const FieldLayout& rhs) const;

  // For printing purposes
  std::string to_string () const;

protected:
  void compute_type ();
  void set_extents ();

  int                       m_rank;
  std::vector<FieldTag>     m_tags;
  std::vector<std::string>  m_names;
  std::vector<int>          m_dims;
  extents_type              m_extents;
  extents_type::HostMirror  m_extents_h;

  LayoutType                m_type;
};

bool operator== (const FieldLayout& fl1, const FieldLayout& fl2);

// ========================== IMPLEMENTATION ======================= //

inline int FieldLayout::dim_idx (const FieldTag t) const {
  // Check exactly one tag (no ambiguity)
  EKAT_REQUIRE_MSG(ekat::count(m_tags,t)==1,
      "Error! FieldTag::dim_idx requires that the tag appears exactly once.\n"
      "  - field tag: " + e2str(t) + "\n"
      "  - tag count: " + std::to_string(ekat::count(m_tags,t)) + "\n");

  return std::distance(m_tags.begin(),ekat::find(m_tags,t));
}

// returns extent
inline int FieldLayout::dim (const FieldTag t) const {
  return m_dims[dim_idx(t)];
}

inline int FieldLayout::dim (const std::string& name) const {
  auto it = ekat::find(m_names,name);

  // Check if found
  EKAT_REQUIRE_MSG(it!=m_names.end(),
      "Error! Dim name '" + name + "' not found in this layout.\n"
      "  - layout dims: " + ekat::join(m_names,",") + "\n");

  // Check only one tag (no ambiguity)
  EKAT_REQUIRE_MSG(ekat::count(m_names,name)==1,
      "Error! Dimension name '" + name + "' appears multiple times.\n"
      "  - layout dims: " + ekat::join(m_names,",") + "\n");

  return m_dims[std::distance(m_names.begin(),it)];
}

inline int FieldLayout::dim (const int idim) const {
  EKAT_REQUIRE_MSG (idim>=0 && idim<m_rank, "Error! Index out of bounds.");
  return m_dims[idim];
}

inline long long FieldLayout::size () const {
  EKAT_REQUIRE_MSG(are_dimensions_set(),
      "Error! Field dimensions not yet set.\n");
  long long prod = 1;
  for (int idim=0; idim<m_rank; ++idim) {
    prod *= m_dims[idim];
  }
  return prod;
}

inline FieldTag FieldLayout::tag (const int idim) const { 
  EKAT_REQUIRE_MSG (idim>=0 && idim<m_rank, "Error! Index out of bounds.");
  return m_tags[idim];
}

inline const std::string& FieldLayout::name (const int idim) const
{
  EKAT_REQUIRE_MSG (idim>=0 && idim<m_rank, "Error! Index out of bounds.");
  return m_names[idim];
}

inline bool FieldLayout::has_tags (const std::vector<FieldTag>& tags) const {
  bool b = true;
  for (auto t : tags) {
    b &= has_tag(t);
  }
  return b;
}

inline bool FieldLayout::are_dimensions_set () const {
  for (int idim=0; idim<m_rank; ++idim) {
    if (m_dims[idim]<0) {
      return false;
    }
  }
  return true;
}

inline bool FieldLayout::congruent (const FieldLayout& rhs) const {
  return rank()==rhs.rank() &&
         tags()==rhs.tags() &&
         dims()==rhs.dims();
}

inline bool operator== (const FieldLayout& fl1, const FieldLayout& fl2) {
  return fl1.congruent(fl2) and fl1.names()==fl2.names();
}

} // namespace scream

#endif // SCREAM_FIELD_LAYOUT_HPP
