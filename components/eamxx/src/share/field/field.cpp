#include "share/field/field.hpp"
#include "share/util/scream_utils.hpp"

namespace scream
{

// ================================= IMPLEMENTATION ================================== //

Field::
Field (const identifier_type& id)
 : m_header     (create_header(id))
{
  // Nothing to do here
}

Field
Field::get_const() const {
  Field f(*this);
  f.m_is_read_only = true;
  return f;
}

Field
Field::clone() const {
  return clone(name());
}

Field
Field::clone(const std::string& name) const {
  // Create new field
  const auto& my_fid = get_header().get_identifier();
  FieldIdentifier fid(name,my_fid.get_layout(),my_fid.get_units(),
                      my_fid.get_grid_name(),my_fid.data_type());
  Field f(fid);

  // Ensure alloc props match
  const auto&  ap = get_header().get_alloc_properties();
        auto& fap = f.get_header().get_alloc_properties();
  fap.request_allocation(ap.get_largest_pack_size());

  // Allocate
  f.allocate_view();

  // Set correct time stamp
  const auto& ts = get_header().get_tracking().get_time_stamp();
  f.get_header().get_tracking().update_time_stamp(ts);

  // Deep copy
  f.deep_copy<Device>(*this);
  f.deep_copy<Host>(*this);

  return f;
}

void Field::
sync_to_host () const {
  // Sanity check
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to sync host and device views.\n");

  Kokkos::deep_copy(m_data.h_view,m_data.d_view);
}

void Field::
sync_to_dev () const {
  // Sanity check
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to sync host and device views.\n");

  // Ensure host view was created (lazy construction)
  Kokkos::deep_copy(m_data.d_view,m_data.h_view);
}

Field Field::
subfield (const std::string& sf_name, const ekat::units::Units& sf_units,
          const int idim, const int index, const bool dynamic) const {

  const auto& id = m_header->get_identifier();
  const auto& lt = id.get_layout();

  // Sanity checks
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to subview it.\n");
  EKAT_REQUIRE_MSG (idim==0 || idim==1,
        "Error! Subview dimension index must be either 0 or 1.\n");

  // Create identifier for subfield
  FieldIdentifier sf_id(sf_name,lt.strip_dim(idim),sf_units,id.get_grid_name());

  // Create empty subfield, then set header and views
  // Note: we can access protected members, since it's the same type
  Field sf;
  sf.m_header = create_subfield_header(sf_id,m_header,idim,index,dynamic);
  sf.m_data = m_data;

  return sf;
}

Field Field::
subfield (const std::string& sf_name, const int idim, const int index, const bool dynamic) const {
  const auto& id = m_header->get_identifier();
  return subfield(sf_name,id.get_units(),idim,index,dynamic);
}

Field Field::
subfield (const int idim, const int index, const bool dynamic) const {
  return subfield(m_header->get_identifier().name(),idim,index,dynamic);
}

Field Field::
get_component (const int i, const bool dynamic) {
  const auto& layout = get_header().get_identifier().get_layout();
  const auto& fname = get_header().get_identifier().name();
  EKAT_REQUIRE_MSG (layout.is_vector_layout(),
      "Error! 'get_component' available only for vector fields.\n"
      "       Layout of '" + fname + "': " + e2str(get_layout_type(layout.tags())) + "\n");

  const int idim = layout.get_vector_dim();
  EKAT_REQUIRE_MSG (i>=0 && i<layout.dim(idim),
      "Error! Component index out of bounds [0," + std::to_string(layout.dim(idim)) + ").\n");

  // Add _$idim to the field name, to avoid issues if the subfield is stored
  // in some structure that requires unique names (e.g., a remapper)
  return subfield (fname + "_" + std::to_string(i),idim,i,dynamic);
}

bool Field::equivalent(const Field& rhs) const
{
  return (m_header==rhs.m_header &&
          is_allocated() &&
          m_data.d_view==rhs.m_data.d_view &&
          m_data.h_view==rhs.m_data.h_view);
}

void Field::allocate_view ()
{
  // Not sure if simply returning would be safe enough. Re-allocating
  // would definitely be error prone (someone may have already gotten
  // a subview of the field). However, it *seems* suspicious to call
  // this method twice, and I think it's more likely than not that
  // such a scenario would indicate a bug. Therefore, I am prohibiting it.
  EKAT_REQUIRE_MSG(!is_allocated(), "Error! View was already allocated.\n");

  // Short names
  const auto& id     = m_header->get_identifier();
  const auto& layout = id.get_layout_ptr();
  auto& alloc_prop   = m_header->get_alloc_properties();

  // Check the identifier has all the dimensions set
  EKAT_REQUIRE_MSG(layout->are_dimensions_set(),
      "Error! Cannot allocate the view until all the field's dimensions are set.\n");

  // Commit the allocation properties
  alloc_prop.commit(layout);

  // Create the view, by quering allocation properties for the allocation size
  const auto view_dim = alloc_prop.get_alloc_size();

  m_data.d_view = decltype(m_data.d_view)(id.name(),view_dim);
  m_data.h_view = Kokkos::create_mirror_view(m_data.d_view);
}

} // namespace scream
