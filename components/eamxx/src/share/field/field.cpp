#include "share/field/field.hpp"
#include "share/util/eamxx_utils.hpp"

namespace scream
{

Field::
Field (const identifier_type& id)
 : m_header (create_header(id))
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
Field::alias (const std::string& name) const {
  Field f;
  f.m_header = get_header().alias(name);
  f.m_data = m_data;
  f.m_is_read_only = m_is_read_only;
  return f;
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
sync_to_host (const bool fence) const {
  // Sanity check
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to sync host and device views.\n");

  // Check for early return if Host and Device are the same memory space
  if (host_and_device_share_memory_space()) return;

  // We allow sync_to_host for constant fields. Temporarily disable read only flag.
  const bool original_read_only = m_is_read_only;
  m_is_read_only = false;

  switch (data_type()) {
    case DataType::IntType:
      sync_views_impl<int, Device, Host>();
      break;
    case DataType::FloatType:
      sync_views_impl<float, Device, Host>();
      break;
    case DataType::DoubleType:
      sync_views_impl<double, Device, Host>();
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized field data type in Field::sync_to_host.\n");
  }

  if (fence) Kokkos::fence();

  // Return field to read-only state
  m_is_read_only = original_read_only;
}

void Field::
sync_to_dev (const bool fence) const {
  // Sanity check
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to sync host and device views.\n");

  // Check for early return if Host and Device are the same memory space
  if (host_and_device_share_memory_space()) return;

  switch (data_type()) {
    case DataType::IntType:
      sync_views_impl<int, Host, Device>();
      break;
    case DataType::FloatType:
      sync_views_impl<float, Host, Device>();
      break;
    case DataType::DoubleType:
      sync_views_impl<double, Host, Device>();
      break;
    default:
      EKAT_ERROR_MSG("Error! Unrecognized field data type in Field::sync_to_dev.\n");
  }

  if (fence) Kokkos::fence();
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
  FieldIdentifier sf_id(sf_name,lt.clone().strip_dim(idim),sf_units,id.get_grid_name(),id.data_type());

  // Create empty subfield, then set header and views
  // Note: we can access protected members, since it's the same type
  Field sf;
  sf.m_header = create_subfield_header(sf_id,m_header,idim,index,dynamic);
  sf.m_data = m_data;
  sf.m_is_read_only = m_is_read_only;

  if (not sf.m_header->get_alloc_properties().contiguous() and
      not sf.host_and_device_share_memory_space()) {
    // If subfield is not contiguous and Host and Device do not
    // share a memory space, we must initialize the helper field
    // for sync_to functions.
    sf.initialize_contiguous_helper_field();
  }

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
subfield (const FieldTag tag, const int index, const bool dynamic) const {
  int idim = get_header().get_identifier().get_layout().dim_idx(tag);
  return subfield(idim,index,dynamic);
}

// slice at index idim, extracting the N = (index_end - index_beg) entries
// written in math notation: [index_beg, index_end)
// or equivalently, subF = F(index_beg, ... , index_beg + N)
Field Field::subfield(const std::string& sf_name,
                      const ekat::units::Units& sf_units, const int idim,
                      const int index_beg, const int index_end) const {

  const auto& id = m_header->get_identifier();
  const auto& lt = id.get_layout();

  // Sanity checks
  EKAT_REQUIRE_MSG(
      is_allocated(),
      "Error! Input field must be allocated in order to subview it.\n");

  auto sf_layout = lt.clone();
  sf_layout.reset_dim(idim, index_end - index_beg);
  // Create identifier for subfield
  FieldIdentifier sf_id(sf_name, sf_layout, sf_units, id.get_grid_name(), id.data_type());

  // Create empty subfield, then set header and views
  // Note: we can access protected members, since it's the same type
  Field sf;
  sf.m_header = create_subfield_header(sf_id, m_header, idim, index_beg,
                                       index_end);
  sf.m_data = m_data;

  if (not sf.m_header->get_alloc_properties().contiguous() and
      not sf.host_and_device_share_memory_space()) {
    // If subfield is not contiguous and Host and Device do not
    // share a memory space, we must initialize the helper field
    // for sync_to functions.
    sf.initialize_contiguous_helper_field();
  }

  return sf;
}

Field Field::subfield(const std::string& sf_name, const int idim,
                      const int index_beg, const int index_end) const {
  const auto& id = m_header->get_identifier();
  return subfield(sf_name, id.get_units(), idim, index_beg, index_end);
}

Field Field::subfield(const int idim, const int index_beg,
                      const int index_end) const {
  return subfield(m_header->get_identifier().name(), idim, index_beg,
                  index_end);
}

Field Field::
get_component (const int i, const bool dynamic) {
  const auto& layout = get_header().get_identifier().get_layout();
  const auto& fname = get_header().get_identifier().name();
  EKAT_REQUIRE_MSG (layout.is_vector_layout(),
      "Error! 'get_component' available only for vector fields.\n"
      "       Layout of '" + fname + "': " + e2str(layout.type()) + "\n");

  const int idim = layout.get_vector_component_idx();
  EKAT_REQUIRE_MSG (i>=0 && i<layout.dim(idim),
      "Error! Component index out of bounds [0," + std::to_string(layout.dim(idim)) + ").\n");

  // Add _$i to the field name, to avoid issues if the subfield is stored
  // in some structure that requires unique names (e.g., a remapper)
  return subfield (fname + "_" + std::to_string(i),idim,i,dynamic);
}

Field Field::get_components(const int beg, const int end) {
  const auto& layout = get_header().get_identifier().get_layout();
  const auto& fname = get_header().get_identifier().name();
  EKAT_REQUIRE_MSG(layout.is_vector_layout(),
                   "Error! 'get_component' available only for vector fields.\n"
                   "       Layout of '" +
                       fname + "': " + e2str(layout.type()) + "\n");

  const int idim = layout.get_vector_component_idx();
  EKAT_REQUIRE_MSG(beg >= 0 && end < layout.dim(idim),
                   "Error! Component index range out of bounds [0," +
                       std::to_string(layout.dim(idim)) + ").\n");
  EKAT_REQUIRE_MSG(beg < end, "Error! Invalid component indices (beg >= end).\n");

  // Add _$beg-$end to the field name, to avoid issues if the subfield is stored
  // in some structure that requires unique names (e.g., a remapper)
  return subfield(fname + "_" + std::to_string(beg) + "-" + std::to_string(end),
                  idim, beg, end);
}

bool Field::is_aliasing(const Field& rhs) const
{
  if (this==&rhs)
    return true;  // Same object

  if (not is_allocated() or not rhs.is_allocated())
    return false; // Once allocated, they will be different

  // NOTE: I'm not sure we NEED to check m_data, but we might as well
  return m_header->is_aliasing(rhs.get_header()) and
         m_data.d_view==rhs.m_data.d_view;
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
  const auto& layout = id.get_layout();
  auto& alloc_prop   = m_header->get_alloc_properties();

  // Commit the allocation properties
  alloc_prop.commit(layout);

  // Create the view, by quering allocation properties for the allocation size
  const auto view_dim = alloc_prop.get_alloc_size();

  m_data.d_view = decltype(m_data.d_view)(id.name(),view_dim);
  m_data.h_view = Kokkos::create_mirror_view(m_data.d_view);
}

} // namespace scream
