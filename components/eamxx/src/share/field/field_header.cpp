#include "share/field/field_header.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

namespace scream {

FieldHeader::FieldHeader (const identifier_type& id)
 : m_identifier (id)
 , m_tracking (create_tracking())
{
  m_alloc_prop = std::make_shared<FieldAllocProp>(get_type_size(id.data_type()));
  m_extra_data = std::make_shared<extra_data_type>();

  // Let's add immediately this att, so that users don't need to check
  // if it already exist before adding string attributes for io.
  using stratts_t = std::map<std::string,std::string>;
  set_extra_data("io: string attributes",stratts_t());
}

void FieldHeader::
set_extra_data (const std::string& key,
                const ekat::any& data,
                const bool throw_if_existing)
{
  if (throw_if_existing) {
    EKAT_REQUIRE_MSG (m_extra_data->find(key)==m_extra_data->end(),
                        "Error! Key '" + key + "' already existing in "
                        "the extra data map of field '" + m_identifier.get_id_string() + "'.\n");
    (*m_extra_data)[key] = data;
  } else {
    (*m_extra_data)[key] = data;
  }
}

std::shared_ptr<FieldHeader> FieldHeader::alias(const std::string& name) const {
  auto fh = create_header(get_identifier().alias(name));
  fh->m_tracking = m_tracking;
  fh->m_alloc_prop = m_alloc_prop;
  fh->m_extra_data = m_extra_data;
  return fh;
}

bool FieldHeader::is_aliasing (const FieldHeader& rhs) const
{
  if (this==&rhs)
    return true;

  if (m_tracking==rhs.m_tracking and
      m_alloc_prop==rhs.m_alloc_prop and
      m_extra_data==rhs.m_extra_data)
    return true;

  auto p = get_parent();
  auto rhs_p = rhs.get_parent();
  if (p!=nullptr and rhs_p!=nullptr) {
    return p->is_aliasing(*rhs_p) and
           m_alloc_prop->get_subview_info()==rhs.m_alloc_prop->get_subview_info();
  }

  return false;
}

// ---------------- Free function -------------------- //

std::shared_ptr<FieldHeader>
create_subfield_header (const FieldIdentifier& id,
                        std::shared_ptr<FieldHeader> parent,
                        const int idim, const int k, const bool dynamic)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (parent!=nullptr,
      "Error! Invalid pointer for parent header.\n");

  // Create header, and set up parent/child
  auto fh = create_header(id);
  fh->create_parent_child_link(parent);

  // Create tracking, and set up parent/child
  fh->m_tracking = create_tracking();
  fh->m_tracking->create_parent_child_link(parent->get_tracking_ptr());
  if (parent->get_tracking().get_time_stamp().is_valid()) {
    fh->m_tracking->update_time_stamp(parent->get_tracking().get_time_stamp());
  }

  // Create alloc props
  fh->m_alloc_prop = std::make_shared<FieldAllocProp>(parent->get_alloc_properties().subview(idim,k,dynamic));

  return fh;
}

// subfield with multiple, contiguous slices
std::shared_ptr<FieldHeader>
create_subfield_header(const FieldIdentifier& id,
                       std::shared_ptr<FieldHeader> parent, const int idim,
                       const int k_beg, const int k_end) {
  // Sanity checks
  EKAT_REQUIRE_MSG(parent != nullptr,
                   "Error! Invalid pointer for parent header.\n");
  EKAT_REQUIRE_MSG(k_end > k_beg,
                   "Error! Slice indices are invalid (non-increasing).\n");

  // Create header, and set up parent/child
  auto fh = create_header(id);
  fh->create_parent_child_link(parent);

  // Create tracking, and set up parent/child
  fh->m_tracking = create_tracking();
  fh->m_tracking->create_parent_child_link(parent->get_tracking_ptr());
  if (parent->get_tracking().get_time_stamp().is_valid()) {
    fh->m_tracking->update_time_stamp(parent->get_tracking().get_time_stamp());
  }

  // Create alloc props
  fh->m_alloc_prop = std::make_shared<FieldAllocProp>(
      parent->get_alloc_properties().subview(idim, k_beg, k_end));

  return fh;
}

} // namespace scream
