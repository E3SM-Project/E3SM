#include "share/field/field_header.hpp"

#include "ekat/std_meta/ekat_std_utils.hpp"

namespace scream
{

FieldHeader::FieldHeader (const identifier_type& id)
 : m_identifier (id)
 , m_tracking (create_tracking())
 , m_alloc_prop (get_type_size(id.data_type()))
{
  // Nothing to be done here
}

void FieldHeader::
set_extra_data (const std::string& key,
                const ekat::any& data,
                const bool throw_if_existing)
{
  if (throw_if_existing) {
    EKAT_REQUIRE_MSG (m_extra_data.find(key)==m_extra_data.end(),
                        "Error! Key '" + key + "' already existing in "
                        "the extra data map of field '" + m_identifier.get_id_string() + "'.\n");
    m_extra_data[key] = data;
  } else {
    m_extra_data[key] = data;
  }
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
  EKAT_REQUIRE_MSG (id.get_layout_ptr()!=nullptr,
      "Error! Input field identifier has an invalid layout pointer.\n");

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
  fh->m_alloc_prop = parent->get_alloc_properties().subview(idim,k,dynamic);
  fh->m_alloc_prop.commit(id.get_layout_ptr());

  return fh;
}

} // namespace scream
