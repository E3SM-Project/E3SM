#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

AbstractRemapper::
AbstractRemapper (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid)
{
  set_grids (src_grid,tgt_grid);
}

void AbstractRemapper::
registration_begins () {
  EKAT_REQUIRE_MSG(m_state==RepoState::Clean,
                       "Error! Cannot start registration on a non-clean repo.\n"
                       "       Did you call 'registration_begins' already?\n");

  do_registration_begins();

  m_state = RepoState::Open;
}

void AbstractRemapper::
register_field (const identifier_type& src, const identifier_type& tgt) {
  EKAT_REQUIRE_MSG(m_state!=RepoState::Clean,
                       "Error! Cannot register fields in the remapper at this time.\n"
                       "       Did you forget to call 'registration_begins' ?");
  EKAT_REQUIRE_MSG(m_state!=RepoState::Closed,
                       "Error! Cannot register fields in the remapper at this time.\n"
                       "       Did you accidentally call 'registration_ends' already?");

  EKAT_REQUIRE_MSG(src.get_grid_name()==m_src_grid->name(),
                       "Error! Source field stores the wrong grid.\n");
  EKAT_REQUIRE_MSG(tgt.get_grid_name()==m_tgt_grid->name(),
                       "Error! Target field stores the wrong grid.\n");

  EKAT_REQUIRE_MSG(compatible_layouts(src.get_layout(),tgt.get_layout()),
                     "Error! Source and target layouts are not compatible.\n");

  do_register_field (src,tgt);

  m_fields_are_bound.push_back(false);
  ++m_num_registered_fields;
}

void AbstractRemapper::
register_field (const field_type& src, const field_type& tgt) {
  register_field(src.get_header().get_identifier(),
                 tgt.get_header().get_identifier());
  bind_field(src,tgt);
}

void AbstractRemapper::
bind_field (const field_type& src, const field_type& tgt) {
  EKAT_REQUIRE_MSG(m_state!=RepoState::Clean,
                     "Error! Cannot bind fields in the remapper at this time.\n"
                     "       Did you forget to call 'registration_begins' ?");

  const auto& src_fid = src.get_header().get_identifier();
  const auto& tgt_fid = tgt.get_header().get_identifier();

  // Try to locate the pair of fields
  const int ifield = find_field(src_fid, tgt_fid);
  EKAT_REQUIRE_MSG(ifield>=0,
                     "Error! The src/tgt field pair\n"
                     "         " + src_fid.get_id_string() + "\n"
                     "         " + tgt_fid.get_id_string() + "\n"
                     "       was not registered. Please, register fields before binding them.\n");

  EKAT_REQUIRE_MSG(src.is_allocated(), "Error! Source field is not yet allocated.\n");
  EKAT_REQUIRE_MSG(tgt.is_allocated(), "Error! Target field is not yet allocated.\n");

  EKAT_REQUIRE_MSG(!m_fields_are_bound[ifield],
                     "Error! Field " + src_fid.get_id_string() + " already bound.\n");

  do_bind_field(ifield,src,tgt);

  m_fields_are_bound[ifield] = true;
  ++m_num_bound_fields;

  // Check if src and/or tgt is read-only, precluding bwd and/or fwd remap
  if (src.is_read_only()) {
    m_bwd_allowed = false;
  }
  if (tgt.is_read_only()) {
    m_fwd_allowed = false;
  }
}

void AbstractRemapper::
registration_ends () {
  EKAT_REQUIRE_MSG(m_state!=RepoState::Closed,
                       "Error! Cannot call registration_ends at this time.\n"
                       "       Did you accidentally call 'registration_ends' already?");

  m_num_fields = m_num_registered_fields;

  do_registration_ends();

  m_state = RepoState::Closed;
}

void AbstractRemapper::remap (const bool forward) {
  EKAT_REQUIRE_MSG(m_state!=RepoState::Open,
                     "Error! Cannot perform remapping at this time.\n"
                     "       Did you forget to call 'registration_ends'?\n");

  EKAT_REQUIRE_MSG(m_num_bound_fields==m_num_fields,
                     "Error! Not all fields have been set in the remapper.\n"
                     "       In particular, field " +
                     std::to_string(std::distance(m_fields_are_bound.begin(),std::find(m_fields_are_bound.begin(),m_fields_are_bound.end(),false))) +
                     " has not been bound.\n");

  if (m_state!=RepoState::Clean) {
    if (forward) {
      EKAT_REQUIRE_MSG (m_fwd_allowed,
                       "Error! Forward remap is not allowed by this remapper.\n"
                       "       This means that some fields on the target grid are read-only.\n");
      do_remap_fwd ();
    } else {
      EKAT_REQUIRE_MSG (m_bwd_allowed,
                       "Error! Backward remap is not allowed by this remapper.\n"
                       "       This means that some fields on the source grid are read-only.\n");
      do_remap_bwd ();
    }
  }
}

void AbstractRemapper::
set_grids (const grid_ptr_type& src_grid,
           const grid_ptr_type& tgt_grid)
{
  EKAT_REQUIRE_MSG(static_cast<bool>(src_grid), "Error! Invalid source grid pointer.\n");
  EKAT_REQUIRE_MSG(static_cast<bool>(tgt_grid), "Error! Invalid target grid pointer.\n");

  m_src_grid = src_grid;
  m_tgt_grid = tgt_grid;
}

} // namespace scream
