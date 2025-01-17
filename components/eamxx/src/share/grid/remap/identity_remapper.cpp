#include "share/grid/remap/identity_remapper.hpp"

namespace scream
{

IdentityRemapper::
IdentityRemapper (const grid_ptr_type grid,
                  const Aliasing aliasing)
 : AbstractRemapper (grid,grid)
{
  set_aliasing(aliasing);
}

void IdentityRemapper::set_aliasing (const Aliasing aliasing) {
  EKAT_REQUIRE_MSG (get_state()==RepoState::Clean,
      "Error! Aliasing in IdentityRemapper must be set *before* registration starts.\n");
  m_aliasing = aliasing;
}

void IdentityRemapper::register_field_from_src (const Field& src)
{
  EKAT_REQUIRE_MSG (m_aliasing!=SrcAliasTgt,
      "Error! Makes no sense to register from src and ask that src alias tgt.\n");
  if (m_aliasing==TgtAliasSrc) {
    register_field(src,src);
  } else {
    AbstractRemapper::register_field_from_src(src);
  }
}
void IdentityRemapper::register_field_from_tgt (const Field& tgt)
{
  EKAT_REQUIRE_MSG (m_aliasing!=TgtAliasSrc,
      "Error! Makes no sense to register from tgt and ask that tgt alias src.\n");
  if (m_aliasing==SrcAliasTgt) {
    register_field(tgt,tgt);
  } else {
    AbstractRemapper::register_field_from_tgt(tgt);
  }
}

void IdentityRemapper::registration_ends_impl ()
{
  // If there is aliasing, src and tgt must be the SAME field
  if (m_aliasing!=Aliasing::NoAliasing) {
    for (int i=0; i<m_num_fields; ++i) {
      const auto& src = m_src_fields[i];
      const auto& tgt = m_tgt_fields[i];
      EKAT_REQUIRE_MSG (src.is_aliasing(tgt),
          "Error! Input fields are not aliasing each other, but aliasing was requested in the IdentityRemapper.\n"
          "       To register field when aliasing is active, use register_field_from_tgt/src.\n");
    }
  }
}

void IdentityRemapper::remap_fwd_impl () {
  // Only do deep copies if there's no aliasing
  if (m_aliasing==Aliasing::NoAliasing) {
    for (int i=0; i<m_num_fields; ++i) {
      m_tgt_fields[i].deep_copy(m_src_fields[i]);
    }
  }
}

void IdentityRemapper::remap_bwd_impl () {
  // Only do deep copies if there's no aliasing
  if (m_aliasing==Aliasing::NoAliasing) {
    for (int i=0; i<m_num_fields; ++i) {
      m_src_fields[i].deep_copy(m_tgt_fields[i]);
    }
  }
}

} // namespace scream
