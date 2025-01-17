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

FieldLayout IdentityRemapper::create_src_layout (const FieldLayout& tgt_layout) const
{
  EKAT_REQUIRE_MSG (is_valid_tgt_layout(tgt_layout),
      "[IdentityRemapper] Error! Input target layout is not valid for this remapper.\n"
      " - input layout: " + tgt_layout.to_string());

  // Src and tgt grids are the same, so return the input
  return tgt_layout;
}

FieldLayout IdentityRemapper::create_tgt_layout (const FieldLayout& src_layout) const
{
  EKAT_REQUIRE_MSG (is_valid_src_layout(src_layout),
      "[IdentityRemapper] Error! Input source layout is not valid for this remapper.\n"
      " - input layout: " + src_layout.to_string());

  // Src and tgt grids are the same, so return the input
  return src_layout;
}

bool IdentityRemapper::
compatible_layouts (const FieldLayout& src, const FieldLayout& tgt) const
{
  return src.congruent(tgt);
}

void IdentityRemapper::registration_ends_impl ()
{
  // If there is aliasing, src and tgt must be the SAME field
  if (m_aliasing!=Aliasing::NoAliasing) {
    for (int i=0; i<m_num_fields; ++i) {
      const auto& src = m_src_fields[i];
      const auto& tgt = m_tgt_fields[i];
      EKAT_REQUIRE_MSG (src.equivalent(tgt),
          "Error! Input fields are not aliasing eahc other, but aliasing was requested.\n"
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
