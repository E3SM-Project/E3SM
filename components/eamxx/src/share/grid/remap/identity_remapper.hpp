#ifndef SCREAM_IDENTITY_REMAPPER_HPP
#define SCREAM_IDENTITY_REMAPPER_HPP

#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

/*
 *  A remapper representing an 'identity'.
 *
 *  This remapper effectively does nothing, since its source and target
 *  grids are the same. There is no *real* need for this routine,
 *  but it makes it easier to 'generically' and 'agnostically' create
 *  remappers. If one hits the case of src_grid=tgt_grid, they can
 *  simply create an identity remapper. When remap methods are called,
 *  this remapper will
 *    - call Field::deep_copy if m_aliasing = NoAliasing
 *    - do nothing if m_aliasing != NoAliasing
 */

class IdentityRemapper : public AbstractRemapper
{
public:
  using base_type = AbstractRemapper;

  // If 
  enum Aliasing {
    SrcAliasTgt,
    TgtAliasSrc,
    NoAliasing
  };

  IdentityRemapper (const grid_ptr_type grid,
                    const Aliasing aliasing = NoAliasing)
   : base_type  (grid,grid)
  {
    set_aliasing(aliasing);
  }

  ~IdentityRemapper () = default;

  void set_aliasing (const Aliasing aliasing) {
    EKAT_REQUIRE_MSG (get_state()!=RepoState::Closed,
        "Error! Aliasing in IdentityRemapper must be set *before* registration ends.\n");
    m_aliasing = aliasing;
  }

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override {
    EKAT_REQUIRE_MSG (is_valid_tgt_layout(tgt_layout),
        "[IdentityRemapper] Error! Input target layout is not valid for this remapper.\n"
        " - input layout: " + to_string(tgt_layout));

    // Src and tgt grids are the same, so return the input
    return tgt_layout;
  }
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override {
    EKAT_REQUIRE_MSG (is_valid_src_layout(src_layout),
        "[IdentityRemapper] Error! Input source layout is not valid for this remapper.\n"
        " - input layout: " + to_string(src_layout));

    // Src and tgt grids are the same, so return the input
    return src_layout;
  }

  void register_field_from_src (const field_type& src) {
    EKAT_REQUIRE_MSG (m_aliasing!=SrcAliasTgt,
        "Error! Makes no sense to register from src and ask that src alias tgt.\n");
    if (m_aliasing==TgtAliasSrc) {
      register_field(src,src);
    } else {
      AbstractRemapper::register_field_from_src(src);
    }
  }
  void register_field_from_tgt (const field_type& tgt) {
    EKAT_REQUIRE_MSG (m_aliasing!=TgtAliasSrc,
        "Error! Makes no sense to register from tgt and ask that tgt alias src.\n");
    if (m_aliasing==SrcAliasTgt) {
      register_field(tgt,tgt);
    } else {
      AbstractRemapper::register_field_from_tgt(tgt);
    }
  }
protected:

  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_fields[ifield].first.get_header().get_identifier();
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_fields[ifield].second.get_header().get_identifier();
  }
  const field_type& do_get_src_field (const int ifield) const override {
    return m_fields[ifield].first;
  }
  const field_type& do_get_tgt_field (const int ifield) const override {
    return m_fields[ifield].second;
  }

  void do_registration_begins () override {
    // Reserve space for the m_fields_ids vector, in case the user set the number
    m_fields.reserve(this->m_num_fields);
  }
  void do_register_field (const identifier_type& src, const identifier_type& tgt) override {
    field_type src_f(src);
    field_type tgt_f(tgt);
    m_fields.emplace_back(src_f,tgt_f);
  }
  void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) override {
    m_fields[ifield].first  = src;
    m_fields[ifield].second = tgt;
  }
  void do_registration_ends () override {
    // If src is an alias of tgt (or viceversa), make the pair of fields store the same field
    if (m_aliasing==Aliasing::SrcAliasTgt) {
      for (auto& it : m_fields) {
        it.first = it.second;
      }
    } else if (m_aliasing==Aliasing::TgtAliasSrc) {
      for (auto& it : m_fields) {
        it.second = it.first;
      }
    }
  }

  void do_remap_fwd () override {
    // If src is an alias of tgt (or viceversa), no need to run the remapper.
    // Otherwise, we can simply run Field::deep_copy.
    if (m_aliasing==Aliasing::NoAliasing) {
      for (auto& it : m_fields) {
        it.second.deep_copy(it.first);
      }
    }
  }
  void do_remap_bwd () override {
    // If src is an alias of tgt (or viceversa), no need to run the remapper.
    // Otherwise, we can simply run Field::deep_copy.
    if (m_aliasing==Aliasing::NoAliasing) {
      for (auto& it : m_fields) {
        it.first.deep_copy(it.second);
      }
    }
  }

  std::vector<std::pair<field_type,field_type>>   m_fields;

  Aliasing m_aliasing;
};

} // namespace scream

#endif // SCREAM_IDENTITY_REMAPPER_HPP
