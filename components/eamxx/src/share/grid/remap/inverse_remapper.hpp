#ifndef SCREAM_INVERSE_REMAPPER_HPP
#define SCREAM_INVERSE_REMAPPER_HPP

#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

// Performs remap by chaining two remap strategies
class InverseRemapper : public AbstractRemapper
{
public:
  using base_type       = AbstractRemapper;

  InverseRemapper (std::shared_ptr<base_type> remapper) :
    base_type(remapper->get_tgt_grid(),remapper->get_src_grid())
  {
    ekat::error::runtime_check(static_cast<bool>(remapper), "Error! Null pointer for inner remapper.\n");

    m_remapper = remapper;
  }

  ~InverseRemapper () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override {
    return m_remapper->create_tgt_layout(tgt_layout);
  }
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override {
    return m_remapper->create_src_layout(src_layout);
  }

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const override {
    return m_remapper->compatible_layouts(tgt,src);
  }

protected:

  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_remapper->get_tgt_field_id(ifield);
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_remapper->get_src_field_id(ifield);
  }
  const field_type& do_get_src_field (const int ifield) const override {
    return m_remapper->get_tgt_field(ifield);
  }
  const field_type& do_get_tgt_field (const int ifield) const override {
    return m_remapper->get_src_field(ifield);
  }

  void do_remap_fwd () override {
    m_remapper->remap(false);
  }
  void do_remap_bwd () override {
    m_remapper->remap(true);
  }

  void do_registration_begins () override {
    m_remapper->registration_begins();
  }
  void do_register_field (const identifier_type& src, const identifier_type& tgt) override {
    m_remapper->register_field(tgt,src);
  }
  void do_bind_field (const int /* ifield */,
                      const field_type& src,
                      const field_type& tgt) override {
    m_remapper->bind_field(tgt,src);
  }
  void do_registration_ends () override {
    m_remapper->registration_ends();
  }

  std::shared_ptr<base_type>  m_remapper;
};

} // namespace scream

#endif // SCREAM_INVERSE_REMAPPER_HPP
