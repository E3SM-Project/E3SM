#ifndef SCREAM_INVERSE_REMAPPER_HPP
#define SCREAM_INVERSE_REMAPPER_HPP

#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

// Performs remap by chaining two remap strategies
class InverseRemapper : public AbstractRemapper
{
public:
  InverseRemapper (std::shared_ptr<AbstractRemapper> remapper)
  {
    EKAT_REQUIRE_MSG (remapper!=nullptr, "Error! Null pointer for inner remapper.\n");

    set_grids(remapper->get_tgt_grid(),remapper->get_src_grid());

    m_remapper = remapper;

    m_fwd_allowed = m_remapper->bwd_allowed();
    m_bwd_allowed = m_remapper->fwd_allowed();
  }

  ~InverseRemapper () = default;

  // Target and src are flipped
  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override {
    return m_remapper->create_tgt_layout(tgt_layout);
  }
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override {
    return m_remapper->create_src_layout(src_layout);
  }

  bool compatible_layouts (const FieldLayout& src,
                           const FieldLayout& tgt) const override {
    return m_remapper->compatible_layouts(tgt,src);
  }
  virtual bool is_valid_src_layout (const FieldLayout& layout) const {
    return m_remapper->is_valid_tgt_layout(layout);
  }
  virtual bool is_valid_tgt_layout (const FieldLayout& layout) const {
    return m_remapper->is_valid_src_layout(layout);
  }

protected:

  void remap_fwd_impl () override {
    m_remapper->remap_bwd();
  }
  void remap_bwd_impl () override {
    m_remapper->remap_fwd();
  }

  void registration_ends_impl () override {
    m_remapper->registration_begins();
    for (int i=0; i<m_num_fields; ++i) {
      m_remapper->register_field(m_tgt_fields[i],m_src_fields[i]);
    }
    m_remapper->registration_ends();
  }

  std::shared_ptr<AbstractRemapper>  m_remapper;
};

} // namespace scream

#endif // SCREAM_INVERSE_REMAPPER_HPP
