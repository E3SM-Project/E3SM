#ifndef SCREAM_COMBO_REMAPPER_HPP
#define SCREAM_COMBO_REMAPPER_HPP

#include "ekat/grid/remap/abstract_remapper.hpp"

namespace scream
{

/*
 *  A remapper representing the chaining of two remappers
 *  
 *  This remapper is rather simple: like InverseRemapper,
 *  it delegates all the work to the underlying remapper(s).
 *  The only additional duty of this remapper is to create
 *  the intermediate fields, corresponding to the target grid
 *  of the first remapper and the source grid of the second
 *  remapper.
 *  NOTE: this class is ASSUMING you won't be calling methods
 *        of the stored remappers once they have been wrapped
 *        in this class. If you do, then you will likely get
 *        errors. Instead, create the two remappers, use them
 *        to create this remapper, and then only interact with
 *        the combo remapper.
 */
template<typename ScalarType, typename DeviceType>
class ComboRemapper : public AbstractRemapper<ScalarType,DeviceType>
{
public:
  using base_type       = AbstractRemapper<ScalarType,DeviceType>;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;
  using layout_type     = typename base_type::layout_type;

  ComboRemapper (std::shared_ptr<base_type> r1, std::shared_ptr<base_type> r2) {
    error::runtime_check(static_cast<bool>(r1), "Error! Null pointer for first remapper.\n");
    error::runtime_check(static_cast<bool>(r2), "Error! Null pointer for second remapper.\n");
    error::runtime_check(r1->get_tgt_grid()->name()==r2->get_src_grid()->name(),
                         "Error! The input remappers are not compatible.\n");

    m_remapper_1 = r1;
    m_remapper_2 = r2;
  }

  ~ComboRemapper () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override {
    return m_remapper_1->create_src_layou(m_remapper_2->create_src_layout(tgt_layout));
  }
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override {
    return m_remapper_2->create_tgt_layou(m_remapper_1->create_tgt_layout(src_layout));
  }

protected:

  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_remapper_1->get_src_field_id(ifield);
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_remapper_2->get_tgt_field_id(ifield);
  }

  void do_remap_fwd () const override {
    m_remapper_1->remap(true);
    m_remapper_2->remap(true);
  }
  void do_remap_bwd () const override {
    m_remapper_2->remap(false);
    m_remapper_1->remap(false);
  }

  void do_registration_begins () override {
    // Reserve space for the fields, in case the user set the number
    m_remapper_1->set_num_fields(this->m_num_fields);
    m_remapper_2->set_num_fields(this->m_num_fields);

    m_remapper_1->registration_begins();
    m_remapper_2->registration_begins();
  }

  void do_register_field (const identifier_type& src, const identifier_type& tgt) override {
    layout_type tmp_layout = m_remapper_1->create_tgt_layout(src.get_layout());
    identifier_type tmp_id(src.name(),tmp_layout,m_remapper_1->get_tgt_grid()->name());

    m_remapper_1->register_field(src,tmp_id);
    m_remapper_2->register_field(tmp_id,tgt);

    m_tmps.emplace_back(field_type(tmp_id));
    m_tmps.back().alloate_view();
  }

  void do_bind_field (const int ifield,
                      const field_type& src,
                      const field_type& tgt) override {
    m_remapper_1->bind_field(ifield,src,m_tmps[ifield]);
    m_remapper_2->bind_field(ifield,m_tmps[ifield],tgt);
  }

  void do_registration_complete () override {
    m_remapper_1->registration_complete();
    m_remapper_2->registration_complete();
  }

  std::vector<field_type> m_tmps;

  std::shared_ptr<base_type>  m_remapper_1;
  std::shared_ptr<base_type>  m_remapper_2;
};

} // namespace scream

#endif // SCREAM_COMBO_REMAPPER_HPP
