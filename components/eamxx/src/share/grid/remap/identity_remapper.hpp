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
 *  remappers. If one hits the case of src_grid=tgt_grid, he/she can
 *  simply create an identity remapper. This remapper is guaranteed
 *  to do absolutely nothing (except for possibly some correctness
 *  checks, mostly through the base class interface) when the remap
 *  method is called. It *does* still store ids for source and target
 *  fields, mostly to allow queries via the base class interface.
 *  However, no field is actually stored (no views, that is), since
 *  there is no need to actually access the data.
 */

class IdentityRemapper : public AbstractRemapper
{
public:
  using base_type       = AbstractRemapper;

  IdentityRemapper (const grid_ptr_type grid)
   : base_type(grid,grid)
  {
    // Nothing to do here
  }

  ~IdentityRemapper () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override {
    // Src and tgt grids are the same, so return the input
    return tgt_layout;
  }
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override {
    // Src and tgt grids are the same, so return the input
    return src_layout;
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
    // Do nothing
  }

  void do_remap_fwd () override {
    // Do nothing
  }
  void do_remap_bwd () override {
    // Do nothing
  }

  std::vector<std::pair<field_type,field_type>>   m_fields;
};

} // namespace scream

#endif // SCREAM_IDENTITY_REMAPPER_HPP
