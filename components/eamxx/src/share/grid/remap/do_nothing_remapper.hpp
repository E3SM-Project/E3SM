#ifndef SCREAM_DO_NOTHING_REMAPPER_HPP
#define SCREAM_DO_NOTHING_REMAPPER_HPP

#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

/*
 *  A placeholder remapper for SE grids
 *
 *  This class is NOT meant to be used outside of unit testing.
 *  This 'remapper' doesn't actually remap anything, and it should only
 *  be used in certain unit tests that require a remapper to be set,
 *  with its status query-able.
 *  You are allowed to use this remapper for things like registering
 *  fields, get number of registered fields and their ids.
 *  You can even ask to create a tgt layout from a src one, but the
 *  result is purely a placeholder result. Do not rely on any
 *  data consistency.
 *  Every call to the actual remap methods doesn't do anything.
 */

class DoNothingRemapper : public AbstractRemapper
{
public:
  using base_type       = AbstractRemapper;

  DoNothingRemapper (const grid_ptr_type from, const grid_ptr_type to)
   : base_type(from,to)
  {
    // Nothing to do here
  }

  ~DoNothingRemapper () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt) const override {
    using namespace ShortFieldTagsNames;

    auto type = get_layout_type(tgt.tags());
    FieldLayout src = {{}};
    switch (type) {
      case LayoutType::Scalar2D:
        src = this->m_src_grid->get_2d_scalar_layout();
        break;
      case LayoutType::Vector2D:
        src = this->m_src_grid->get_2d_vector_layout(CMP,tgt.dim(CMP));
        break;
      case LayoutType::Scalar3D:
        src = this->m_src_grid->get_3d_scalar_layout(tgt.has_tag(LEV));
        break;
      case LayoutType::Vector3D:
        src = this->m_src_grid->get_3d_vector_layout(tgt.has_tag(LEV),CMP,tgt.dim(CMP));
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field layout.\n");
    }
    return src;
  }

  FieldLayout create_tgt_layout (const FieldLayout& src) const override {
    using namespace ShortFieldTagsNames;

    auto type = get_layout_type(src.tags());
    FieldLayout tgt = {{}};
    switch (type) {
      case LayoutType::Scalar2D:
        tgt = this->m_tgt_grid->get_2d_scalar_layout();
        break;
      case LayoutType::Vector2D:
        tgt = this->m_tgt_grid->get_2d_vector_layout(CMP,tgt.dim(CMP));
        break;
      case LayoutType::Scalar3D:
        tgt = this->m_tgt_grid->get_3d_scalar_layout(tgt.has_tag(LEV));
        break;
      case LayoutType::Vector3D:
        tgt = this->m_tgt_grid->get_3d_vector_layout(tgt.has_tag(LEV),CMP,tgt.dim(CMP));
        break;
      default:
        EKAT_ERROR_MSG ("Error! Unsupported field layout.\n");
    }
    return tgt;
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

  // We need to store fields in case the remapper is queried for them
  std::vector<std::pair<field_type,field_type>>   m_fields;
};

} // namespace scream

#endif // SCREAM_DO_NOTHING_REMAPPER_HPP
