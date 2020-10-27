#ifndef SCREAM_DUMMY_SE_POINT_REMAPPER_HPP
#define SCREAM_DUMMY_SE_POINT_REMAPPER_HPP

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

template<typename RealType>
class DummySEPointRemapper : public AbstractRemapper<RealType>
{
public:
  using base_type       = AbstractRemapper<RealType>;
  using field_type      = typename base_type::field_type;
  using identifier_type = typename base_type::identifier_type;
  using layout_type     = typename base_type::layout_type;
  using grid_ptr_type   = typename base_type::grid_ptr_type;

  DummySEPointRemapper (const grid_ptr_type point_grid, const grid_ptr_type se_grid)
   : base_type(point_grid,se_grid)
  {
    EKAT_REQUIRE_MSG(se_grid->type()==GridType::SE,
                       "Error in DummySEPointRemapper! Invalid input se grid.\n");
    EKAT_REQUIRE_MSG(point_grid->type()==GridType::Point,
                       "Error in DummySEPointRemapper! Invalid input point grid.\n");
  }

  ~DummySEPointRemapper () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt) const override {
    using namespace ShortFieldTagsNames;
    auto ncol = this->get_src_grid()->get_native_dof_layout().dim(0);

    if (tgt.rank()==3) {
      FieldLayout src({COL},{ncol});
      return src;
    } else if (tgt.rank()==4) {
      FieldLayout src({COL,VL},{ncol,tgt.dim(3)});
      return src;
    } else if (tgt.rank()==5) {
      FieldLayout src({COL,CMP,VL},{ncol,tgt.dim(1),tgt.dim(4)});
      return src;
    } else {
      EKAT_ERROR_MSG ("Error! Target layout not supported. Remember that this class has limited support.\n");
    }
    FieldLayout src({},{});
    return src;
  }
  FieldLayout create_tgt_layout (const FieldLayout& src) const override {
    using namespace ShortFieldTagsNames;
    auto nele = this->get_tgt_grid()->get_native_dof_layout().dim(0);
    if (src.rank()==1) {
      FieldLayout tgt({EL,GP,GP},{nele,4,4});
      return tgt;
    } else if (src.rank()==2) {
      auto tag = src.tag(1);
      auto dim = src.dim(1);
      FieldLayout tgt({EL,GP,GP,tag},{nele,4,4,dim});
      return tgt;
    } else if (src.rank()==3) {
      FieldLayout tgt({EL,CMP,GP,GP,VL},{nele,4,4,src.dim(1),src.dim(2)});
      return tgt;
    } else {
      EKAT_ERROR_MSG ("Error! Target layout not supported. Remember that this class has limited support.\n");
    }
    FieldLayout tgt ({},{});
    return tgt;
  }

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const {
    return get_layout_type(src.tags())==get_layout_type(tgt.tags());
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
  void do_unregister_field (const int ifield) override {
    m_fields.erase(m_fields.begin()+ifield);
  }
  void do_registration_ends () override {
    // Do nothing
  }

  void do_remap_fwd () const override {
    // Do nothing
  }
  void do_remap_bwd () const override {
    // Do nothing
  }

  std::vector<std::pair<field_type,field_type>>   m_fields;
};

} // namespace scream

#endif // SCREAM_DUMMY_SE_POINT_REMAPPER_HPP
