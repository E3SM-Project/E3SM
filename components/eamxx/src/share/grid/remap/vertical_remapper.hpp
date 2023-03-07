#ifndef EAMXX_VERTICAL_REMAPPER_HPP
#define EAMXX_VERTICAL_REMAPPER_HPP

#include "share/field/field_tag.hpp"
#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

/*
 * A remapper to interpolate fields on a separate vertical grid
 */

class VerticalRemapper : public AbstractRemapper
{
public:

  VerticalRemapper (const grid_ptr_type& src_grid,
                    const std::string& map_file,
                    const Field& lev_prof,
                    const Field& ilev_prof,
                    const Real mask_val);

  VerticalRemapper (const grid_ptr_type& src_grid,
                    const std::string& map_file,
                    const Field& lev_prof,
                    const Field& ilev_prof);

  ~VerticalRemapper () = default;

  FieldLayout create_src_layout (const FieldLayout& tgt_layout) const override;
  FieldLayout create_tgt_layout (const FieldLayout& src_layout) const override;

  bool compatible_layouts (const layout_type& src,
                           const layout_type& tgt) const override {
    // Same type of layout, and same sizes except for possibly the first one
    // Note: we can't do tgt.size()/tgt.dim(0), since there may be 0 tgt gids
    //       on some ranks, which means tgt.dim(0)=0.
    // Note: for vertical remapping we strip out the LEV or ILEV dimension when
    //       calculating the size.
    auto src_dims = src.dims();
    auto tgt_dims = tgt.dims();
    auto src_size = src.rank();
    auto tgt_size = tgt.rank();

    using namespace ShortFieldTagsNames;
    if (src.has_tag(LEV) || src.has_tag(ILEV)) {
      // Then we ignore the last dimension:
      src_size -= 1; 
    } 
    if (tgt.has_tag(LEV) || tgt.has_tag(ILEV)) {
      // Then we ignore the last dimension:
      tgt_size -= 1; 
    } 

    int src_col_size = 1;
    for (int i=0; i<src_size; ++i) {
      src_col_size *= src_dims[i];
    }
    int tgt_col_size = 1;
    for (int i=0; i<tgt_size; ++i) {
      tgt_col_size *= tgt_dims[i];
    }
    return get_layout_type(src.tags())==get_layout_type(tgt.tags()) &&
           src_col_size == tgt_col_size;
  }


protected:

  void register_vertical_source_field(const Field& src, const std::string& mode);

  const identifier_type& do_get_src_field_id (const int ifield) const override {
    return m_src_fields[ifield].get_header().get_identifier();
  }
  const identifier_type& do_get_tgt_field_id (const int ifield) const override {
    return m_tgt_fields[ifield].get_header().get_identifier();
  }
  const field_type& do_get_src_field (const int ifield) const override {
    return m_src_fields[ifield];
  }
  const field_type& do_get_tgt_field (const int ifield) const override {
    return m_tgt_fields[ifield];
  }

  void do_registration_begins () override { /* Nothing to do here */ }

  void do_register_field (const identifier_type& src, const identifier_type& tgt) override;

  void do_bind_field (const int ifield, const field_type& src, const field_type& tgt) override;

  void do_registration_ends () override;

  void do_remap_fwd () override;

  void do_remap_bwd () override {
    EKAT_ERROR_MSG ("VerticalRemapper only supports fwd remapping.\n");
  }

  void set_pressure_levels (const std::string& map_file);
  void do_print();

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int N>
  void apply_vertical_interpolation (const Field& f_src, const Field& f_tgt, const bool mask_interp=false) const;
protected:

  using KT = KokkosTypes<DefaultDevice>;
  using gid_t = AbstractGrid::gid_type;

  template<int N>
  using RPack = ekat::Pack<Real,N>;

  using mPack = RPack<SCREAM_PACK_SIZE>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;
  template<typename T>
  using view_2d = typename KT::template view_2d<T>;

  ekat::Comm            m_comm;

  // Source and target fields
  std::vector<Field>    m_src_fields;
  std::vector<Field>    m_tgt_fields;
  std::vector<Field>    m_tgt_masks;
  std::vector<Field>    m_src_masks;

  // Vertical profile fields, both for source and target
  int                   m_num_remap_levs;
  Real                  m_mask_val;
  Field                 m_remap_pres;
  Field                 m_src_mid;  // Src vertical profile for LEV layouts
  Field                 m_src_int;  // Src vertical profile for ILEV layouts
  bool                  m_mid_set = false;
  bool                  m_int_set = false;
};

} // namespace scream

#endif // EAMXX_VERTICAL_REMAPPER_HPP
