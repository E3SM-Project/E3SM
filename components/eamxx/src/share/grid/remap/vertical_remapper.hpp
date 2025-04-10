#ifndef EAMXX_VERTICAL_REMAPPER_HPP
#define EAMXX_VERTICAL_REMAPPER_HPP

#include "share/grid/remap/abstract_remapper.hpp"

#include <ekat/util/ekat_lin_interp.hpp>

namespace scream
{

/*
 * A remapper to interpolate fields on a separate vertical grid
 */

class VerticalRemapper : public AbstractRemapper
{
public:
  enum ExtrapType {
    Mask,  // Use fixed value
    P0     // Constant extrapolation
  };

  enum TopBot {
    Top = 1,
    Bot = 2,
    TopAndBot = Top | Bot
  };

  // When setting src/tgt pressure profiles, we may not care to distinguish
  // between midpoints/interface. E.g., in output, we may remap both midpoints
  // and interface quantities to the SAME set of pressure levels.
  enum ProfileType {
    Midpoints,
    Interfaces,
    Both
  };

  // Use fixed value as mask value
  VerticalRemapper (const grid_ptr_type& src_grid,
                    const std::string& map_file);

  VerticalRemapper (const grid_ptr_type& src_grid,
                    const grid_ptr_type& tgt_grid);

  ~VerticalRemapper () = default;

  void set_extrapolation_type (const ExtrapType etype, const TopBot where = TopAndBot);
  void set_mask_value (const Real mask_val);

  void set_source_pressure (const Field& p, const ProfileType ptype);
  void set_target_pressure (const Field& p, const ProfileType ptype);

  void set_source_pressure (const Field& pmid, const Field& pint) {
    set_source_pressure (pmid, Midpoints);
    set_source_pressure (pint, Interfaces);
  }
  void set_target_pressure (const Field& pmid, const Field& pint) {
    set_target_pressure (pmid, Midpoints);
    set_target_pressure (pint, Interfaces);
  }

  Field get_source_pressure (bool midpoints) const {
    return midpoints ? m_src_pmid : m_src_pint;
  }
  Field get_target_pressure (bool midpoints) const {
    return midpoints ? m_tgt_pmid : m_tgt_pint;
  }

  // This method simply creates the tgt grid from a map file
  static std::shared_ptr<AbstractGrid>
  create_tgt_grid (const grid_ptr_type& src_grid, const std::string& map_file);

  bool compatible_layouts (const FieldLayout& src, const FieldLayout& tgt) const override;

  bool is_valid_tgt_layout (const FieldLayout& layout) const override;
  bool is_valid_src_layout (const FieldLayout& layout) const override;
protected:

  void set_pressure (const Field& p, const std::string& src_or_tgt, const ProfileType ptype);
  FieldLayout create_layout (const FieldLayout& from_layout,
                             const std::shared_ptr<const AbstractGrid>& to_grid) const override;

  void registration_ends_impl () override;

  void remap_fwd_impl () override;

#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  template<int N>
  void apply_vertical_interpolation (const ekat::LinInterp<Real,N>& lin_interp,
                                     const Field& f_src, const Field& f_tgt,
                                     const Field& p_src, const Field& p_tgt) const;

  void extrapolate (const Field& f_src, const Field& f_tgt,
                    const Field& p_src, const Field& p_tgt,
                    const Real mask_val) const;

  template<int N>
  void setup_lin_interp (const ekat::LinInterp<Real,N>& lin_interp,
                         const Field& p_src, const Field& p_tgt) const;
protected:

  void create_lin_interp ();
  
  using KT = KokkosTypes<DefaultDevice>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;
  template<typename T>
  using view_2d = typename KT::template view_2d<T>;

  ekat::Comm            m_comm;

  // Source and target fields
  std::vector<Field>    m_src_masks;
  std::vector<Field>    m_tgt_masks;

  // Vertical profile fields, both for source and target
  Field                 m_src_pmid;
  Field                 m_src_pint;
  Field                 m_tgt_pmid;
  Field                 m_tgt_pint;

  bool m_src_int_same_as_mid = false;
  bool m_tgt_int_same_as_mid = false;

  // Extrapolation settings at top/bottom. Default to P0 extrapolation
  ExtrapType            m_etype_top = P0;
  ExtrapType            m_etype_bot = P0;
  Real                  m_mask_val = std::numeric_limits<Real>::quiet_NaN();

  // We need to remap mid/int fields separately, and we want to use packs if possible,
  // so we need to divide input fields into 4 separate categories

  // Map field id to whether it's packed/scalar and midpoint/interface
  struct FType {
    bool packed = false;
    bool midpoints = true;
  };
  std::map<std::string,FType> m_field2type;

  std::shared_ptr<ekat::LinInterp<Real,SCREAM_PACK_SIZE>> m_lin_interp_mid_packed;
  std::shared_ptr<ekat::LinInterp<Real,SCREAM_PACK_SIZE>> m_lin_interp_int_packed;
  std::shared_ptr<ekat::LinInterp<Real,1>>                m_lin_interp_mid_scalar;
  std::shared_ptr<ekat::LinInterp<Real,1>>                m_lin_interp_int_scalar;
};

} // namespace scream

#endif // EAMXX_VERTICAL_REMAPPER_HPP
