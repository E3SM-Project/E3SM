#ifndef EAMXX_VERTICAL_REMAPPER_HPP
#define EAMXX_VERTICAL_REMAPPER_HPP

#include "share/remap/abstract_remapper.hpp"

#include <ekat_lin_interp.hpp>

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

  enum InterpType {
    Linear,     // Standard linear interpolation
    LogLinear   // Log-linear interpolation (pressure coordinates are log-transformed)
  };

  VerticalRemapper (const grid_ptr_type& src_grid,
                    const std::string& map_file);

  VerticalRemapper (const grid_ptr_type& src_grid,
                    const grid_ptr_type& tgt_grid);

  ~VerticalRemapper () = default;

  void set_extrapolation_type (const ExtrapType etype, const TopBot where = TopAndBot);

  // Set the interpolation type. If pressure fields are already set, this method
  // will re-transform them: Linear->LogLinear applies log to stored pressures,
  // LogLinear->Linear applies exp to recover raw pressure values.
  void set_interp_type (const InterpType itype);

  // Declare that the source pressure passed to set_source_pressure is already
  // in log-space (i.e., values are log(Pa), not Pa). When true, no further log
  // transformation is applied to the source pressure at remap time. Useful when
  // the source grid has a logarithmic pressure coordinate.
  // Should be called before set_interp_type and set_source_pressure.
  void set_src_log_space (bool v = true);

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
  // Apply log element-wise in-place to f (rank 1 or 2). Caller must clone first if a copy is needed.
  void log_pressure (Field& f) const;
  // Apply exp element-wise in-place to f (rank 1 or 2). Caller must clone first if a copy is needed.
  void exp_pressure (Field& f) const;

  template<int N>
  void apply_vertical_interpolation (const ekat::LinInterp<Real,N>& lin_interp,
                                     const Field& f_src, const Field& f_tgt,
                                     const Field& p_src, const Field& p_tgt) const;

  void extrapolate (const Field& f_src, const Field& f_tgt,
                    const Field& p_src, const Field& p_tgt) const;

  template<int N>
  void setup_lin_interp (const ekat::LinInterp<Real,N>& lin_interp,
                         const Field& p_src, const Field& p_tgt) const;
protected:

  void create_lin_interp ();

  // Allocate and initialise m_src_pmid_log / m_src_pint_log from the current
  // raw source pressure fields.  Called whenever LogLinear mode is enabled and
  // m_src_log_space is false.
  void create_src_log_workspaces ();
  
  using KT = KokkosTypes<DefaultDevice>;

  template<typename T>
  using view_1d = typename KT::template view_1d<T>;
  template<typename T>
  using view_2d = typename KT::template view_2d<T>;

  ekat::Comm            m_comm;

  // Tgt grid masks (in case extrap type at top or bot is Mask)
  std::map<std::string,Field>    m_masks;

  // Vertical profile fields, both for source and target.
  // NOTE: m_src_pmid/m_src_pint ALWAYS hold the raw (non-log) pressure values,
  //       regardless of the interpolation type. Log transformation for the source
  //       is handled at remap time via the m_src_pmid_log/m_src_pint_log workspaces.
  Field                 m_src_pmid;
  Field                 m_src_pint;
  Field                 m_tgt_pmid;
  Field                 m_tgt_pint;

  // For log-linear mode: pre-allocated workspaces that hold log(p_src) at
  // remap time. The source pressure fields (p_mid, p_int) change every model
  // timestep, so we cannot freeze a log-clone at construction time.  Instead
  // we keep the raw shared reference in m_src_pmid / m_src_pint and refresh
  // these workspaces (copy raw values then apply log) at the start of every
  // remap_fwd_impl call. The initial values set during construction are
  // overwritten at the first remap call.
  Field                 m_src_pmid_log;
  Field                 m_src_pint_log;

  // If user provides pressure profiles that are NOT compatible with SCREAM_PACK_SIZE,
  // we will set these booleans to false, and use ONLY the "scalar" LinInterp structures
  bool m_int_packs_supported = true;
  bool m_mid_packs_supported = true;

  // Extrapolation settings at top/bottom. Default to P0 extrapolation
  ExtrapType            m_etype_top = P0;
  ExtrapType            m_etype_bot = P0;

  // Interpolation type: linear or log-linear
  InterpType            m_interp_type = Linear;

  // When true, the source pressure provided by the user is already in log-space
  // (log(Pa)), and no further log transformation is applied in remap_fwd_impl.
  // This supports source grids whose vertical coordinate is logarithmic pressure.
  bool                  m_src_log_space = false;

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
