#ifndef EAMXX_HEIGHT_LEVEL_INDEX_HPP
#define EAMXX_HEIGHT_LEVEL_INDEX_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic computes, for a given target height and a given
 * vertical layer (mid or int), the index of the first level whose height
 * is strictly below the target, in each column (nlevs if there is none).
 * Heights decrease with the level index, so this is the counterpart of
 * numpy's searchsorted (with side='right') for a descending array.
 * The output field has layout (ncol) and IntType data, with values in [0,nlevs].
 *
 * Whether a target height outside the column's height range is
 * considered an error depends on which end, and (for the bottom end)
 * on the surface reference:
 *  - Above the top entry (k==0) is always flagged as invalid (mask set to 0),
 *    regardless of the surface reference: "how far above the model top" is not
 *    something we can meaningfully extrapolate.
 *  - At or below the bottom entry (k==nlevs):
 *     - "above surface" heights are always well defined (every column
 *       has a surface), so the mask is 1, and consumers are expected to extrapolate
 *       using the closest (bottom) level. This is what makes something like
 *       wind_speed.interp(z_mid=2,reference='surface') well-defined even when
 *       the lowest model level is well above 2m.
 *     - "above sealevel" heights are NOT always well defined (e.g. a
 *       target elevation below a mountain's surface), so in that case
 *       the mask is set to 0, just like at the top.
 *
 * The output is masked: its valid mask is 1
 * where the target is valid (including when we extrapolate at the bottom),
 * and 0 where it is flagged as invalid. The mask only depends on the height
 * profile, so it has no dimension besides the columns, and it can be broadcast
 * to the layout of whichever field is interpolated.
 *
 * To linearly interpolate in a column where the mask is 1 and k<nlevs, use the
 * levels (k-1,k). Where the mask is 1 and k==nlevs, extrapolate with level nlevs-1.
 * Clamping k to [1,nlevs-1] gives levels one can always access without out-of-bounds
 * accesses (and ignore, where the mask is 0).
 *
 * This diagnostic exists so that multiple diagnostics/fields requesting
 * the same "value at height" (e.g. via FieldAtHeight) can share a single
 * search, rather than each one performing its own binary search, and can
 * share the same mask.
 */

class HeightLevelIndex : public AbstractDiagnostic
{
public:

  // Constructors
  HeightLevelIndex (const ekat::Comm& comm, const ekat::ParameterList& params,
                    const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "HeightLevelIndex"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl ();
protected:
  void initialize_impl ();

  std::string    m_diag_name;
  std::string    m_z_field_name;
  std::string    m_mask_name;

  // Whether a target below the bottom entry is valid, so consumers extrapolate
  // to the nearest (bottom) level (true, for "above surface"), or flagged
  // as invalid via the mask (false, for "above sealevel"). The top
  // entry never extrapolates, regardless of this flag. See class doc.
  bool           m_extrapolate_bottom;

  Real           m_z;
};

} //namespace scream

#endif // EAMXX_HEIGHT_LEVEL_INDEX_HPP
