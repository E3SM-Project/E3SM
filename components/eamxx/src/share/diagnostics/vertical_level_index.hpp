#ifndef EAMXX_VERTICAL_LEVEL_INDEX_HPP
#define EAMXX_VERTICAL_LEVEL_INDEX_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"
#include "share/diagnostics/level_index_utils.hpp"

namespace scream
{

/*
 * This diagnostic computes, for a given vertical target (a pressure, a height above
 * sea level, or a height above the surface) and a given vertical layer (mid or int),
 * the index of the first level that is past the target, in each column. This is numpy's
 * searchsorted (with side='right') along the level dimension, as std::upper_bound:
 *  - for the pressure, which increases with the level index, it is the first level with
 *    pressure strictly greater than the target;
 *  - for the heights, which decrease with the level index, it is the first level with
 *    height strictly below the target.
 * The output field has layout (ncol) and IntType data, with values in [0,nlevs].
 *
 * Whether a target outside the column's range is valid depends on which end, and
 * (for the bottom end) on the coordinate:
 *  - A target before the top entry (k==0) is always flagged as invalid, regardless
 *    of the coordinate: "how far above the model top" is not something we can
 *    meaningfully extrapolate.
 *  - A target beyond the bottom entry (k==nlevs):
 *     - for the height above the surface, it is valid, and consumers are expected to
 *       extrapolate using the closest (bottom) level. This is what makes something like
 *       wind_speed.interp(z_mid=2,reference='surface') well-defined even when the lowest
 *       model level is well above 2m. Every column has a surface, after all.
 *     - for the pressure and the height above sea level, it is flagged as invalid
 *       (e.g., a target elevation below a mountain's surface is not well defined).
 *
 * The output is masked: its valid mask is 1 where the target is valid, and 0 where it is
 * flagged as invalid. The index alone cannot tell, since nlevs is returned both when the
 * target is beyond the last level and when it is exactly equal to it. The mask only depends
 * on the vertical coordinate profile, so it has no dimension besides the columns, and it can
 * be broadcast to the layout of whichever field is interpolated.
 *
 * To linearly interpolate in a column where the mask is 1 and k<nlevs, use the levels
 * (k-1,k), with k clamped to [1,nlevs-1] (which also covers a target equal to the last entry).
 * Where the mask is 1 and k==nlevs, a column where we extrapolate, use level nlevs-1.
 * Where the mask is 0, the clamped levels are still valid, so one can access data without
 * out-of-bounds accesses, and then ignore it.
 *
 * This diagnostic exists so that multiple diagnostics/fields requesting the same "value at
 * vertical target" (e.g. via FieldVerticalInterp) can share a single search, rather than
 * each one performing its own binary search, and can share the same mask.
 */

class VerticalLevelIndex : public AbstractDiagnostic
{
public:

  // Constructors
  VerticalLevelIndex (const ekat::Comm& comm, const ekat::ParameterList& params,
                      const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "VerticalLevelIndex"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl ();
protected:
  void initialize_impl ();

  std::string    m_diag_name;
  std::string    m_mask_name;
  std::string    m_coord_name;

  // The vertical target, and what we need to know about its coordinate: if the coordinate
  // decreases with the level index, and if a target past the bottom entry is valid
  level_index_utils::VerticalTarget m_target;
};

} //namespace scream

#endif // EAMXX_VERTICAL_LEVEL_INDEX_HPP
