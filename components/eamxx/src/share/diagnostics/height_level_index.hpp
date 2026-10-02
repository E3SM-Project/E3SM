#ifndef EAMXX_HEIGHT_LEVEL_INDEX_HPP
#define EAMXX_HEIGHT_LEVEL_INDEX_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream
{

/*
 * This diagnostic computes, for a given target height and a given
 * vertical layer (mid or int), the pair of level indices that bracket
 * the target height in each column, ready to be used for a linear
 * interpolation. The output field has layout (ncol,2) and IntType data.
 *
 * Whether a target height outside the column's height range is
 * considered an error depends on which end, and (for the bottom end)
 * on the surface reference:
 *  - Above the top entry (i.e. above the model top) is always flagged
 *    as invalid (both indices set to -1), regardless of the surface
 *    reference: this matches PressureLevelIndex's behavior at the
 *    model top, where "how far above the model top" is not something
 *    we can meaningfully extrapolate.
 *  - Below the bottom entry:
 *     - "above surface" heights are always well defined (every column
 *       has a surface), so we extrapolate using the closest (bottom)
 *       level, and both indices are set equal to nlevs-1. This is what
 *       makes something like wind_speed_at_2m_above_surface
 *       well-defined even when the lowest model level is well above 2m.
 *     - "above sealevel" heights are NOT always well defined (e.g. a
 *       target elevation below a mountain's surface), so in that case
 *       both indices are set to -1, just like at the top.
 *
 * This diagnostic exists so that multiple diagnostics/fields requesting
 * the same "value at height" (e.g. via FieldAtHeight) can share a single
 * bracket-index computation, rather than each one performing its own
 * binary search.
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

  // Whether a target below the bottom entry should extrapolate to the
  // nearest (bottom) level (true, for "above surface"), or be flagged
  // as invalid via a -1 index (false, for "above sealevel"). The top
  // entry never extrapolates, regardless of this flag. See class doc.
  bool           m_extrapolate_bottom;

  Real           m_z;
};

} //namespace scream

#endif // EAMXX_HEIGHT_LEVEL_INDEX_HPP
