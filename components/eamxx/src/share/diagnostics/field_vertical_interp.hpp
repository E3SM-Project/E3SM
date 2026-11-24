#ifndef EAMXX_FIELD_VERTICAL_INTERP_HPP
#define EAMXX_FIELD_VERTICAL_INTERP_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"
#include "share/diagnostics/level_index_utils.hpp"

namespace scream
{

/*
 * This diagnostic will produce a slice of a field at a given vertical target, by linearly
 * interpolating the field, along the vertical dimension, in the vertical coordinate:
 *  - a pressure,
 *  - a height above sea level, or
 *  - a height above the surface.
 *
 * The level index of the target, and the mask of the columns where the target is valid,
 * are computed by VerticalLevelIndex (see there for the details), which can be shared with
 * all the diagnostics interpolating at the same target.
 *
 * Parameters:
 *  - field_name: the field to interpolate. It must have layout (ncol,[dim,]lev), with lev
 *    either LEV or ILEV.
 *  - vertical_coordinate: 'p' (pressure), 'z' (height above sea level) or 'height'
 *    (height above surface)
 *  - value: the value of the target
 *  - units: the units of the value ('mb', 'hPa' or 'Pa' for 'p', and 'm' otherwise)
 */

class FieldVerticalInterp : public AbstractDiagnostic
{
public:

  // Constructors
  FieldVerticalInterp (const ekat::Comm& comm, const ekat::ParameterList& params,
                       const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name () const { return "FieldVerticalInterp"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl ();
protected:
  void initialize_impl ();

  std::string    m_diag_name;
  std::string    m_field_name;

  // The vertical target. It tells the names of the (mid/int) coordinate and level-index fields
  // we depend on, and whether a valid column past the bottom entry extrapolates (see
  // VerticalLevelIndex).
  level_index_utils::VerticalTarget m_target;

  // The coordinate and level-index fields we use. They are the mid or int ones, which we only
  // know in initialize_impl, once we know whether the input field lives on LEV or ILEV.
  std::string    m_coord_name;
  std::string    m_index_name;
};

} //namespace scream

#endif // EAMXX_FIELD_VERTICAL_INTERP_HPP
