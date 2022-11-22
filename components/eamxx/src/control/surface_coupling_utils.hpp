#ifndef SCREAM_SURFACE_COUPLING_UTILS_HPP
#define SCREAM_SURFACE_COUPLING_UTILS_HPP

#include "share/scream_types.hpp"
#include "share/field/field.hpp"

namespace scream {

// Enum for distiguishing between an import or export 
enum class SurfaceCouplingTransferType {
  Import,
  Export
};

// A device-friendly helper struct, storing column information about the import/export.
struct SurfaceCouplingColumnInfo {
  // Set to invalid, for ease of checking
  KOKKOS_INLINE_FUNCTION
  SurfaceCouplingColumnInfo () : data(nullptr) {}

  KOKKOS_INLINE_FUNCTION
  SurfaceCouplingColumnInfo& operator= (const SurfaceCouplingColumnInfo&) = default;

  // Index of this column in cpl data.
  int cpl_indx;
 
  // Stride between the 1st entry of two consecutive columns to be imported.
  // Note: this is >= that number of scalars in a column. E.g., for a vector field layout like
  //       (ncols,2,nlevs), where we import only the 1st vector component, the stride
  //       is 2*nlevs
  int col_stride;

  // Offset to surface field from the column start. Should be 0 for scalar fields, but
  // may be non-zero for vector quantities for which we import the 2nd (or larger)
  // component (the layout would be something like (num_cols,2,num_levs), so the 1st
  // entry to import would be at index num_levs.
  int col_offset;

  // Constant multple applied to import data.
  // An example where this is useful is for fluxes, where cpl and SCREAM have different interpretations
  // of the positive direction. In this case, the constant multiple should be set to -1.
  Real constant_multiple;

  // Boolean that dictates if the field can be imported/exported if do_import/export() is called during
  // initialization. An example of when this useful inside SCREAM is that some exported fields require computation
  // done in run_impl() of various processes, and therefore do not have valid entries during the do_export() call
  // in initialization.
  bool transfer_during_initialization;

  // Pointer to the scream field device memory
  Real* data;
};

// For a given field and vector component (set vecComp=-1 for scalar fields),
// this function calculates the col_offset and col_stride for iterating through
// surface values. In this case, col_offset is the distance to the surface
// value within a mesh point (column) i, and col_stride is the distance in between
// surface values of adjacent columns. To access the surface value on column i,
// simply calculate surface_value_index=col_offset + i*col_stride.
//
// Some examples:
// 2D scalar field (ncols):
//     col_offset=0,                                     col_stride=1
// 2D vector field (ncols, num_comp):
//     col_offset=vec_comp,                              col_stride=num_comp
// 3D scalar field (ncols, nlevs):
//     col_offset=nlevs-1,                               col_stride=nlevs_with_padding
// 3D vector field (ncols, num_comp, nlevs):
//     col_offset=vec_comp*nlevs_with_padding + nlevs-1, col_stride=nlevs_with_padding*num_comp
//
// Field with parents are allowed if the parent has layout 3D vector (ncols, num_comp, nlevs),
// and the field was sliced along num_comp. See comments below for details.
void get_col_info_for_surface_values(const std::shared_ptr<const FieldHeader>& fh,
                                     int vecComp, int& col_offset, int& col_stride);

} // namespace scream

#endif // SCREAM_SURFACE_COUPLING_UTILS_HPP
