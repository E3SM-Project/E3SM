#ifndef SCREAM_SURFACE_COUPLING_HPP
#define SCREAM_SURFACE_COUPLING_HPP

#include "share/field/field.hpp"
#include "share/field/field_repository.hpp"
#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"

#include "ekat/kokkos/ekat_kokkos_meta.hpp"

#include <set>

namespace scream {
namespace control {

// This class is responsible to import/export fields from/to the rest of
// E3SM, via the mct coupler.
class SurfaceCoupling {
public:

  using host_device_type = HostDevice;
  using device_type      = DefaultDevice;

  template<typename DevT, typename DataT>
  using view_1d = typename KokkosTypes<DevT>::template view_1d<DataT>;
  template<typename DevT, typename DataT>
  using view_2d = typename KokkosTypes<DevT>::template view_2d<DataT>;


  // Use const type for exports, to protect against bugs
  using import_field_type  = Field<      Real>;
  using export_field_type  = Field<const Real>;

  // MCT always passes data as double. Also, we can't enforce const correctness there
  using cpl_data_ptr_type = double*;

  // The input grid is the one where import/export happens
  SurfaceCoupling (const std::shared_ptr<const AbstractGrid>& grid,
                   const FieldRepository<Real>& repo);

  // This allocates some service views
  void set_num_fields (const int num_imports, const int num_exports);

  // Register import/export scream fields
  void register_import (const std::string& fname,
                        const int cpl_idx,
                        const int vecComp = -1);
  void register_export (const std::string& fname,
                        const int cpl_idx,
                        const int vecComp = -1);

  // Marks the end of the registration phase. Here, we check that there are no
  // import/export Info struct with only partial information
  void registration_ends (cpl_data_ptr_type cpl_imports_2d_array_ptr,
                          cpl_data_ptr_type cpl_exports_2d_array_ptr);

  // Import host fields from the component coupler to device fields in the AD
  void do_import ();

  // Export device fields from the AD to host fields in the component coupler
  void do_export ();

  // Getters
  RepoState get_repo_state () const { return m_state; }

  const std::set<FieldIdentifier>& get_import_fids () const { return m_imports_fids; }
  const std::set<FieldIdentifier>& get_export_fids () const { return m_exports_fids; }

#ifndef KOKKOS_ENABLE_CUDA
protected:
#endif

  // A device-friendly helper struct, storing information about the import/export.
  // Templated on the value type of the scream field associated with the import/export.
  // This can be either 'Real' or 'const Real'
  template<typename ValueType>
  struct Info {
    // Set to invalid, for ease of checking
    KOKKOS_INLINE_FUNCTION
    Info () : data(nullptr) {}

    KOKKOS_INLINE_FUNCTION
    Info& operator= (const Info&) = default;

    // Index of the field inside the coupler 2d array
    int cpl_idx;

    // Stride between the 1st entry of two consecutive columns to be imported/exported.
    // Note: this is >= that number of scalars in a column. E.g., for a vector field layout like
    //       (ncols,2,nlevs), where we import/export only the 1st vector component, the stride
    //       is 2*nlevs
    int col_stride;

    // Offset to surface field from the column start. Should be 0 for scalar fields, but
    // may be non-zero for vector quantities for which we import/export the 2nd (or larger)
    // component (the layout would be something like (num_cols,2,num_levs), so the 1st
    // entry to import export would be at index num_levs.
    int col_offset;

    // Pointer to the scream field device memory
    ValueType*  data;
  };

#ifdef KOKKOS_ENABLE_CUDA
protected:
#endif

  using import_value_type  = import_field_type::value_type;
  using export_value_type  = export_field_type::value_type;

  // Packed, device-friendly verson of the helper structures above
  view_1d<device_type,Info<import_value_type>>  m_scream_imports_dev;
  view_1d<device_type,Info<export_value_type>>  m_scream_exports_dev;
  decltype(m_scream_imports_dev)::HostMirror    m_scream_imports_host;
  decltype(m_scream_exports_dev)::HostMirror    m_scream_exports_host;

  // The type for device and host view storing a 2d array with dims (num_cols,num_fields).
  // The field idx strides faster, since that's what mct does (so we can "view" the
  // pointer to the whole x2a and a2x arrays from Fortran)
  using cpl_dview_type = view_2d<device_type,double>;
  using cpl_hview_type = ekat::Unmanaged<view_2d<host_device_type,double>>;

  cpl_dview_type       m_cpl_imports_view_d;
  cpl_hview_type       m_cpl_imports_view_h;

  cpl_dview_type       m_cpl_exports_view_d;
  cpl_hview_type       m_cpl_exports_view_h;

  // The following is only stored for debug/inspection routines
  std::set<FieldIdentifier>  m_imports_fids;
  std::set<FieldIdentifier>  m_exports_fids;

  const FieldRepository<Real>& m_field_repo;

  int           m_num_imports;
  int           m_num_exports;

  int           m_num_cols;
  std::string   m_grid_name;

  RepoState     m_state;
};

} // namespace control
} // namespace scream

#endif // SCREAM_SURFACE_COUPLING_HPP
