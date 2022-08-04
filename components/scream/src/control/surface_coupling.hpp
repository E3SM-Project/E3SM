#ifndef SCREAM_SURFACE_COUPLING_HPP
#define SCREAM_SURFACE_COUPLING_HPP

#include "share/field/field.hpp"
#include "share/field/field_manager.hpp"
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

  using field_mgr_type = const FieldManager;
  using field_mgr_ptr  = std::shared_ptr<field_mgr_type>;
  using grid_type       = const AbstractGrid;
  using grid_ptr        = std::shared_ptr<grid_type>;

  // MCT always passes data as double. Also, we can't enforce const correctness there
  using cpl_data_ptr_type = double*;

  // The input grid is the one where import/export happens
  explicit SurfaceCoupling (const field_mgr_ptr& field_mgr);

  // This allocates some service views. Since not all imported
  // data is used by SCREAM, we distinguish between cpl
  // imports and SCREAM imports for book keeping.
  void set_num_fields (const int num_cpl_imports, const int num_scream_imports,
                       const int num_cpl_exports);
  // Version of the above function when num_cpl_imports = num_scream_imports
  void set_num_fields (const int num_imports, const int num_exports)
  { set_num_fields(num_imports, num_imports, num_exports); }

  // Register import scream fields.
  void register_import (const std::string& fname,
                        const int cpl_idx,
                        const int vecComp = -1);

  // Register export scream fields. Since do_export() can be called during initialization,
  // some scream fields may not have valid entries (i.e., computed fields with no IC).
  // For these fields, set export_during_init=false and they will be skipped.
  void register_export (const std::string& fname,
                        const int cpl_idx,
                        const int vecComp = -1,
                        const bool export_during_init = true);

  // Marks the end of the registration phase. Here, we check that there are no
  // import/export Info struct with only partial information
  void registration_ends (cpl_data_ptr_type cpl_imports_2d_array_ptr,
                          cpl_data_ptr_type cpl_exports_2d_array_ptr);

  // Import host fields from the component coupler to device fields in the AD
  void do_import ();

  // Export device fields from the AD to host fields in the component coupler.
  // If this export is called during the init phase, set init_phase=true
  // so that fields which are computed inside SCREAM during the run phase are skipped.
  void do_export (const int dt, const bool init_phase = false);

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

    // Boolean that dictates if the field can be exported if do_export() is called during
    // initialization. This is useful since inside SCREAM some fields require computation
    // internally.
    bool do_initial_export;

    // Pointer to the scream field device memory
    Real*  data;
  };

#ifdef KOKKOS_ENABLE_CUDA
protected:
#endif

  void get_col_info (const std::shared_ptr<const FieldHeader>& fh,
                     int vecComp, int& col_offset, int& col_stride) const;

  using value_type = Real;

  // Packed, device-friendly verson of the helper structures above
  view_1d<device_type,Info>  m_scream_imports_dev;
  view_1d<device_type,Info>  m_scream_exports_dev;
  decltype(m_scream_imports_dev)::HostMirror    m_scream_imports_host;
  decltype(m_scream_exports_dev)::HostMirror    m_scream_exports_host;

  // Views needed for export computations
  view_2d<device_type,Real> dz, z_int, z_mid;

  // Views for storing export values for various fields that need to be computed
  view_1d<device_type,Real> Sa_z;
  view_1d<device_type,Real> Sa_ptem;
  view_1d<device_type,Real> Sa_dens;
  view_1d<device_type,Real> Sa_pslv;
  view_1d<device_type,Real> Faxa_rainl;
  view_1d<device_type,Real> Faxa_snowl;
  view_1d<device_type,Real> zero_view;

  // Dummy field to allow for the storage of the FieldIdentifier for debugging (not currently used)
  // and the FieldHeader needed to set up column info for export fields which must be computed.
  Field dummy_field;

  // The type for device and host view storing a 2d array with dims (num_cols,num_fields).
  // The field idx strides faster, since that's what mct does (so we can "view" the
  // pointer to the whole x2a and a2x arrays from Fortran)
  using cpl_dview_type = view_2d<device_type,double>;
  using cpl_hview_type = ekat::Unmanaged<view_2d<host_device_type,double>>;

  cpl_dview_type       m_cpl_imports_view_d;
  cpl_hview_type       m_cpl_imports_view_h;

  cpl_dview_type       m_cpl_exports_view_d;
  cpl_hview_type       m_cpl_exports_view_h;

  // Views which store 1 for fields whose convention with regard to sign is the same between cpl
  // and SCREAM, and -1 for fields that differ.
  // Currently, this is needed for fluxes as surface at atmosphere have different interpretations
  // of the positive direction.
  // TODO: This should be imporved in the future. Currently, if conventions change wrt. the flux
  //       direction in either SCREAM or cpl, a bug related to the sign could be hard
  //       to track down.
  view_1d<device_type,      int>  m_cpl_scream_sign_change_dev;
  view_1d<host_device_type, int>  m_cpl_scream_sign_change_host;

  // The following is only stored for debug/inspection routines
  std::set<FieldIdentifier>  m_imports_fids;
  std::set<FieldIdentifier>  m_exports_fids;

  field_mgr_ptr m_field_mgr;

  int           m_num_cpl_imports;
  int           m_num_scream_imports;
  int           m_num_cpl_exports;
  int           m_num_scream_exports;

  int           m_num_cols;
  int           m_num_levs;

  RepoState     m_state;
};

} // namespace control
} // namespace scream

#endif // SCREAM_SURFACE_COUPLING_HPP
