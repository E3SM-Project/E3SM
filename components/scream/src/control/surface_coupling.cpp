#include "control/surface_coupling.hpp"

#include "share/field/field_utils.hpp"

namespace scream {
namespace control {

SurfaceCoupling::
SurfaceCoupling (const field_mgr_ptr& field_mgr)
 : m_field_mgr (field_mgr)
 , m_state (RepoState::Clean)
{
  auto grid = m_field_mgr->get_grid();

  m_num_cols = grid->get_num_local_dofs();
  m_num_levs = grid->get_num_vertical_levels();

  using namespace ShortFieldTagsNames;
  using namespace ekat::units;

  FieldLayout layout { {COL}, {m_num_cols} };
  const auto nondim = Units::nondimensional();
  auto grid_name = grid->name();
  FieldIdentifier id ("dummy_field", layout, nondim, grid_name);
  dummy_field = Field(id);

  EKAT_REQUIRE_MSG(grid->type()==GridType::Point,
      "Error! Surface coupling only implemented for 'Point' grids.\n"
      "       Input grid type: " + e2str(grid->type()) + "\n");
}

void SurfaceCoupling::
set_num_fields (const int num_cpl_imports, const int num_scream_imports,
                const int num_cpl_exports)
{
  EKAT_REQUIRE_MSG(m_state==RepoState::Clean,
                     "Error! You can only set the number of fields in SurfaceCoupling once.\n");

  m_scream_imports_dev = decltype(m_scream_imports_dev)("",num_scream_imports);
  m_scream_exports_dev = decltype(m_scream_exports_dev)("",num_cpl_exports);

  m_scream_imports_host = Kokkos::create_mirror_view(m_scream_imports_dev);
  m_scream_exports_host = Kokkos::create_mirror_view(m_scream_exports_dev);

  m_cpl_scream_sign_change_dev  = decltype(m_cpl_scream_sign_change_dev)("",num_scream_imports);
  m_cpl_scream_sign_change_host = Kokkos::create_mirror_view(m_cpl_scream_sign_change_dev);

  // These fields contain computation needed for some export fields
  dz    = decltype(dz)    ("", m_num_cols, m_num_levs);
  z_int = decltype(z_int) ("", m_num_cols, m_num_levs+1);
  z_mid = decltype(z_mid) ("", m_num_cols, m_num_levs);

  // These fields contain export data
  Sa_z       = decltype(Sa_z)       ("", m_num_cols);
  Sa_ptem    = decltype(Sa_ptem)    ("", m_num_cols);
  Sa_dens    = decltype(Sa_dens)    ("", m_num_cols);
  Sa_pslv    = decltype(Sa_pslv)    ("", m_num_cols);
  Faxa_rainl = decltype(Faxa_rainl) ("", m_num_cols);
  Faxa_snowl = decltype(Faxa_snowl) ("", m_num_cols);
  zero_view  = decltype(zero_view)  ("", m_num_cols);
  Kokkos::deep_copy(zero_view, 0.0);

  // These will be incremented every time we register an import/export.
  // At registration end, they should match the length of the corresponding view above.
  m_num_scream_imports = 0;
  m_num_scream_exports = 0;

  // We keep up with the total number of cpl imports/exports for indexing input data
  m_num_cpl_imports = num_cpl_imports;
  m_num_cpl_exports = num_cpl_exports;

  m_state = RepoState::Open;
}

void SurfaceCoupling::
register_import(const std::string& fname,
                const int cpl_idx,
                const int vecComp)
{
  // Two separate checks rather than state==Open, so we can print more specific error messages
  EKAT_REQUIRE_MSG (m_state!=RepoState::Clean,
                      "Error! Registration phase has not started yet.\n"
                      "       Did you forget to call set_num_fields(..) ?\n");
  EKAT_REQUIRE_MSG (m_state!=RepoState::Closed, "Error! Registration phase has already ended.\n");

  if (fname == "unused")
  {
    // Do nothing
  } else {
    // Check that we still have room
    EKAT_REQUIRE_MSG(m_num_scream_imports<m_scream_imports_host.extent_int(0),
                       "Error! Imports view is already full. Did you call 'set_num_fields' with the wrong arguments?\n");

    // Get the field, and check that is valid
    Field field = m_field_mgr->get_field(fname);

    EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Import field view has not been allocated yet.\n");

    EKAT_REQUIRE_MSG(cpl_idx>=0, "Error! Input cpl_idx is negative.\n");

    // Check that this cpl_idx wasn't already registered
    for (int i=0; i<m_num_scream_imports; ++i) {
      EKAT_REQUIRE_MSG(cpl_idx!=m_scream_imports_host(i).cpl_idx,
                       "Error! An import with cpl_idx " + std::to_string(cpl_idx) + " was already registered.\n");
    }

    auto& info = m_scream_imports_host(m_num_scream_imports);

    // Set view data ptr
    info.data = field.get_internal_view_data<Real>();

    // Set cpl index
    info.cpl_idx = cpl_idx;

    // For import fluxes, we must change the sign as cpl and atm interprete the direction differently.
    if (fname == "surf_mom_flux" || fname == "surf_sens_flux" ||
        fname == "surf_evap"     || fname == "surf_lw_flux_up") {
      m_cpl_scream_sign_change_host(m_num_scream_imports) = -1;
    } else {
      m_cpl_scream_sign_change_host(m_num_scream_imports) = 1;
    }

    // Get column offset and stride
    get_col_info (field.get_header_ptr(), vecComp, info.col_offset, info.col_stride);

    // Store the identifier of this field, for debug purposes
    m_imports_fids.insert(field.get_header().get_identifier());

    // Update number of imports stored
    ++m_num_scream_imports;
  }
}

void SurfaceCoupling::
register_export (const std::string& fname,
                 const int cpl_idx,
                 const int vecComp,
                 const bool export_during_init)
{
  // Two separate checks rather than state==Open, so we can print more specific error messages
  EKAT_REQUIRE_MSG (m_state!=RepoState::Clean,
                      "Error! Registration phase has not started yet.\n"
                      "       Did you forget to call set_num_fields(..) ?\n");
  EKAT_REQUIRE_MSG (m_state!=RepoState::Closed, "Error! Registration phase has already ended.\n");

  // Check that we still have room
  EKAT_REQUIRE_MSG(m_num_scream_exports<m_scream_exports_host.extent_int(0),
                     "Error! Exports view is already full. Did you call 'set_num_fields' with the wrong arguments?\n");

  EKAT_REQUIRE_MSG(cpl_idx>=0, "Error! Input cpl_idx is negative.\n");

  // Check that this cpl_idx wasn't already registered
  for (int i=0; i<m_num_scream_exports; ++i) {
    EKAT_REQUIRE_MSG(cpl_idx!=m_scream_exports_host(i).cpl_idx,
                       "Error! An export with cpl_idx " + std::to_string(cpl_idx) + " was already registered.\n");
  }

  auto& info = m_scream_exports_host(m_num_scream_exports);

  // 2 cases for setting export info:
  //   1. General case: fname describes a field in the FieldManager which has no parent.
  //      In this case, gather info directly from the field.
  //   2. Special cases
  //      a. fname == set_zero, which indicates that this field is not used in SCREAM.
  //      c. fname corresponeds to a field which is a combination of SCREAM fields.
  //      For these special cases, we have member variable 1d views which we will store the correct
  //      values at the surface for each field during do_export(). Here, just set the data
  //      to point to these views, and fill in stride/offset info using a field setup in
  //      the constructor of the class
  if (m_field_mgr->has_field(fname)) {

    // If fname is a field in the the field manager, then set the data and column info based on this field.
    Field field = m_field_mgr->get_field(fname);

    EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Export field view has not been allocated yet.\n");

    // Set view data ptr
    info.data = field.get_internal_view_data<Real>();

    // Get column offset and stride
    get_col_info (field.get_header_ptr(), vecComp, info.col_offset, info.col_stride);

    // Store the identifier of this field, for debug purposes
    m_exports_fids.insert(field.get_header().get_identifier());

  } else {

    // Set view data ptr
    if (fname == "set_zero")        info.data = zero_view.data();
    else if (fname == "Sa_z")       info.data = Sa_z.data();
    else if (fname == "Sa_ptem")    info.data = Sa_ptem.data();
    else if (fname == "Sa_dens")    info.data = Sa_dens.data();
    else if (fname == "Sa_pslv")    info.data = Sa_pslv.data();
    else if (fname == "Faxa_rainl") info.data = Faxa_rainl.data();
    else if (fname == "Faxa_snowl") info.data = Faxa_snowl.data();
    else                            EKAT_ERROR_MSG("Error! Unrecognized export field name \"" + fname + "\".");

    // Get column offset and stride
    get_col_info (dummy_field.get_header_ptr(), vecComp, info.col_offset, info.col_stride);

    // Store the identifier of this field, for debug purposes
    m_exports_fids.insert(dummy_field.get_header().get_identifier());
  }

  // Fields which are computed inside SCREAM should skip the initial export
  info.do_initial_export = export_during_init;

  // Set cpl index
  info.cpl_idx = cpl_idx;

  // Update number of exports stored
  ++m_num_scream_exports;
}

void SurfaceCoupling::
registration_ends (cpl_data_ptr_type cpl_imports_ptr,
                   cpl_data_ptr_type cpl_exports_ptr)
{
  // Two separate checks rather than state==Open, so we can print more specific error messages
  EKAT_REQUIRE_MSG (m_state!=RepoState::Clean,  "Error! Registration phase hasn't started yet.\n");
  EKAT_REQUIRE_MSG (m_state!=RepoState::Closed, "Error! Registration phase has already ended.\n");

  // Check that we registered all the imports/exports we said we would
  // TODO: this is not really needed for correctness. It only means we over-allocated.
  //       Should we allow it?
  EKAT_REQUIRE_MSG (m_num_scream_imports==m_scream_imports_host.extent_int(0),
                      "Error! You registered less imports than you said you would.\n"
                      "       Imports registered: " + std::to_string(m_num_scream_imports) + "\n"
                      "       Imports declared:   " + std::to_string(m_scream_imports_host.extent_int(0)) + "\n");
  EKAT_REQUIRE_MSG (m_num_scream_exports==m_scream_exports_host.extent_int(0),
                      "Error! You registered less exports than you said you would.\n"
                      "       Exports registered: " + std::to_string(m_num_scream_exports) + "\n"
                      "       Exports declared:   " + std::to_string(m_scream_exports_host.extent_int(0)) + "\n");

  // Loop over import/exports; make sure both data and field are set,
  // and that cpl indices form a range [0,N]
  for (int i=0; i<m_num_scream_imports; ++i) {
    EKAT_REQUIRE_MSG (m_scream_imports_host(i).data!=nullptr,
                        "Error! No field set for import index " + std::to_string(i) + ".\n");
  }
  for (int i=0; i<m_num_scream_exports; ++i) {
    EKAT_REQUIRE_MSG (m_scream_exports_host(i).data!=nullptr,
                        "Error! No field set for import index " + std::to_string(i) + ".\n");
  }

  // Deep copy the Info host view to device
  Kokkos::deep_copy(m_scream_imports_dev, m_scream_imports_host);
  Kokkos::deep_copy(m_scream_exports_dev, m_scream_exports_host);

  // Deep copy sign change view to device
  Kokkos::deep_copy(m_cpl_scream_sign_change_dev, m_cpl_scream_sign_change_host);

  if (m_num_scream_imports>0) {
    // Check input pointer
    EKAT_REQUIRE_MSG(cpl_imports_ptr!=nullptr, "Error! Data pointer for imports is null.\n");

    // Setup the host and device 2d views
    m_cpl_imports_view_h = decltype(m_cpl_imports_view_h)(cpl_imports_ptr,m_num_cols,m_num_cpl_imports);
    m_cpl_imports_view_d = Kokkos::create_mirror_view(device_type(),m_cpl_imports_view_h);
  }

  if (m_num_scream_exports>0) {
    // Check input pointer
    EKAT_REQUIRE_MSG(cpl_exports_ptr!=nullptr, "Error! Data pointer for exports is null.\n");

    // Setup the host and device 2d views
    m_cpl_exports_view_h = decltype(m_cpl_exports_view_h)(cpl_exports_ptr,m_num_cols,m_num_scream_exports);
    m_cpl_exports_view_d = Kokkos::create_mirror_view(device_type(),m_cpl_exports_view_h);

    // Deep copy to preserve any existing data in cpl_exports_ptr
    Kokkos::deep_copy(m_cpl_exports_view_d,m_cpl_exports_view_h);
  }

  // Finally, mark registration as completed.
  m_state = RepoState::Closed;
}

void SurfaceCoupling::do_import ()
{
  if (m_num_scream_imports==0) {
    return;
  }

  using policy_type = KokkosTypes<device_type>::RangePolicy;

  // Local copies, to deal with CUDA's handling of *this.
  const auto scream_imports = m_scream_imports_dev;
  const auto cpl_imports_view_d = m_cpl_imports_view_d;
  const auto cpl_scream_sign_change = m_cpl_scream_sign_change_dev;
  const int num_cols = m_num_cols;

  // Deep copy cpl host array to device
  Kokkos::deep_copy(m_cpl_imports_view_d,m_cpl_imports_view_h);

  // Unpack the fields
  auto unpack_policy = policy_type(0,m_num_scream_imports*num_cols);
  Kokkos::parallel_for(unpack_policy, KOKKOS_LAMBDA(const int& i) {
    const int ifield = i / num_cols;
    const int icol   = i % num_cols;

    const auto& info = scream_imports(ifield);

    auto offset = icol*info.col_stride + info.col_offset;
    info.data[offset] = cpl_imports_view_d(icol,info.cpl_idx)*cpl_scream_sign_change(ifield);
  });
}

void SurfaceCoupling::
get_col_info(const std::shared_ptr<const FieldHeader>& fh,
             int vecComp, int& col_offset, int& col_stride) const
{
  // Each field is seen as a certain number of columns worth of data (CWD)
  // at every column point. E.g., a field with layout (ncols,nlevs) is 1 CWD,
  // while a field with layout (ncols,3,nlevs) is 3 CWD.
  // The col_offset encodes the offset of this colum *within a single mesh
  // point icol". So for a vector layout (ncols,3,nlevs), to import/export
  // the icomp-th component of the field, we'd have col_offset=icomp. For scalar
  // layouts, the offset is always 0 (only 1 CWD per mesh point).
  // HOWEVER: if this field was a subview of another field, things are a bit
  // more involved. Say the field F1 (ncols,3,nlevs) is a subview of a larger
  // field F2 with layout (ncols,nvars,3,nlevs), corresponding to ivar=2.
  // Then, F2 has 3*nvars CWD per mesh point. Since the underlying data ptr
  // is the same, to get the correct col_offset, we need to figure out the
  // offset inside F2. In this case, it would be col_offset = ivar*3+icomp.

  using namespace ShortFieldTagsNames;

  // Get this field layout info
  const auto& layout = fh->get_identifier().get_layout();
  const auto& dims = layout.dims();

  auto lt = get_layout_type(layout.tags());
  const bool scalar   = lt==LayoutType::Scalar2D || lt==LayoutType::Scalar3D;
  const bool vector   = lt==LayoutType::Vector2D || lt==LayoutType::Vector3D;
  const bool layout3d = lt==LayoutType::Scalar3D || lt==LayoutType::Vector3D;
  EKAT_REQUIRE_MSG(scalar || vector,
      "Error! Support for tensor fields not yet implemented in surface coupling.\n");
  EKAT_REQUIRE_MSG( (vecComp<0) == scalar,
      "Error! You can and must specify a vector component only for vector fields.\n");

  vecComp = std::max(0,vecComp);

  // Compute initial offset. We use get_last_extent() as it
  // accounts for padding (if any).
  col_offset = vecComp;
  if (layout3d) {
    if (lt==LayoutType::Vector3D) {
      col_offset *= fh->get_alloc_properties().get_last_extent();
    }
    col_offset += dims.back()-1;
  }

  // If rank>1, there always is some stride
  col_stride = 1;
  if (layout.rank()>1) {
    col_stride = fh->get_alloc_properties().get_last_extent();
    if (layout.rank()>2) {
      col_stride *= dims[1];
    }
  }

  // If this field has a parent, then the underlying data includes all entires
  // in the larger field. We then must treat the child field as having col_stride
  // and col_offset of a vector field with component as its subview_idx.
  std::shared_ptr<const FieldHeader> parent = fh->get_parent().lock();
  if (parent != nullptr) {

    EKAT_REQUIRE_MSG(parent->get_parent().lock() == nullptr,
                     "Error! Currently support isn't added for fields with grandparents.\n");

    const auto parent_lt = get_layout_type(parent->get_identifier().get_layout().tags());

    EKAT_REQUIRE_MSG(parent_lt==LayoutType::Vector3D,
                     "Error! SurfaceCoupling expects all subfields to have parents "
                     "with LayoutType::Vector3D.\n");

    const auto& sv_info = fh->get_alloc_properties().get_subview_info();

    // Recall: idx = (idim,k) = (dimension where slice happened, index along said dimension).
    // Field class only allows idim=0,1. But we should never be in the case of idim=0, here.
    // If we have idim=0, it means that the parent field did not have COL as tag[0].
    EKAT_REQUIRE_MSG(sv_info.dim_idx==1, "Error! Bizarre scenario discovered. Contact developers.\n");

    // Additional col_offset
    col_offset += sv_info.slice_idx*parent->get_alloc_properties().get_last_extent();

    // Additional product for col_stride
    col_stride *= parent->get_identifier().get_layout().dim(1);
  }
}

} // namespace control
} // namespace scream
