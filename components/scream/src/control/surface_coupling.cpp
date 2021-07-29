#include "control/surface_coupling.hpp"

#include "share/field/field_utils.hpp"
#include "share/util/scream_common_physics_functions.hpp"

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
  dummy_field = Field<const Real>(id);

  EKAT_REQUIRE_MSG(grid->type()==GridType::Point,
      "Error! Surface coupling only implemented for 'Point' grids.\n"
      "       Input grid type: " + e2str(grid->type()) + "\n");
}

void SurfaceCoupling::
set_num_fields (const int num_cpl_imports, const int num_scream_imports, const int num_exports)
{
  EKAT_REQUIRE_MSG(m_state==RepoState::Clean,
                     "Error! You can only set the number of fields in SurfaceCoupling once.\n");

  m_scream_imports_dev = decltype(m_scream_imports_dev)("",num_scream_imports);
  m_scream_exports_dev = decltype(m_scream_exports_dev)("",num_exports);

  m_scream_imports_host = Kokkos::create_mirror_view(m_scream_imports_dev);
  m_scream_exports_host = Kokkos::create_mirror_view(m_scream_exports_dev);

  // These fields contain computation needed for some export fields
  Sa_ptem      = decltype(Sa_ptem)     ("", m_num_cols);
  Sa_u         = decltype(Sa_u)        ("", m_num_cols);
  Sa_v         = decltype(Sa_v)        ("", m_num_cols);
  Sa_dens      = decltype(Sa_dens)     ("", m_num_cols);
  zero_view    = decltype(zero_view)   ("", m_num_cols);
  Kokkos::deep_copy(zero_view, 0.0);

  // These will be incremented every time we register an import/export.
  // At registration end, they should match the length of the corresponding view above.
  m_num_scream_imports = 0;
  m_num_exports        = 0;

  // We keep up with the total number of imports for indexing input data
  m_num_cpl_imports = num_cpl_imports;

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

  if (fname == "unused" || fname == "RRTMGP")
  {
    // Do nothing
  } else {
    // Check that we still have room
    EKAT_REQUIRE_MSG(m_num_scream_imports<m_scream_imports_host.extent_int(0),
                       "Error! Imports view is already full. Did you call 'set_num_fields' with the wrong arguments?\n");

    // Get the field, and check that is valid
    import_field_type field = m_field_mgr->get_field(fname);

    EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Import field view has not been allocated yet.\n");

    EKAT_REQUIRE_MSG(cpl_idx>=0, "Error! Input cpl_idx is negative.\n");

    // Check that this cpl_idx wasn't already registered
    for (int i=0; i<m_num_scream_imports; ++i) {
      EKAT_REQUIRE_MSG(cpl_idx!=m_scream_imports_host(i).cpl_idx,
                         "Error! An import with cpl_idx " + std::to_string(cpl_idx) + " was already registered.\n");
    }

    auto& info = m_scream_imports_host(m_num_scream_imports);

    // Set view data ptr
    info.data = field.get_view().data();

    // Set cpl index
    info.cpl_idx = cpl_idx;

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
                 const int vecComp)
{
  // Two separate checks rather than state==Open, so we can print more specific error messages
  EKAT_REQUIRE_MSG (m_state!=RepoState::Clean,
                      "Error! Registration phase has not started yet.\n"
                      "       Did you forget to call set_num_fields(..) ?\n");
  EKAT_REQUIRE_MSG (m_state!=RepoState::Closed, "Error! Registration phase has already ended.\n");

  // Check that we still have room
  EKAT_REQUIRE_MSG(m_num_exports<m_scream_exports_host.extent_int(0),
                     "Error! Exports view is already full. Did you call 'set_num_fields' with the wrong arguments?\n");

  EKAT_REQUIRE_MSG(cpl_idx>=0, "Error! Input cpl_idx is negative.\n");

  // Check that this cpl_idx wasn't already registered
  for (int i=0; i<m_num_exports; ++i) {
    EKAT_REQUIRE_MSG(cpl_idx!=m_scream_exports_host(i).cpl_idx,
                       "Error! An export with cpl_idx " + std::to_string(cpl_idx) + " was already registered.\n");
  }

  auto& info = m_scream_exports_host(m_num_exports);

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
    // Note that the field must not have a parent as get_view() would be unsafe.
    export_field_type field = m_field_mgr->get_field(fname);

    EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Export field view has not been allocated yet.\n");

    // Set view data ptr
    info.data = field.get_view().data();

    // Get column offset and stride
    get_col_info (field.get_header_ptr(), vecComp, info.col_offset, info.col_stride);

    // Store the identifier of this field, for debug purposes
    m_exports_fids.insert(field.get_header().get_identifier());

  } else {

    // Set view data ptr
    if (fname == "set_zero")     info.data = zero_view.data();
    else if (fname == "Sa_ptem") info.data = Sa_ptem.data();
    else if (fname == "Sa_u")    info.data = Sa_u.data();
    else if (fname == "Sa_v")    info.data = Sa_v.data();
    else if (fname == "Sa_dens") info.data = Sa_dens.data();
    else                         EKAT_ERROR_MSG("Error! Unrecognized export field name \"" + fname + "\".");

    // Get column offset and stride
    get_col_info (dummy_field.get_header_ptr(), vecComp, info.col_offset, info.col_stride);

    // Store the identifier of this field, for debug purposes
    m_exports_fids.insert(dummy_field.get_header().get_identifier());
  }

  // Set cpl index
  info.cpl_idx = cpl_idx;

  // Update number of exports stored
  ++m_num_exports;
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
  EKAT_REQUIRE_MSG (m_num_exports==m_scream_exports_host.extent_int(0),
                      "Error! You registered less exports than you said you would.\n"
                      "       Exports registered: " + std::to_string(m_num_exports) + "\n"
                      "       Exports declared:   " + std::to_string(m_scream_exports_host.extent_int(0)) + "\n");

  // Loop over import/exports; make sure both data and field are set,
  // and that cpl indices form a range [0,N]
  for (int i=0; i<m_num_scream_imports; ++i) {
    EKAT_REQUIRE_MSG (m_scream_imports_host(i).data!=nullptr,
                        "Error! No field set for import index " + std::to_string(i) + ".\n");
  }
  for (int i=0; i<m_num_exports; ++i) {
    EKAT_REQUIRE_MSG (m_scream_exports_host(i).data!=nullptr,
                        "Error! No field set for import index " + std::to_string(i) + ".\n");
  }

  // Deep copy the Info host view to device
  Kokkos::deep_copy(m_scream_imports_dev, m_scream_imports_host);
  Kokkos::deep_copy(m_scream_exports_dev, m_scream_exports_host);

  if (m_num_scream_imports>0) {
    // Check input pointer
    EKAT_REQUIRE_MSG(cpl_imports_ptr!=nullptr, "Error! Data pointer for imports is null.\n");

    // Setup the host and device 2d views
    m_cpl_imports_view_h = decltype(m_cpl_imports_view_h)(cpl_imports_ptr,m_num_cols,m_num_cpl_imports);
    m_cpl_imports_view_d = Kokkos::create_mirror_view(device_type(),m_cpl_imports_view_h);
  }

  if (m_num_exports>0) {
    // Check input pointer
    EKAT_REQUIRE_MSG(cpl_exports_ptr!=nullptr, "Error! Data pointer for exports is null.\n");

    // Setup the host and device 2d views
    m_cpl_exports_view_h = decltype(m_cpl_exports_view_h)(cpl_exports_ptr,m_num_cols,m_num_exports);
    m_cpl_exports_view_d = Kokkos::create_mirror_view(device_type(),m_cpl_exports_view_h);
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
    info.data[offset] = cpl_imports_view_d(icol,info.cpl_idx);
  });
}

void SurfaceCoupling::do_export ()
{
  if (m_num_exports==0) {
    return;
  }

  using policy_type = KokkosTypes<device_type>::RangePolicy;
  using PF = PhysicsFunctions<device_type>;

  // For each export fields that is not trivially exist in the field
  // manager (see Case 2 in register_export()), calculate correct data
  // values.
  const bool scream_ad_run =
      (m_field_mgr->has_field("qv") && m_field_mgr->has_field("T_mid") && m_field_mgr->has_field("p_mid") &&
       m_field_mgr->has_field("horiz_winds") && m_field_mgr->has_field("pseudo_density"));
  if (scream_ad_run) {
    const int last_entry = m_num_levs-1;
    const auto qv             = m_field_mgr->get_field("qv").get_reshaped_view<const Real**>();
    const auto T_mid          = m_field_mgr->get_field("T_mid").get_reshaped_view<const Real**>();
    const auto p_mid          = m_field_mgr->get_field("p_mid").get_reshaped_view<const Real**>();
    const auto horiz_winds    = m_field_mgr->get_field("horiz_winds").get_reshaped_view<Real***>();
    const auto pseudo_density = m_field_mgr->get_field("pseudo_density").get_reshaped_view<const Real**>();

    const auto policy = policy_type (0, m_num_cols);
    Kokkos::parallel_for(policy, KOKKOS_LAMBDA(const int& i) {
      const auto dz = PF::calculate_dz(pseudo_density(i, last_entry), p_mid(i, last_entry),
                                       T_mid(i, last_entry), qv(i, last_entry));

      Sa_ptem(i) = PF::calculate_theta_from_T(T_mid(i, last_entry), p_mid(i, last_entry));
      Sa_u(i)    = horiz_winds(i, 0, last_entry);
      Sa_v(i)    = horiz_winds(i, 1, last_entry);
      Sa_dens(i) = PF::calculate_density(pseudo_density(i, last_entry), dz);
    });
  }

  // Local copies, to deal with CUDA's handling of *this.
  const auto scream_exports = m_scream_exports_dev;
  const auto cpl_exports_view_d = m_cpl_exports_view_d;
  const int num_cols = m_num_cols;

  // Pack the fields
  auto pack_policy   = policy_type (0,m_num_exports*num_cols);
  Kokkos::parallel_for(pack_policy, KOKKOS_LAMBDA(const int& i) {
    const int ifield = i / num_cols;
    const int icol   = i % num_cols;
    const auto& info = scream_exports(ifield);

    const auto offset = icol*info.col_stride + info.col_offset;
    cpl_exports_view_d(icol,info.cpl_idx) = info.data[offset];
  });

  // Deep copy fields from device to cpl host array
  Kokkos::deep_copy(m_cpl_exports_view_h,m_cpl_exports_view_d);
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

    // If the field has a parent field, it should be the case that the
    // parent field is a vectored field. Currently we conly expect
    EKAT_REQUIRE_MSG(parent_lt==LayoutType::Vector3D,
                     "Error! SurfaceCoupling expects all subfields to have parents"
                     " with LayoutType::Vector3D.\n");

    const auto& idx = fh->get_alloc_properties().get_subview_idx();

    // Recall: idx = (idim,k) = (dimension where slice happened, index along said dimension).
    // Field class only allows idim=0,1. But we should never be in the case of idim=0, here.
    // If we have idim=0, it means that the parent field did not have CL has tag[0]...
    EKAT_REQUIRE_MSG(idx.first==1, "Error! Bizarre scenario discovered. Contact developers.\n");

    // Additional col_offset
    col_offset += idx.second*parent->get_alloc_properties().get_last_extent();

    // Additional product for col_stride
    col_stride *= parent->get_identifier().get_layout().dim(1);
  }
}

} // namespace control
} // namespace scream
