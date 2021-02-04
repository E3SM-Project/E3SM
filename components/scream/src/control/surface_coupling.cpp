#include "control/surface_coupling.hpp"

#include "share/field/field_utils.hpp"

namespace scream {
namespace control {

SurfaceCoupling::
SurfaceCoupling (const std::shared_ptr<const AbstractGrid>& grid,
                 const FieldRepository<Real>& repo)
 : m_field_repo (repo)
 , m_state (RepoState::Clean)
{
  EKAT_REQUIRE_MSG (grid!=nullptr, "Error! Invalid grid pointer.\n");

  m_num_cols  = grid->get_num_local_dofs();
  m_grid_name = grid->name();

  EKAT_REQUIRE_MSG(grid->type()==GridType::Point,
                     "Error! Surface coupling only implemented for 'Point' grids.\n"
                     "       Input grid type: " + e2str(grid->type()) + "\n");
}

void SurfaceCoupling::
set_num_fields (const int num_imports, const int num_exports)
{
  EKAT_REQUIRE_MSG(m_state==RepoState::Clean,
                     "Error! You can only set the number of fields in SurfaceCoupling once.\n");

  m_scream_imports_dev = decltype(m_scream_imports_dev)("",num_imports);
  m_scream_exports_dev = decltype(m_scream_exports_dev)("",num_exports);

  m_scream_imports_host = Kokkos::create_mirror_view(m_scream_imports_dev);
  m_scream_exports_host = Kokkos::create_mirror_view(m_scream_exports_dev);

  // These will be incremented every time we register an import/export.
  // At registration end, they should match the length of the corresponding view above.
  m_num_imports = 0;
  m_num_exports = 0;

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

  // Check that we still have room
  EKAT_REQUIRE_MSG(m_num_imports<m_scream_imports_host.extent_int(0),
                     "Error! Imports view is already full. Did you call 'set_num_fields' with the wrong arguments?\n");

  // Get the field, and check that is valid
  import_field_type field = m_field_repo.get_field(fname,m_grid_name);
  EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Import field view has not been allocated yet.\n");
  
  EKAT_REQUIRE_MSG(cpl_idx>=0, "Error! Input cpl_idx is negative.\n");

  // Check that this cpl_idx wasn't already registered
  for (int i=0; i<m_num_imports; ++i) {
    EKAT_REQUIRE_MSG(cpl_idx!=m_scream_imports_host(i).cpl_idx,
                       "Error! An import with cpl_idx " + std::to_string(cpl_idx) + " was already registered.\n");
  }

  auto& info = m_scream_imports_host(m_num_imports);

  // Set view data ptr
  info.data = field.get_view().data();

  // Set cpl index
  info.cpl_idx = cpl_idx;

  // Get column offset and stride
  get_col_info (field.get_header_ptr(), vecComp, info.col_offset, info.col_stride);

  // Store the identifier of this field, for debug purposes
  m_imports_fids.insert(field.get_header().get_identifier());

  // Update number of imports stored
  ++m_num_imports;
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

  // Get the field, and check that is valid
  export_field_type field = m_field_repo.get_field(fname,m_grid_name);
  EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Export field view has not been allocated yet.\n");

  const std::string& fgn = field.get_header().get_identifier().get_grid_name();
  EKAT_REQUIRE_MSG (fgn==m_grid_name,
                      "Error! Input field '" + fname + "' is defined on the wrong grid:\n"
                      "       expected grid name: " + m_grid_name + "\n"
                      "       field grid name: " + fgn + "\n");

  EKAT_REQUIRE_MSG(cpl_idx>=0, "Error! Input cpl_idx is negative.\n");

  // Check that this cpl_idx wasn't already registered
  for (int i=0; i<m_num_exports; ++i) {
    EKAT_REQUIRE_MSG(cpl_idx!=m_scream_exports_host(i).cpl_idx,
                       "Error! An export with cpl_idx " + std::to_string(cpl_idx) + " was already registered.\n");
  }

  auto& info = m_scream_exports_host(m_num_exports);

  // Set view data ptr
  info.data = field.get_view().data();

  // Set cpl index
  info.cpl_idx = cpl_idx;

  // Get column offset and stride
  get_col_info (field.get_header_ptr(), vecComp, info.col_offset, info.col_stride);

  // Store the identifier of this field, for debug purposes
  m_exports_fids.insert(field.get_header().get_identifier());

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
  EKAT_REQUIRE_MSG (m_num_imports==m_scream_imports_host.extent_int(0),
                      "Error! You registered less imports than you said you would.\n"
                      "       Imports registered: " + std::to_string(m_num_imports) + "\n"
                      "       Imports declared:   " + std::to_string(m_scream_imports_host.extent_int(0)) + "\n");
  EKAT_REQUIRE_MSG (m_num_exports==m_scream_exports_host.extent_int(0),
                      "Error! You registered less exports than you said you would.\n"
                      "       Exports registered: " + std::to_string(m_num_exports) + "\n"
                      "       Exports declared:   " + std::to_string(m_scream_exports_host.extent_int(0)) + "\n");

  // Loop over import/exports; make sure both data and field are set,
  // and that cpl indices form a range [0,N]
  for (int i=0; i<m_num_imports; ++i) {
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

  if (m_num_imports>0) {
    // Check input pointer
    EKAT_REQUIRE_MSG(cpl_imports_ptr!=nullptr, "Error! Data pointer for imports is null.\n");

    // Setup the host and device 2d views
    m_cpl_imports_view_h = decltype(m_cpl_imports_view_h)(cpl_imports_ptr,m_num_cols,m_num_imports);
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
  if (m_num_imports==0) {
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
  auto unpack_policy = policy_type(0,m_num_imports*num_cols);
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

    auto offset = icol*info.col_stride + info.col_offset;
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
  const bool scalar = lt==LayoutType::Scalar2D || lt==LayoutType::Scalar3D;
  const bool vector = lt==LayoutType::Vector2D || lt==LayoutType::Vector3D;
  EKAT_REQUIRE_MSG(scalar || vector,
      "Error! Support for tensor fields not yet implemented in surface coupling.\n");
  EKAT_REQUIRE_MSG( (vecComp<0) == scalar,
      "Error! You can and must specify a vector component only for vector fields.\n");

  vecComp = std::max(0,vecComp);

  // Product of dimensions of field, without counting space dims, that is,
  // the product of dimensions of the "analytical" field.
  int dims_prod = 1;

  // Compute initial offset based on this field rank and requested component
  col_offset = vecComp;
  if (vector) {
    dims_prod = dims[1];
  }

  // If rank>1, there always is some stride
  col_stride = 1;
  if (layout.rank()>1) {
    // Note: use alloc prop's last extent, cause it includes padding (if any);
    col_stride = fh->get_alloc_properties().get_last_extent();
    if (layout.rank()>2) {
      col_stride *= dims[1];
    }
  }

  // Recursively go through parent field, adding offset depending
  // on the location of the slice.
  // Note: as of now, I don't expect there to be "a grandparent" for this field,
  //       so we could avoid the while loop. However, the logic is so simple, that
  //       I may as well add it, and be safe down the road if it happens.
  std::shared_ptr<const FieldHeader> me = fh;
  std::shared_ptr<const FieldHeader> p  = me->get_parent().lock();
  while (p!=nullptr) {
    const auto& idx = me->get_alloc_properties().get_subview_idx();

    // Recall: idx = (idim,k) = (dimension where slice happened, index along said dimension).
    // Field class only allows idim=0,1. But we should never be in the case of idim=0, here.
    // If we have idim=0, it means that the parent field did not have CL has tag[0]...
    EKAT_REQUIRE_MSG(idx.first==1, "Error! Bizarre scenario discovered. Contact developers.\n");

    // Update the offset
    col_offset += idx.second*dims_prod;

    // Update stride.
    col_stride *= p->get_identifier().get_layout().dim(1);

    // Update trail_dims, in case there's another parent.
    dims_prod *= p->get_identifier().get_layout().dim(1);

    me = p;
    p = me->get_parent().lock();
  }
}

} // namespace control
} // namespace scream
