#include "control/surface_coupling.hpp"

namespace scream {
namespace control {

SurfaceCoupling::
SurfaceCoupling (const std::shared_ptr<const AbstractGrid>& grid)
 : m_state (RepoState::Clean)
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
                const import_field_type& field,
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

  // Check that input field is valid
  EKAT_REQUIRE_MSG (field.is_allocated(), "Error! Import field view has not been allocated yet.\n");
  
  const std::string& fgn = field.get_header().get_identifier().get_grid_name();
  EKAT_REQUIRE_MSG (fgn==m_grid_name,
                      "Error! Input field is defined on the wrong grid:\n"
                      "         expected grid name: " + m_grid_name + "\n"
                      "         field grid name: " + fgn + "\n");

  EKAT_REQUIRE_MSG(cpl_idx>=0, "Error! Input cpl_idx is negative.\n");

  // Check that this cpl_idx wasn't already registered
  for (int i=0; i<m_num_imports; ++i) {
    EKAT_REQUIRE_MSG(cpl_idx!=m_scream_imports_host(i).cpl_idx,
                       "Error! An import with cpl_idx " + std::to_string(cpl_idx) + " was already registered.\n");
  }

  auto& info = m_scream_imports_host(m_num_imports);

  // Set view data ptr
  info.data = field.get_view().data();

  const auto& layout = field.get_header().get_identifier().get_layout();
  const auto& alloc_prop = field.get_header().get_alloc_properties();
  const auto& tags = layout.tags();

  using namespace ShortFieldTagsNames;

  switch (layout.rank()) {
    case 1:
      EKAT_REQUIRE_MSG(vecComp==-1, "Error! Vector component specified, but field '" + fname + "' is a 2d scalar.\n");
      info.col_size = 1;
      break;
    case 2:
      EKAT_REQUIRE_MSG(tags[1]==VL || tags[1]==VAR || tags[1]==CMP,
                         "Error! Unexpected tag '" + e2str(tags[1]) + "' for second dimension of export field '" + fname + "'.\n");
      EKAT_REQUIRE_MSG(tags.back()!=VL || vecComp==-1, "Error! Vector component specified, but field '" + fname + "' is a 3d scalar.\n");

      info.col_size = alloc_prop.get_last_dim_extent<Real>();
      info.col_offset = vecComp;
      break;
    case 3:
      EKAT_REQUIRE_MSG(tags[1]==VAR || tags[1]==CMP,
                         "Error! Unexpected tag '" + e2str(tags[1]) + "' for second dimension of export field '" + fname + "'.\n");
      EKAT_REQUIRE_MSG(tags[2]==VL,
                         "Error! Unexpected tag '" + e2str(tags[2]) + "' for third dimension of export field '" + fname + "'.\n");
      EKAT_REQUIRE_MSG(vecComp>=0, "Error! Vector component not specified for 3d vector field '" + fname + "/.\n");

      info.col_size = layout.dim(1)*alloc_prop.get_last_dim_extent<Real>();
      info.col_offset = vecComp;
      break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank for import field '" + fname + "'.\n");
  }

  m_imports_fids.insert(field.get_header().get_identifier());

  // Update number of imports stored
  ++m_num_imports;
}

void SurfaceCoupling::
register_export (const std::string& fname,
                 const int cpl_idx,
                 const import_field_type& field,
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

  // Check that input field is valid
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

  // To figure out col_size and col_offset, we need to inspect the field layout,
  // as well as its allocation properties (to handle the case od padding)
  const auto& layout = field.get_header().get_identifier().get_layout();
  const auto& alloc_prop = field.get_header().get_alloc_properties();
  const auto& tags = layout.tags();

  using namespace ShortFieldTagsNames;

  switch (layout.rank()) {
    case 1:
      EKAT_REQUIRE_MSG(vecComp==-1, "Error! Vector component specified, but field '" + fname + "' is a 2d scalar.\n");
      info.col_size = 1;
      break;
    case 2:
      EKAT_REQUIRE_MSG(tags[1]==VL || tags[1]==VAR || tags[1]==CMP,
                         "Error! Unexpected tag '" + e2str(tags[1]) + "' for second dimension of import field '" + fname + "'.\n");
      EKAT_REQUIRE_MSG(tags.back()!=VL || vecComp==-1, "Error! Vector component specified, but field '" + fname + "' is a 3d scalar.\n");

      info.col_size = alloc_prop.get_last_dim_extent<Real>();
      info.col_offset = vecComp;
      break;
    case 3:
      EKAT_REQUIRE_MSG(tags[1]==VAR || tags[1]==CMP,
                         "Error! Unexpected tag '" + e2str(tags[1]) + "' for second dimension of import field '" + fname + "'.\n");
      EKAT_REQUIRE_MSG(tags[2]==VL,
                         "Error! Unexpected tag '" + e2str(tags[2]) + "' for third dimension of import field '" + fname + "'.\n");
      EKAT_REQUIRE_MSG(vecComp>=0, "Error! Vector component not specified for 3d vector field '" + fname + "/.\n");

      info.col_size = layout.dim(1)*alloc_prop.get_last_dim_extent<Real>();
      info.col_offset = vecComp;
      break;
    default:
      EKAT_ERROR_MSG("Error! Unsupported field rank for import field '" + fname + "'.\n");
  }

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

  // Check input pointers
  EKAT_REQUIRE_MSG(cpl_exports_ptr!=nullptr, "Error! Data pointer for exports is null.\n");
  EKAT_REQUIRE_MSG(cpl_imports_ptr!=nullptr, "Error! Data pointer for imports is null.\n");

  // Setup the host and device 2d views
  m_cpl_imports_view_h = decltype(m_cpl_imports_view_h)(cpl_imports_ptr,m_num_cols,m_num_imports);
  m_cpl_exports_view_h = decltype(m_cpl_exports_view_h)(cpl_exports_ptr,m_num_cols,m_num_exports);
  m_cpl_imports_view_d = Kokkos::create_mirror_view(device_type(),m_cpl_imports_view_h);
  m_cpl_exports_view_d = Kokkos::create_mirror_view(device_type(),m_cpl_exports_view_h);

  // Finally, mark registration as completed.
  m_state = RepoState::Closed;
}

void SurfaceCoupling::do_import ()
{
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

    auto offset = num_cols*info.col_size + info.col_offset;
    info.data[offset] = cpl_imports_view_d(icol,info.cpl_idx);
  });
}

void SurfaceCoupling::do_export ()
{
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

    auto offset = num_cols*info.col_size + info.col_offset;
    cpl_exports_view_d(icol,info.cpl_idx) = info.data[offset];
  });

  // Deep copy fields from device to cpl host array
  Kokkos::deep_copy(m_cpl_exports_view_h,m_cpl_exports_view_d);
}

} // namespace control
} // namespace scream
