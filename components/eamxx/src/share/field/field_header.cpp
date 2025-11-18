#include "share/field/field_header.hpp"

#include <ekat_std_utils.hpp>

namespace {
std::vector<std::int64_t> get_strides (const scream::FieldHeader& fh)
{
  const auto& fid = fh.get_identifier();
  const auto& fl  = fid.get_layout();
  const int rank = fl.rank();

  std::vector<std::int64_t> strides(rank);
  if (rank==0) {
    return strides;
  }

  const auto& fap = fh.get_alloc_properties ();
  auto p = fh.get_parent();
  if (p) {
    strides = get_strides(*p);
    const auto& sv_info = fap.get_subview_info();
    strides.erase(strides.begin()+sv_info.dim_idx);
  } else {
    auto dims = fl.dims();
    dims.back() = fap.get_last_extent();
    strides.back() = 1;
    for (int i=rank-1; i>=1; --i) {
      strides[i-1] = strides[i]*dims[i];
    }
  }

  return strides;
}

// Get offset of field actual data from the start of the internal data pointer
std::uint64_t get_offset (const scream::FieldHeader& fh)
{
  auto p = fh.get_parent();
  if (not p) {
    return 0;
  }

  std::uint64_t p_offset = get_offset(*p);

  const auto& fap = fh.get_alloc_properties ();
  const auto& sv_info = fap.get_subview_info();
  if (sv_info.dim_idx>0) {
    return p_offset;
  }

  return p_offset*get_strides(fh)[0]*sv_info.slice_idx;
}

} // anonymous namespace

namespace scream {

FieldHeader::FieldHeader (const identifier_type& id)
 : m_identifier (id)
{
  m_tracking   = std::make_shared<FieldTracking>();
  m_alloc_prop = std::make_shared<FieldAllocProp>(get_type_size(id.data_type()));
  m_extra_data = std::make_shared<extra_data_type>();

  // Let's add immediately this att, so that users don't need to check
  // if it already exist before adding string attributes for io.
  using stratts_t = std::map<std::string,std::string>;
  set_extra_data("io: string attributes",stratts_t());
}

void FieldHeader::
set_extra_data (const std::string& key,
                const std::any& data,
                const bool throw_if_existing)
{
  if (throw_if_existing) {
    EKAT_REQUIRE_MSG (m_extra_data->find(key)==m_extra_data->end(),
                        "Error! Key '" + key + "' already existing in "
                        "the extra data map of field '" + m_identifier.get_id_string() + "'.\n");
    (*m_extra_data)[key] = data;
  } else {
    (*m_extra_data)[key] = data;
  }
}

std::shared_ptr<FieldHeader> FieldHeader::alias(const std::string& name) const {
  auto fh = std::make_shared<FieldHeader>(get_identifier().alias(name));
  if (get_parent() != nullptr) {
    // If we're aliasing, we MUST keep track of the parent
    fh->create_parent_child_link(get_parent());
  }
  fh->m_tracking = m_tracking;
  fh->m_alloc_prop = m_alloc_prop;
  fh->m_extra_data = m_extra_data;
  return fh;
}

bool FieldHeader::is_aliasing (const FieldHeader& rhs) const
{
  if (this==&rhs)
    return true;

  if (m_tracking==rhs.m_tracking and
      m_alloc_prop==rhs.m_alloc_prop and
      m_extra_data==rhs.m_extra_data)
    return true;

  auto p = get_parent();
  auto rhs_p = rhs.get_parent();
  if (p!=nullptr and rhs_p!=nullptr) {
    return p->is_aliasing(*rhs_p) and
           m_alloc_prop->get_subview_info()==rhs.m_alloc_prop->get_subview_info();
  }

  return false;
}

#ifdef EAMXX_HAS_PYTHON
void FieldHeader::create_dldensor ()
{
  // FieldHeader is non-const, as we may add extra data
  if (has_extra_data("dltensor")) {
    return;
  }

  EKAT_REQUIRE_MSG (m_alloc_prop.is_committed(),
      "Error! Cannot crate dlpack data until field alloc props are committed.\n"
      " - field name: " + get_identifier().name() + "\n");
  EKAT_REQUIRE_MSG (not m_alloc_prop.get_subview_info().dynamic,
      "Error! We cannot create dlpack data for a field that is a dynamic subfield of another.\n"
      " - field name: " + get_identifier().name() + "\n");
  
  DLTensor tensor;
  tensor.ndim = f.rank();

  const auto& fid = get_identifier();
  const auto& fl  = fid.get_layout();

  // We don't have strides, and even shape cannot be taken from the layout dims,
  // since DLTensor expects int64_t*, while layout stores int* (inside a vector)
  if (not has_extra_data("dlpack_shape")) {
    std::vector<std::int64_t> shape;
    for (auto d : fl.dims()) {
      shape.push_back(d);
    }
    set_extra_data("dlpack_shape",shape);
  }
  if (not has_extra_data("dlpack_strides")) {
    std::vector<std::int64_t> strides = get_strides(*this);
    set_extra_data("dlpack_strides",strides);
  }

  tensor.shape   = get_extra_data<std::vector<std::int64_t>>("dlpack_shape").data();
  tensor.strides = get_extra_data<std::vector<std::int64_t>>("dlpack_strides").data();
  tensor.byte_offset = get_offset(*this);

  switch (fid.data_type()) {
    case DataType::IntType:
      tensor.dtype.code = kDLInt;
      tensor.dtype.bits = 8*sizeof(int);
      tensor.data = f.get_internal_view_data<int,HD>();
      break;
    case DataType::FloatType:
      tensor.dtype.code = kDLFloat;
      tensor.dtype.bits = 8*sizeof(float);
      tensor.data = f.get_internal_view_data<float,HD>();
      break;
    case DataType::DoubleType:
      tensor.dtype.code = kDLFloat;
      tensor.dtype.bits = 8*sizeof(double);
      tensor.data = f.get_internal_view_data<double,HD>();
      break;
    default:
      EKAT_ERROR_MSG ("Unrecognized/unsupported data type.\n");
  }

  return tensor;

}
void FieldHeader::create_dldevice ()
{
  if (has_extra_data("dldevice"))
    return;

  DLDevice device;
  device.device_type = PySession::dl_device_type();
  if constexpr (ekat::OnGpu<typename DefaultDevice::execution_space>::value) {
    device.device_id = Kokkos::device_id();
  } else {
    device.device_id = 0;
  }
}
#endif
// ---------------- Free function -------------------- //

std::shared_ptr<FieldHeader>
create_subfield_header (const FieldIdentifier& id,
                        std::shared_ptr<FieldHeader> parent,
                        const int idim, const int k, const bool dynamic)
{
  // Sanity checks
  EKAT_REQUIRE_MSG (parent!=nullptr,
      "Error! Invalid pointer for parent header.\n");

  // Create header, and set up parent/child
  auto fh = std::make_shared<FieldHeader>(id);
  fh->create_parent_child_link(parent);

  // Create tracking, and set up parent/child
  fh->m_tracking = std::make_shared<FieldTracking>();
  fh->m_tracking->create_parent_child_link(parent->get_tracking_ptr());
  if (parent->get_tracking().get_time_stamp().is_valid()) {
    fh->m_tracking->update_time_stamp(parent->get_tracking().get_time_stamp());
  }

  // Create alloc props
  fh->m_alloc_prop = std::make_shared<FieldAllocProp>(parent->get_alloc_properties().subview(idim,k,dynamic));

  return fh;
}

// subfield with multiple, contiguous slices
std::shared_ptr<FieldHeader>
create_subfield_header(const FieldIdentifier& id,
                       std::shared_ptr<FieldHeader> parent, const int idim,
                       const int k_beg, const int k_end) {
  // Sanity checks
  EKAT_REQUIRE_MSG(parent != nullptr,
                   "Error! Invalid pointer for parent header.\n");
  EKAT_REQUIRE_MSG(k_end > k_beg,
                   "Error! Slice indices are invalid (non-increasing).\n");

  // Create header, and set up parent/child
  auto fh = std::make_shared<FieldHeader>(id);
  fh->create_parent_child_link(parent);

  // Create tracking, and set up parent/child
  fh->m_tracking = std::make_shared<FieldTracking>();
  fh->m_tracking->create_parent_child_link(parent->get_tracking_ptr());
  if (parent->get_tracking().get_time_stamp().is_valid()) {
    fh->m_tracking->update_time_stamp(parent->get_tracking().get_time_stamp());
  }

  // Create alloc props
  fh->m_alloc_prop = std::make_shared<FieldAllocProp>(
      parent->get_alloc_properties().subview(idim, k_beg, k_end));

  return fh;
}

} // namespace scream
