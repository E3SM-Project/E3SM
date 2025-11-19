#include "share/field/field.hpp"
#include "share/util/eamxx_utils.hpp"

#include "share/field/field_impl_details.hpp"

namespace scream
{

Field::
Field (const identifier_type& id, bool allocate)
{
  m_header = std::make_shared<FieldHeader>(id);
  if (allocate) {
    allocate_view();
  }
}

Field
Field::get_const() const {
  Field f(*this);
  f.m_is_read_only = true;
  return f;
}

Field
Field::clone() const {
  return clone(name());
}

Field
Field::alias (const std::string& name) const {
  Field f;
  f.m_header = get_header().alias(name);
  f.m_data = m_data;
  f.m_is_read_only = m_is_read_only;
  return f;
}

Field
Field::clone(const std::string& name) const {
  return clone(name, get_header().get_identifier().get_grid_name());
}

Field
Field::clone(const std::string& name, const std::string& grid_name) const {
  // Create new field
  const auto& my_fid = get_header().get_identifier();
  FieldIdentifier fid(name,my_fid.get_layout(),my_fid.get_units(),
                      grid_name,my_fid.data_type());
  Field f(fid);

  // Ensure alloc props match
  const auto&  ap = get_header().get_alloc_properties();
        auto& fap = f.get_header().get_alloc_properties();
  fap.request_allocation(ap.get_largest_pack_size());

  // Allocate
  f.allocate_view();

  // Set correct time stamp
  const auto& ts = get_header().get_tracking().get_time_stamp();
  f.get_header().get_tracking().update_time_stamp(ts);

  // Deep copy
  f.deep_copy(*this);
  f.sync_to_host();

  return f;
}

Field Field::
subfield (const std::string& sf_name, const ekat::units::Units& sf_units,
          const int idim, const int index, const bool dynamic) const {

  const auto& id = m_header->get_identifier();
  const auto& lt = id.get_layout();

  // Sanity checks
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Input field must be allocated in order to subview it.\n");
  EKAT_REQUIRE_MSG (idim==0 || idim==1,
        "Error! Subview dimension index must be either 0 or 1.\n");

  // Create identifier for subfield
  FieldIdentifier sf_id(sf_name,lt.clone().strip_dim(idim),sf_units,id.get_grid_name(),id.data_type());

  // Create empty subfield, then set header and views
  // Note: we can access protected members, since it's the same type
  Field sf;
  sf.m_header = create_subfield_header(sf_id,m_header,idim,index,dynamic);
  sf.m_data = m_data;
  sf.m_is_read_only = m_is_read_only;

  if (not sf.m_header->get_alloc_properties().contiguous() and
      not sf.host_and_device_share_memory_space()) {
    // If subfield is not contiguous and Host and Device do not
    // share a memory space, we must initialize the helper field
    // for sync_to functions.
    sf.initialize_contiguous_helper_field();
  }

  return sf;
}

Field Field::
subfield (const std::string& sf_name, const int idim, const int index, const bool dynamic) const {
  const auto& id = m_header->get_identifier();
  return subfield(sf_name,id.get_units(),idim,index,dynamic);
}

Field Field::
subfield (const int idim, const int index, const bool dynamic) const {
  return subfield(m_header->get_identifier().name(),idim,index,dynamic);
}

Field Field::
subfield (const FieldTag tag, const int index, const bool dynamic) const {
  int idim = get_header().get_identifier().get_layout().dim_idx(tag);
  return subfield(idim,index,dynamic);
}

// slice at index idim, extracting the N = (index_end - index_beg) entries
// written in math notation: [index_beg, index_end)
// or equivalently, subF = F(index_beg, ... , index_beg + N)
Field Field::subfield(const std::string& sf_name,
                      const ekat::units::Units& sf_units, const int idim,
                      const int index_beg, const int index_end) const {

  const auto& id = m_header->get_identifier();
  const auto& lt = id.get_layout();

  // Sanity checks
  EKAT_REQUIRE_MSG(
      is_allocated(),
      "Error! Input field must be allocated in order to subview it.\n");

  auto sf_layout = lt.clone();
  sf_layout.reset_dim(idim, index_end - index_beg);

  // Create identifier for subfield
  FieldIdentifier sf_id(sf_name, sf_layout, sf_units, id.get_grid_name(), id.data_type());

  // Create empty subfield, then set header and views
  // Note: we can access protected members, since it's the same type
  Field sf;
  sf.m_header = create_subfield_header(sf_id, m_header, idim, index_beg,
                                       index_end);
  sf.m_data = m_data;

  if (not sf.m_header->get_alloc_properties().contiguous() and
      not sf.host_and_device_share_memory_space()) {
    // If subfield is not contiguous and Host and Device do not
    // share a memory space, we must initialize the helper field
    // for sync_to functions.
    sf.initialize_contiguous_helper_field();
  }

  return sf;
}

Field Field::subfield(const std::string& sf_name, const int idim,
                      const int index_beg, const int index_end) const {
  const auto& id = m_header->get_identifier();
  return subfield(sf_name, id.get_units(), idim, index_beg, index_end);
}

Field Field::subfield(const int idim, const int index_beg,
                      const int index_end) const {
  return subfield(m_header->get_identifier().name(), idim, index_beg,
                  index_end);
}

Field Field::
get_component (const int i, const bool dynamic) const {
  const auto& layout = get_header().get_identifier().get_layout();
  const auto& fname = get_header().get_identifier().name();
  EKAT_REQUIRE_MSG (layout.is_vector_layout(),
      "Error! 'get_component' available only for vector fields.\n"
      "       Layout of '" + fname + "': " + e2str(layout.type()) + "\n");

  const int idim = layout.get_vector_component_idx();
  EKAT_REQUIRE_MSG (i>=0 && i<layout.dim(idim),
      "Error! Component index out of bounds [0," + std::to_string(layout.dim(idim)) + ").\n");

  // Add _$i to the field name, to avoid issues if the subfield is stored
  // in some structure that requires unique names (e.g., a remapper)
  return subfield (fname + "_" + std::to_string(i),idim,i,dynamic);
}

Field Field::get_components(const int beg, const int end) const {
  const auto& layout = get_header().get_identifier().get_layout();
  const auto& fname = get_header().get_identifier().name();
  EKAT_REQUIRE_MSG(layout.is_vector_layout(),
                   "Error! 'get_component' available only for vector fields.\n"
                   "       Layout of '" +
                       fname + "': " + e2str(layout.type()) + "\n");

  const int idim = layout.get_vector_component_idx();
  EKAT_REQUIRE_MSG(beg >= 0 && end < layout.dim(idim),
                   "Error! Component index range out of bounds [0," +
                       std::to_string(layout.dim(idim)) + ").\n");
  EKAT_REQUIRE_MSG(beg < end, "Error! Invalid component indices (beg >= end).\n");

  // Add _$beg-$end to the field name, to avoid issues if the subfield is stored
  // in some structure that requires unique names (e.g., a remapper)
  return subfield(fname + "_" + std::to_string(beg) + "-" + std::to_string(end),
                  idim, beg, end);
}

bool Field::is_aliasing(const Field& rhs) const
{
  if (this==&rhs)
    return true;  // Same object

  if (not is_allocated() or not rhs.is_allocated())
    return false; // Once allocated, they will be different

  // NOTE: I'm not sure we NEED to check m_data, but we might as well
  return m_header->is_aliasing(rhs.get_header()) and
         m_data.d_view==rhs.m_data.d_view;
}

void Field::allocate_view ()
{
  // Not sure if simply returning would be safe enough. Re-allocating
  // would definitely be error prone (someone may have already gotten
  // a subview of the field). However, it *seems* suspicious to call
  // this method twice, and I think it's more likely than not that
  // such a scenario would indicate a bug. Therefore, I am prohibiting it.
  EKAT_REQUIRE_MSG(!is_allocated(), "Error! View was already allocated.\n");

  // Short names
  const auto& id     = m_header->get_identifier();
  const auto& layout = id.get_layout();
  auto& alloc_prop   = m_header->get_alloc_properties();

  // Commit the allocation properties
  alloc_prop.commit(layout);

  // Create the view, by quering allocation properties for the allocation size
  const auto view_dim = alloc_prop.get_alloc_size();

  m_data.d_view = decltype(m_data.d_view)(id.name(),view_dim);
  m_data.h_view = Kokkos::create_mirror_view(m_data.d_view);
}


void Field::deep_copy (const Field& src)
{
  update<CombineMode::Replace>(src,1,0);
}

void Field::deep_copy (const ScalarWrapper value) {
  EKAT_REQUIRE_MSG (not m_is_read_only,
      "Error! Cannot call deep_copy on read-only fields.\n");

  const auto my_data_type = data_type();
  switch (my_data_type) {
    case DataType::IntType:
      deep_copy_impl<false>(value.as<int>(),*this); // 2nd arg unused
      break;
    case DataType::FloatType:
      deep_copy_impl<false>(value.as<float>(),*this); // 2nd arg unused
      break;
    case DataType::DoubleType:
      deep_copy_impl<false>(value.as<double>(),*this); // 2nd arg unused
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized field data type in Field::deep_copy.\n");
  }
}

void Field::deep_copy (const ScalarWrapper value, const Field& mask)
{
  EKAT_REQUIRE_MSG (not m_is_read_only,
      "Error! Cannot call deep_copy on read-only fields.\n");

  const auto my_data_type = data_type();
  switch (my_data_type) {
    case DataType::IntType:
      deep_copy_impl<true>(value.as<int>(),mask);
      break;
    case DataType::FloatType:
      deep_copy_impl<true>(value.as<float>(),mask);
      break;
    case DataType::DoubleType:
      deep_copy_impl<true>(value.as<double>(),mask);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized field data type in Field::deep_copy.\n");
  }
}

void Field::scale (const ScalarWrapper beta)
{
  auto zero = ScalarWrapper::zero();
  update<CombineMode::Update>(*this,zero,beta);
}

void Field::scale (const Field& x)
{
  auto one = ScalarWrapper::one();
  update<CombineMode::Multiply>(x,one,one);
}

void Field::scale_inv (const Field& x)
{
  auto one = ScalarWrapper::one();
  update<CombineMode::Divide>(x,one,one);
}

void Field::max (const Field& x)
{
  auto one = ScalarWrapper::one();
  update<CombineMode::Max>(x,one,one);
}

void Field::min (const Field& x)
{
  auto one = ScalarWrapper::one();
  update<CombineMode::Min>(x,one,one);
}

template<CombineMode CM>
void Field::
update (const Field& x, const ScalarWrapper alpha, const ScalarWrapper beta)
{
  // Check this field is writable
  EKAT_REQUIRE_MSG (not is_read_only(),
      "Error! Cannot update field, as it is read-only.\n"
      " - field name: " + name() + "\n");

  const auto& dt = data_type();
  const auto& rhs_dt = x.data_type();

  // If user passes, say, double alpha/beta for an int field, we should error out, warning about
  // a potential narrowing rounding. The other way around, otoh, is allowed (even though
  // there's an upper limit to the int values that a double can store, it is unlikely the user
  // will use such large factors).
  // Similarly, we allow updating a field Y with another X as long as converting the data type of X
  // to the data type of Y does not require narrowing
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(alpha.type,dt),
      "Error! Coefficient alpha may be narrowed when converted to x/y data type.\n"
      " - x/y data type  : " + e2str(dt) + "\n"
      " - alpha data type: " + e2str(alpha.type) + "\n");
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(beta.type,dt),
      "Error! Coefficient beta may be narrowed when converted to x/y data type.\n"
      " - x/y data type  : " + e2str(dt) + "\n"
      " - beta data type: " + e2str(beta.type) + "\n");
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(rhs_dt,dt),
      "Error! Right hand side data type may be narrowed when converted to x data type.\n"
      " - rhs data type: " + e2str(rhs_dt) + "\n"
      " - lhs data type: " + e2str(dt) + "\n");

  // Check x/y are allocated
  EKAT_REQUIRE_MSG (is_allocated(),
      "Error! Cannot update field, since it is not allocated.\n"
      " - field name: " + name() + "\n");
  EKAT_REQUIRE_MSG (x.is_allocated(),
      "Error! Cannot update field, since source field is not allocated.\n"
      " - field name: " + x.name() + "\n");

  const auto& y_l = get_header().get_identifier().get_layout();
  const auto& x_l = x.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (y_l==x_l,
      "Error! Incompatible layouts for update_field.\n"
      " - x name: " + x.name() + "\n"
      " - y name: " + name() + "\n"
      " - x layout: " + x_l.to_string() + "\n"
      " - y layout: " + y_l.to_string() + "\n");

  // Determine if the RHS can contain fill_value entries
  bool fill_aware = x.get_header().may_be_filled();

  if (dt==DataType::IntType) {
    auto a = alpha.as<int>();
    auto b = beta.as<int>();
    if (fill_aware) {
      return update_impl<CM,true,int,int>(x,a,b);
    } else {
      return update_impl<CM,false,int,int>(x,a,b);
    }
  } else if (dt==DataType::FloatType) {
    auto a = alpha.as<float>();
    auto b = beta.as<float>();
    if (fill_aware) {
      if (rhs_dt==DataType::FloatType)
        return update_impl<CM,true,float,float>(x,a,b);
      else
        return update_impl<CM,true,float,int>(x,a,b);
    } else {
      if (rhs_dt==DataType::FloatType)
        return update_impl<CM,false,float,float>(x,a,b);
      else
        return update_impl<CM,false,float,int>(x,a,b);
    }
  } else if (dt==DataType::DoubleType) {
    auto a = alpha.as<double>();
    auto b = beta.as<double>();
    if (fill_aware) {
      if (rhs_dt==DataType::DoubleType)
        return update_impl<CM,true,double,double>(x,a,b);
      else if (rhs_dt==DataType::FloatType)
        return update_impl<CM,true,double,float>(x,a,b);
      else
        return update_impl<CM,true,double,int>(x,a,b);
    } else {
      if (rhs_dt==DataType::DoubleType)
        return update_impl<CM,false,double,double>(x,a,b);
      else if (rhs_dt==DataType::FloatType)
        return update_impl<CM,false,double,float>(x,a,b);
      else
        return update_impl<CM,false,double,int>(x,a,b);
    }
  } else {
    EKAT_ERROR_MSG ("Error! Unrecognized/unsupported field data type in Field::update.\n");
  }
}

template<CombineMode CM, bool FillAware, typename ST, typename XST>
void Field::
update_impl (const Field& x, const ST alpha, const ST beta)
{
  const auto& layout = x.get_header().get_identifier().get_layout();
  const auto& dims = layout.dims();

  // Must handle the case where one of the two views is strided (or both)
  const auto x_contig = x.get_header().get_alloc_properties().contiguous();
  const auto y_contig = get_header().get_alloc_properties().contiguous();
  switch (layout.rank()) {
    case 0:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST>(),
                               x.get_view<const XST>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST>(),
                               x.get_strided_view<const XST>(),
                               alpha,beta,dims);
      break;
    case 1:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST*>(),
                               x.get_view<const XST*>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST*>(),
                               x.get_strided_view<const XST*>(),
                               alpha,beta,dims);
      break;
    case 2:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST**>(),
                               x.get_view<const XST**>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST**>(),
                               x.get_strided_view<const XST**>(),
                               alpha,beta,dims);
      break;
    case 3:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST***>(),
                               x.get_view<const XST***>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST***>(),
                               x.get_strided_view<const XST***>(),
                               alpha,beta,dims);
      break;
    case 4:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST****>(),
                               x.get_view<const XST****>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST****>(),
                               x.get_strided_view<const XST****>(),
                               alpha,beta,dims);
      break;
    case 5:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST*****>(),
                               x.get_view<const XST*****>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST*****>(),
                               x.get_strided_view<const XST*****>(),
                               alpha,beta,dims);
      break;
    case 6:
      if (x_contig and y_contig)
        details::cvh<CM,FillAware>(get_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims);
      else if (x_contig)
        details::cvh<CM,FillAware>(get_strided_view<ST******>(),
                               x.get_view<const XST******>(),
                               alpha,beta,dims);
      else if (y_contig)
        details::cvh<CM,FillAware>(get_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims);
      else
        details::cvh<CM,FillAware>(get_strided_view<ST******>(),
                               x.get_strided_view<const XST******>(),
                               alpha,beta,dims);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Rank not supported in update_field.\n"
          " - x name: " + x.name() + "\n"
          " - y name: " + name() + "\n");
  }
  Kokkos::fence();
}

template<bool use_mask, typename ST>
void Field::deep_copy_impl (const ST value, const Field& mask)
{
  const auto& layout = get_header().get_identifier().get_layout();
  const auto  rank   = layout.rank();
  const auto& dims   = layout.dims();
  const auto contig = get_header().get_alloc_properties().contiguous();

  switch (rank) {
    case 0:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST>(),value,dims,
                                 mask.get_view<const int>());
        else
          details::svm<use_mask>(get_strided_view<ST>(),value,dims,
                                 mask.get_view<const int>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST>(),value,dims);
      }
      break;
    case 1:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST*>(),value,dims,
                                 mask.get_view<const int*>());
        else
          details::svm<use_mask>(get_strided_view<ST*>(),value,dims,
                                 mask.get_view<const int*>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST*>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST*>(),value,dims);
      }
      break;
    case 2:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST**>(),value,dims,
                                 mask.get_view<const int**>());
        else
          details::svm<use_mask>(get_strided_view<ST**>(),value,dims,
                                 mask.get_view<const int**>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST**>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST**>(),value,dims);
      }
      break;
    case 3:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST***>(),value,dims,
                                 mask.get_view<const int***>());
        else
          details::svm<use_mask>(get_strided_view<ST***>(),value,dims,
                                 mask.get_view<const int***>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST***>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST***>(),value,dims);
      }
      break;
    case 4:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST****>(),value,dims,
                                 mask.get_view<const int****>());
        else
          details::svm<use_mask>(get_strided_view<ST****>(),value,dims,
                                 mask.get_view<const int****>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST****>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST****>(),value,dims);
      }
      break;
    case 5:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST*****>(),value,dims,
                                 mask.get_view<const int*****>());
        else
          details::svm<use_mask>(get_strided_view<ST*****>(),value,dims,
                                 mask.get_view<const int*****>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST*****>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST*****>(),value,dims);
      }
      break;
    case 6:
      if constexpr (use_mask) {
        if (contig)
          details::svm<use_mask>(get_view<ST******>(),value,dims,
                                 mask.get_view<const int******>());
        else
          details::svm<use_mask>(get_strided_view<ST******>(),value,dims,
                                 mask.get_view<const int******>());
      } else {
        if (contig)
          details::svm<use_mask>(get_view<ST******>(),value,dims);
        else
          details::svm<use_mask>(get_strided_view<ST******>(),value,dims);
      }
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unsupported field rank in 'deep_copy'.\n");
  }
}

} // namespace scream
