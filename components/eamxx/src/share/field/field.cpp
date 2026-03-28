#include "share/field/field.hpp"
#include "share/util/eamxx_utils.hpp"

namespace scream
{

namespace
{
void update_checks (const std::string& caller,
                    const Field& y, const Field& x,
                    const ScalarWrapper& alpha,
                    const ScalarWrapper& beta,
                    const Field* mask = nullptr)
{
  // Check output field is writable
  EKAT_REQUIRE_MSG (not y.is_read_only(),
      "[" + caller + "] Error! Cannot modify field, as it is read-only.\n"
      " - field name: " + y.name() + "\n");

  const auto& y_dt = y.data_type();
  const auto& x_dt = x.data_type();
  const auto& a_dt = alpha.type;
  const auto& b_dt = beta.type;

  // Ensure data types are not somehow invalid
  EKAT_REQUIRE_MSG (y_dt!=DataType::Invalid,
      "[" + caller + "] Error! Lhs data type is invalid.\n"
      " - lhs name: " + y.name() + "\n");
  EKAT_REQUIRE_MSG (x_dt!=DataType::Invalid,
      "[" + caller + "] Error! Rhs data type is invalid.\n"
      " - rhs name: " + x.name() + "\n");
  EKAT_REQUIRE_MSG (a_dt!=DataType::Invalid,
      "[" + caller + "] Error! Alpha coeff data type is invalid.\n"
      " - lhs name: " + y.name() + "\n");
  EKAT_REQUIRE_MSG (b_dt!=DataType::Invalid,
      "[" + caller + "] Error! Beta coeff data type is invalid.\n"
      " - lhs name: " + y.name() + "\n");

  // If user passes, say, double alpha/beta for an int field, we should error out, warning about
  // a potential narrowing rounding. The other way around, otoh, is allowed (even though
  // there's an upper limit to the int values that a double can store, it is unlikely the user
  // will use such large factors).
  // Similarly, we allow updating a field Y with another X as long as converting the data type of X
  // to the data type of Y does not require narrowing
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(x_dt,y_dt),
      "[" + caller + "] Error! Rhs data type may be narrowed when converted to lhs data type.\n"
      " - lhs name: " + y.name() + "\n"
      " - rhs name: " + x.name() + "\n"
      " - rhs data type: " + e2str(x_dt) + "\n"
      " - lhs data type: " + e2str(y_dt) + "\n");
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(a_dt,y_dt),
      "[" + caller + "] Error! Coefficient alpha may be narrowed when converted to lhs data type.\n"
      " - lhs name: " + y.name() + "\n"
      " - lhs data type  : " + e2str(y_dt) + "\n"
      " - alpha data type: " + e2str(a_dt) + "\n");
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(b_dt,y_dt),
      "[" + caller + "] Error! Coefficient beta may be narrowed when converted to lhs data type.\n"
      " - lhs name: " + y.name() + "\n"
      " - lhs data type  : " + e2str(y_dt) + "\n"
      " - beta data type: " + e2str(b_dt) + "\n");

  // Check x/y are allocated
  EKAT_REQUIRE_MSG (y.is_allocated(),
      "[" + caller + "] Error! Lhs field is not yet allocated.\n"
      " - field name: " + y.name() + "\n");
  EKAT_REQUIRE_MSG (x.is_allocated(),
      "[" + caller + "] Error! Rhs field is not yet allocated.\n"
      " - field name: " + x.name() + "\n");

  const auto& y_l = y.get_header().get_identifier().get_layout();
  const auto& x_l = x.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG (y_l==x_l,
      "[" + caller + "] Error! Incompatible fields layouts.\n"
      " - rhs name: " + x.name() + "\n"
      " - lhs name: " + y.name() + "\n"
      " - rhs layout: " + x_l.to_string() + "\n"
      " - lhs layout: " + y_l.to_string() + "\n");

  // Now check mask
  if (mask) {
    EKAT_REQUIRE_MSG (mask->is_allocated(),
        "[" + caller + "] Error! Mask field was not yet allocated.\n"
        " - mask name: " + mask->name() + "\n");

    EKAT_REQUIRE_MSG (mask->data_type()==DataType::IntType,
        "[" + caller + "] Error! Mask data type MUST be IntType.\n"
        " - mask name: " + mask->name() + "\n"
        " - mask data type: " + e2str(mask->data_type()) + "\n");

    const auto& m_l = mask->get_header().get_identifier().get_layout();
    EKAT_REQUIRE_MSG (y_l.congruent(m_l),
      "[" + caller + "] Error! Incompatible mask layout.\n"
      " - lhs name: " + y.name() + "\n"
      " - mask name: " + mask->name() + "\n"
      " - lhs layout: " + y_l.to_string() + "\n"
      " - mask layout: " + m_l.to_string() + "\n");
  }
}
} // anonymous namespace


Field::
Field ()
 : Field (identifier_type("UNSET",
                          FieldLayout::invalid(),
                          ekat::units::Units::invalid(),
                          "UNKNOWN",
                          DataType::Invalid))
{
  // Nothing to do here
}

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
  auto fid = my_fid.clone(name).reset_grid(grid_name);
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

  if (m_header->has_extra_data("valid_mask")) {
    const auto& mask = m_header->get_extra_data<Field>("valid_mask");
    const auto& mfid = mask.get_header().get_identifier();
    sf.m_header->set_extra_data("valid_mask",mask.subfield(mask.name(),mfid.get_units(),idim,index,dynamic));
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

void Field::deep_copy (const ScalarWrapper value)
{
  // Check consistency of inputs
  update_checks("Field::deep_copy (scalar)",*this,*this,value,value);

  const auto my_data_type = data_type();
  switch (my_data_type) {
    case DataType::IntType:
      deep_copy_impl(value.as<int>());
      break;
    case DataType::FloatType:
      deep_copy_impl(value.as<float>());
      break;
    case DataType::DoubleType:
      deep_copy_impl(value.as<double>());
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized field data type in Field::deep_copy.\n");
  }
}

void Field::deep_copy (const ScalarWrapper value, const Field& mask, const bool negate_mask)
{
  update_checks("Field::deep_copy (scalar, masked)",*this,*this,value,value,&mask);

  const auto my_data_type = data_type();
  switch (my_data_type) {
    case DataType::IntType:
      negate_mask ? deep_copy_masked<true>(value.as<int>(),mask)
                  : deep_copy_masked<false>(value.as<int>(),mask);
      break;
    case DataType::FloatType:
      negate_mask ? deep_copy_masked<true>(value.as<float>(),mask)
                  : deep_copy_masked<false>(value.as<float>(),mask);
      break;
    case DataType::DoubleType:
      negate_mask ? deep_copy_masked<true>(value.as<double>(),mask)
                  : deep_copy_masked<false>(value.as<double>(),mask);
      break;
    default:
      EKAT_ERROR_MSG ("Error! Unrecognized field data type in Field::deep_copy.\n");
  }
}

void Field::deep_copy (const Field& x)
{
  constexpr auto CM = CombineMode::Replace;
  update_cm<CM>("Field::deep_copy",x,1,0);
}

void Field::deep_copy (const Field& x, const Field& mask)
{
  constexpr auto CM = CombineMode::Replace;
  update_cm<CM>("Field::deep_copy (masked)",x,1,0,mask);
}

void Field::scale (const ScalarWrapper beta)
{
  constexpr auto CM = CombineMode::Update;
  update_cm<CM>("Field::scale (scalar)",*this,0,beta);
}

void Field::scale (const ScalarWrapper beta, const Field& mask)
{
  constexpr auto CM = CombineMode::Update;
  update_cm<CM>("Field::scale (scalar, masked)",*this,0,beta,mask);
}

void Field::scale (const Field& x)
{
  constexpr auto CM = CombineMode::Multiply;
  update_cm<CM>("Field::scale",x,1,1);
}

void Field::scale (const Field& x, const Field& mask)
{
  constexpr auto CM = CombineMode::Multiply;
  update_cm<CM>("Field::scale (masked)",x,1,1,mask);
}

void Field::scale_inv (const Field& x)
{
  constexpr auto CM = CombineMode::Divide;
  update_cm<CM>("Field::scale_inv",x,1,1);
}

void Field::scale_inv (const Field& x, const Field& mask)
{
  constexpr auto CM = CombineMode::Divide;
  update_cm<CM>("Field::scale_inv (masked)",x,1,1,mask);
}

void Field::max (const Field& x)
{
  constexpr auto CM = CombineMode::Max;
  update_cm<CM>("Field::max",x,1,1);
}

void Field::max (const Field& x, const Field& mask)
{
  constexpr auto CM = CombineMode::Max;
  update_cm<CM>("Field::max (masked)",x,1,1,mask);
}

void Field::min (const Field& x)
{
  constexpr auto CM = CombineMode::Min;
  update_cm<CM>("Field::min",x,1,1);
}

void Field::min (const Field& x, const Field& mask)
{
  constexpr auto CM = CombineMode::Min;
  update_cm<CM>("Field::min (masked)",x,1,1,mask);
}

void Field::
update (const Field& x, const ScalarWrapper alpha, const ScalarWrapper beta)
{
  constexpr auto CM = CombineMode::Update;
  update_cm<CM>("Field::update",x,alpha,beta);
}

void Field::
update (const Field& x, const ScalarWrapper alpha, const ScalarWrapper beta, const Field& mask)
{
  constexpr auto CM = CombineMode::Update;
  update_cm<CM>("Field::update (masked)",x,alpha,beta,mask);
}

template<CombineMode CM>
void Field::
update_cm (const std::string& caller, const Field& x, const ScalarWrapper alpha, const ScalarWrapper beta)
{
  // Determine if the RHS can contain fill_value entries
  // Check consistency of inputs
  update_checks(caller,*this,x,alpha,beta);

  bool fill_aware = x.get_header().may_be_filled();

  if (data_type()==DataType::IntType) {
    return fill_aware ? update_fill_aware<CM,int,int>(x,alpha.as<int>(),beta.as<int>())
                      : update_impl<CM,int,int>(x,alpha.as<int>(),beta.as<int>());
  } else if (data_type()==DataType::FloatType) {
    if (x.data_type()==DataType::FloatType)
      return fill_aware ? update_fill_aware<CM,float,float>(x,alpha.as<float>(),beta.as<float>())
                        : update_impl<CM,float,float>(x,alpha.as<float>(),beta.as<float>());
    else
      return fill_aware ? update_fill_aware<CM,float,int>(x,alpha.as<float>(),beta.as<float>())
                        : update_impl<CM,float,int>(x,alpha.as<float>(),beta.as<float>());
  } else if (data_type()==DataType::DoubleType) {
    if (x.data_type()==DataType::DoubleType)
      return fill_aware ? update_fill_aware<CM,double,double>(x,alpha.as<double>(),beta.as<double>())
                        : update_impl<CM,double,double>(x,alpha.as<double>(),beta.as<double>());
    else if (x.data_type()==DataType::FloatType)
      return fill_aware ? update_fill_aware<CM,double,float>(x,alpha.as<double>(),beta.as<double>())
                        : update_impl<CM,double,float>(x,alpha.as<double>(),beta.as<double>());
    else
      return fill_aware ? update_fill_aware<CM,double,int>(x,alpha.as<double>(),beta.as<double>())
                        : update_impl<CM,double,int>(x,alpha.as<double>(),beta.as<double>());
  }
}

template<CombineMode CM>
void Field::
update_cm (const std::string& caller, const Field& x, const ScalarWrapper alpha, const ScalarWrapper beta, const Field& mask)
{
  // Determine if the RHS can contain fill_value entries
  // Check consistency of inputs
  update_checks(caller,*this,x,alpha,beta,&mask);

  if (data_type()==DataType::IntType) {
    return update_masked<CM,int,int>(x,alpha.as<int>(),beta.as<int>(),mask);
  } else if (data_type()==DataType::FloatType) {
    if (x.data_type()==DataType::FloatType)
      return update_masked<CM,float,float>(x,alpha.as<float>(),beta.as<float>(),mask);
    else
      return update_masked<CM,float,int>(x,alpha.as<float>(),beta.as<float>(),mask);
  } else if (data_type()==DataType::DoubleType) {
    if (x.data_type()==DataType::DoubleType)
      return update_masked<CM,double,double>(x,alpha.as<double>(),beta.as<double>(),mask);
    else if (x.data_type()==DataType::FloatType)
      return update_masked<CM,double,float>(x,alpha.as<double>(),beta.as<double>(),mask);
    else
      return update_masked<CM,double,int>(x,alpha.as<double>(),beta.as<double>(),mask);
  }
}

} // namespace scream
