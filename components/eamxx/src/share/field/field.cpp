#include "share/field/field.hpp"
#include "share/util/eamxx_utils.hpp"
#include "share/field/eamxx_fill_aware_binary_op_expression.hpp"

#include <ekat_expression_eval.hpp>
#include <ekat_expression_binary_op.hpp>
#include <ekat_expression_view.hpp>

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
  f.deep_copy<Device>(*this);
  f.deep_copy<Host>(*this);

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

template<HostOrDevice HD>
void Field::scale (const ScalarWrapper& beta)
{
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(beta.type,data_type()),
      "[Field::scale] Error! Coefficient beta may be narrowed when converted to the field data type.\n"
      " - field name     : " + name() + "\n"
      " - field data type: " + e2str(data_type()) + "\n"
      " - beta data type : " + e2str(beta.type) + "\n");

  auto run = [&] (auto coeff) {
    using ST = decltype(coeff);
    auto do_scale = [&](auto view) {
      auto ve = ekat::view_expression(view);
      auto expr = ve*coeff;
      ekat::evaluate(expr,view);
    };
    switch(rank()) {
      case 0: do_scale(get_view<ST>()); break;
      case 1: do_scale(get_view<ST*>()); break;
      case 2: do_scale(get_view<ST**>()); break;
      case 3: do_scale(get_view<ST***>()); break;
      case 4: do_scale(get_view<ST****>()); break;
      case 5: do_scale(get_view<ST*****>()); break;
      case 6: do_scale(get_view<ST******>()); break;
      default:
        EKAT_ERROR_MSG (
          "[Field::scale] Error! Unsupported field rank.\n"
          " - field name: " + name() + "\n"
          " - field rank: " + std::to_string(rank()) + "\n");
    };
  };

  switch (data_type()) {
    case DataType::IntType:
      run (beta.as<int>());
      break;
    case DataType::FloatType:
      run (beta.as<float>());
      break;
    case DataType::DoubleType:
      run (beta.as<double>());
      break;
    default:
      EKAT_ERROR_MSG (
          "[Field::scale] Error! Unrecognized field data type.\n"
          " - field name: " + name() + "\n");
  }
}

template<HostOrDevice HD>
void Field::scale (const Field& x)
{
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(x.data_type(),data_type()),
      "[Field::scale] Error! Coefficient field scalars may be narrowed when converted to the field data type.\n"
      " - field name           : " + name() + "\n"
      " - coeff field name     : " + x.name() + "\n"
      " - field data type      : " + e2str(data_type()) + "\n"
      " - coeff field data type: " + e2str(x.data_type()) + "\n");

  auto run = [&] (auto y_scalar, auto x_scalar) {
    using YST = decltype(y_scalar);
    using XST = decltype(x_scalar);
    auto do_scale = [&](auto yview, auto xview) {
      auto ye = ekat::view_expression(yview);
      auto xe = ekat::view_expression(xview);
      auto expr = ye*xe;
      ekat::evaluate(expr,yview);
    };
    switch(rank()) {
      case 0: do_scale(get_view<YST>(),       x.get_view<const XST>()      ); break;
      case 1: do_scale(get_view<YST*>(),      x.get_view<const XST*>()     ); break;
      case 2: do_scale(get_view<YST**>(),     x.get_view<const XST**>()    ); break;
      case 3: do_scale(get_view<YST***>(),    x.get_view<const XST***>()   ); break;
      case 4: do_scale(get_view<YST****>(),   x.get_view<const XST****>()  ); break;
      case 5: do_scale(get_view<YST*****>(),  x.get_view<const XST*****>() ); break;
      case 6: do_scale(get_view<YST******>(), x.get_view<const XST******>()); break;
      default:
        EKAT_ERROR_MSG (
          "[Field::scale] Error! Unsupported field rank.\n"
          " - field name: " + name() + "\n"
          " - field rank: " + std::to_string(rank()) + "\n");
    };
  };

  switch (data_type()) {
    case DataType::IntType:
      run(int(0),int(0));
      break;
    case DataType::FloatType:
      if (x.data_type()==DataType::FloatType)
        run(float(0),float(0));
      else
        run(float(0),int(0));
      break;
    case DataType::DoubleType:
      if (x.data_type()==DataType::DoubleType)
        run(double(0),double(0));
      else if (x.data_type()==DataType::FloatType)
        run(double(0),float(0));
      else
        run(double(0),int(0));
      break;
    default:
      EKAT_ERROR_MSG (
          "[Field::scale] Error! Unrecognized field data type.\n"
          " - field name: " + name() + "\n");
  }
}

template<HostOrDevice HD>
void Field::scale_inv (const Field& x)
{
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(x.data_type(),data_type()),
      "[Field::scale_inv] Error! Coefficient field scalars may be narrowed when converted to the field data type.\n"
      " - field name           : " + name() + "\n"
      " - coeff field name     : " + x.name() + "\n"
      " - field data type      : " + e2str(data_type()) + "\n"
      " - coeff field data type: " + e2str(x.data_type()) + "\n");

  auto run = [&] (auto y_scalar, auto x_scalar) {
    using YST = decltype(y_scalar);
    using XST = decltype(x_scalar);
    auto do_scale = [&](auto yview, auto xview) {
      auto ye = ekat::view_expression(yview);
      auto xe = ekat::view_expression(xview);
      auto expr = ye/xe;
      ekat::evaluate(expr,yview);
    };
    switch(rank()) {
      case 0: do_scale(get_view<YST>(),       x.get_view<const XST>()      ); break;
      case 1: do_scale(get_view<YST*>(),      x.get_view<const XST*>()     ); break;
      case 2: do_scale(get_view<YST**>(),     x.get_view<const XST**>()    ); break;
      case 3: do_scale(get_view<YST***>(),    x.get_view<const XST***>()   ); break;
      case 4: do_scale(get_view<YST****>(),   x.get_view<const XST****>()  ); break;
      case 5: do_scale(get_view<YST*****>(),  x.get_view<const XST*****>() ); break;
      case 6: do_scale(get_view<YST******>(), x.get_view<const XST******>()); break;
      default:
        EKAT_ERROR_MSG (
          "[Field::scale] Error! Unsupported field rank.\n"
          " - field name: " + name() + "\n"
          " - field rank: " + std::to_string(rank()) + "\n");
    };
  };

  switch (data_type()) {
    case DataType::IntType:
      run(int(0),int(0));
      break;
    case DataType::FloatType:
      if (x.data_type()==DataType::FloatType)
        run(float(0),float(0));
      else
        run(float(0),int(0));
      break;
    case DataType::DoubleType:
      if (x.data_type()==DataType::DoubleType)
        run(double(0),double(0));
      else if (x.data_type()==DataType::FloatType)
        run(double(0),float(0));
      else
        run(double(0),int(0));
      break;
    default:
      EKAT_ERROR_MSG (
          "[Field::scale] Error! Unrecognized field data type.\n"
          " - field name: " + name() + "\n");
  }
}

template<HostOrDevice HD>
void Field::max (const Field& x)
{
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(x.data_type(),data_type()),
      "[Field::max] Error! RHS field scalars may be narrowed when converted to the field data type.\n"
      " - field name           : " + name() + "\n"
      " - coeff field name     : " + x.name() + "\n"
      " - field data type      : " + e2str(data_type()) + "\n"
      " - coeff field data type: " + e2str(x.data_type()) + "\n");

  bool fill_aware = get_header().may_be_filled() or x.get_header().may_be_filled();

  auto run = [&] (auto y_scalar, auto x_scalar) {
    using YST = decltype(y_scalar);
    using XST = decltype(x_scalar);
    auto do_scale = [&](auto yview, auto xview) {
      auto ye = ekat::view_expression(yview);
      auto xe = ekat::view_expression(xview);
      if (fill_aware) {
        auto expr = ekat::fa_max(ye,xe);
        ekat::evaluate(expr,yview);
      } else {
        auto expr = ekat::max(ye,xe);
        ekat::evaluate(expr,yview);
      }
    };
    switch(rank()) {
      case 0: do_scale(get_view<YST>(),       x.get_view<const XST>()      ); break;
      case 1: do_scale(get_view<YST*>(),      x.get_view<const XST*>()     ); break;
      case 2: do_scale(get_view<YST**>(),     x.get_view<const XST**>()    ); break;
      case 3: do_scale(get_view<YST***>(),    x.get_view<const XST***>()   ); break;
      case 4: do_scale(get_view<YST****>(),   x.get_view<const XST****>()  ); break;
      case 5: do_scale(get_view<YST*****>(),  x.get_view<const XST*****>() ); break;
      case 6: do_scale(get_view<YST******>(), x.get_view<const XST******>()); break;
      default:
        EKAT_ERROR_MSG (
          "[Field::max] Error! Unsupported field rank.\n"
          " - field name: " + name() + "\n"
          " - field rank: " + std::to_string(rank()) + "\n");
    };
  };

  switch (data_type()) {
    case DataType::IntType:
      run(int(0),int(0));
      break;
    case DataType::FloatType:
      if (x.data_type()==DataType::FloatType)
        run(float(0),float(0));
      else
        run(float(0),int(0));
      break;
    case DataType::DoubleType:
      if (x.data_type()==DataType::DoubleType)
        run(double(0),double(0));
      else if (x.data_type()==DataType::FloatType)
        run(double(0),float(0));
      else
        run(double(0),int(0));
      break;
    default:
      EKAT_ERROR_MSG (
          "[Field::max] Error! Unrecognized field data type.\n"
          " - field name: " + name() + "\n");
  }
}

template<HostOrDevice HD>
void Field::min (const Field& x)
{
  EKAT_REQUIRE_MSG (not is_narrowing_conversion(x.data_type(),data_type()),
      "[Field::min] Error! RHS field scalars may be narrowed when converted to the field data type.\n"
      " - field name           : " + name() + "\n"
      " - coeff field name     : " + x.name() + "\n"
      " - field data type      : " + e2str(data_type()) + "\n"
      " - coeff field data type: " + e2str(x.data_type()) + "\n");

  bool fill_aware = get_header().may_be_filled() or x.get_header().may_be_filled();

  auto run = [&] (auto y_scalar, auto x_scalar) {
    using YST = decltype(y_scalar);
    using XST = decltype(x_scalar);
    auto do_scale = [&](auto yview, auto xview) {
      auto ye = ekat::view_expression(yview);
      auto xe = ekat::view_expression(xview);
      if (fill_aware) {
        auto expr = ekat::fa_min(ye,xe);
        ekat::evaluate(expr,yview);
      } else {
        auto expr = ekat::min(ye,xe);
        ekat::evaluate(expr,yview);
      }
    };
    switch(rank()) {
      case 0: do_scale(get_view<YST>(),       x.get_view<const XST>()      ); break;
      case 1: do_scale(get_view<YST*>(),      x.get_view<const XST*>()     ); break;
      case 2: do_scale(get_view<YST**>(),     x.get_view<const XST**>()    ); break;
      case 3: do_scale(get_view<YST***>(),    x.get_view<const XST***>()   ); break;
      case 4: do_scale(get_view<YST****>(),   x.get_view<const XST****>()  ); break;
      case 5: do_scale(get_view<YST*****>(),  x.get_view<const XST*****>() ); break;
      case 6: do_scale(get_view<YST******>(), x.get_view<const XST******>()); break;
      default:
        EKAT_ERROR_MSG (
          "[Field::min] Error! Unsupported field rank.\n"
          " - field name: " + name() + "\n"
          " - field rank: " + std::to_string(rank()) + "\n");
    };
  };

  switch (data_type()) {
    case DataType::IntType:
      run(int(0),int(0));
      break;
    case DataType::FloatType:
      if (x.data_type()==DataType::FloatType)
        run(float(0),float(0));
      else
        run(float(0),int(0));
      break;
    case DataType::DoubleType:
      if (x.data_type()==DataType::DoubleType)
        run(double(0),double(0));
      else if (x.data_type()==DataType::FloatType)
        run(double(0),float(0));
      else
        run(double(0),int(0));
      break;
    default:
      EKAT_ERROR_MSG (
          "[Field::min] Error! Unrecognized field data type.\n"
          " - field name: " + name() + "\n");
  }
}

template void Field::scale<Device> (const ScalarWrapper&);
template void Field::scale<Device> (const Field&);
template void Field::scale_inv<Device> (const Field&);
template void Field::max<Device> (const Field&);
template void Field::min<Device> (const Field&);
#ifdef EAMXX_ENABLE_GPU
template void Field::scale<Host>   (const ScalarWrapper&);
template void Field::scale<Host> (const Field&);
template void Field::scale_inv<Host> (const Field&);
template void Field::max<Host> (const Field&);
template void Field::min<Host> (const Field&);
#endif

} // namespace scream
