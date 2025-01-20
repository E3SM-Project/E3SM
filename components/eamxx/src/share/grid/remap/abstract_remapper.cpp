#include "share/grid/remap/abstract_remapper.hpp"

namespace scream
{

AbstractRemapper::
AbstractRemapper (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid)
{
  set_grids (src_grid,tgt_grid);
}

void AbstractRemapper::
registration_begins () {
  EKAT_REQUIRE_MSG(m_state==RepoState::Clean,
      "Error! Cannot start registration on a non-clean repo.\n"
      "       Did you call 'registration_begins' already?\n");

  m_state = RepoState::Open;
}

void AbstractRemapper::
register_field (const Field& src, const Field& tgt)
{
  EKAT_REQUIRE_MSG(m_state==RepoState::Open,
      "Error! Cannot register fields in the remapper at this time.\n"
      "       Did you forget to call 'registration_begins' or called 'registeration_ends' already?");

  EKAT_REQUIRE_MSG(src.is_allocated(), "Error! Source field is not yet allocated.\n");
  EKAT_REQUIRE_MSG(tgt.is_allocated(), "Error! Target field is not yet allocated.\n");

  const auto& src_layout = src.get_header().get_identifier().get_layout();
  const auto& tgt_layout = tgt.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG(is_valid_src_layout(src_layout),
      "Error! Source field has an invalid layout.\n"
      " - field name  : " + src.name() + "\n"
      " - field layout: " + src_layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(is_valid_tgt_layout(tgt_layout),
      "Error! Source field has an invalid layout.\n"
      " - field name  : " + tgt.name() + "\n"
      " - field layout: " + tgt_layout.to_string() + "\n");
  EKAT_REQUIRE_MSG(compatible_layouts(src_layout,tgt_layout),
      "Error! Source and target layouts are not compatible.\n"
      " - src name: " + src.name() + "\n"
      " - tgt name: " + tgt.name() + "\n"
      " - src layout: " + src_layout.to_string() + "\n"
      " - tgt layout: " + tgt_layout.to_string() + "\n");

  if (src.is_read_only()) {
    m_has_read_only_src_fields = true;
  }
  if (tgt.is_read_only()) {
    m_has_read_only_tgt_fields = true;
  }

  m_src_fields.emplace_back(src);
  m_tgt_fields.emplace_back(tgt);

  ++m_num_fields;
}

void AbstractRemapper::
register_field_from_src (const Field& src) {
  const auto& src_fid = src.get_header().get_identifier();
  const auto& tgt_fid = create_tgt_fid(src_fid);

  Field tgt(tgt_fid);
  const auto& src_ap = src.get_header().get_alloc_properties();
        auto& tgt_ap = tgt.get_header().get_alloc_properties();
  tgt_ap.request_allocation(src_ap.get_largest_pack_size());
  tgt.allocate_view();

  register_field(src,tgt);
}

void AbstractRemapper::
register_field_from_tgt (const Field& tgt) {
  const auto& tgt_fid = tgt.get_header().get_identifier();
  const auto& src_fid = create_src_fid(tgt_fid);

  Field src(src_fid);
  const auto& tgt_ap = tgt.get_header().get_alloc_properties();
        auto& src_ap = src.get_header().get_alloc_properties();
  src_ap.request_allocation(tgt_ap.get_largest_pack_size());
  src.allocate_view();

  register_field(src,tgt);
}

void AbstractRemapper::registration_ends ()
{
  EKAT_REQUIRE_MSG(m_state!=RepoState::Closed,
      "Error! Cannot call registration_ends at this time.\n"
      "       Did you accidentally call 'registration_ends' already?");

  // Call derived class impl first. They may register extra/internal fields,
  // so we must keep the repo OPEN until they're done.
  registration_ends_impl();

  m_state = RepoState::Closed;
}

void AbstractRemapper::remap_fwd ()
{
  EKAT_REQUIRE_MSG(m_state!=RepoState::Open,
      "Error! Cannot perform remapping at this time.\n"
      "       Did you forget to call 'registration_ends'?\n");

  if (m_state!=RepoState::Clean) {
    EKAT_REQUIRE_MSG (m_fwd_allowed,
        "Error! Forward remap is not allowed by this remapper.\n");
    EKAT_REQUIRE_MSG (not m_has_read_only_tgt_fields,
        "Error! Forward remap IS allowed by this remapper, but some of the tgt fields are read-only\n");
    remap_fwd_impl ();
  }
}

void AbstractRemapper::remap_bwd ()
{
  EKAT_REQUIRE_MSG(m_state!=RepoState::Open,
      "Error! Cannot perform remapping at this time.\n"
      "       Did you forget to call 'registration_ends'?\n");

  if (m_state!=RepoState::Clean) {
    EKAT_REQUIRE_MSG (m_bwd_allowed,
        "Error! Backward remap is not allowed by this remapper.\n");
    EKAT_REQUIRE_MSG (not m_has_read_only_src_fields,
        "Error! Backward remap IS allowed by this remapper, but some of the src fields are read-only\n");
    remap_bwd_impl ();
  }
}

void AbstractRemapper::
set_grids (const grid_ptr_type& src_grid,
           const grid_ptr_type& tgt_grid)
{
  EKAT_REQUIRE_MSG (src_grid!=nullptr, "Error! Invalid source grid pointer.\n");
  EKAT_REQUIRE_MSG (tgt_grid!=nullptr, "Error! Invalid target grid pointer.\n");

  m_src_grid = src_grid;
  m_tgt_grid = tgt_grid;
}

// Returns the source field for the given field index.
const Field& AbstractRemapper::get_src_field (const int i) const
{
  EKAT_REQUIRE_MSG(i>=0 && i<m_num_fields, "Error! Field index out of bounds.\n");
  return m_src_fields[i];
}

// Returns the target field for the given field index.
const Field& AbstractRemapper::get_tgt_field (const int i) const
{
  EKAT_REQUIRE_MSG (i>=0 && i<m_num_fields, "Error! Field index out of bounds.\n");
  return m_tgt_fields[i];
}

FieldIdentifier AbstractRemapper::create_src_fid (const FieldIdentifier& tgt_fid) const
{
  const auto& name = tgt_fid.name();
  const auto& layout = create_src_layout(tgt_fid.get_layout());
  const auto& units = tgt_fid.get_units();

  return FieldIdentifier(name,layout,units,m_src_grid->name());
}

FieldIdentifier AbstractRemapper::create_tgt_fid (const FieldIdentifier& src_fid) const
{
  const auto& name = src_fid.name();
  const auto& layout = create_tgt_layout(src_fid.get_layout());
  const auto& units = src_fid.get_units();

  return FieldIdentifier(name,layout,units,m_tgt_grid->name());
}

FieldLayout AbstractRemapper::
create_src_layout (const FieldLayout& tgt_layout) const
{
  EKAT_REQUIRE_MSG (m_src_grid!=nullptr,
      "Error! Cannot create source layout until the source grid has been set.\n");

  EKAT_REQUIRE_MSG (is_valid_tgt_layout(tgt_layout),
      "[HorizInterpRemapperBase] Error! Input target layout is not valid for this remapper.\n"
      " - input layout: " + tgt_layout.to_string());

  return create_layout(tgt_layout,m_src_grid);
}

FieldLayout AbstractRemapper::
create_tgt_layout (const FieldLayout& src_layout) const
{
  EKAT_REQUIRE_MSG (m_tgt_grid!=nullptr,
      "Error! Cannot create target layout until the target grid has been set.\n");

  EKAT_REQUIRE_MSG (is_valid_src_layout(src_layout),
      "[HorizInterpRemapperBase] Error! Input source layout is not valid for this remapper.\n"
      " - input layout: " + src_layout.to_string());

  return create_layout(src_layout,m_tgt_grid);
}

bool AbstractRemapper::
compatible_layouts (const FieldLayout& src, const FieldLayout& tgt) const
{
  if (src.type()!=tgt.type())
    return false;

  if (src.is_vector_layout()) {
    return src.get_vector_dim()==tgt.get_vector_dim();
  } else if (src.is_tensor_layout()) {
    return src.get_tensor_dims()==tgt.get_tensor_dims();
  }

  return true;
}

FieldLayout AbstractRemapper::
create_layout (const FieldLayout& from_layout,
               const grid_ptr_type& to_grid) const
{
  using namespace ShortFieldTagsNames;


  const bool midpoints = from_layout.has_tag(LEV);

  auto fl_out = FieldLayout::invalid();
  switch (from_layout.type()) {
    case LayoutType::Scalar0D: [[ fallthrough ]]; 
    case LayoutType::Vector0D: [[ fallthrough ]]; 
    case LayoutType::Tensor0D:
      // 0d layouts are the same on all grids
      fl_out = from_layout;
      break;
    case LayoutType::Scalar1D:
      // 1d layouts require the grid correct number of levs
      fl_out = to_grid->get_vertical_layout(midpoints);
      break;
    case LayoutType::Vector1D:
    {
      auto vdim_idx  = from_layout.get_vector_component_idx();
      auto vdim_name = from_layout.names()[vdim_idx];
      auto vdim_len  = from_layout.get_vector_dim();
      fl_out = to_grid->get_vertical_layout(midpoints,vdim_len,vdim_name);
      break;
    }
    case LayoutType::Scalar2D:
      fl_out = to_grid->get_2d_scalar_layout();
      break;
    case LayoutType::Vector2D:
    {
      auto vdim_idx  = from_layout.get_vector_component_idx();
      auto vdim_name = from_layout.names()[vdim_idx];
      auto vdim_len  = from_layout.get_vector_dim();
      fl_out = to_grid->get_2d_vector_layout(vdim_len,vdim_name);
      break;
    }
    case LayoutType::Tensor2D:
    {
      std::vector<std::string> tdims_names;
      for (auto idx  : from_layout.get_tensor_components_ids()) {
        tdims_names.push_back(from_layout.names()[idx]);
      }   
      fl_out = to_grid->get_2d_tensor_layout(from_layout.get_tensor_dims(),tdims_names);
      break;
    }
    case LayoutType::Scalar3D:
      fl_out = to_grid->get_3d_scalar_layout(midpoints);
      break;
    case LayoutType::Vector3D:
    {
      auto vdim_idx  = from_layout.get_vector_component_idx();
      auto vdim_name = from_layout.names()[vdim_idx];
      auto vdim_len  = from_layout.get_vector_dim();
      fl_out = to_grid->get_3d_vector_layout(midpoints,vdim_len,vdim_name);
      break;
    }
    case LayoutType::Tensor3D:
    {
      std::vector<std::string> tdims_names;
      for (auto idx  : from_layout.get_tensor_components_ids()) {
        tdims_names.push_back(from_layout.names()[idx]);
      }   
      fl_out = to_grid->get_3d_tensor_layout(midpoints,from_layout.get_tensor_dims(),tdims_names);
    }
    default:
      EKAT_ERROR_MSG ("Layout not supported by this remapper.\n"
                      " - layout: " + from_layout.to_string() + "\n");
  }
  return fl_out;
}

} // namespace scream
