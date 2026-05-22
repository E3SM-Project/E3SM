#include "vertical_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/field/field_tag.hpp"
#include "share/field/field_identifier.hpp"
#include "share/util/eamxx_timing.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"

#include <ekat_units.hpp>
#include <ekat_team_policy_utils.hpp>
#include <ekat_pack_utils.hpp>
#include <ekat_pack_kokkos.hpp>

#include <numeric>
#include <filesystem>

namespace scream
{

std::shared_ptr<AbstractGrid>
VerticalRemapper::
create_tgt_grid (const grid_ptr_type& src_grid,
                 const std::string& map_file)
{
  // Create tgt_grid as a clone of src_grid with different nlevs
  scorpio::register_file(map_file,scorpio::FileMode::Read);
  auto nlevs_tgt = scorpio::get_dimlen(map_file,"lev");

  auto tgt_grid = src_grid->clone(src_grid->name()+"_vremap_tgt",true);
  tgt_grid->reset_vertical_configuration(nlevs_tgt, AbstractGrid::VKind::Pressure);

  // Gather the pressure level data for vertical remapping
  using namespace ShortFieldTagsNames;
  auto layout = tgt_grid->get_vertical_layout(LEVP);
  // Add tgt pressure levels to the tgt grid
  auto p_tgt = tgt_grid->create_geometry_data<Real>("p_levs",layout,ekat::units::Pa,SCREAM_PACK_SIZE);
  scorpio::read_var(map_file,"p_levs",p_tgt.get_view<Real*,Host>().data());
  p_tgt.sync_to_dev();

  scorpio::release_file(map_file);

  return tgt_grid;
}

VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const std::string& map_file)
 : VerticalRemapper(src_grid,create_tgt_grid(src_grid,map_file))
{
  set_target_pressure (m_tgt_grid->get_geometry_data("p_levs"));
  std::filesystem::path p(map_file);

  set_name("VRemap " + p.filename().string());
}

VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid)
{
  set_name("VRemap " + tgt_grid->name());

  // We only go in one direction for simplicity, since we need to setup some
  // infrsatructures, and we don't want to setup 2x as many "just in case".
  // If you need to remap bwd, just create another remapper with src/tgt grids swapped.
  m_bwd_allowed = false;

  EKAT_REQUIRE_MSG (src_grid->get_2d_scalar_layout().congruent(tgt_grid->get_2d_scalar_layout()),
      "Error! Source and target grid can only differ for their number of level.\n");

  this->set_grids (src_grid,tgt_grid);
}

void VerticalRemapper::
set_extrapolation_type (const ExtrapType etype, const TopBot where)
{
  if (where & Top) {
    m_etype_top = etype;
  }
  if (where & Bot) {
    m_etype_bot = etype;
  }
}

void VerticalRemapper::
set_source_pressure (const Field& p)
{
  set_pressure (p, "source");
}

void VerticalRemapper::
set_target_pressure (const Field& p)
{
  set_pressure (p, "target");
}

void VerticalRemapper::
set_pressure (const Field& p, const std::string& src_or_tgt)
{
  using namespace ShortFieldTagsNames;
  using PackT = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  bool src = src_or_tgt=="source";
  auto grid = src ? m_src_grid : m_tgt_grid;

  std::string msg_prefix = "[VerticalRemapper::set_" + src_or_tgt + "_pressure] ";

  EKAT_REQUIRE_MSG(p.is_allocated(),
      msg_prefix + "Field is not yet allocated.\n"
      " - field name: " + p.name() + "\n");

  // Check layout:
  //  1. valid for this grid (e.g., no LEV/ILEV for a Pressure grid)
  //  2. either a Scalar1D (vert profile) or Scalar3D (horiz+vert dims)
  const auto& p_layout = p.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG(grid->is_valid_layout(p_layout) and
                   (p_layout.type()==LayoutType::Scalar1D or p_layout.type()==LayoutType::Scalar3D),
      "[VerticalRemapper::set_pressure] Unexpected/unsupported pressure layout.\n"
      " - pressure name  : " + p.name() + "\n"
      " - pressure layout: " + p_layout.to_string() + "\n");

  const auto vtag = p_layout.tags().back();
  if (src)
    m_src_pressure[vtag] = p;
  else
    m_tgt_pressure[vtag] = p;

  // If already in the map, do &=, otherwise set.
  // Equivalent to inserting 'true' (if missing) and doing &= after
  bool pack_compatible = p.get_header().get_alloc_properties().is_compatible<PackT>();
  auto res = m_packs_supported.emplace(vtag,true);
  res.first->second &= pack_compatible;
}

void VerticalRemapper::
registration_ends_impl ()
{
  using namespace ShortFieldTagsNames;
  using PackT = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  for (int i=0; i<m_num_fields; ++i) {
    const auto& src = m_src_fields[i];
          auto& tgt = m_tgt_fields[i];

    auto& ft = m_field2type[src.name()];

    const auto& src_layout = src.get_header().get_identifier().get_layout();
    const auto& tgt_layout = tgt.get_header().get_identifier().get_layout();

    ft.src_vtag = src.rank()>0 ? src_layout.tags().back() : INV;
    ft.tgt_vtag = tgt.rank()>0 ? tgt_layout.tags().back() : INV;
    if (ft.src_vtag==LEV or ft.src_vtag==ILEV or ft.src_vtag==LEVP) {
      // Sanity check: pressure fields MUST be set by now
      EKAT_REQUIRE_MSG (m_src_pressure.count(ft.src_vtag)>0,
          "[VerticalRemapper::registration_ends_impl] Error! Missing source pressure field.\n"
          " - src field name: " + src.name() + "\n"
          " - src layout    : " + src_layout.to_string() + "\n");
      EKAT_REQUIRE_MSG (m_tgt_pressure.count(ft.tgt_vtag)>0,
          "[VerticalRemapper::registration_ends_impl] Error! Missing target pressure field.\n"
          " - tgt field name: " + tgt.name() + "\n"
          " - tgt layout    : " + tgt_layout.to_string() + "\n");

      // Determine if this field can be handled with packs, and what's the vertical tag to use
      // for LinInterp lookup. If the src vtag is LEVP, we grab the vtag from the tgt layout.
      // This way, if the tags differ we are GUARANTEED we'll be using the grid with Model vkind.
      // Also, add mask tracking to the target field. The mask tracks tgt levels that are outside the
      // bounds of the src pressure field, and hence cannot be recovered by interpolation.

      ft.li_vtag = ft.src_vtag == LEVP ? ft.tgt_vtag : ft.src_vtag;

      // Packs are supported if both src and tgt fields do, and if src/tgt pressures also do
      ft.packs_supported = src.get_header().get_alloc_properties().is_compatible<PackT>() and
                           tgt.get_header().get_alloc_properties().is_compatible<PackT>() and
                           m_packs_supported[ft.src_vtag] and m_packs_supported[ft.tgt_vtag];

      if (m_etype_top==Mask or m_etype_bot==Mask) {
        // NOTE: for now we assume that masking is determined only by the COL,LEV location in space
        //       and that fields with multiple components will have the same masking for each component
        //       at a specific COL,LEV

        auto tgt_layout = create_tgt_layout(src_layout);

        // I this mask has already been created, retrieve it, otherwise create it
        // CAVEATS:
        //  1. We need different masks for fields with different vertical tags (LEV, ILEV, LEVP),
        //     so use src_layout's vtag to craft a unique mask name.
        //  2. for vector dimensions, we must include the vector dim length, as there may be
        //     2+ vector fields with different vector length, which need 2 different masks
        std::vector<std::string> tagdim_names;
        for (int i=0; i<src_layout.rank(); ++i) {
          tagdim_names.push_back(src_layout.names()[i]);
          if (src_layout.tags()[i]==CMP) {
            tagdim_names.back() += std::to_string(src_layout.dims()[i]);
          }
        }
        const auto mask_name = m_tgt_grid->name() + "_" + ekat::join(tagdim_names,"_") + "_mask";
        auto& mask = m_masks[mask_name];
        if (not mask.is_allocated()) {
          // Create this src/tgt mask fields, and assign them to these src/tgt fields extra data

          FieldIdentifier mask_fid (mask_name, tgt_layout, ekat::units::none, m_tgt_grid->name(), DataType::IntType );
          mask  = Field (mask_fid);
          if (ft.packs_supported)
            mask.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
          mask.allocate_view();
        }

        EKAT_REQUIRE_MSG(not tgt.has_valid_mask(),
            "[VerticalRemapper::registration_ends_impl] Error! Target field already has mask data assigned.\n"
            " - tgt field name: " + tgt.name() + "\n");

        tgt.set_valid_mask(mask);

        // Since we do mask (at top and/or bot), the tgt field MAY be contain fill_value entries
        tgt.get_header().set_may_be_filled(true);
      }
    } else {
      // If a field does not have any vertical tag (LEV, ILEV, or LEVP) it may still have
      // a mask assigned from somewhere else.
      // For instance, this could be a 2d field computed by FieldAtPressureLevel diagnostic.
      // In those cases we want to copy that mask to the target field.
      if (src.has_valid_mask()) {
        EKAT_REQUIRE_MSG(not tgt.has_valid_mask(),
            "[VerticalRemapper::registration_ends_impl] Error! Target field already has mask data assigned.\n"
            " - tgt field name: " + tgt.name() + "\n");
        auto src_mask = src.get_valid_mask();
        tgt.set_valid_mask(src_mask.alias(src_mask.name(),m_tgt_grid->name()));
      }

      // TODO: remove when we get rid of fill-aware Field manipulation methods
      if (src.get_header().may_be_filled()) {
        tgt.get_header().set_may_be_filled(true);
      }
    }
  }
  create_lin_interp ();
}

void VerticalRemapper::create_lin_interp()
{
  const int ncols = m_src_grid->get_num_local_dofs();
  for (auto& [name, ft] : m_field2type) {
    if (ft.li_vtag == FieldTag::Invalid)
      continue; // Does not have the vertical dimension

    EKAT_REQUIRE_MSG (m_src_pressure.count(ft.src_vtag)>0,
        "[VerticalRemapper::create_lin_interp] Error! Required src pressure field was not set.\n"
        " - source vert tag: " + e2str(ft.src_vtag) + "\n");
    EKAT_REQUIRE_MSG (m_tgt_pressure.count(ft.tgt_vtag)>0,
        "[VerticalRemapper::create_lin_interp] Error! Required tgt pressure field was not set.\n"
        " - target vert tag: " + e2str(ft.tgt_vtag) + "\n");

    const auto& src_p = m_src_pressure.at(ft.src_vtag);
    const auto& tgt_p = m_tgt_pressure.at(ft.tgt_vtag);

    const int src_nlevs = src_p.get_header().get_identifier().get_layout().dims().back();
    const int tgt_nlevs = tgt_p.get_header().get_identifier().get_layout().dims().back();
    if (ft.packs_supported) {
      m_lin_interp_packed.try_emplace(ft.li_vtag,ncols,src_nlevs,tgt_nlevs);
    } else {
      m_lin_interp_scalar.try_emplace(ft.li_vtag,ncols,src_nlevs,tgt_nlevs);
    }
  }
}

bool VerticalRemapper::
is_valid_tgt_layout (const FieldLayout& layout) const {
  using namespace ShortFieldTagsNames;
  const auto vkind = m_tgt_grid->get_vkind();
  const bool has_model_vtag    = layout.has_tag(LEV) or layout.has_tag(ILEV);
  const bool has_pressure_vtag = layout.has_tag(LEVP);
  if (vkind==AbstractGrid::VKind::Pressure and has_model_vtag) return false;
  if (vkind==AbstractGrid::VKind::Model    and has_pressure_vtag) return false;
  return AbstractRemapper::is_valid_tgt_layout(layout);
}

bool VerticalRemapper::
is_valid_src_layout (const FieldLayout& layout) const {
  using namespace ShortFieldTagsNames;
  const auto vkind = m_src_grid->get_vkind();
  const bool has_model_vtag    = layout.has_tag(LEV) or layout.has_tag(ILEV);
  const bool has_pressure_vtag = layout.has_tag(LEVP);
  if (vkind==AbstractGrid::VKind::Pressure and has_model_vtag) return false;
  if (vkind==AbstractGrid::VKind::Model    and has_pressure_vtag) return false;
  return AbstractRemapper::is_valid_src_layout(layout);
}

bool VerticalRemapper::
compatible_layouts (const FieldLayout& src,
                    const FieldLayout& tgt) const
{
  // Layouts are compatible if their non-vertical parts are congruent and their
  // vertical tags (if any) are compatible.
  // Rules for the vertical tag pair:
  //   - LEV <-> LEV or LEVP is allowed
  //   - ILEV <-> ILEV or LEVP is allowed
  //   - LEV <-> ILEV is NOT allowed (midpoints and interfaces are incompatible)

  using namespace ShortFieldTagsNames;

  auto src_vtag = src.rank()>0 ? src.tags().back() : INV;
  auto tgt_vtag = tgt.rank()>0 ? tgt.tags().back() : INV;

  auto src_stripped = src.clone().strip_dims({LEV,ILEV,LEVP});
  auto tgt_stripped = tgt.clone().strip_dims({LEV,ILEV,LEVP});

  return src.rank()==tgt.rank() and
         src_stripped.congruent(tgt_stripped) and
         (src_vtag==tgt_vtag or src_vtag==LEVP or tgt_vtag==LEVP);
}

FieldLayout VerticalRemapper::
create_layout (const FieldLayout& from_layout,
               const std::shared_ptr<const AbstractGrid>& to_grid) const
{
  using namespace ShortFieldTagsNames;

  auto from_grid = to_grid==m_src_grid ? m_tgt_grid : m_src_grid;

  auto check_vkind = [&]() {
    // If the from_grid is a Pressure grid, its layout uses LEVP.
    // We cannot map LEVP to LEV or ILEV without additional information.
    EKAT_REQUIRE_MSG (from_grid->get_vkind()!=AbstractGrid::VKind::Pressure or
                      to_grid->get_vkind()==AbstractGrid::VKind::Pressure,
        "[VerticalRemapper::create_layout] Error! Starting layout uses LEVP which cannot be mapped to LEV/ILEV.\n"
        "  - from grid: " + from_grid->name() + "\n"
        "  - to grid  : " + to_grid->name() + "\n");
  };

  auto to_layout = FieldLayout::invalid();
  const bool to_grid_is_pressure = to_grid->get_vkind()==AbstractGrid::VKind::Pressure;
  FieldTag vtag;
  std::string vdim_name;
  switch (from_layout.type()) {
    case LayoutType::Scalar0D: [[ fallthrough ]];
    case LayoutType::Vector0D: [[ fallthrough ]];
    case LayoutType::Scalar2D: [[ fallthrough ]];
    case LayoutType::Vector2D: [[ fallthrough ]];
    case LayoutType::Tensor2D:
      // These layouts do not have vertical dim tags, so no change
      to_layout = from_layout;
      break;
    case LayoutType::Scalar1D:
      vtag = to_grid_is_pressure ? LEVP
           : (from_layout.tags().back()==LEV ? LEV : ILEV);
      to_layout = to_grid->get_vertical_layout(vtag);
      break;
    case LayoutType::Scalar3D:
      check_vkind();
      vtag = to_grid_is_pressure ? LEVP
           : (from_layout.tags().back()==LEV ? LEV : ILEV);
      to_layout = to_grid->get_3d_scalar_layout(vtag);
      break;
    case LayoutType::Vector3D:
      check_vkind();
      vdim_name = from_layout.name(from_layout.get_vector_component_idx());
      vtag = to_grid_is_pressure ? LEVP
           : (from_layout.tags().back()==LEV ? LEV : ILEV);
      to_layout = to_grid->get_3d_vector_layout(vtag,from_layout.get_vector_dim(),vdim_name);
      break;
    default:
      // NOTE: this also include Tensor3D. We don't really have any atm proc
      //       that needs to handle a tensor3d quantity, so no need to add it
      EKAT_ERROR_MSG (
        "[VerticalRemapper] Error! Layout not supported by VerticalRemapper.\n"
        " - input layout: " + from_layout.to_string() + "\n");
  }
  return to_layout;
}

void VerticalRemapper::remap_fwd_impl ()
{
  using namespace ShortFieldTagsNames;

  // 1. Setup any interp object that was created (if nullptr, no fields need it)
  if (m_timers_enabled)
    start_timer(name() + " setup LI");

  bool src_grid_levp = m_src_grid->get_vkind()==AbstractGrid::VKind::Pressure;
  bool tgt_grid_levp = m_tgt_grid->get_vkind()==AbstractGrid::VKind::Pressure;

  // For a Pressure grid, there is a single pressure (keyed LEVP) used for all fields.
  // For a Model grid, the pressure key matches the field's vtag (LEV or ILEV).
  auto src_pressure = [&](FieldTag vtag) -> const Field& {
    return src_grid_levp ? m_src_pressure.at(LEVP) : m_src_pressure.at(vtag);
  };
  auto tgt_pressure = [&](FieldTag vtag) -> const Field& {
    return tgt_grid_levp ? m_tgt_pressure.at(LEVP) : m_tgt_pressure.at(vtag);
  };

  for (auto& [vtag, li] : m_lin_interp_packed) {
    setup_lin_interp(li, src_pressure(vtag), tgt_pressure(vtag));
  }
  for (auto& [vtag, li] : m_lin_interp_scalar) {
    setup_lin_interp(li, src_pressure(vtag), tgt_pressure(vtag));
  }
  if (m_timers_enabled)
    stop_timer(name() + " setup LI");

  // 2. Init all masks fields (if any) to 1 (signaling no masked entries)
  for (auto& [name, mask] : m_masks) {
    mask.deep_copy(1);
  }

  // 3. Interpolate (and extrapolate) the fields
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src    = m_src_fields[i];
          auto& f_tgt    = m_tgt_fields[i];
    const auto& type = m_field2type.at(f_src.name());
    if (type.li_vtag!=INV) {
      const auto& x_src = src_pressure(type.src_vtag);
      const auto& x_tgt = tgt_pressure(type.tgt_vtag);
      if (type.packs_supported) {
        apply_vertical_interpolation(m_lin_interp_packed.at(type.li_vtag),f_src,f_tgt,x_src,x_tgt);
      } else {
        apply_vertical_interpolation(m_lin_interp_scalar.at(type.li_vtag),f_src,f_tgt,x_src,x_tgt);
      }
      extrapolate(f_src,f_tgt,x_src,x_tgt);
    } else {
      // There is nothing to do, this field does not need vertical interpolation,
      // so just copy it over.  Note, if this field has its own mask data make
      // sure that is copied too.
      f_tgt.deep_copy(f_src);
      if (f_tgt.has_valid_mask()) {
        auto& f_tgt_mask = f_tgt.get_valid_mask();
        auto& f_src_mask = f_src.get_valid_mask();
        f_tgt_mask.deep_copy(f_src_mask);
      }
    }
  }
}

template<int Packsize>
void VerticalRemapper::
setup_lin_interp (const ekat::LinInterp<Real,Packsize>& lin_interp,
                  const Field& x_src, const Field& x_tgt) const
{
  using LI_t   = ekat::LinInterp<Real,Packsize>;
  using TPF    = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;
  using PackT  = ekat::Pack<Real,Packsize>;
  using view2d = typename KokkosTypes<DefaultDevice>::view<const PackT**>;
  using view1d = typename KokkosTypes<DefaultDevice>::view<const PackT*>;

  auto src1d = x_src.rank()==1;
  auto tgt1d = x_tgt.rank()==1;

  view2d p_src2d_v, p_tgt2d_v;
  view1d p_src1d_v, p_tgt1d_v;
  if (src1d) {
    p_src1d_v = x_src.get_view<const PackT*>();
  } else {
    p_src2d_v = x_src.get_view<const PackT**>();
  }
  if (tgt1d) {
    p_tgt1d_v = x_tgt.get_view<const PackT*>();
  } else {
    p_tgt2d_v = x_tgt.get_view<const PackT**>();
  }

  auto lambda = KOKKOS_LAMBDA(typename LI_t::MemberType const& team) {
    const int icol = team.league_rank();
    // Extract subviews if src/tgt were not 1d to start with
    auto x_src = p_src1d_v;
    if (not src1d)
      x_src = ekat::subview(p_src2d_v,icol);
    auto x_tgt = p_tgt1d_v;
    if (not tgt1d)
      x_tgt = ekat::subview(p_tgt2d_v,icol);

    lin_interp.setup(team,x_src,x_tgt);
  };
  const int ncols = m_src_grid->get_num_local_dofs();
  const int nlevs_tgt = x_tgt.get_header().get_identifier().get_layout().dims().back();
  const int npacks_tgt = ekat::PackInfo<Packsize>::num_packs(nlevs_tgt);
  auto policy = TPF::get_default_team_policy(ncols,npacks_tgt);
  Kokkos::parallel_for("VerticalRemapper::interp_setup",policy,lambda);
  Kokkos::fence();
}

template<int Packsize>
void VerticalRemapper::
apply_vertical_interpolation(const ekat::LinInterp<Real,Packsize>& lin_interp,
                             const Field& f_src, const Field& f_tgt,
                             const Field& p_src, const Field& p_tgt) const
{
  if (m_timers_enabled)
    start_timer(name() + " run LI");

  // Note: if Packsize==1, we grab packs of size 1, which are for sure
  //       compatible with the allocation
  using LI_t   = ekat::LinInterp<Real,Packsize>;
  using PackT  = ekat::Pack<Real,Packsize>;
  using TPF    = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  using view2d = typename KokkosTypes<DefaultDevice>::view<const PackT**>;
  using view1d = typename KokkosTypes<DefaultDevice>::view<const PackT*>;

  auto src1d = p_src.rank()==1;
  auto tgt1d = p_tgt.rank()==1;

  view2d p_src2d_v, p_tgt2d_v;
  view1d p_src1d_v, p_tgt1d_v;
  if (src1d) {
    p_src1d_v = p_src.get_view<const PackT*>();
  } else {
    p_src2d_v = p_src.get_view<const PackT**>();
  }
  if (tgt1d) {
    p_tgt1d_v = p_tgt.get_view<const PackT*>();
  } else {
    p_tgt2d_v = p_tgt.get_view<const PackT**>();
  }

  const auto& f_tgt_l = f_tgt.get_header().get_identifier().get_layout();
  const int ncols = m_src_grid->get_num_local_dofs();
  const int nlevs_tgt = f_tgt_l.dims().back();
  const int npacks_tgt = ekat::PackInfo<Packsize>::num_packs(nlevs_tgt);

  switch(f_src.rank()) {
    case 2:
    {
      auto f_src_v = f_src.get_view<const PackT**>();
      auto f_tgt_v = f_tgt.get_view<      PackT**>();
      auto policy = TPF::get_default_team_policy(ncols,npacks_tgt);
      auto lambda = KOKKOS_LAMBDA(typename LI_t::MemberType const& team)
      {
        const int icol = team.league_rank();

        // Extract subviews if src/tgt pressures were not 1d to start with
        auto x_src = p_src1d_v;
        auto x_tgt = p_tgt1d_v;
        if (not src1d)
          x_src = ekat::subview(p_src2d_v,icol);
        if (not tgt1d)
          x_tgt = ekat::subview(p_tgt2d_v,icol);

        auto y_src = ekat::subview(f_src_v,icol);
        auto y_tgt = ekat::subview(f_tgt_v,icol);
        lin_interp.lin_interp(team,x_src,x_tgt,y_src,y_tgt,icol);
      };
      Kokkos::parallel_for("VerticalRemapper::apply_vertical_interpolation",policy,lambda);
      break;
    }
    case 3:
    {
      auto f_src_v = f_src.get_view<const PackT***>();
      auto f_tgt_v = f_tgt.get_view<      PackT***>();
      const int ncomps = f_tgt_l.get_vector_dim();
      auto policy = TPF::get_default_team_policy(ncols*ncomps,npacks_tgt);

      auto lambda = KOKKOS_LAMBDA(typename LI_t::MemberType const& team)
      {
        const int icol = team.league_rank() / ncomps;
        const int icmp = team.league_rank() % ncomps;

        // Extract subviews if src/tgt pressures were not 1d to start with
        auto x_src = p_src1d_v;
        auto x_tgt = p_tgt1d_v;
        if (not src1d)
          x_src = ekat::subview(p_src2d_v,icol);
        if (not tgt1d)
          x_tgt = ekat::subview(p_tgt2d_v,icol);

        auto y_src = ekat::subview(f_src_v,icol,icmp);
        auto y_tgt = ekat::subview(f_tgt_v,icol,icmp);
        lin_interp.lin_interp(team,x_src,x_tgt,y_src,y_tgt,icol);
      };
      Kokkos::parallel_for("VerticalRemapper::apply_vertical_interpolation",policy,lambda);
      break;
    }
    default:
      EKAT_ERROR_MSG (
          "[VerticalRemapper::apply_vertical_interpolation] Error! Unsupported field rank.\n"
          " - src field name: " + f_src.name() + "\n"
          " - src field rank: " + std::to_string(f_src.rank()) + "\n");
  }
  if (m_timers_enabled)
    stop_timer(name() + " run LI");
}

void VerticalRemapper::
extrapolate (const Field& f_src,
             const Field& f_tgt,
             const Field& p_src,
             const Field& p_tgt) const
{
  if (m_timers_enabled)
    start_timer(name() + " extrapolate");

  using TPF = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  using view2d = typename KokkosTypes<DefaultDevice>::view<const Real**>;
  using view1d = typename KokkosTypes<DefaultDevice>::view<const Real*>;

  constexpr auto fill_val = constants::fill_value<Real>;

  auto src1d = p_src.rank()==1;
  auto tgt1d = p_tgt.rank()==1;

  view2d p_src2d_v, p_tgt2d_v;
  view1d p_src1d_v, p_tgt1d_v;
  if (src1d) {
    p_src1d_v = p_src.get_view<const Real*>();
  } else {
    p_src2d_v = p_src.get_view<const Real**>();
  }
  if (tgt1d) {
    p_tgt1d_v = p_tgt.get_view<const Real*>();
  } else {
    p_tgt2d_v = p_tgt.get_view<const Real**>();
  }

  const auto& f_tgt_l = f_tgt.get_header().get_identifier().get_layout();
  const auto& f_src_l = f_src.get_header().get_identifier().get_layout();
  const int ncols = m_src_grid->get_num_local_dofs();
  const int nlevs_tgt = f_tgt_l.dims().back();
  const int nlevs_src = f_src_l.dims().back();

  auto etop = m_etype_top;
  auto ebot = m_etype_bot;
  auto mid = nlevs_tgt / 2;
  auto do_mask = etop==Mask or ebot==Mask;

  switch(f_src.rank()) {
    case 2:
    {
      auto f_src_v = f_src.get_view<const Real**>();
      auto f_tgt_v = f_tgt.get_view<      Real**>();
      auto mask_v = do_mask ? f_tgt.get_valid_mask().get_view<int**>()
                            : typename Field::view_dev_t<int**>{};

      auto policy = TPF::get_default_team_policy(ncols,nlevs_tgt);

      using MemberType = typename decltype(policy)::member_type;
      auto lambda = KOKKOS_LAMBDA(const MemberType& team)
      {
        const int icol = team.league_rank();

        // Extract subviews if src/tgt pressures were not 1d to start with
        auto x_src = p_src1d_v;
        auto x_tgt = p_tgt1d_v;
        if (not src1d)
          x_src = ekat::subview(p_src2d_v,icol);
        if (not tgt1d)
          x_tgt = ekat::subview(p_tgt2d_v,icol);

        auto y_src = ekat::subview(f_src_v,icol);
        auto y_tgt = ekat::subview(f_tgt_v,icol);

        auto x_min = x_src[0];
        auto x_max = x_src[nlevs_src-1];
        auto extrapolate = [&](const int ilev) {
          if (ilev>=mid) {
            // Near surface
            if (x_tgt[ilev]>x_max) {
              if (ebot==P0) {
                y_tgt[ilev] = y_src[nlevs_src-1];
              } else {
                y_tgt[ilev] = fill_val;
                mask_v(icol,ilev) = 0;
              }
            }
          } else {
            // Near top
            if (x_tgt[ilev]<x_min) {
              if (etop==P0) {
                y_tgt[ilev] = y_src[0];
              } else {
                y_tgt[ilev] = fill_val;
                mask_v(icol,ilev) = 0;
              }
            }
          }
        };
        Kokkos::parallel_for (Kokkos::TeamVectorRange(team,nlevs_tgt), extrapolate);
      };
      Kokkos::parallel_for("VerticalRemapper::extrapolate",policy,lambda);
      break;
    }
    case 3:
    {
      auto f_src_v = f_src.get_view<const Real***>();
      auto f_tgt_v = f_tgt.get_view<      Real***>();
      auto mask_v = do_mask ? f_tgt.get_valid_mask().get_view<int***>()
                            : typename Field::view_dev_t<int***>{};
      const int ncomps = f_tgt_l.get_vector_dim();
      auto policy = TPF::get_default_team_policy(ncols*ncomps,nlevs_tgt);

      using MemberType = typename decltype(policy)::member_type;
      auto lambda = KOKKOS_LAMBDA(const MemberType& team)
      {
        const int icol = team.league_rank() / ncomps;
        const int icmp = team.league_rank() % ncomps;

        // Extract subviews if src/tgt pressures were not 1d to start with
        auto x_src = p_src1d_v;
        auto x_tgt = p_tgt1d_v;
        if (not src1d)
          x_src = ekat::subview(p_src2d_v,icol);
        if (not tgt1d)
          x_tgt = ekat::subview(p_tgt2d_v,icol);

        auto y_src = ekat::subview(f_src_v,icol,icmp);
        auto y_tgt = ekat::subview(f_tgt_v,icol,icmp);
        auto x_min = x_src[0];
        auto x_max = x_src[nlevs_src-1];
        auto extrapolate = [&](const int ilev) {
          if (ilev>=mid) {
            // Near surface
            if (x_tgt[ilev]>x_max) {
              if (ebot==P0) {
                y_tgt[ilev] = y_src[nlevs_src-1];
              } else {
                y_tgt[ilev] = fill_val;
                mask_v(icol,icmp,ilev) = 0;
              }
            }
          } else {
            // Near top
            if (x_tgt[ilev]<x_min) {
              if (etop==P0) {
                y_tgt[ilev] = y_src[0];
              } else {
                y_tgt[ilev] = fill_val;
                mask_v(icol,icmp,ilev) = 0;
              }
            }
          }
        };
        Kokkos::parallel_for (Kokkos::TeamVectorRange(team,nlevs_tgt), extrapolate);
      };
      Kokkos::parallel_for("VerticalRemapper::extrapolate",policy,lambda);
      break;
    }
    default:
      EKAT_ERROR_MSG (
          "[VerticalRemapper::extrapolate] Error! Unsupported field rank.\n"
          " - src field name: " + f_src.name() + "\n"
          " - src field rank: " + std::to_string(f_src.rank()) + "\n");
  }
  if (m_timers_enabled)
    stop_timer(name() + " extrapolate");
}

} // namespace scream
