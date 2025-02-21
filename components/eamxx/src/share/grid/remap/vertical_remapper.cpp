#include "vertical_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/field/field_tag.hpp"
#include "share/field/field_identifier.hpp"
#include "share/util/eamxx_universal_constants.hpp"
#include "share/io/eamxx_scorpio_interface.hpp"

#include <ekat/util/ekat_units.hpp>
#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>
#include <ekat/ekat_pack_kokkos.hpp>

#include <numeric>

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

  auto tgt_grid = src_grid->clone("vertical_remap_tgt_grid",true);
  tgt_grid->reset_num_vertical_lev(nlevs_tgt);

  // Gather the pressure level data for vertical remapping
  auto layout = tgt_grid->get_vertical_layout(true);
  Field p_tgt(FieldIdentifier("p_levs",layout,ekat::units::Pa,tgt_grid->name()));
  p_tgt.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  p_tgt.allocate_view();
  scorpio::read_var(map_file,"p_levs",p_tgt.get_view<Real*,Host>().data());
  p_tgt.sync_to_dev();

  // Add tgt pressure levels to the tgt grid
  tgt_grid->set_geometry_data(p_tgt);

  scorpio::release_file(map_file);

  return tgt_grid;
}

VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const std::string& map_file)
 : VerticalRemapper(src_grid,create_tgt_grid(src_grid,map_file))
{
  set_target_pressure (m_tgt_grid->get_geometry_data("p_levs"),Both);
}

VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const grid_ptr_type& tgt_grid)
{
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
set_mask_value (const Real mask_val)
{
  EKAT_REQUIRE_MSG (not ekat::is_invalid(mask_val),
      "[VerticalRemapper::set_mask_value] Error! Input mask value must be a valid number.\n");

  m_mask_val = mask_val;
}

void VerticalRemapper::
set_source_pressure (const Field& p, const ProfileType ptype)
{
  set_pressure (p, "source", ptype);
}

void VerticalRemapper::
set_target_pressure (const Field& p, const ProfileType ptype)
{
  set_pressure (p, "target", ptype);
}

void VerticalRemapper::
set_pressure (const Field& p, const std::string& src_or_tgt, const ProfileType ptype)
{
  using namespace ShortFieldTagsNames;
  using PackT = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  bool src = src_or_tgt=="source";

  std::string msg_prefix = "[VerticalRemapper::set_" + src_or_tgt + "_pressure] ";

  EKAT_REQUIRE_MSG(p.is_allocated(),
      msg_prefix + "Field is not yet allocated.\n"
      " - field name: " + p.name() + "\n");

  EKAT_REQUIRE_MSG(p.get_header().get_alloc_properties().is_compatible<PackT>(),
      msg_prefix + "Field not compatible with default pack size.\n"
      " - pack size: " + std::to_string(SCREAM_PACK_SIZE) + "\n");

  const int nlevs = src ? m_src_grid->get_num_vertical_levels()
                        : m_tgt_grid->get_num_vertical_levels();
  const auto& p_layout = p.get_header().get_identifier().get_layout();
  const auto vtag = p_layout.tags().back();
  const auto vdim = p_layout.dims().back();

  FieldTag expected_tag;
  int      expected_dim;
  switch (ptype) {
    case Midpoints:
      expected_tag = LEV;
      expected_dim = nlevs;
      if (src) {
        m_src_pmid = p;
      } else {
        m_tgt_pmid = p;
      }
      break;
    case Interfaces:
      expected_tag = ILEV;
      expected_dim = nlevs+1;
      if (src) {
        m_src_pint = p;
      } else {
        m_tgt_pint = p;
      }
      break;
    case Both:
      expected_tag = LEV;
      expected_dim = nlevs;
      if (src) {
        m_src_pint = p;
        m_src_pmid = p;
        m_src_int_same_as_mid = true;
      } else {
        m_tgt_pint = p;
        m_tgt_pmid = p;
        m_tgt_int_same_as_mid = true;
      }
      break;
    default:
      EKAT_ERROR_MSG ("[VerticalRemapper::set_source_pressure] Error! Unrecognized value for 'ptype'.\n");
  }
  EKAT_REQUIRE_MSG (vtag==expected_tag and vdim==expected_dim,
      msg_prefix + "Invalid pressure layout.\n"
      "  - layout: " + p_layout.to_string() + "\n"
      "  - expected last layout tag: " + e2str(expected_tag) + "\n"
      "  - expected last layout dim: " + std::to_string(expected_dim) + "\n");
}

void VerticalRemapper::
registration_ends_impl ()
{
  using namespace ShortFieldTagsNames;
  using PackT = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  for (int i=0; i<m_num_fields; ++i) {
    const auto& src = m_src_fields[i];
          auto& tgt = m_tgt_fields[i];

    // Clone src layout, since we may strip dims later for mask creation
    auto src_layout = src.get_header().get_identifier().get_layout().clone();

    if (src_layout.has_tag(LEV) or src_layout.has_tag(ILEV)) {
      // Determine if this field can be handled with packs, and whether it's at midpoints
      // NOTE: we don't know if mid==int on src or tgt. If it is, we use the other to determine mid-vs-int
      // Add mask tracking to the target field. The mask tracks location of tgt pressure levs that are outside the
      // bounds of the src pressure field, and hence cannot be recovered by interpolation
      auto& ft = m_field2type[src.name()];
      ft.midpoints = m_src_int_same_as_mid
                   ? tgt.get_header().get_identifier().get_layout().has_tag(LEV)
                   : src.get_header().get_identifier().get_layout().has_tag(LEV);
      ft.packed    = src.get_header().get_alloc_properties().is_compatible<PackT>() and
                     tgt.get_header().get_alloc_properties().is_compatible<PackT>();

      if (m_etype_top==Mask or m_etype_bot==Mask) {
        // NOTE: for now we assume that masking is determined only by the COL,LEV location in space
        //       and that fields with multiple components will have the same masking for each component
        //       at a specific COL,LEV
        src_layout.strip_dims({CMP});

        // I this mask has already been created, retrieve it, otherwise create it
        const auto mask_name = m_tgt_grid->name() + "_" + ekat::join(src_layout.names(),"_") + "_mask";
        Field tgt_mask;
        if (m_field2type.count(mask_name)==0) {
          auto nondim = ekat::units::Units::nondimensional();
          // Create this src/tgt mask fields, and assign them to these src/tgt fields extra data

          FieldIdentifier src_mask_fid (mask_name, src_layout, nondim, m_src_grid->name() );
          FieldIdentifier tgt_mask_fid = create_tgt_fid(src_mask_fid);

          Field src_mask (src_mask_fid);
          src_mask.allocate_view();

          tgt_mask  = Field (tgt_mask_fid);
          tgt_mask.allocate_view();

          // Initialize the src mask values to 1.0
          src_mask.deep_copy(1.0);

          m_src_masks.push_back(src_mask);
          m_tgt_masks.push_back(tgt_mask);

          auto& mt = m_field2type[src_mask_fid.name()];
          mt.packed = false;
          mt.midpoints = src_layout.has_tag(LEV);
        } else {
          for (size_t i=0; i<m_tgt_masks.size(); ++i) {
            if (m_tgt_masks[i].name()==mask_name) {
              tgt_mask = m_tgt_masks[i];
              break;
            }
          }
        }

        EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_data"),
            "[VerticalRemapper::registration_ends_impl] Error! Target field already has mask data assigned.\n"
            " - tgt field name: " + tgt.name() + "\n");
        EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_value"),
            "[VerticalRemapper::registration_ends_impl] Error! Target field already has mask value assigned.\n"
            " - tgt field name: " + tgt.name() + "\n");

        tgt.get_header().set_extra_data("mask_data",tgt_mask);
        tgt.get_header().set_extra_data("mask_value",m_mask_val);
      }
    } else {
      // If a field does not have LEV or ILEV it may still have mask tracking assigned from somewhere else.
      // For instance, this could be a 2d field computed by FieldAtPressureLevel diagnostic.
      // In those cases we want to copy that mask tracking to the target field.
      if (src.get_header().has_extra_data("mask_data")) {
        EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_data"),
            "[VerticalRemapper::registration_ends_impl] Error! Target field already has mask data assigned.\n"
            " - tgt field name: " + tgt.name() + "\n");
        auto src_mask = src.get_header().get_extra_data<Field>("mask_data");
        tgt.get_header().set_extra_data("mask_data",src_mask);
      }
      if (src.get_header().has_extra_data("mask_value")) {
        EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_value"),
            "[VerticalRemapper::registration_ends_impl] Error! Target field already has mask value assigned.\n"
            " - tgt field name: " + tgt.name() + "\n");
        auto src_mask_val = src.get_header().get_extra_data<Real>("mask_value");
        tgt.get_header().set_extra_data("mask_value",src_mask_val);
      }
    }
  }
  create_lin_interp ();
}

void VerticalRemapper::create_lin_interp()
{
  // Count number fields for each packed-midpoints value
  auto beg = m_field2type.begin();
  auto end = m_field2type.end();

  int num_packed_mid =
    std::count_if(beg,end,[](const std::pair<std::string,FType>& it)
        {
          return it.second.midpoints and it.second.packed;
        });
  int num_packed_int =
    std::count_if(beg,end,[](const std::pair<std::string,FType>& it)
        {
          return not it.second.midpoints and it.second.packed;
        });
  int num_scalar_mid =
    std::count_if(beg,end,[](const std::pair<std::string,FType>& it)
        {
          return it.second.midpoints and not it.second.packed;
        });
  int num_scalar_int =
    std::count_if(beg,end,[](const std::pair<std::string,FType>& it)
        {
          return not it.second.midpoints and not it.second.packed;
        });

  // Create the linear interpolation object
  const auto ncols     = m_src_grid->get_num_local_dofs();
  const auto nlevs_src = m_src_grid->get_num_vertical_levels();
  const auto nlevs_tgt = m_tgt_grid->get_num_vertical_levels();

  if (num_packed_mid>0) {
    m_lin_interp_mid_packed =
      std::make_shared<ekat::LinInterp<Real,SCREAM_PACK_SIZE>>(ncols,nlevs_src,nlevs_tgt);
  }
  if (num_scalar_mid>0) {
    m_lin_interp_mid_scalar =
      std::make_shared<ekat::LinInterp<Real,1>>(ncols,nlevs_src,nlevs_tgt);
  }
  if (num_packed_int>0) {
    m_lin_interp_int_packed =
      std::make_shared<ekat::LinInterp<Real,SCREAM_PACK_SIZE>>(ncols,nlevs_src,nlevs_tgt);
  }
  if (num_scalar_int>0) {
    m_lin_interp_int_scalar =
      std::make_shared<ekat::LinInterp<Real,1>>(ncols,nlevs_src,nlevs_tgt);
  }
}

bool VerticalRemapper::
is_valid_tgt_layout (const FieldLayout& layout) const {
  using namespace ShortFieldTagsNames;
  return !(m_tgt_int_same_as_mid and layout.has_tag(ILEV))
         and AbstractRemapper::is_valid_tgt_layout(layout);
}

bool VerticalRemapper::
is_valid_src_layout (const FieldLayout& layout) const {
  using namespace ShortFieldTagsNames;
  return !(m_src_int_same_as_mid and layout.has_tag(ILEV))
         and AbstractRemapper::is_valid_src_layout(layout);
}

bool VerticalRemapper::
compatible_layouts (const FieldLayout& src,
                    const FieldLayout& tgt) const {
  // Strip the LEV/ILEV tags, and check if they are the same
  // Also, check rank compatibility, in case one has LEV/ILEV and the other doesn't
  // NOTE: tgt layouts always use LEV (not ILEV), while src can have ILEV or LEV.

  using namespace ShortFieldTagsNames;
  auto src_stripped = src.clone().strip_dims({LEV,ILEV});
  auto tgt_stripped = tgt.clone().strip_dims({LEV,ILEV});

  return src.rank()==tgt.rank() and
         src_stripped.congruent(tgt_stripped);
}

FieldLayout VerticalRemapper::
create_layout (const FieldLayout& from_layout,
               const std::shared_ptr<const AbstractGrid>& to_grid) const
{
  using namespace ShortFieldTagsNames;

  // Detect if for the output grid we distinguish between midpoints and interfaces or not
  // If we don't distinguish, we just use the LEV tag (for layout with the vertical dim)
  auto from_grid = to_grid==m_src_grid ? m_tgt_grid : m_src_grid;
  bool output_int_same_as_mid = to_grid==m_src_grid ? m_src_int_same_as_mid : m_tgt_int_same_as_mid;
  bool input_int_same_as_mid  = from_grid==m_src_grid ? m_src_int_same_as_mid : m_tgt_int_same_as_mid;

  // If the input layout does not distinguish between LEV/ILEV, we cannot deduce the output layout
  EKAT_REQUIRE_MSG (not input_int_same_as_mid,
      "[VerticalRemapper::create_layout] Error! Starting layout does not distinguish between LEV and ILEV.\n"
      "  - from grid: " + from_grid->name() + "\n"
      "  - to grid  : " + to_grid->name() + "\n");

  auto to_layout = FieldLayout::invalid();
  bool midpoints;
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
      midpoints = output_int_same_as_mid || from_layout.tags().back()==LEV;
      to_layout = to_grid->get_vertical_layout(midpoints);
      break;
    case LayoutType::Scalar3D:
      midpoints = output_int_same_as_mid || from_layout.tags().back()==LEV;
      to_layout = to_grid->get_3d_scalar_layout(midpoints);
      break;
    case LayoutType::Vector3D:
      vdim_name = from_layout.name(from_layout.get_vector_component_idx());
      midpoints = output_int_same_as_mid || from_layout.tags().back()==LEV;
      to_layout = to_grid->get_3d_vector_layout(midpoints,from_layout.get_vector_dim(),vdim_name);
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
  // 1. Setup any interp object that was created (if nullptr, no fields need it)
  if (m_lin_interp_mid_packed) {
    setup_lin_interp(*m_lin_interp_mid_packed,m_src_pmid,m_tgt_pmid);
  }
  if (m_lin_interp_int_packed) {
    setup_lin_interp(*m_lin_interp_int_packed,m_src_pint,m_tgt_pint);
  }
  if (m_lin_interp_mid_scalar) {
    setup_lin_interp(*m_lin_interp_mid_scalar,m_src_pmid,m_tgt_pmid);
  }
  if (m_lin_interp_int_scalar) {
    setup_lin_interp(*m_lin_interp_int_scalar,m_src_pint,m_tgt_pint);
  }

  using namespace ShortFieldTagsNames;

  // 2. Interpolate the fields
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src    = m_src_fields[i];
          auto& f_tgt    = m_tgt_fields[i];
    const auto& tgt_layout   = f_tgt.get_header().get_identifier().get_layout();
    if (tgt_layout.has_tag(LEV) or tgt_layout.has_tag(ILEV)) {
      const auto& type = m_field2type.at(f_src.name());
      // Dispatch interpolation to the proper lin interp object
      if (type.midpoints) {
        if (type.packed) {
          apply_vertical_interpolation(*m_lin_interp_mid_packed,f_src,f_tgt,m_src_pmid,m_tgt_pmid);
        } else {
          apply_vertical_interpolation(*m_lin_interp_mid_scalar,f_src,f_tgt,m_src_pmid,m_tgt_pmid);
        }
        extrapolate(f_src,f_tgt,m_src_pmid,m_tgt_pmid,m_mask_val);
      } else {
        if (type.packed) {
          apply_vertical_interpolation(*m_lin_interp_int_packed,f_src,f_tgt,m_src_pint,m_tgt_pint);
        } else {
          apply_vertical_interpolation(*m_lin_interp_int_scalar,f_src,f_tgt,m_src_pint,m_tgt_pint);
        }
        extrapolate(f_src,f_tgt,m_src_pint,m_tgt_pint,m_mask_val);
      }
    } else {
      // There is nothing to do, this field does not need vertical interpolation,
      // so just copy it over.  Note, if this field has its own mask data make
      // sure that is copied too.
      f_tgt.deep_copy(f_src);
      if (f_tgt.get_header().has_extra_data("mask_data")) {
        auto f_tgt_mask = f_tgt.get_header().get_extra_data<Field>("mask_data");
        auto f_src_mask = f_src.get_header().get_extra_data<Field>("mask_data");
        f_tgt_mask.deep_copy(f_src_mask);
      }
    }
  }

  // 3. Interpolate the mask fields
  for (unsigned i=0; i<m_tgt_masks.size(); ++i) {
          auto& f_src = m_src_masks[i];
          auto& f_tgt = m_tgt_masks[i];
    const auto& type = m_field2type.at(f_src.name());

    // Dispatch interpolation to the proper lin interp object
    if (type.midpoints) {
      if (type.packed) {
        apply_vertical_interpolation(*m_lin_interp_mid_packed,f_src,f_tgt,m_src_pmid,m_tgt_pmid);
      } else {
        apply_vertical_interpolation(*m_lin_interp_mid_scalar,f_src,f_tgt,m_src_pmid,m_tgt_pmid);
      }
      extrapolate(f_src,f_tgt,m_src_pmid,m_tgt_pmid,0);
    } else {
      if (type.packed) {
        apply_vertical_interpolation(*m_lin_interp_int_packed,f_src,f_tgt,m_src_pint,m_tgt_pint);
      } else {
        apply_vertical_interpolation(*m_lin_interp_int_scalar,f_src,f_tgt,m_src_pint,m_tgt_pint);
      }
      extrapolate(f_src,f_tgt,m_src_pint,m_tgt_pint,0);
    }
  }
}

template<int Packsize>
void VerticalRemapper::
setup_lin_interp (const ekat::LinInterp<Real,Packsize>& lin_interp,
                  const Field& p_src, const Field& p_tgt) const
{
  using LI_t = ekat::LinInterp<Real,Packsize>;
  using ESU = ekat::ExeSpaceUtils<DefaultDevice::execution_space>;
  using PackT = ekat::Pack<Real,Packsize>;
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
  const int nlevs_tgt = m_tgt_grid->get_num_vertical_levels();
  const int npacks_tgt = ekat::PackInfo<Packsize>::num_packs(nlevs_tgt);
  auto policy = ESU::get_default_team_policy(ncols,npacks_tgt);
  Kokkos::parallel_for("VerticalRemapper::interp_setup",policy,lambda);
  Kokkos::fence();
}

template<int Packsize>
void VerticalRemapper::
apply_vertical_interpolation(const ekat::LinInterp<Real,Packsize>& lin_interp,
                             const Field& f_src, const Field& f_tgt,
                             const Field& p_src, const Field& p_tgt) const
{
  // Note: if Packsize==1, we grab packs of size 1, which are for sure
  //       compatible with the allocation
  using LI_t = ekat::LinInterp<Real,Packsize>;
  using PackT = ekat::Pack<Real,Packsize>;
  using ESU = ekat::ExeSpaceUtils<DefaultDevice::execution_space>;

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
      auto policy = ESU::get_default_team_policy(ncols,npacks_tgt);
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
      auto policy = ESU::get_default_team_policy(ncols*ncomps,npacks_tgt);

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
}

void VerticalRemapper::
extrapolate (const Field& f_src,
             const Field& f_tgt,
             const Field& p_src,
             const Field& p_tgt,
             const Real mask_val) const
{
  using ESU = ekat::ExeSpaceUtils<DefaultDevice::execution_space>;

  using view2d = typename KokkosTypes<DefaultDevice>::view<const Real**>;
  using view1d = typename KokkosTypes<DefaultDevice>::view<const Real*>;

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
  switch(f_src.rank()) {
    case 2:
    {
      auto f_src_v = f_src.get_view<const Real**>();
      auto f_tgt_v = f_tgt.get_view<      Real**>();
      auto policy = ESU::get_default_team_policy(ncols,nlevs_tgt);

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
                y_tgt[ilev] = mask_val;
              }
            }
          } else {
            // Near top
            if (x_tgt[ilev]<x_min) {
              if (etop==P0) {
                y_tgt[ilev] = y_src[0];
              } else {
                y_tgt[ilev] = mask_val;
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
      const int ncomps = f_tgt_l.get_vector_dim();
      auto policy = ESU::get_default_team_policy(ncols*ncomps,nlevs_tgt);

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
                y_tgt[ilev] = mask_val;
              }
            }
          } else {
            // Near top
            if (x_tgt[ilev]<x_min) {
              if (etop==P0) {
                y_tgt[ilev] = y_src[0];
              } else {
                y_tgt[ilev] = mask_val;
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
}

} // namespace scream
