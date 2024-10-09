#include "vertical_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/field/field_tag.hpp"
#include "share/field/field_identifier.hpp"
#include "share/util/scream_universal_constants.hpp"
#include "share/io/scream_scorpio_interface.hpp"

#include <ekat/util/ekat_units.hpp>
#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>
#include <ekat/ekat_pack_kokkos.hpp>

#include <numeric>

namespace scream
{

VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const std::string& map_file,
                  const Field& pmid_src,
                  const Field& pint_src)
  : VerticalRemapper(src_grid,map_file,pmid_src,pint_src,constants::DefaultFillValue<float>::value)
{
  // Nothing to do here
}

VerticalRemapper::
VerticalRemapper (const grid_ptr_type& src_grid,
                  const std::string& map_file,
                  const Field& pmid_src,
                  const Field& pint_src,
                  const Real mask_val)
 : AbstractRemapper()
 , m_comm (src_grid->get_comm())
 , m_mask_val(mask_val)
{
  using namespace ShortFieldTagsNames;

  // Sanity checks
  EKAT_REQUIRE_MSG (src_grid->type()==GridType::Point,
      "Error! VerticalRemapper only works on PointGrid grids.\n"
      "  - src grid name: " + src_grid->name() + "\n"
      "  - src_grid_type: " + e2str(src_grid->type()) + "\n");
  EKAT_REQUIRE_MSG (src_grid->is_unique(),
      "Error! VerticalRemapper requires a unique source grid.\n");

  // This is a vertical remapper. We only go in one direction
  m_bwd_allowed = false;

  // Create tgt_grid that is a clone of the src grid but with
  // the correct number of levels.  Note that when vertically
  // remapping the target field will be defined on the same DOFs
  // as the source field, but will have a different number of 
  // vertical levels.
  scorpio::register_file(map_file,scorpio::FileMode::Read);
  auto nlevs_tgt = scorpio::get_dimlen(map_file,"lev");

  auto tgt_grid = src_grid->clone("vertical_remap_tgt_grid",true);
  tgt_grid->reset_num_vertical_lev(nlevs_tgt);
  this->set_grids(src_grid,tgt_grid);

  // Set the LEV and ILEV vertical profiles for interpolation from
  set_source_pressure_fields(pmid_src,pint_src);

  // Gather the pressure level data for vertical remapping
  set_pressure_levels(map_file);

  // Add tgt pressure levels to the tgt grid
  tgt_grid->set_geometry_data(m_tgt_pressure);

  scorpio::release_file(map_file);
}

FieldLayout VerticalRemapper::
create_src_layout (const FieldLayout& tgt_layout) const
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG (is_valid_tgt_layout(tgt_layout),
      "[VerticalRemapper] Error! Input target layout is not valid for this remapper.\n"
      " - input layout: " + tgt_layout.to_string());

  return create_layout(tgt_layout,m_src_grid);
}

FieldLayout VerticalRemapper::
create_tgt_layout (const FieldLayout& src_layout) const
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG (is_valid_src_layout(src_layout),
      "[VerticalRemapper] Error! Input source layout is not valid for this remapper.\n"
      " - input layout: " + src_layout.to_string());

  return create_layout(src_layout,m_tgt_grid);
}

FieldLayout VerticalRemapper::
create_layout (const FieldLayout& fl_in,
               const grid_ptr_type& grid_out) const
{
  // NOTE: for the vert remapper, it doesn't really make sense to distinguish
  //       between midpoints and interfaces: we're simply asking for a quantity
  //       at a given set of pressure levels. So we choose to have fl_out
  //       to *always* have LEV as vertical tag.
        auto fl_out = FieldLayout::invalid();
  switch (fl_in.type()) {
    case LayoutType::Scalar0D: [[ fallthrough ]];
    case LayoutType::Vector0D: [[ fallthrough ]];
    case LayoutType::Scalar2D: [[ fallthrough ]];
    case LayoutType::Vector2D: [[ fallthrough ]];
    case LayoutType::Tensor2D:
      // These layouts do not have vertical dim tags, so no change
      fl_out = fl_in;
      break;
    case LayoutType::Scalar1D:
      fl_out = grid_out->get_vertical_layout(true);
      break;
    case LayoutType::Scalar3D:
      fl_out = grid_out->get_3d_scalar_layout(true);
      break;
    case LayoutType::Vector3D:
      fl_out = grid_out->get_3d_vector_layout(true,fl_in.get_vector_dim());
      break;
    default:
      // NOTE: this also include Tensor3D. We don't really have any atm proc
      //       that needs to handle a tensor3d quantity, so no need to add it
      EKAT_ERROR_MSG (
        "[VerticalRemapper] Error! Layout not supported by VerticalRemapper.\n"
        " - input layout: " + fl_in.to_string() + "\n");
  }
  return fl_out;
}

void VerticalRemapper::
set_pressure_levels(const std::string& map_file)
{
  // Ensure each map file gets a different decomp name
  static std::map<std::string,int> file2idx;
  if (file2idx.find(map_file)==file2idx.end()) {
    file2idx[map_file] = file2idx.size();
  }

  using namespace ShortFieldTagsNames;
  auto layout = m_tgt_grid->get_vertical_layout(true);
  FieldIdentifier fid("p_levs",layout,ekat::units::Pa,m_tgt_grid->name());
  m_tgt_pressure = Field(fid);
  // Just in case input fields are packed
  m_tgt_pressure.get_header().get_alloc_properties().request_allocation(SCREAM_PACK_SIZE);
  m_tgt_pressure.allocate_view();

  auto remap_pres_data = m_tgt_pressure.get_view<Real*,Host>().data();
  scorpio::read_var(map_file,"p_levs",remap_pres_data);

  m_tgt_pressure.sync_to_dev();
}

void VerticalRemapper::
set_source_pressure_fields(const Field& pmid, const Field& pint)
{
  using namespace ShortFieldTagsNames;

  EKAT_REQUIRE_MSG(pmid.is_allocated(),
      "Error! Source midpoint pressure field is not yet allocated.\n"
      " - field name: " + pmid.name() + "\n");

  EKAT_REQUIRE_MSG(pint.is_allocated(),
      "Error! Source interface pressure field is not yet allocated.\n"
      " - field name: " + pint.name() + "\n");

  const auto& pmid_layout = pmid.get_header().get_identifier().get_layout();
  const auto& pint_layout = pint.get_header().get_identifier().get_layout();
  EKAT_REQUIRE_MSG(pmid_layout.congruent(m_src_grid->get_3d_scalar_layout(true)),
      "Error! Source midpoint pressure field has the wrong layout.\n"
      " - field name: " + pmid.name() + "\n"
      " - field layout: " + pmid_layout.to_string() + "\n"
      " - expected layout: " + m_src_grid->get_3d_scalar_layout(true).to_string() + "\n");
  EKAT_REQUIRE_MSG(pint_layout.congruent(m_src_grid->get_3d_scalar_layout(false)),
      "Error! Source interface pressure field has the wrong layout.\n"
      " - field name: " + pint.name() + "\n"
      " - field layout: " + pint_layout.to_string() + "\n"
      " - expected layout: " + m_src_grid->get_3d_scalar_layout(false).to_string() + "\n");

  m_src_pmid = pmid;
  m_src_pint = pint;
}

void VerticalRemapper::
do_register_field (const identifier_type& src, const identifier_type& tgt)
{
  using namespace ShortFieldTagsNames;

  // Note, for vertical remapper we set all target fields as having LEV as the vertical dimension.
  // So we check that all other tags between source and target match, but skip vert tag (since we
  // could have src with ILEV and tgt with LEV)
  auto src_layout = src.get_layout().clone();
  auto tgt_layout = tgt.get_layout().clone();
  EKAT_REQUIRE_MSG(src_layout.strip_dims({ILEV,LEV}).congruent(tgt_layout.strip_dims({LEV})),
    "[VerticalRemapper] Error! Once vertical level tag is stripped, src/tgt layouts are incompatible.\n"
    "  - src field name: " + src.name() + "\n"
    "  - tgt field name: " + tgt.name() + "\n"
    "  - src field layout: " + src_layout.to_string() + "\n"
    "  - tgt field layout: " + tgt_layout.to_string() + "\n");

  m_src_fields.emplace_back(src);
  m_tgt_fields.emplace_back(tgt);
}

void VerticalRemapper::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  using namespace ShortFieldTagsNames;
  using PackT = ekat::Pack<Real,SCREAM_PACK_SIZE>;

  m_src_fields[ifield] = src;
  m_tgt_fields[ifield] = tgt;

  // Clone src layout, since we may strip dims later for mask creation
  auto src_layout = src.get_header().get_identifier().get_layout().clone();

  auto& f_tgt = m_tgt_fields[ifield]; // Nonconst, since we need to set extra data in the header
  if (src_layout.has_tag(LEV) or src_layout.has_tag(ILEV)) {
    // Determine if this field can be handled with packs, and whether it's at midpoints
    // Add mask tracking to the target field. The mask tracks location of tgt pressure levs that are outside the
    // bounds of the src pressure field, and hence cannot be recovered by interpolation
    auto& ft = m_field2type[src.name()];
    ft.midpoints = src.get_header().get_identifier().get_layout().has_tag(LEV);
    ft.packed    = src.get_header().get_alloc_properties().is_compatible<PackT>() and
                   tgt.get_header().get_alloc_properties().is_compatible<PackT>();

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
        "[VerticalRemapper::do_bind_field] Error! Target field already has mask data assigned.\n"
        " - tgt field name: " + tgt.name() + "\n");
    EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_value"),
        "[VerticalRemapper::do_bind_field] Error! Target field already has mask value assigned.\n"
        " - tgt field name: " + tgt.name() + "\n");

    f_tgt.get_header().set_extra_data("mask_data",tgt_mask);
    f_tgt.get_header().set_extra_data("mask_value",m_mask_val);
  } else {
    // If a field does not have LEV or ILEV it may still have mask tracking assigned from somewhere else.
    // For instance, this could be a 2d field computed by FieldAtPressureLevel diagnostic.
    // In those cases we want to copy that mask tracking to the target field.
    if (src.get_header().has_extra_data("mask_data")) {
      EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_data"),
          "[VerticalRemapper::do_bind_field] Error! Target field already has mask data assigned.\n"
          " - tgt field name: " + tgt.name() + "\n");
      auto src_mask = src.get_header().get_extra_data<Field>("mask_data");
      f_tgt.get_header().set_extra_data("mask_data",src_mask);
    }
    if (src.get_header().has_extra_data("mask_value")) {
      EKAT_REQUIRE_MSG(not tgt.get_header().has_extra_data("mask_value"),
          "[VerticalRemapper::do_bind_field] Error! Target field already has mask value assigned.\n"
          " - tgt field name: " + tgt.name() + "\n");
      auto src_mask_val = src.get_header().get_extra_data<Real>("mask_value");
      f_tgt.get_header().set_extra_data("mask_value",src_mask_val);
    }
  }

  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    create_lin_interp ();
  }
}

void VerticalRemapper::do_registration_ends ()
{
  if (this->m_num_bound_fields==this->m_num_registered_fields) {
    create_lin_interp ();
  }
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

void VerticalRemapper::do_remap_fwd ()
{
  // 1. Setup any interp object that was created (if nullptr, no fields need it)
  if (m_lin_interp_mid_packed) {
    setup_lin_interp(*m_lin_interp_mid_packed,m_src_pmid);
  }
  if (m_lin_interp_int_packed) {
    setup_lin_interp(*m_lin_interp_int_packed,m_src_pint);
  }
  if (m_lin_interp_mid_scalar) {
    setup_lin_interp(*m_lin_interp_mid_scalar,m_src_pmid);
  }
  if (m_lin_interp_int_scalar) {
    setup_lin_interp(*m_lin_interp_int_scalar,m_src_pint);
  }

  using namespace ShortFieldTagsNames;

  // 2. Interpolate the fields
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f_src    = m_src_fields[i];
          auto& f_tgt    = m_tgt_fields[i];
    const auto& tgt_layout   = f_tgt.get_header().get_identifier().get_layout();
    if (tgt_layout.has_tag(LEV)) {
      const auto& type = m_field2type.at(f_src.name());
      // Dispatch interpolation to the proper lin interp object
      if (type.midpoints) {
        if (type.packed) {
          apply_vertical_interpolation(*m_lin_interp_mid_packed,f_src,f_tgt,m_src_pmid,m_mask_val);
        } else {
          apply_vertical_interpolation(*m_lin_interp_mid_scalar,f_src,f_tgt,m_src_pmid,m_mask_val);
        }
      } else {
        if (type.packed) {
          apply_vertical_interpolation(*m_lin_interp_int_packed,f_src,f_tgt,m_src_pint,m_mask_val);
        } else {
          apply_vertical_interpolation(*m_lin_interp_int_scalar,f_src,f_tgt,m_src_pint,m_mask_val);
        }
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
        apply_vertical_interpolation(*m_lin_interp_mid_packed,f_src,f_tgt,m_src_pmid,0);
      } else {
        apply_vertical_interpolation(*m_lin_interp_mid_scalar,f_src,f_tgt,m_src_pmid,0);
      }
    } else {
      if (type.packed) {
        apply_vertical_interpolation(*m_lin_interp_int_packed,f_src,f_tgt,m_src_pint,0);
      } else {
        apply_vertical_interpolation(*m_lin_interp_int_scalar,f_src,f_tgt,m_src_pint,0);
      }
    }
  }
}

template<int Packsize>
void VerticalRemapper::
setup_lin_interp (const ekat::LinInterp<Real,Packsize>& lin_interp,
                  const Field& p_src) const
{
  using LI_t = ekat::LinInterp<Real,Packsize>;
  using ESU = ekat::ExeSpaceUtils<DefaultDevice::execution_space>;
  using PackT = ekat::Pack<Real,Packsize>;
  auto p_src_v = p_src.get_view<const PackT**>();
  auto p_tgt_v = m_tgt_pressure.get_view<const PackT*>();

  auto lambda = KOKKOS_LAMBDA(typename LI_t::MemberType const& team) {
    const int icol = team.league_rank();
    lin_interp.setup(team,ekat::subview(p_src_v,icol),
                          p_tgt_v);
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
                             const Field& p_src,
                             const Real mask_val) const
{
  // Note: if Packsize==1, we grab packs of size 1, which are for sure
  //       compatible with the allocation
  using LI_t = ekat::LinInterp<Real,Packsize>;
  using PackT = ekat::Pack<Real,Packsize>;
  using ESU = ekat::ExeSpaceUtils<DefaultDevice::execution_space>;

  auto p_src_v = p_src.get_view<const PackT**>();
  auto x_tgt = m_tgt_pressure.get_view<const PackT*>();
  const auto& f_src_l = f_src.get_header().get_identifier().get_layout();
  const int ncols = m_src_grid->get_num_local_dofs();
  const int nlevs_tgt = m_tgt_grid->get_num_vertical_levels();
  const int nlevs_src = f_src_l.dims().back();
  const int npacks_tgt = ekat::PackInfo<Packsize>::num_packs(nlevs_tgt);

  const int last_src_pack_idx = ekat::PackInfo<Packsize>::last_pack_idx(nlevs_src);
  const int last_src_pack_end = ekat::PackInfo<Packsize>::last_vec_end(nlevs_src);
  
  switch(f_src.rank()) {
    case 2:
    {
      auto f_src_v = f_src.get_view<const PackT**>();
      auto f_tgt_v = f_tgt.get_view<      PackT**>();
      auto policy = ESU::get_default_team_policy(ncols,npacks_tgt);
      auto lambda = KOKKOS_LAMBDA(typename LI_t::MemberType const& team) {

        // Interpolate
        const int icol = team.league_rank();
        auto x_src = ekat::subview(p_src_v,icol);
        auto y_src = ekat::subview(f_src_v,icol);
        auto y_tgt = ekat::subview(f_tgt_v,icol);
        lin_interp.lin_interp(team,x_src,x_tgt,y_src,y_tgt,icol);
        team.team_barrier();

        // If x_tgt is extrapolated, set to mask_val
        auto x_min = x_src[0][0];
        auto x_max = x_src[last_src_pack_idx][last_src_pack_end-1];
        auto set_mask = [&](const int ipack) {
          auto in_range = ekat::range<PackT>(ipack*Packsize) < nlevs_tgt;
          auto oob = (x_tgt[ipack]<x_min or x_tgt[ipack]>x_max) and in_range;
          if (oob.any()) {
            y_tgt[ipack].set(oob,mask_val);
          }
        };
        Kokkos::parallel_for (Kokkos::TeamThreadRange(team,npacks_tgt), set_mask);
      };
      Kokkos::parallel_for("VerticalRemapper::apply_vertical_interpolation",policy,lambda);
      break;
    }
    case 3:
    {
      auto f_src_v = f_src.get_view<const PackT***>();
      auto f_tgt_v = f_tgt.get_view<      PackT***>();
      const auto& layout = f_src.get_header().get_identifier().get_layout();
      const int ncomps = layout.get_vector_dim();
      auto policy = ESU::get_default_team_policy(ncols*ncomps,npacks_tgt);

      auto lambda = KOKKOS_LAMBDA(typename LI_t::MemberType const& team)
      {
        // Interpolate
        const int icol = team.league_rank() / ncomps;
        const int icmp = team.league_rank() % ncomps;
        auto x_src = ekat::subview(p_src_v,icol);
        auto y_src = ekat::subview(f_src_v,icol,icmp);
        auto y_tgt = ekat::subview(f_tgt_v,icol,icmp);
        lin_interp.lin_interp(team,x_src,x_tgt,y_src,y_tgt,icol);
        team.team_barrier();

        // If x_tgt is extrapolated, set to mask_val
        auto x_min = x_src[0][0];
        auto x_max = x_src[last_src_pack_idx][last_src_pack_end-1];
        auto set_mask = [&](const int ipack) {
          auto oob = x_tgt[ipack]<x_min or x_tgt[ipack]>x_max;
          if (oob.any()) {
            y_tgt[ipack].set(oob,mask_val);
          }
        };
        Kokkos::parallel_for (Kokkos::TeamThreadRange(team,npacks_tgt), set_mask);
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

} // namespace scream
