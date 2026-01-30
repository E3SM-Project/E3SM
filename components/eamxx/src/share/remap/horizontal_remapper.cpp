#include "horizontal_remapper.hpp"

#include "share/scorpio_interface/eamxx_scorpio_interface.hpp"
#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/field/field.hpp"

#include <ekat_team_policy_utils.hpp>
#include <ekat_pack_utils.hpp>

#include <numeric>

namespace scream
{

HorizontalRemapper::
HorizontalRemapper (const grid_ptr_type& grid,
                    const std::string& map_file,
                    const bool track_mask)
 : m_track_mask (track_mask)
{
  // Horiz remappers are built from a map file, which only goes in one direction
  m_bwd_allowed = false;

  // Sanity checks
  EKAT_REQUIRE_MSG (grid->type()==GridType::Point,
      "Error! Horizontal interpolatory remap only works on PointGrid grids.\n"
      "  - grid name: " + grid->name() + "\n"
      "  - grid_type: " + e2str(grid->type()) + "\n");
  EKAT_REQUIRE_MSG (grid->is_unique(),
      "Error! HorizInterpRemapperBase requires a unique grid.\n");

  // Get the remap data (if not already present, it will be built)
  m_remap_data = HorizRemapperDataRepo::instance().get_data(grid,map_file);

  // The grids really only matter for the horiz part. We may have 2+ remappers with
  // grids that only differ in terms of number of levs. Such remappers cannot
  // store the same generated grid.
  // So we soft-clone the generated grid, and reset the number of levels.
  // This requires to also delete any geo data that has the lev dim, as we cannot map it
  auto gen_grid = m_remap_data->m_generated_grid->clone(map_file,true);
  gen_grid->reset_num_vertical_lev(grid->get_num_vertical_levels());
  using namespace ShortFieldTagsNames;
  for (const auto& name : gen_grid->get_geometry_data_names()) {
    const auto& f = gen_grid->get_geometry_data(name);
    const auto& fl = f.get_header().get_identifier().get_layout();
    if (fl.has_tag(LEV) or fl.has_tag(ILEV)) {
      gen_grid->delete_geometry_data(name);
    }
  }

  if (m_remap_data->m_built_from_src) {
    set_grids(grid,gen_grid);
  } else {
    set_grids(gen_grid,grid);
  }

  if (m_remap_data->m_built_from_src) {
    // If the src grid contains geo data that is not in the tgt grid, then transfer it.
    // If the geo data has the COL dimension, then remap it.
    const auto& src_geo_data_names = m_src_grid->get_geometry_data_names();
    for (const auto& name : src_geo_data_names) {
      // Since different remappers may share the same data (if the map file is the same)
      // the coarse grid may already have the geo data.
      if (m_tgt_grid->has_geometry_data(name)) {
        continue;
      }
      const auto& src_data = m_src_grid->get_geometry_data(name);
      auto tgt_data = register_field_from_src(src_data);
      m_tgt_grid->set_geometry_data(tgt_data);
    }

    // It is possible that the fields we registered (if any) did not have the COL tag.
    registration_ends();
    if (m_num_fields>0) {
      remap_fwd();
      // The remap phase only alters the fields on device.
      // We need to sync them to host as well
      for (auto& f : m_tgt_fields) {
        f.sync_to_host();
      }
    }
    clean_up();
  }
}

HorizontalRemapper::
~HorizontalRemapper ()
{
  clean_up();
}

void HorizontalRemapper::
registration_ends_impl ()
{
  if (m_track_mask) {
    // We store masks (src-tgt) here, and register them AFTER we parse all currently registered fields.
    // That makes the for loop below easier, since we can take references without worrring that they
    // would get invalidated. In fact, if you call register_field inside the loop, the src/tgt fields
    // vectors will grow, which may cause reallocation and references invalidation
    std::vector<std::pair<Field,Field>> masks;

    auto get_mask_idx = [&](const FieldIdentifier& src_mask_fid) {

      // Masks will be registered AFTER all fields, so the 1st mask will
      // be right after the last registered "regular" field.
      int idx = 0;
      for (const auto& it : masks) {
        if (it.first.get_header().get_identifier()==src_mask_fid) {
          return idx;
        }
        ++idx;
      }
      return -1;
    };

    for (int i=0; i<m_num_fields; ++i) {
      const auto& src = m_src_fields[i];
            auto& tgt = m_tgt_fields[i];
      if (not src.get_header().has_extra_data("mask_field"))
        continue;

      const auto& src_mask = src.get_header().get_extra_data<Field>("mask_field");

      // Make sure fields representing masks are not themselves meant to be masked.
      EKAT_REQUIRE_MSG(not src_mask.get_header().has_extra_data("mask_field"),
          "Error! A mask field cannot be itself masked.\n"
          "  - field name: " + src.name() + "\n"
          "  - mask field name: " + src_mask.name() + "\n");

      // Check that the mask field has the correct layout
      const auto& f_lt = src.get_header().get_identifier().get_layout();
      const auto& m_lt = src_mask.get_header().get_identifier().get_layout();
      using namespace ShortFieldTagsNames;
      EKAT_REQUIRE_MSG(f_lt.has_tag(COL) == m_lt.has_tag(COL),
          "Error! Incompatible field and mask layouts.\n"
          "  - field name: " + src.name() + "\n"
          "  - field layout: " + f_lt.to_string() + "\n"
          "  - mask layout: " + m_lt.to_string() + "\n");
      EKAT_REQUIRE_MSG(f_lt.has_tag(LEV) == m_lt.has_tag(LEV),
          "Error! Incompatible field and mask layouts.\n"
          "  - field name: " + src.name() + "\n"
          "  - field layout: " + f_lt.to_string() + "\n"
          "  - mask layout: " + m_lt.to_string() + "\n");

      // If it's the first time we find this mask, store it, so we can register later
      const auto& src_mask_fid = src_mask.get_header().get_identifier();
      int mask_idx = get_mask_idx(src_mask_fid);
      if (mask_idx==-1) {
        Field tgt_mask(create_tgt_fid(src_mask_fid));
        auto src_pack_size = src_mask.get_header().get_alloc_properties().get_largest_pack_size();
        tgt_mask.get_header().get_alloc_properties().request_allocation(src_pack_size);
        tgt_mask.allocate_view();

        masks.push_back(std::make_pair(src_mask,tgt_mask));
        mask_idx = masks.size()-1;
      }
      tgt.get_header().set_extra_data("mask_field",masks[mask_idx].second);
    }

    // Add all masks to the fields to remap
    for (const auto& it : masks) {
      register_field(it.first,it.second);
    }
  }

  using namespace ShortFieldTagsNames;

  m_needs_remap.resize(m_num_fields,1);
  for (int i=0; i<m_num_fields; ++i) {
    const auto& src_dt = m_src_fields[i].get_header().get_identifier().data_type();
    const auto& tgt_dt = m_tgt_fields[i].get_header().get_identifier().data_type();
    EKAT_REQUIRE_MSG (src_dt==DataType::RealType and tgt_dt==DataType::RealType,
        "Error! HorizInterpRmapperBase requires src/tgt fields to have Real data type.\n"
        "  - src field name: " + m_src_fields[i].name() + "\n"
        "  - tgt field name: " + m_tgt_fields[i].name() + "\n"
        "  - src data type : " + e2str(src_dt) + "\n"
        "  - tgt data type : " + e2str(tgt_dt) + "\n");

    const auto& src_fl = m_src_fields[i].get_header().get_identifier().get_layout();
    if (not src_fl.has_tag(COL)) {
      // This field will be skipped in several of the remap steps
      m_needs_remap[i] = 0;
    } else {
      EKAT_REQUIRE_MSG (src_fl.tag(0)==COL,
          "[HorizInterpRemapperBase::registration_ends_impl] Error! If present, the COL dimension MUST be the first one.\n"
          " - field name: " + m_src_fields[i].name() + "\n"
          " - field layout: " + src_fl.to_string() + "\n");
    }
  }

  create_ov_fields ();
  setup_mpi_data_structures ();
}

void HorizontalRemapper::create_ov_fields ()
{
  using namespace ShortFieldTagsNames;

  m_ov_fields.reserve(m_num_fields);
  const auto num_ov_gids = m_remap_data->m_overlap_grid->get_num_local_dofs();
  const auto ov_gn = m_remap_data->m_overlap_grid->name();
  const auto dt = DataType::RealType;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_remap_data->m_coarsening ? m_src_fields[i] : m_tgt_fields[i];
    const auto& fid = f.get_header().get_identifier();
    if (m_needs_remap[i]==0) {
      // This field won't be remapped. We can simply emplace an empty field (which won't be used),
      // to make sure m_ov_fields[i] always returns the ov field for the i-th field
      m_ov_fields.emplace_back();
      continue;
    }

    const auto layout = fid.get_layout().clone().reset_dim(0,num_ov_gids);
    FieldIdentifier ov_fid (fid.name(),layout,fid.get_units(),ov_gn,dt);

    auto& ov_f = m_ov_fields.emplace_back(ov_fid);

    // Use same alloc props as fine fields, to allow packing in local_mat_vec
    const auto pack_size = f.get_header().get_alloc_properties().get_largest_pack_size();
    ov_f.get_header().get_alloc_properties().request_allocation(pack_size);
    ov_f.allocate_view();
  }
}

void HorizontalRemapper::remap_fwd_impl ()
{
  const auto& comm = m_src_grid->get_comm();

  // Fire the recv requests right away, so that if some other ranks
  // is done packing before us, we can start receiving their data
  if (not m_recv_req.empty()) {
    int ierr = MPI_Startall(m_recv_req.size(),m_recv_req.data());
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something went wrong while starting persistent recv requests.\n"
        "  - recv rank: " + std::to_string(comm.rank()) + "\n");
  }

  // TODO: Add check that if there are mask values they are either 1's or 0's for unmasked/masked.

  // Helper function, to establish if a field can be handled with packs
  auto can_pack_field = [](const Field& f) {
    const auto& ap = f.get_header().get_alloc_properties();
    return (ap.get_last_extent() % SCREAM_PACK_SIZE) == 0;
  };

  bool coarsen = m_remap_data->m_coarsening;

  if (not coarsen) {
    // For refining, MPI happens on the src grid
    pack_and_send();
    recv_and_unpack();
  }

  // Perform the local mat-vec using the proper fields depending on coarsen
  for (int i=0; i<m_num_fields; ++i) {
    if (m_needs_remap[i]==0) {
      // No need to do a mat-vec here. Just deep copy and move on
      m_tgt_fields[i].deep_copy(m_src_fields[i]);
      continue;
    }

    const auto& x = coarsen ? m_src_fields[i] : m_ov_fields[i];
    const auto& y = coarsen ? m_ov_fields[i] : m_tgt_fields[i];

    const bool masked = m_track_mask and x.get_header().has_extra_data("mask_field");
    if (masked) {
      // Pass the mask to the local_mat_vec routine
      const auto& mask = x.get_header().get_extra_data<Field>("mask_field");

      // If possible, dispatch kernel with SCREAM_PACK_SIZE
      if (can_pack_field(x) and can_pack_field(y) and can_pack_field(mask)) {
        local_mat_vec<SCREAM_PACK_SIZE>(x,y,mask);
      } else {
        local_mat_vec<1>(x,y,mask);
      }
    } else {
      // If possible, dispatch kernel with SCREAM_PACK_SIZE
      if (can_pack_field(x) and can_pack_field(y)) {
        local_mat_vec<SCREAM_PACK_SIZE>(x,y);
      } else {
        local_mat_vec<1>(x,y);
      }
    }
  }

  if (coarsen) {
    // For coarsening, MPI happens on the tgt grid
    pack_and_send ();
    recv_and_unpack ();
  }

  // Wait for all sends to be completed
  if (not m_send_req.empty()) {
    int ierr = MPI_Waitall(m_send_req.size(),m_send_req.data(), MPI_STATUSES_IGNORE);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something went wrong while waiting on persistent send requests.\n"
        "  - send rank: " + std::to_string(comm.rank()) + "\n");
  }

  // Rescale any fields that had the mask applied.
  if (m_track_mask) {
    for (int i=0; i<m_num_fields; ++i) {
      const auto& f_tgt = m_tgt_fields[i];
      if (f_tgt.get_header().has_extra_data("mask_field")) {
        // Then this field did use a mask
        const auto& mask = f_tgt.get_header().get_extra_data<Field>("mask_field");
        if (can_pack_field(f_tgt) and can_pack_field(mask)) {
          rescale_masked_fields<SCREAM_PACK_SIZE>(f_tgt,mask);
        } else {
          rescale_masked_fields<1>(f_tgt,mask);
        }
      }
    }
  }
}

template<int PackSize>
void HorizontalRemapper::
local_mat_vec (const Field& x, const Field& y) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using TPF         = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto row_grid = m_remap_data->m_coarsening ? m_remap_data->m_overlap_grid : m_tgt_grid;
  const int  nrows    = row_grid->get_num_local_dofs();

  const auto& src_layout = x.get_header().get_identifier().get_layout();
  const int   rank       = src_layout.rank();

  auto row_offsets = m_remap_data->m_row_offsets;
  auto col_lids    = m_remap_data->m_col_lids;
  auto weights     = m_remap_data->m_weights;

  switch (rank) {
    // Note: in each case, handle 1st contribution to each row separately,
    //       using = instead of +=. This allows to avoid doing an extra
    //       loop to zero out y before the mat-vec.
    case 1:
    {
      // Unlike get_view, get_strided_view returns a LayoutStride view,
      // therefore allowing the 1d field to be a subfield of a 2d field
      // along the 2nd dimension.
      auto x_view = x.get_strided_view<const Real*>();
      auto y_view = y.get_strided_view<      Real*>();
      Kokkos::parallel_for(RangePolicy(0,nrows),
                           KOKKOS_LAMBDA(const int& row) {
        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        y_view(row) = weights(beg)*x_view(col_lids(beg));
        for (int icol=beg+1; icol<end; ++icol) {
          y_view(row) += weights(icol)*x_view(col_lids(icol));
        }
      });
      break;
    }
    case 2:
    {
      auto x_view = x.get_view<const Pack**>();
      auto y_view = y.get_view<      Pack**>();
      const int dim1 = PackInfo::num_packs(src_layout.dim(1));
      auto policy = TPF::get_default_team_policy(nrows,dim1);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                            [&](const int j){
          y_view(row,j) = weights(beg)*x_view(col_lids(beg),j);
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j) += weights(icol)*x_view(col_lids(icol),j);
          }
        });
      });
      break;
    }
    case 3:
    {
      auto x_view = x.get_view<const Pack***>();
      auto y_view = y.get_view<      Pack***>();
      const int dim1 = src_layout.dim(1);
      const int dim2 = PackInfo::num_packs(src_layout.dim(2));
      auto policy = TPF::get_default_team_policy(nrows,dim1*dim2);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2),
                            [&](const int idx){
          const int j = idx / dim2;
          const int k = idx % dim2;
          y_view(row,j,k) = weights(beg)*x_view(col_lids(beg),j,k);
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j,k) += weights(icol)*x_view(col_lids(icol),j,k);
          }
        });
      });
      break;
    }
    case 4:
    {
      auto x_view = x.get_view<const Pack****>();
      auto y_view = y.get_view<      Pack****>();
      const int dim1 = src_layout.dim(1);
      const int dim2 = src_layout.dim(2);
      const int dim3 = PackInfo::num_packs(src_layout.dim(3));
      auto policy = TPF::get_default_team_policy(nrows,dim1*dim2*dim3);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2*dim3),
                            [&](const int idx){
          const int j = (idx / dim3) / dim2;
          const int k = (idx / dim3) % dim2;
          const int l =  idx % dim3;
          y_view(row,j,k,l) = weights(beg)*x_view(col_lids(beg),j,k,l);
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j,k,l) += weights(icol)*x_view(col_lids(icol),j,k,l);
          }
        });
      });
      break;
    }
    default:
      EKAT_ERROR_MSG("[HorizInterpRemapperBase::local_mat_vec] Error! Fields of rank 4 or greater are not supported.\n");
  }
}

template<int PackSize>
void HorizontalRemapper::
rescale_masked_fields (const Field& x, const Field& mask) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using TPF         = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  constexpr auto fill_val = constants::fill_value<Real>;

  const auto& layout = x.get_header().get_identifier().get_layout();
  const int rank = layout.rank();
  const int ncols = m_tgt_grid->get_num_local_dofs();
  const Real mask_threshold = std::numeric_limits<Real>::epsilon();  // TODO: Should we not hardcode the threshold for simply masking out the column.

  switch (rank) {
    case 1:
    {
      // Unlike get_view, get_strided_view returns a LayoutStride view,
      // therefore allowing the 1d field to be a subfield of a 2d field
      // along the 2nd dimension.
      auto x_view =    x.get_strided_view<      Real*>();
      auto m_view = mask.get_strided_view<const Real*>();
      Kokkos::parallel_for(RangePolicy(0,ncols),
                           KOKKOS_LAMBDA(const int& icol) {
        if (m_view(icol)>mask_threshold) {
          x_view(icol) /= m_view(icol);
        } else {
          x_view(icol) = fill_val;
        }
      });
      break;
    }
    case 2:
    {
      auto x_view =    x.get_view<      Pack**>();
      bool mask1d = mask.rank()==1;
      view_1d<const Real> mask_1d;
      view_2d<const Pack> mask_2d;
      // If the mask comes from FieldAtLevel, it's only defined on columns (rank=1)
      // If the mask comes from vert interpolation remapper, it is defined on ncols x nlevs (rank=2)
      if (mask.rank()==1) {
        mask_1d = mask.get_view<const Real*>();
      } else {
        mask_2d = mask.get_view<const Pack**>();
      }
      const int dim1 = PackInfo::num_packs(layout.dim(1));
      auto policy = TPF::get_default_team_policy(ncols,dim1);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto icol = team.league_rank();
        auto x_sub = ekat::subview(x_view,icol);
        if (mask1d) {
          auto mask = mask_1d(icol);
          if (mask>mask_threshold) {
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                                [&](const int j){
                x_sub(j) /= mask;
            });
          }
        } else {
          auto m_sub = ekat::subview(mask_2d,icol);
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                              [&](const int j){
            auto masked = m_sub(j) > mask_threshold;
            if (masked.any()) {
              x_sub(j).set(masked,x_sub(j)/m_sub(j));
            }
            x_sub(j).set(!masked,fill_val);
          });
        }
      });
      break;
    }
    case 3:
    {
      auto x_view =    x.get_view<      Pack***>();
      bool mask1d = mask.rank()==1;
      view_1d<const Real> mask_1d;
      view_2d<const Pack> mask_2d;
      // If the mask comes from FieldAtLevel, it's only defined on columns (rank=1)
      // If the mask comes from vert interpolation remapper, it is defined on ncols x nlevs (rank=2)
      if (mask.rank()==1) {
        mask_1d = mask.get_view<const Real*>();
      } else {
        mask_2d = mask.get_view<const Pack**>();
      }
      const int dim1 = layout.dim(1);
      const int dim2 = PackInfo::num_packs(layout.dim(2));
      auto policy = TPF::get_default_team_policy(ncols,dim1*dim2);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto icol = team.league_rank();
        if (mask1d) {
          auto mask = mask_1d(icol);
          if (mask>mask_threshold) {
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2),
                                [&](const int idx){
              const int j = idx / dim2;
              const int k = idx % dim2;
              auto x_sub = ekat::subview(x_view,icol,j);
              x_sub(k) /= mask;
            });
          }
        } else {
          auto m_sub      = ekat::subview(mask_2d,icol);
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2),
                              [&](const int idx){
            const int j = idx / dim2;
            const int k = idx % dim2;
            auto x_sub = ekat::subview(x_view,icol,j);
            auto masked = m_sub(k) > mask_threshold;

            if (masked.any()) {
              x_sub(k).set(masked,x_sub(k)/m_sub(k));
            }
            x_sub(k).set(!masked,fill_val);
          });
        }
      });
      break;
    }
    case 4:
    {
      auto x_view = x.get_view<Pack****>();
      bool mask1d = mask.rank()==1;
      view_1d<const Real> mask_1d;
      view_2d<const Pack> mask_2d;
      // If the mask comes from FieldAtLevel, it's only defined on columns (rank=1)
      // If the mask comes from vert interpolation remapper, it is defined on ncols x nlevs (rank=2)
      if (mask.rank()==1) {
        mask_1d = mask.get_view<const Real*>();
      } else {
        mask_2d = mask.get_view<const Pack**>();
      }
      const int dim1 = layout.dim(1);
      const int dim2 = layout.dim(2);
      const int dim3 = PackInfo::num_packs(layout.dim(3));
      auto policy = TPF::get_default_team_policy(ncols,dim1*dim2*dim3);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto icol = team.league_rank();
        if (mask1d) {
          auto mask = mask_1d(icol);
          if (mask>mask_threshold) {
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2*dim3),
                                [&](const int idx){
              const int j = (idx / dim3) / dim2;
              const int k = (idx / dim3) % dim2;
              const int l =  idx % dim3;
              auto x_sub = ekat::subview(x_view,icol,j,k);
              x_sub(l) /= mask;
            });
          }
        } else {
          auto m_sub      = ekat::subview(mask_2d,icol);
          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2*dim3),
                              [&](const int idx){
            const int j = (idx / dim3) / dim2;
            const int k = (idx / dim3) % dim2;
            const int l =  idx % dim3;
            auto x_sub = ekat::subview(x_view,icol,j,k);
            auto masked = m_sub(l) > mask_threshold;

            if (masked.any()) {
              x_sub(l).set(masked,x_sub(l)/m_sub(l));
            }
            x_sub(l).set(!masked,fill_val);
          });
        }
      });
      break;
    }
  }
}

template<int PackSize>
void HorizontalRemapper::
local_mat_vec (const Field& x, const Field& y, const Field& mask) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using TPF         = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto& src_layout = x.get_header().get_identifier().get_layout();
  const int rank = src_layout.rank();
  const int nrows  = m_remap_data->m_overlap_grid->get_num_local_dofs();
  auto row_offsets = m_remap_data->m_row_offsets;
  auto col_lids    = m_remap_data->m_col_lids;
  auto weights     = m_remap_data->m_weights;
  switch (rank) {
    // Note: in each case, handle 1st contribution to each row separately,
    //       using = instead of +=. This allows to avoid doing an extra
    //       loop to zero out y before the mat-vec.
    case 1:
    {
      // Unlike get_view, get_strided_view returns a LayoutStride view,
      // therefore allowing the 1d field to be a subfield of a 2d field
      // along the 2nd dimension.
      auto x_view    =    x.get_strided_view<const Real*>();
      auto y_view    =    y.get_strided_view<      Real*>();
      auto mask_view = mask.get_strided_view<const Real*>();
      Kokkos::parallel_for(RangePolicy(0,nrows),
                           KOKKOS_LAMBDA(const int& row) {
        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        y_view(row) = weights(beg)*x_view(col_lids(beg))*mask_view(col_lids(beg));
        for (int icol=beg+1; icol<end; ++icol) {
          y_view(row) += weights(icol)*x_view(col_lids(icol))*mask_view(col_lids(icol));
        }
      });
      break;
    }
    case 2:
    {
      auto x_view = x.get_view<const Pack**>();
      auto y_view = y.get_view<      Pack**>();
      view_1d<const Real> mask_1d;
      view_2d<const Pack> mask_2d;
      // If the mask comes from FieldAtLevel, it's only defined on columns (rank=1)
      // If the mask comes from vert interpolation remapper, it is defined on ncols x nlevs (rank=2)
      bool mask1d = mask.rank()==1;
      if (mask1d) {
        mask_1d = mask.get_view<const Real*>();
      } else {
        mask_2d = mask.get_view<const Pack**>();
      }
      const int dim1 = PackInfo::num_packs(src_layout.dim(1));
      auto policy = TPF::get_default_team_policy(nrows,dim1);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                            [&](const int j){
          y_view(row,j) = weights(beg)*x_view(col_lids(beg),j) *
                          (mask1d ? mask_1d (col_lids(beg)) : mask_2d(col_lids(beg),j));
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j) += weights(icol)*x_view(col_lids(icol),j) *
                          (mask1d ? mask_1d (col_lids(icol)) : mask_2d(col_lids(icol),j));
          }
        });
      });
      break;
    }
    case 3:
    {
      auto x_view = x.get_view<const Pack***>();
      auto y_view = y.get_view<      Pack***>();
      // Note, the mask is still assumed to be defined on COLxLEV so still only 2D for case 3.
      view_1d<const Real> mask_1d;
      view_2d<const Pack> mask_2d;
      bool mask1d = mask.rank()==1;
      // If the mask comes from FieldAtLevel, it's only defined on columns (rank=1)
      // If the mask comes from vert interpolation remapper, it is defined on ncols x nlevs (rank=2)
      if (mask1d) {
        mask_1d = mask.get_view<const Real*>();
      } else {
        mask_2d = mask.get_view<const Pack**>();
      }
      const int dim1 = src_layout.dim(1);
      const int dim2 = PackInfo::num_packs(src_layout.dim(2));
      auto policy = TPF::get_default_team_policy(nrows,dim1*dim2);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2),
                            [&](const int idx){
          const int j = idx / dim2;
          const int k = idx % dim2;
          y_view(row,j,k) = weights(beg)*x_view(col_lids(beg),j,k) * 
                          (mask1d ? mask_1d (col_lids(beg)) : mask_2d(col_lids(beg),k));
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j,k) += weights(icol)*x_view(col_lids(icol),j,k) *
                          (mask1d ? mask_1d (col_lids(icol)) : mask_2d(col_lids(icol),k));
          }
        });
      });
      break;
    }
    case 4:
    {
      auto x_view = x.get_view<const Pack****>();
      auto y_view = y.get_view<      Pack****>();
      // Note, the mask is still assumed to be defined on COLxLEV so still only 2D for case 3.
      view_1d<const Real> mask_1d;
      view_2d<const Pack> mask_2d;
      bool mask1d = mask.rank()==1;
      // If the mask comes from FieldAtLevel, it's only defined on columns (rank=1)
      // If the mask comes from vert interpolation remapper, it is defined on ncols x nlevs (rank=2)
      if (mask1d) {
        mask_1d = mask.get_view<const Real*>();
      } else {
        mask_2d = mask.get_view<const Pack**>();
      }
      const int dim1 = src_layout.dim(1);
      const int dim2 = src_layout.dim(2);
      const int dim3 = PackInfo::num_packs(src_layout.dim(3));
      auto policy = TPF::get_default_team_policy(nrows,dim1*dim2*dim3);
      Kokkos::parallel_for(policy,
                           KOKKOS_LAMBDA(const MemberType& team) {
        const auto row = team.league_rank();

        const auto beg = row_offsets(row);
        const auto end = row_offsets(row+1);
        Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2*dim3),
                            [&](const int idx){
          const int j = (idx / dim3) / dim2;
          const int k = (idx / dim3) % dim2;
          const int l =  idx % dim3;
          y_view(row,j,k,l) = weights(beg)*x_view(col_lids(beg),j,k,l) * 
                          (mask1d ? mask_1d (col_lids(beg)) : mask_2d(col_lids(beg),l));
          for (int icol=beg+1; icol<end; ++icol) {
            y_view(row,j,k,l) += weights(icol)*x_view(col_lids(icol),j,k,l) *
                          (mask1d ? mask_1d (col_lids(icol)) : mask_2d(col_lids(icol),l));
          }
        });
      });
      break;
    }
    default:
    {
      EKAT_ERROR_MSG("Error::coarsening_remapper::local_mat_vec doesn't support fields of rank 4 or greater");
    }
  }
}

void HorizontalRemapper::pack_and_send ()
{
  using RangePolicy = typename KT::RangePolicy;
  using TeamMember  = typename KT::MemberType;
  using TPF         = ekat::TeamPolicyFactory<DefaultDevice::execution_space>;

  bool coarsen = m_remap_data->m_coarsening;

  auto pids = coarsen ? m_imp_exp->import_pids() : m_imp_exp->export_pids();
  auto lids = coarsen ? m_imp_exp->import_lids() : m_imp_exp->export_lids();
  auto ncols_send  = coarsen ? m_imp_exp->num_imports_per_pid() : m_imp_exp->num_exports_per_pid();
  auto pids_send_offsets = m_pids_send_offsets;
  auto send_buf = m_send_buffer;
  const int num_iters = pids.size();
  const int total_col_size = m_field_offset.back();
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    if (m_needs_remap[ifield]==0)
      // No need to process this field. We'll deep copy src->tgt later
      continue;

    const auto& f = coarsen ? m_ov_fields[ifield] : m_src_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto field_offset = m_field_offset[ifield];
    switch (fl.rank()) {
      case 1: 
      { 
        const auto v = f.get_strided_view<const Real*>();
        auto pack = KOKKOS_LAMBDA(const int idx) {
          auto pid = pids(idx);
          auto icol = lids(idx);
          auto pid_offset = pids_send_offsets(pid); 
          auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_send(pid)*field_offset
                      + pos_within_pid;
          send_buf(offset) = v(icol);
        };
        Kokkos::parallel_for(RangePolicy(0,num_iters),pack);
        break;        
      }   
      case 2:
      {     
        const auto v = f.get_view<const Real**>();
        const int dim1 = fl.dim(1);
        auto policy = TPF::get_default_team_policy(num_iters,dim1);
        auto pack = KOKKOS_LAMBDA(const TeamMember& team) {
          const int idx = team.league_rank();
          const int icol = lids(idx);
          const int pid  = pids(idx);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_send(pid)*field_offset
                      + pos_within_pid*dim1;
          auto col_pack = [&](const int& k) {
            send_buf(offset+k) = v(icol,k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,dim1);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,pack);
        break;
      }
      case 3:
      {
        const auto v = f.get_view<const Real***>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int f_col_size = dim1*dim2;
        auto policy = TPF::get_default_team_policy(num_iters,dim1*dim2);
        auto pack = KOKKOS_LAMBDA(const TeamMember& team) {
          const int idx = team.league_rank();
          const int icol = lids(idx);
          const int pid  = pids(idx);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_send(pid)*field_offset
                      + pos_within_pid*f_col_size;
          auto col_pack = [&](const int& idx) {
            const int j = idx / dim2;
            const int k = idx % dim2;
            send_buf(offset+idx) = v(icol,j,k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,f_col_size);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,pack);
        break;
      }
      case 4:
      {
        const auto v = f.get_view<const Real****>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int dim3 = fl.dim(3);
        const int f_col_size = dim1*dim2*dim3;
        auto policy = TPF::get_default_team_policy(num_iters,dim1*dim2*dim3);
        auto pack = KOKKOS_LAMBDA(const TeamMember& team) {
          const int idx = team.league_rank();
          const int icol = lids(idx);
          const int pid  = pids(idx);
          auto pid_offset = pids_send_offsets(pid);
          auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_send(pid)*field_offset
                      + pos_within_pid*f_col_size;
          auto col_pack = [&](const int& idx) {
            const int j = (idx / dim3) / dim2;
            const int k = (idx / dim3) % dim2;
            const int l =  idx % dim3;
            send_buf(offset+idx) = v(icol,j,k,l);
          };
          auto tvr = Kokkos::TeamVectorRange(team,f_col_size);
          Kokkos::parallel_for(tvr,col_pack);
        };
        Kokkos::parallel_for(policy,pack);
        break;
      }
      default:
        EKAT_ERROR_MSG ("Unexpected field rank in RefiningRemapperP2P::pack.\n"
            "  - MPI rank  : " + std::to_string(m_src_grid->get_comm().rank()) + "\n"
            "  - field name: " + f.name() + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }

  // Ensure all threads are done packing before firing off the sends
  Kokkos::fence();

  // If MPI does not use dev pointers, we need to deep copy from dev to host
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_mpi_send_buffer,m_send_buffer);
  }

  if (not m_send_req.empty()) {
    int ierr = MPI_Startall(m_send_req.size(),m_send_req.data());
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something went wrong while starting persistent send requests.\n"
        "  - send rank: " + std::to_string(m_src_grid->get_comm().rank()) + "\n");
  }
}

void HorizontalRemapper::recv_and_unpack ()
{
  if (not m_recv_req.empty()) {
    check_mpi_call(MPI_Waitall(m_recv_req.size(),m_recv_req.data(), MPI_STATUSES_IGNORE),
                   "[RefiningRemapperP2P] waiting on persistent recv requests.\n");
  }
  // If MPI does not use dev pointers, we need to deep copy from host to dev
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_recv_buffer,m_mpi_recv_buffer);
  }

  using RangePolicy = typename KT::RangePolicy;
  using TeamMember  = typename KT::MemberType;
  using TPF         = ekat::TeamPolicyFactory<typename KT::ExeSpace>;

  bool coarsen = m_remap_data->m_coarsening;

  auto pids = coarsen ? m_imp_exp->export_pids(): m_imp_exp->import_pids();
  auto lids = coarsen ? m_imp_exp->export_lids(): m_imp_exp->import_lids();
  auto ncols_recv  = coarsen ? m_imp_exp->num_exports_per_pid() : m_imp_exp->num_imports_per_pid();
  auto pids_recv_offsets = m_pids_recv_offsets;
  auto recv_buf = m_recv_buffer;
  const int num_iters = pids.size();
  const int total_col_size = m_field_offset.back();
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    if (m_needs_remap[ifield]==0)
      // No need to process this field. We'll deep copy src->tgt later
      continue;

          auto& f  = coarsen ? m_tgt_fields[ifield] : m_ov_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto field_offset = m_field_offset[ifield];

    // We accummulate contributions, so init to 0
    f.deep_copy(0);

    switch (fl.rank()) {
      case 1:
      {
        auto v = f.get_view<Real*>();
        auto unpack = KOKKOS_LAMBDA (const int idx) {
          const int pid  = pids(idx);
          const int icol = lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_recv(pid)*field_offset
                      + pos_within_pid;
          v(icol) += recv_buf(offset);
        };
        Kokkos::parallel_for(RangePolicy(0,num_iters),unpack);
        break;
      }
      case 2:
      {
        auto v = f.get_view<Real**>();
        const int dim1 = fl.dim(1);
        auto policy = TPF::get_default_team_policy(num_iters,dim1);
        auto unpack = KOKKOS_LAMBDA (const TeamMember& team) {
          const int idx  = team.league_rank();
          const int pid  = pids(idx);
          const int icol = lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_recv(pid)*field_offset
                      + pos_within_pid*dim1;
          auto col_unpack = [&](const int& k) {
            v(icol,k) += recv_buf(offset+k);
          };
          auto tvr = Kokkos::TeamVectorRange(team,dim1);
          Kokkos::parallel_for(tvr,col_unpack);
        };
        Kokkos::parallel_for(policy,unpack);
        break;
      }
      case 3:
      {
        auto v = f.get_view<Real***>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int f_col_size = dim1*dim2;
        auto policy = TPF::get_default_team_policy(num_iters,dim1*dim2);
        auto unpack = KOKKOS_LAMBDA (const TeamMember& team) {
          const int idx  = team.league_rank();
          const int pid  = pids(idx);
          const int icol = lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_recv(pid)*field_offset
                      + pos_within_pid*f_col_size;
          auto col_unpack = [&](const int& idx) {
            const int j = idx / dim2;
            const int k = idx % dim2;
            v(icol,j,k) += recv_buf(offset+idx);
          };
          auto tvr = Kokkos::TeamVectorRange(team,f_col_size);
          Kokkos::parallel_for(tvr,col_unpack);
        };
        Kokkos::parallel_for(policy,unpack);
        break;
      }
      case 4:
      {
        auto v = f.get_view<Real****>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int dim3 = fl.dim(3);
        const int f_col_size = dim1*dim2*dim3;
        auto policy = TPF::get_default_team_policy(num_iters,dim1*dim2*dim3);
        auto unpack = KOKKOS_LAMBDA (const TeamMember& team) {
          const int idx  = team.league_rank();
          const int pid  = pids(idx);
          const int icol = lids(idx);
          const auto pid_offset = pids_recv_offsets(pid);
          const auto pos_within_pid = idx - pid_offset;
          auto offset = pid_offset*total_col_size
                      + ncols_recv(pid)*field_offset
                      + pos_within_pid*f_col_size;
          auto col_unpack = [&](const int& idx) {
            const int j = (idx / dim3) / dim2;
            const int k = (idx / dim3) % dim2;
            const int l =  idx % dim3;
            v(icol,j,k,l) += recv_buf(offset+idx);
          };
          auto tvr = Kokkos::TeamVectorRange(team,f_col_size);
          Kokkos::parallel_for(tvr,col_unpack);
        };
        Kokkos::parallel_for(policy,unpack);
        break;
      }
      default:
        EKAT_ERROR_MSG ("Unexpected field rank in RefiningRemapperP2P::unpack.\n"
            "  - MPI rank  : " + std::to_string(m_src_grid->get_comm().rank()) + "\n"
            "  - field name: " + f.name() + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }
}

void HorizontalRemapper::setup_mpi_data_structures ()
{
  using namespace ShortFieldTagsNames;

  const int nranks = m_src_grid->get_comm().size();
  const bool coarsen = m_remap_data->m_coarsening;

  // Compute offset of each field when we splice together a col for each
  m_field_offset.resize(m_num_fields+1,0);
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f = m_src_fields[i];
    const auto& fl = f.get_header().get_identifier().get_layout();

    // Fields without COL tag are nor remapped, so consider their col size as 0
    auto col_size = m_needs_remap[i] ? fl.clone().strip_dim(COL).size() : 0;
    m_field_offset[i+1] = m_field_offset[i] + col_size;
  }
  auto total_col_size = m_field_offset.back();

  auto ov_grid = m_remap_data->m_overlap_grid;
  if (coarsen) {
    m_imp_exp = std::make_shared<GridImportExport>(m_tgt_grid,ov_grid);
  } else {
    m_imp_exp = std::make_shared<GridImportExport>(m_src_grid,ov_grid);
  }
  const int ncols_ov = ov_grid->get_num_local_dofs();

  // ----------- Compute RECV metadata -------------- //

  // We can now compute the offset of each pid in the recv buffer
  m_pids_recv_offsets = view_1d<int>("",nranks+1);
  auto ncols_recv_h = coarsen ? m_imp_exp->num_exports_per_pid_h()
                              : m_imp_exp->num_imports_per_pid_h();
  auto pids_recv_offsets_h = Kokkos::create_mirror_view(m_pids_recv_offsets);
  pids_recv_offsets_h[0] = 0;
  for (int pid=0; pid<nranks; ++pid) {
    pids_recv_offsets_h(pid+1) = pids_recv_offsets_h(pid)
                               + ncols_recv_h(pid);
  }
  Kokkos::deep_copy(m_pids_recv_offsets,pids_recv_offsets_h);

  // Create the recv buffer(s)
  auto recv_buf_size = coarsen ? pids_recv_offsets_h(nranks)*total_col_size
                               : ncols_ov*total_col_size;
  m_recv_buffer = decltype(m_recv_buffer)("HorizontalRemapper::recv_buf",recv_buf_size);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // ----------- Compute SEND metadata -------------- //
  
  m_pids_send_offsets = view_1d<int>("",nranks+1);
  auto ncols_send_h = coarsen ? m_imp_exp->num_imports_per_pid_h()
                              : m_imp_exp->num_exports_per_pid_h();
  auto pids_send_offsets_h = Kokkos::create_mirror_view(m_pids_send_offsets);
  pids_send_offsets_h[0] = 0;
  for (int pid=0; pid<nranks; ++pid) {
    pids_send_offsets_h(pid+1) = pids_send_offsets_h(pid)
                               + ncols_send_h(pid);
  }
  Kokkos::deep_copy(m_pids_send_offsets,pids_send_offsets_h);
  
  // Create the send buffer(s)
  auto send_buf_size = coarsen ? ncols_ov*total_col_size
                               : pids_send_offsets_h(nranks)*total_col_size;
  m_send_buffer = decltype(m_send_buffer)("HorizontalRemapper::send_buf",send_buf_size);
  m_mpi_send_buffer = Kokkos::create_mirror_view(decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);

  // ----------- Create Requests ------------ //

  const auto mpi_comm = m_src_grid->get_comm().mpi_comm();
  const auto mpi_real  = ekat::get_mpi_type<Real>();
  for (int pid=0; pid<nranks; ++pid) {
    // Send request
    if (ncols_send_h(pid)>0) {
      auto send_ptr = m_mpi_send_buffer.data() + pids_send_offsets_h(pid)*total_col_size;
      auto send_count = ncols_send_h(pid)*total_col_size;
      auto& req = m_send_req.emplace_back();
      MPI_Send_init (send_ptr, send_count, mpi_real, pid,
                     0, mpi_comm, &req);
    }
    // Recv request
    if (ncols_recv_h(pid)>0) {
      auto recv_ptr = m_mpi_recv_buffer.data() + pids_recv_offsets_h(pid)*total_col_size;
      auto recv_count = ncols_recv_h(pid)*total_col_size;
      auto& req = m_recv_req.emplace_back();
      MPI_Recv_init (recv_ptr, recv_count, mpi_real, pid,
                     0, mpi_comm, &req);
    }
  }
}

void HorizontalRemapper::clean_up ()
{
  // Free MPI requests
  for (auto& req : m_send_req)
    MPI_Request_free(&req);
  for (auto& req : m_recv_req)
    MPI_Request_free(&req);

  // Clear all MPI related structures
  m_send_buffer     = view_1d<Real>();
  m_recv_buffer     = view_1d<Real>();
  m_mpi_send_buffer = mpi_view_1d<Real>();
  m_mpi_recv_buffer = mpi_view_1d<Real>();
  m_send_req.clear();
  m_recv_req.clear();
  m_imp_exp = nullptr;

  // Clear all fields
  m_src_fields.clear();
  m_tgt_fields.clear();
  m_ov_fields.clear();
  m_needs_remap.clear();

  // Reset the state of the base class
  m_state = RepoState::Clean;
  m_num_fields = 0;
}

} // namespace scream
