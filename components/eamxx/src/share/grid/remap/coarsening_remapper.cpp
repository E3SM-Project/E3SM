#include "coarsening_remapper.hpp"

#include "share/grid/point_grid.hpp"
#include "share/grid/grid_import_export.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/field/field.hpp"

#include <ekat/kokkos/ekat_kokkos_utils.hpp>
#include <ekat/ekat_pack_utils.hpp>

#include <numeric>

namespace scream
{

CoarseningRemapper::
CoarseningRemapper (const grid_ptr_type& src_grid,
                    const std::string& map_file,
                    const bool track_mask,
                    const bool populate_tgt_grid_geo_data)
 : HorizInterpRemapperBase (src_grid,map_file,InterpType::Coarsen)
 , m_track_mask (track_mask)
{
  using namespace ShortFieldTagsNames;

  if (populate_tgt_grid_geo_data) {
    // Replicate the src grid geo data in the tgt grid. We use this remapper to do
    // the remapping (if needed), and clean it up afterwards.
    const auto& src_geo_data_names = src_grid->get_geometry_data_names();
    registration_begins();
    for (const auto& name : src_geo_data_names) {
      // Since different remappers may share the same data (if the map file is the same)
      // the coarse grid may already have the geo data.
      if (m_coarse_grid->has_geometry_data(name)) {
        continue;
      }
      const auto& src_data = src_grid->get_geometry_data(name);
      const auto& src_data_fid = src_data.get_header().get_identifier();
      const auto& layout = src_data_fid.get_layout();
      if (layout.tags().empty()) {
        // This is a scalar field, so won't be coarsened.
        // Simply copy it in the tgt grid, but we still need to assign the new grid name.
        FieldIdentifier tgt_data_fid(src_data_fid.name(),src_data_fid.get_layout(),src_data_fid.get_units(),m_tgt_grid->name());
        auto tgt_data = m_coarse_grid->create_geometry_data(tgt_data_fid);
        tgt_data.deep_copy(src_data);
      } else if (layout.tags()[0]!=COL) {
        // Not a field to be coarsened (perhaps a vertical coordinate field).
        // Simply copy it in the tgt grid, but we still need to assign the new grid name.
        FieldIdentifier tgt_data_fid(src_data_fid.name(),src_data_fid.get_layout(),src_data_fid.get_units(),m_tgt_grid->name());
        auto tgt_data = m_coarse_grid->create_geometry_data(tgt_data_fid);
        tgt_data.deep_copy(src_data);
      } else {
        // This field needs to be remapped
        auto tgt_data_fid = create_tgt_fid(src_data_fid);
        auto tgt_data = m_coarse_grid->create_geometry_data(tgt_data_fid);
        register_field(src_data,tgt_data);
      }
    }
    registration_ends();
    if (get_num_fields()>0) {
      remap(true);

      // The remap phase only alters the fields on device.
      // We need to sync them to host as well
      for (int i=0; i<get_num_fields(); ++i) {
        auto tgt_data = get_tgt_field(i);
        tgt_data.sync_to_host();
      }
    }
    clean_up();
  }
}

CoarseningRemapper::
~CoarseningRemapper ()
{
  // We need to free MPI requests
  for (size_t i=0; i<m_send_req.size(); ++i) {
    MPI_Request_free(&m_send_req[i]);
  }
  for (size_t i=0; i<m_recv_req.size(); ++i) {
    MPI_Request_free(&m_recv_req[i]);
  }
}

void CoarseningRemapper::
do_bind_field (const int ifield, const field_type& src, const field_type& tgt)
{
  // Assume no mask tracking for this field. Can correct below
  m_field_idx_to_mask_idx[ifield] = -1;

  if (m_track_mask) {
    if (src.get_header().has_extra_data("mask_data")) {
      // First, check that we also have the mask value, to be used if mask_data is too small
      EKAT_REQUIRE_MSG (src.get_header().has_extra_data("mask_value"),
          "Error! Field " + src.name() + " stores a mask field but not a mask value.\n");
      const auto& src_mask_val = src.get_header().get_extra_data<Real>("mask_value");

      Field tgt_copy = tgt;

      auto& tgt_hdr = tgt_copy.get_header();
      if (tgt_hdr.has_extra_data("mask_value")) {
        const auto& tgt_mask_val = tgt_hdr.get_extra_data<Real>("mask_value");

        EKAT_REQUIRE_MSG (tgt_mask_val==src_mask_val,
            "Error! Target field stores a mask data different from the src field.\n"
            "  - src field name: " + src.name() + "\n"
            "  - tgt field name: " + tgt.name() + "\n"
            "  - src mask value: " << src_mask_val << "\n"
            "  - tgt mask value: " << tgt_mask_val << "\n");
      } else {
        tgt_hdr.set_extra_data("mask_value",src_mask_val);
      }

      // Then, register the mask field, if not yet registered
      const auto& src_mask = src.get_header().get_extra_data<Field>("mask_data");
      const auto& src_mask_fid = src_mask.get_header().get_identifier();
      // Make sure fields representing masks are not themselves meant to be masked.
      EKAT_REQUIRE_MSG(not src_mask.get_header().has_extra_data("mask_data"),
          "Error! A mask field cannot be itself masked.\n"
          "  - field name: " + src.name() + "\n"
          "  - mask field name: " + src_mask.name() + "\n");

      const auto tgt_mask_fid = create_tgt_fid(src_mask_fid);
      if (not has_src_field (src_mask_fid)) {
        Field tgt_mask(tgt_mask_fid);
        tgt_mask.get_header().get_alloc_properties().request_allocation(src_mask.get_header().get_alloc_properties().get_largest_pack_size());
        tgt_mask.allocate_view();
        register_field(src_mask,tgt_mask);
      }

      // Store position of the mask for this field
      const int mask_idx = find_field(src_mask_fid,tgt_mask_fid);
      m_field_idx_to_mask_idx[ifield] = mask_idx;

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
    }
  }
  HorizInterpRemapperBase::do_bind_field(ifield,src,tgt);
}

void CoarseningRemapper::do_remap_fwd ()
{
  // Fire the recv requests right away, so that if some other ranks
  // is done packing before us, we can start receiving their data
  if (not m_recv_req.empty()) {
    int ierr = MPI_Startall(m_recv_req.size(),m_recv_req.data());
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while starting persistent recv requests.\n"
        "  - recv rank: " + std::to_string(m_comm.rank()) + "\n");
  }

  // TODO: Add check that if there are mask values they are either 1's or 0's for unmasked/masked.

  // Helpef function, to establish if a field can be handled with packs
  auto can_pack_field = [](const Field& f) {
    const auto& ap = f.get_header().get_alloc_properties();
    return (ap.get_last_extent() % SCREAM_PACK_SIZE) == 0;
  };

  // Loop over each field
  for (int i=0; i<m_num_fields; ++i) {
    // First, perform the local mat-vec. Recall that in these y=Ax products,
    // x is the src field, and y is the overlapped tgt field.
    const auto& f_src = m_src_fields[i];
    const auto& f_ov  = m_ov_fields[i];

    const int mask_idx = m_field_idx_to_mask_idx[i];
    if (mask_idx>0) {
      // Pass the mask to the local_mat_vec routine
      const auto& mask = m_src_fields[mask_idx];

      // If possible, dispatch kernel with SCREAM_PACK_SIZE
      if (can_pack_field(f_src) and can_pack_field(f_ov) and can_pack_field(mask)) {
        local_mat_vec<SCREAM_PACK_SIZE>(f_src,f_ov,mask);
      } else {
        local_mat_vec<1>(f_src,f_ov,mask);
      }
    } else {
      // If possible, dispatch kernel with SCREAM_PACK_SIZE
      if (can_pack_field(f_src) and can_pack_field(f_ov)) {
        local_mat_vec<SCREAM_PACK_SIZE>(f_src,f_ov);
      } else {
        local_mat_vec<1>(f_src,f_ov);
      }
    }
  }

  // Pack, then fire off the sends
  pack_and_send ();

  // Wait for all data to be received, then unpack
  recv_and_unpack ();

  // Wait for all sends to be completed
  if (not m_send_req.empty()) {
    int ierr = MPI_Waitall(m_send_req.size(),m_send_req.data(), MPI_STATUSES_IGNORE);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while waiting on persistent send requests.\n"
        "  - send rank: " + std::to_string(m_comm.rank()) + "\n");
  }

  // Rescale any fields that had the mask applied.
  if (m_track_mask) {
    for (int i=0; i<m_num_fields; ++i) {
      const auto& f_tgt = m_tgt_fields[i];
      const int mask_idx = m_field_idx_to_mask_idx[i];
      if (mask_idx>0) {
        // Then this field did use a mask
        const auto& mask = m_tgt_fields[mask_idx];
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
void CoarseningRemapper::
rescale_masked_fields (const Field& x, const Field& mask) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto& layout = x.get_header().get_identifier().get_layout();
  const int rank = layout.rank();
  const int ncols = m_tgt_grid->get_num_local_dofs();
  Real mask_val = std::numeric_limits<float>::max()/10.0;
  if (x.get_header().has_extra_data("mask_value")) {
    mask_val = x.get_header().get_extra_data<Real>("mask_value");
  } else {
    EKAT_ERROR_MSG ("ERROR! Field " + x.name() + " is masked, but stores no mask_value extra data.\n");
  }
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
          x_view(icol) = mask_val;
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
      auto policy = ESU::get_default_team_policy(ncols,dim1);
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
            x_sub(j).set(!masked,mask_val);
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
      auto policy = ESU::get_default_team_policy(ncols,dim1*dim2);
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
            x_sub(k).set(!masked,mask_val);
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
      auto policy = ESU::get_default_team_policy(ncols,dim1*dim2*dim3);
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
            x_sub(l).set(!masked,mask_val);
          });
        }
      });
      break;
    }
  }
}

template<int PackSize>
void CoarseningRemapper::
local_mat_vec (const Field& x, const Field& y, const Field& mask) const
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;
  using Pack        = ekat::Pack<Real,PackSize>;
  using PackInfo    = ekat::PackInfo<PackSize>;

  const auto& src_layout = x.get_header().get_identifier().get_layout();
  const int rank = src_layout.rank();
  const int nrows = m_ov_coarse_grid->get_num_local_dofs();
  auto row_offsets = m_row_offsets;
  auto col_lids = m_col_lids;
  auto weights = m_weights;
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
      auto mask_view = mask.get_strided_view<Real*>();
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
      auto policy = ESU::get_default_team_policy(nrows,dim1);
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
      auto policy = ESU::get_default_team_policy(nrows,dim1*dim2);
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
      auto policy = ESU::get_default_team_policy(nrows,dim1*dim2*dim3);
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

void CoarseningRemapper::pack_and_send ()
{
  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const int num_send_gids = m_ov_coarse_grid->get_num_local_dofs();
  const auto pid_lid_start = m_send_pid_lids_start;
  const auto lids_pids = m_send_lids_pids;
  const auto buf = m_send_buffer;

  for (int ifield=0; ifield<m_num_fields; ++ifield) {
    const auto& f  = m_ov_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto f_pid_offsets = ekat::subview(m_send_f_pid_offsets,ifield);

    switch (fl.rank()) {
      case 1:
      {
        // Unlike get_view, get_strided_view returns a LayoutStride view,
        // therefore allowing the 1d field to be a subfield of a 2d field
        // along the 2nd dimension.
        auto v = f.get_strided_view<const Real*>();
        Kokkos::parallel_for(RangePolicy(0,num_send_gids),
                             KOKKOS_LAMBDA(const int& i){
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          buf (offset + lidpos) = v(lid);
        });
      } break;
      case 2:
      {
        auto v = f.get_view<const Real**>();
        const int dim1 = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_send_gids,dim1);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                               [&](const int idim) {
            buf(offset + lidpos*dim1 + idim) = v(lid,idim);
          });
        });
      } break;
      case 3:
      {
        auto v = f.get_view<const Real***>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        auto policy = ESU::get_default_team_policy(num_send_gids,dim1*dim2);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2),
                               [&](const int idx) {
            const int idim = idx / dim2;
            const int ilev = idx % dim2;
            buf(offset + lidpos*dim1*dim2 + idim*dim2 + ilev) = v(lid,idim,ilev);
          });
        });
      } break;
      case 4:
      {
        auto v = f.get_view<const Real****>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int dim3 = fl.dim(3);
        auto policy = ESU::get_default_team_policy(num_send_gids,dim1*dim2*dim3);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int i = team.league_rank();
          const int lid = lids_pids(i,0);
          const int pid = lids_pids(i,1);
          const int lidpos = i - pid_lid_start(pid);
          const int offset = f_pid_offsets(pid);

          Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2*dim3),
                               [&](const int idx) {
            const int idim = (idx / dim3) / dim2;
            const int jdim = (idx / dim3) % dim2;
            const int ilev =  idx % dim3;
            buf(offset + lidpos*dim1*dim2*dim3 + idim*dim2*dim3 + jdim*dim3 + ilev) = v(lid,idim,jdim,ilev);
          });
        });
      } break;

      default:
        EKAT_ERROR_MSG ("Unexpected field rank in CoarseningRemapper::pack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
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
        "Error! Something whent wrong while starting persistent send requests.\n"
        "  - send rank: " + std::to_string(m_comm.rank()) + "\n");
  }
}

void CoarseningRemapper::recv_and_unpack ()
{
  if (not m_recv_req.empty()) {
    int ierr = MPI_Waitall(m_recv_req.size(),m_recv_req.data(), MPI_STATUSES_IGNORE);
    EKAT_REQUIRE_MSG (ierr==MPI_SUCCESS,
        "Error! Something whent wrong while waiting on persistent recv requests.\n"
        "  - recv rank: " + std::to_string(m_comm.rank()) + "\n");
  }
  // If MPI does not use dev pointers, we need to deep copy from host to dev
  if (not MpiOnDev) {
    Kokkos::deep_copy (m_recv_buffer,m_mpi_recv_buffer);
  }

  using RangePolicy = typename KT::RangePolicy;
  using MemberType  = typename KT::MemberType;
  using ESU         = ekat::ExeSpaceUtils<typename KT::ExeSpace>;

  const int num_tgt_dofs = m_tgt_grid->get_num_local_dofs();

  const auto buf = m_recv_buffer;
  const auto recv_lids_beg = m_recv_lids_beg;
  const auto recv_lids_end = m_recv_lids_end;
  const auto recv_lids_pidpos = m_recv_lids_pidpos;
  for (int ifield=0; ifield<m_num_fields; ++ifield) {
          auto& f  = m_tgt_fields[ifield];
    const auto& fl = f.get_header().get_identifier().get_layout();
    const auto f_pid_offsets = ekat::subview(m_recv_f_pid_offsets,ifield);

    f.deep_copy(0);
    switch (fl.rank()) {
      case 1:
      {
        // Unlike get_view, get_strided_view returns a LayoutStride view,
        // therefore allowing the 1d field to be a subfield of a 2d field
        // along the 2nd dimension.
        auto v = f.get_strided_view<Real*>();
        Kokkos::parallel_for(RangePolicy(0,num_tgt_dofs),
                             KOKKOS_LAMBDA(const int& lid){
          const int recv_beg = recv_lids_beg(lid);
          const int recv_end = recv_lids_end(lid);
          for (int irecv=recv_beg; irecv<recv_end; ++irecv) {
            const int pid = recv_lids_pidpos(irecv,0);
            const int lidpos = recv_lids_pidpos(irecv,1);
            const int offset = f_pid_offsets(pid) + lidpos;
            v(lid) += buf (offset);
          }
        });
      } break;
      case 2:
      {
        auto v = f.get_view<Real**>();
        const int dim1 = fl.dim(1);
        auto policy = ESU::get_default_team_policy(num_tgt_dofs,dim1);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int lid = team.league_rank();
          const int recv_beg = recv_lids_beg(lid);
          const int recv_end = recv_lids_end(lid);
          for (int irecv=recv_beg; irecv<recv_end; ++irecv) {
            const int pid = recv_lids_pidpos(irecv,0);
            const int lidpos = recv_lids_pidpos(irecv,1);
            const int offset = f_pid_offsets(pid)+lidpos*dim1;
            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1),
                                 [&](const int idim) {
              v(lid,idim) += buf (offset + idim);
            });
          }
        });
      } break;
      case 3:
      {
        auto v = f.get_view<Real***>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dims().back();
        auto policy = ESU::get_default_team_policy(num_tgt_dofs,dim2*dim1);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int lid = team.league_rank();
          const int recv_beg = recv_lids_beg(lid);
          const int recv_end = recv_lids_end(lid);
          for (int irecv=recv_beg; irecv<recv_end; ++irecv) {
            const int pid = recv_lids_pidpos(irecv,0);
            const int lidpos = recv_lids_pidpos(irecv,1);
            const int offset = f_pid_offsets(pid) + lidpos*dim1*dim2;

            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim2*dim1),
                                 [&](const int idx) {
              const int idim = idx / dim2;
              const int ilev = idx % dim2;
              v(lid,idim,ilev) += buf (offset + idim*dim2 + ilev);
            });
          }
        });
      } break;

      case 4:
      {
        auto v = f.get_view<Real****>();
        const int dim1 = fl.dim(1);
        const int dim2 = fl.dim(2);
        const int dim3 = fl.dim(3);
        auto policy = ESU::get_default_team_policy(num_tgt_dofs,dim1*dim2*dim3);
        Kokkos::parallel_for(policy,
                             KOKKOS_LAMBDA(const MemberType& team){
          const int lid = team.league_rank();
          const int recv_beg = recv_lids_beg(lid);
          const int recv_end = recv_lids_end(lid);
          for (int irecv=recv_beg; irecv<recv_end; ++irecv) {
            const int pid = recv_lids_pidpos(irecv,0);
            const int lidpos = recv_lids_pidpos(irecv,1);
            const int offset = f_pid_offsets(pid) + lidpos*dim1*dim2*dim3;

            Kokkos::parallel_for(Kokkos::TeamVectorRange(team,dim1*dim2*dim3),
                                 [&](const int idx) {
              const int idim = (idx / dim3) / dim2;
              const int jdim = (idx / dim3) % dim2;
              const int ilev =  idx % dim3;
              v(lid,idim,jdim,ilev) += buf (offset + idim*dim2*dim3 + jdim*dim3 + ilev);
            });
          }
        });
      } break;

      default:
        EKAT_ERROR_MSG ("Unexpected field rank in CoarseningRemapper::pack.\n"
            "  - MPI rank  : " + std::to_string(m_comm.rank()) + "\n"
            "  - field rank: " + std::to_string(fl.rank()) + "\n");
    }
  }
}

std::vector<int>
CoarseningRemapper::get_pids_for_recv (const std::vector<int>& send_to_pids) const
{
  // Figure out how many sends each PID is doing
  std::vector<int> num_sends(m_comm.size(),0);
  num_sends[m_comm.rank()] = send_to_pids.size();
  m_comm.all_gather(num_sends.data(),1);

  // Offsets for send_pids coming from each pid
  // NOTE: the extra entry at the end if for ease of use later
  std::vector<int> sends_offsets (m_comm.size()+1,0);
  for (int pid=1; pid<=m_comm.size(); ++pid) {
    sends_offsets[pid] = sends_offsets[pid-1] + num_sends[pid-1];
  }

  // Gather all the pids that each rank is sending data to
  auto nglobal_sends = std::accumulate(num_sends.begin(),num_sends.end(),0);
  std::vector<int> global_send_pids(nglobal_sends,-1);
  MPI_Allgatherv (send_to_pids.data(),send_to_pids.size(),MPI_INT,
                  global_send_pids.data(),num_sends.data(),sends_offsets.data(),
                  MPI_INT, m_comm.mpi_comm());

  // Loop over all the ranks, and all the pids they send to, and look for my pid
  std::vector<int> recv_from_pids;
  for (int pid=0; pid<m_comm.size(); ++pid) {
    const int beg = sends_offsets[pid];
    const int end = sends_offsets[pid+1];

    for (int i=beg; i<end; ++i) {
      if (global_send_pids[i]==m_comm.rank()) {
        recv_from_pids.push_back(pid);
        break;
      }
    }
  }

  return recv_from_pids;
}

std::map<int,std::vector<int>>
CoarseningRemapper::
recv_gids_from_pids (const std::map<int,std::vector<int>>& pid2gids_send) const
{
  const auto comm = m_comm.mpi_comm();

  // First, figure out which PIDs I need to recv from
  std::vector<int> send_to;
  for (const auto& it : pid2gids_send) {
    send_to.push_back(it.first);
  }
  auto recv_from = get_pids_for_recv (send_to);

  // Send num of lids send to each send_pid, and recv num of lids
  // recv from each recv_pid
  std::map<int,int> nsends, nrecvs;
  std::vector<MPI_Request> send_req, recv_req;
  for (const auto& it : pid2gids_send) {
    const int pid = it.first;
    nsends[pid] = it.second.size();
    send_req.emplace_back();
    MPI_Isend (&nsends[pid],1,MPI_INT,pid,0,comm,&send_req.back());
  }
  for (auto pid : recv_from) {
    recv_req.emplace_back();
    MPI_Irecv(&nrecvs[pid],1,MPI_INT,pid,0,comm,&recv_req.back());
  }
  MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE);

  send_req.resize(0);
  recv_req.resize(0);

  // Now that we know how many gids we will get from each pid, we can send/recv them
  for (const auto& it : pid2gids_send) {
    send_req.emplace_back();
    MPI_Isend(it.second.data(),it.second.size(),MPI_INT,
              it.first,0,comm,&send_req.back());
  }
  std::map<int,std::vector<int>> pid2gids_recv;
  for (const auto& it : nrecvs) {
    const int pid = it.first;
    const int n   = it.second;
    pid2gids_recv[pid].resize(n);
    recv_req.emplace_back();
    MPI_Irecv(pid2gids_recv[pid].data(),n,MPI_INT,
              pid,0,comm,&recv_req.back());
  }
  MPI_Waitall(send_req.size(),send_req.data(),MPI_STATUSES_IGNORE);
  MPI_Waitall(recv_req.size(),recv_req.data(),MPI_STATUSES_IGNORE);

  return pid2gids_recv;
}

void CoarseningRemapper::setup_mpi_data_structures ()
{
  using namespace ShortFieldTagsNames;
  using gid_type = AbstractGrid::gid_type;

  const auto mpi_comm  = m_comm.mpi_comm();
  const auto mpi_real  = ekat::get_mpi_type<Real>();

  const int last_rank = m_comm.size()-1;

  // Pre-compute the amount of data stored in each field on each dof
  std::vector<int> field_col_size (m_num_fields);
  int sum_fields_col_sizes = 0;
  for (int i=0; i<m_num_fields; ++i) {
    const auto& f  = m_src_fields[i];
    const auto& fl = f.get_header().get_identifier().get_layout();
    field_col_size[i] = fl.clone().strip_dim(COL).size();
    sum_fields_col_sizes += field_col_size[i];
  }

  // --------------------------------------------------------- //
  //                   Setup SEND structures                   //
  // --------------------------------------------------------- //

  // 1. Retrieve pid (and associated lid) of all ov gids
  //    on the tgt grid
  const auto ov_gids = m_ov_coarse_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  auto gids_owners = m_tgt_grid->get_owners (ov_gids);

  // 2. Group dofs to send by remote pid
  const int num_ov_gids = ov_gids.size();
  std::map<int,std::vector<int>> pid2lids_send;
  std::map<int,std::vector<gid_type>> pid2gids_send;
  for (int i=0; i<num_ov_gids; ++i) {
    const int pid = gids_owners[i];
    pid2lids_send[pid].push_back(i);
    pid2gids_send[pid].push_back(ov_gids(i));
  }
  const int num_send_pids = pid2lids_send.size();
  m_send_lids_pids = view_2d<int>("",num_ov_gids,2);
  m_send_pid_lids_start = view_1d<int>("",m_comm.size());
  auto send_lids_pids_h = Kokkos::create_mirror_view(m_send_lids_pids);
  auto send_pid_lids_start_h = Kokkos::create_mirror_view(m_send_pid_lids_start);
  for (int pid=0,pos=0; pid<m_comm.size(); ++pid) {
    send_pid_lids_start_h(pid) = pos;
    for (auto lid : pid2lids_send[pid]) {
      send_lids_pids_h(pos,0) = lid;
      send_lids_pids_h(pos++,1) = pid;
    }
  }
  Kokkos::deep_copy(m_send_lids_pids,send_lids_pids_h);
  Kokkos::deep_copy(m_send_pid_lids_start,send_pid_lids_start_h);

  // 3. Compute offsets in send buffer for each pid/field pair
  m_send_f_pid_offsets = view_2d<int>("",m_num_fields,m_comm.size());
  auto send_f_pid_offsets_h = Kokkos::create_mirror_view(m_send_f_pid_offsets);
  std::vector<int> send_pid_offsets(m_comm.size());
  for (int pid=0,pos=0; pid<m_comm.size(); ++pid) {
    send_pid_offsets[pid] = pos;
    for (int i=0; i<m_num_fields; ++i) {
      send_f_pid_offsets_h(i,pid) = pos;
      pos += field_col_size[i]*pid2lids_send[pid].size();
    }

    // At the end, pos must match the total amount of data in the overlapped fields
    if (pid==last_rank) {
      EKAT_REQUIRE_MSG (pos==num_ov_gids*sum_fields_col_sizes,
          "Error! Something went wrong in CoarseningRemapper::setup_mpi_structures.\n");
    }
  }
  Kokkos::deep_copy (m_send_f_pid_offsets,send_f_pid_offsets_h);

  // 4. Allocate send buffers
  m_send_buffer = view_1d<Real>("",sum_fields_col_sizes*num_ov_gids);
  m_mpi_send_buffer = Kokkos::create_mirror_view(decltype(m_mpi_send_buffer)::execution_space(),m_send_buffer);

  // 5. Setup send requests
  m_send_req.reserve(num_send_pids);
  for (const auto& it : pid2lids_send) {
    const int n = it.second.size()*sum_fields_col_sizes;
    if (n==0) {
      continue;
    }

    const int pid = it.first;
    const auto send_ptr = m_mpi_send_buffer.data() + send_pid_offsets[pid];

    m_send_req.emplace_back();
    auto& req = m_send_req.back();
    MPI_Send_init (send_ptr, n, mpi_real, pid,
                   0, mpi_comm, &req);
  }

  // --------------------------------------------------------- //
  //                   Setup RECV structures                   //
  // --------------------------------------------------------- //

  // 1. Obtain the dual map of send_gids: a list of gids we need to
  //    receive, grouped by the pid we recv them from
  const int num_tgt_dofs = m_tgt_grid->get_num_local_dofs();
  auto pid2gids_recv = recv_gids_from_pids(pid2gids_send);
  const int num_recv_pids = pid2gids_recv.size();

  // 2. Convert the gids to lids, and arrange them by lid
  std::vector<std::vector<int>> lid2pids_recv(num_tgt_dofs);
  int num_total_recv_gids = 0;
  auto tgt_gid2lid = m_tgt_grid->get_gid2lid_map();
  for (const auto& it : pid2gids_recv) {
    const int pid = it.first;
    for (auto gid : it.second) {
      const int lid = tgt_gid2lid[gid];
      lid2pids_recv[lid].push_back(pid);
    }
    num_total_recv_gids += it.second.size();
  }

  // 3. Splice the vector-of-vectors above in a 1d view,
  //    keeping track of where each lid starts/ends
  m_recv_lids_pidpos = view_2d<int>("",num_total_recv_gids,2);
  m_recv_lids_beg = view_1d<int>("",num_tgt_dofs);
  m_recv_lids_end = view_1d<int>("",num_tgt_dofs);
  auto recv_lids_pidpos_h = Kokkos::create_mirror_view(m_recv_lids_pidpos);
  auto recv_lids_beg_h  = Kokkos::create_mirror_view(m_recv_lids_beg);
  auto recv_lids_end_h  = Kokkos::create_mirror_view(m_recv_lids_end);

  auto tgt_dofs_h = m_tgt_grid->get_dofs_gids().get_view<const gid_type*,Host>();
  for (int i=0,pos=0; i<num_tgt_dofs; ++i) {
    recv_lids_beg_h(i) = pos;
    const int gid = tgt_dofs_h[i];
    for (auto pid : lid2pids_recv[i]) {
      auto it = std::find(pid2gids_recv.at(pid).begin(),pid2gids_recv.at(pid).end(),gid);
      EKAT_REQUIRE_MSG (it!=pid2gids_recv.at(pid).end(),
          "Error! Something went wrong in CoarseningRemapper::setup_mpi_structures.\n");
      recv_lids_pidpos_h(pos,0) = pid;
      recv_lids_pidpos_h(pos++,1) = std::distance(pid2gids_recv.at(pid).begin(),it);
    }
    recv_lids_end_h(i) = pos;
  }
  Kokkos::deep_copy(m_recv_lids_pidpos,recv_lids_pidpos_h);
  Kokkos::deep_copy(m_recv_lids_beg,recv_lids_beg_h);
  Kokkos::deep_copy(m_recv_lids_end,recv_lids_end_h);

  // 3. Splice gids from all pids, and convert to lid, and compute the offset of each pid
  //    in the spliced list of gids.
  std::vector<int> recv_pid_start(m_comm.size()+1);
  for (int pid=0,pos=0; pid<=m_comm.size(); ++pid) {
    recv_pid_start[pid] = pos;
    pos += pid2gids_recv[pid].size();
  }

  // 4. Compute offsets in recv buffer for each pid/field pair
  m_recv_f_pid_offsets = view_2d<int>("",m_num_fields,m_comm.size());
  auto recv_f_pid_offsets_h = Kokkos::create_mirror_view(m_recv_f_pid_offsets);
  std::vector<int> recv_pid_offsets(m_comm.size());
  for (int pid=0,pos=0; pid<m_comm.size(); ++pid) {
    recv_pid_offsets[pid] = pos;
    const int num_recv_gids = recv_pid_start[pid+1] - recv_pid_start[pid];
    for (int i=0; i<m_num_fields; ++i) {
      recv_f_pid_offsets_h(i,pid) = pos;
      pos += field_col_size[i]*num_recv_gids;
    }

    // // At the end, pos must match the total amount of data received
    if (pid==last_rank) {
      EKAT_REQUIRE_MSG (pos==num_total_recv_gids*sum_fields_col_sizes,
          "Error! Something went wrong in CoarseningRemapper::setup_mpi_structures.\n");
    }
  }
  Kokkos::deep_copy (m_recv_f_pid_offsets,recv_f_pid_offsets_h);

  // 5. Allocate recv buffers
  m_recv_buffer = view_1d<Real>("",sum_fields_col_sizes*num_total_recv_gids);
  m_mpi_recv_buffer = Kokkos::create_mirror_view(decltype(m_mpi_recv_buffer)::execution_space(),m_recv_buffer);

  // 6. Setup recv requests
  m_recv_req.reserve(num_recv_pids);
  for (int pid=0; pid<m_comm.size(); ++pid) {
    const int num_recv_gids = recv_pid_start[pid+1] - recv_pid_start[pid];
    const int n = num_recv_gids*sum_fields_col_sizes;
    if (n==0) {
      continue;
    }

    const auto recv_ptr = m_mpi_recv_buffer.data() + recv_pid_offsets[pid];

    m_recv_req.emplace_back();
    auto& req = m_recv_req.back();
    MPI_Recv_init (recv_ptr, n, mpi_real, pid,
                   0, mpi_comm, &req);
  }
}

void CoarseningRemapper::clean_up ()
{
  // Clear all MPI related structures
  m_send_buffer         = view_1d<Real>();
  m_recv_buffer         = view_1d<Real>();
  m_mpi_send_buffer     = mpi_view_1d<Real>();
  m_mpi_recv_buffer     = mpi_view_1d<Real>();
  m_send_f_pid_offsets  = view_2d<int>();
  m_recv_f_pid_offsets  = view_2d<int>();
  m_send_lids_pids      = view_2d<int>();
  m_send_pid_lids_start = view_1d<int>();
  m_recv_lids_pidpos    = view_2d<int>();
  m_recv_lids_beg       = view_1d<int>();
  m_recv_lids_end       = view_1d<int>();
  m_send_req.clear();
  m_recv_req.clear();

  HorizInterpRemapperBase::clean_up();
}

} // namespace scream
