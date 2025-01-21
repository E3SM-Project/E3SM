#include "share/grid/remap/iop_remapper.hpp"

// Anonymous namespace, in case other TU's define the same type
namespace {
struct RealInt
{
  scream::Real val;
  int  idx;
};
}
// Specialize ekat's get_mpi_type for RealInt (needed for the all_reduce call)
namespace ekat {
  template<>
  MPI_Datatype get_mpi_type<RealInt> () {
#ifdef SCREAM_DOUBLE_PRECISION
    return MPI_DOUBLE_INT;
#else
    return MPI_FLOAT_INT;
#endif
  }
} // namespace ekat

namespace scream
{

IOPRemapper::
IOPRemapper (const grid_ptr_type src_grid,
             const grid_ptr_type tgt_grid,
             const Real lat, const Real lon)
 : AbstractRemapper (src_grid,tgt_grid)
{
  m_bwd_allowed = false;

  EKAT_REQUIRE_MSG (src_grid->has_geometry_data("lat") and src_grid->has_geometry_data("lon"),
      "Error! IOP remapper requires lat/lon geometry data in the source grid.\n");
  EKAT_REQUIRE_MSG (src_grid->get_num_vertical_levels()==tgt_grid->get_num_vertical_levels(),
      "Error! IOP remapper requires src/tgt grid to have the same number of vertical levels.\n");
  EKAT_REQUIRE_MSG (src_grid->type()==GridType::Point and tgt_grid->type()==GridType::Point,
      "Error! IOP remapper requires src/tgt grid to be PointGrid instances.\n");
  
  m_comm = src_grid->get_comm();

  setup_closest_col_info (lat,lon);
}

void IOPRemapper::
setup_closest_col_info (const Real lat, const Real lon)
{
  auto lat_f = m_src_grid->get_geometry_data("lat");
  auto lon_f = m_src_grid->get_geometry_data("lon");
  EKAT_REQUIRE_MSG (lat_f.get_header().get_identifier().get_layout().type()==LayoutType::Scalar2D,
      "Error! Source grid 'lat' field has the wrong layout.\n"
      " - expected layout: " + e2str(LayoutType::Scalar2D) + "\n"
      " - actual layout  : " + e2str(lat_f.get_header().get_identifier().get_layout().type()) + "\n");
  EKAT_REQUIRE_MSG (lon_f.get_header().get_identifier().get_layout().type()==LayoutType::Scalar2D,
      "Error! Source grid 'lon' field has the wrong layout.\n"
      " - expected layout: " + e2str(LayoutType::Scalar2D) + "\n"
      " - actual layout  : " + e2str(lon_f.get_header().get_identifier().get_layout().type()) + "\n");
  EKAT_REQUIRE_MSG(-90 <= lat and lat <= 90,
      "Error! IOPRemapper lat="+std::to_string(lat) +" outside of expected range [-90, 90].\n");
  EKAT_REQUIRE_MSG(0 <= lon and lon <= 360,
      "Error! IOPRemapper lon="+std::to_string(lon) +" outside of expected range [0, 360].\n");

  auto lat_v = lat_f.get_view<const Real*>();
  auto lon_v = lon_f.get_view<const Real*>();
  using minloc_t = Kokkos::MinLoc<Real,int>;
  using minloc_value_t = typename minloc_t::value_type;
  minloc_value_t minloc;

  auto lambda = KOKKOS_LAMBDA (int icol, minloc_value_t& result) {
    auto dist = Kokkos::abs(lat_v(icol)-lat)+std::abs(lon_v(icol)-lon);
    if(dist<result.val) {
      result.val = dist;
      result.loc = icol;
    }   
  };
  const int ncols = m_src_grid->get_num_local_dofs();
  Kokkos::parallel_reduce(ncols, lambda, minloc_t(minloc));

  RealInt dist_rank;
  dist_rank.val = minloc.val;
  dist_rank.idx = m_comm.rank();

  m_comm.all_reduce(&dist_rank,1,MPI_MINLOC);
  m_closest_col_info.mpi_rank = dist_rank.idx;
  if (dist_rank.idx==m_comm.rank()) {
    m_closest_col_info.col_lid = minloc.loc;
  }
}

void IOPRemapper::registration_ends_impl ()
{
  // NOTE: one could think to make the single-col fields to NOT be clones
  //       on the rank that owns the closest col. However, src fields MAY
  //       be read-only, and MPI_Bcast needs a pointer to non-const (since
  //       the receiving ranks must write on it).
  //       MOREOVER, by cloning, we do away with padding shenaningans, which
  //       may require a bit more attention for bcast operations.
  using namespace ShortFieldTagsNames;
  for (const auto& f : m_src_fields) {
    auto col = f.subfield(COL,0).clone();
    m_single_col_fields.push_back(col);
  }
}

void IOPRemapper::remap_fwd_impl ()
{
  using namespace ShortFieldTagsNames;

  const auto root_id  = m_closest_col_info.mpi_rank;
  const auto iam_root = root_id==m_comm.rank();

  // 1. Rank root_id extracts the closest col for all fields
  if (iam_root) {
    for (int i=0; i<m_num_fields; ++i) {
      auto& dst = m_single_col_fields[i];
      auto  src = m_src_fields[i].subfield(COL,m_closest_col_info.col_lid);
      dst.deep_copy(src);
    }
  }

  // 2. Rank root_id broadcasts the single-col fields
  for (int i=0; i<m_num_fields; ++i) {
    auto& f = m_single_col_fields[i];
    int col_size = f.get_header().get_identifier().get_layout().size();
#if SCREAM_MPI_ON_DEVICE
    m_comm.broadcast(f.get_internal_view_data<Real>(),col_size,root_id);
#else
    if (iam_root) {
      f.sync_to_host();
    }
    m_comm.broadcast(f.get_internal_view_data<Real,Host>(),col_size,root_id);
    if (not iam_root) {
      f.sync_to_dev();
    }
#endif
  }

  // 3. Every rank copies the single col into ALL cols of the tgt fields
  int ncols = m_tgt_grid->get_num_local_dofs();
  for (int i=0; i<m_num_fields; ++i) {
    auto& col = m_single_col_fields[i];
    auto& tgt = m_tgt_fields[i];

    // TODO: one may think of dispatching a TP kernel, to reduce latency.
    //       That's fine, but you make the remapper code longer, since you need
    //       to handle different ranks separately.
    for (int icol=0; icol<ncols; ++icol) {
      tgt.subfield(COL,icol).deep_copy(col);
    }
  }
}

} // namespace scream
