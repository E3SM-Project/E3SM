#include "share/grid/remap/horizontal_remap_utility.hpp"
#include "share/util/scream_timing.hpp"

namespace scream {

/*-----------------------------------------------------------------------------------------------*/
HorizontalMap::HorizontalMap(const ekat::Comm& comm)
  : m_comm (comm)
{
  m_dofs_set = false;
}
/*-----------------------------------------------------------------------------------------------*/
HorizontalMap::HorizontalMap(const ekat::Comm& comm, const std::string& map_name)
  : m_name (map_name)
  , m_comm (comm)
{
  m_dofs_set = false;
}
/*-----------------------------------------------------------------------------------------------*/
HorizontalMap::HorizontalMap(const ekat::Comm& comm, const std::string& map_name, const view_1d<const gid_type>& dofs_gids, const gid_type min_dof)
  : m_name (map_name)
  , m_comm (comm)
{
  set_dof_gids(dofs_gids,min_dof);
}
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// This function sets the remap segments for this map given a remap file that was created offline.
// Note: This function assumes that the remap file follows the convention:
//   col - This variable in the file represents the set of source dofs that map to a specific target.
//   row - This variable represents the corresponding list of target columns mapped to
//   S   - This variable stores the corresponding weights for each col -> row pair
//   n_s - Is the integer number of col -> row mappings.
//   for reference there may also be
//   n_a - Is the size of the source grid.
//   n_b - Is the size of the target grid.
// Following these conventions we assume that 
//   col - will be used to populate a segment's "source_dofs"
//   row - will be used to populate a segment's "m_dof"
//   S   - will be used to populate a segment's "weights"
void HorizontalMap::set_remap_segments_from_file(const std::string& remap_filename)
{
  start_timer("EAMxx::HorizontalMap::set_remap_segments_from_file");
  // Open remap file and determine the amount of data to be read
  scorpio::register_file(remap_filename,scorpio::Read);
  const auto remap_size = scorpio::get_dimlen_c2f(remap_filename.c_str(),"n_s"); // Note, here we assume a standard format of col, row, S
  // Step 1: Read in the "row" data from the file to figure out which mpi ranks care about which
  //         chunk of the remap data.  This step reduces the memory footprint of reading in the
  //         map data, which can be rather large.
  // Distribute responsibility for reading remap data over all ranks
  const int my_rank   = m_comm.rank();
  const int num_ranks = m_comm.size();
  // my_chunk will represent the chunk of data this rank will read from file.
  int my_chunk        = remap_size/num_ranks;
  int remainder       = remap_size - (my_chunk*num_ranks);
  if (remainder != 0) {
    my_chunk += my_rank<remainder ? 1 : 0;
  }
  // now determine where this rank start reading the data.
  std::vector<int> chunks_glob(num_ranks);
  m_comm.all_gather(&my_chunk,chunks_glob.data(),1);
  int my_start = 0;
  for (int ii=0; ii<my_rank;ii++) {
    my_start += chunks_glob[ii];
  }
  // Check that the total set of chunks covers all the data
  {
    int chunk_check = 0;
    for (int ii=0; ii<num_ranks; ii++) {
      chunk_check += chunks_glob[ii];
    }
    EKAT_REQUIRE_MSG(chunk_check==remap_size,"ERROR: HorizontalMap " + m_name +" get_remap_indices - Something went wrong distributing remap data among the MPI ranks");
  }
  // Using scream input routines, read remap data from file by chunk
  view_1d<int>  tgt_col("row",my_chunk); 
  auto tgt_col_h = Kokkos::create_mirror_view(tgt_col);
  std::vector<std::string> vec_of_dims = {"n_s"};
  std::string i_decomp = std::string("int-row-n_s-") + std::to_string(my_chunk);
  scorpio::get_variable(remap_filename, "row", "row", vec_of_dims, "int", i_decomp);
  std::vector<int64_t> var_dof(my_chunk);
  std::iota(var_dof.begin(),var_dof.end(),my_start);
  scorpio::set_dof(remap_filename,"row",var_dof.size(),var_dof.data());
  scorpio::set_decomp(remap_filename);
  scorpio::grid_read_data_array(remap_filename,"row",0,tgt_col_h.data(),tgt_col_h.size()); 
  scorpio::eam_pio_closefile(remap_filename);
  // Step 2: Now that we have the data distributed among all ranks we organize the data
  //         into sets of target column, start location in data and length of data.
  //         At the same time, determine the min_dof for remap column indices.
  std::vector<int> chunk_dof, chunk_start, chunk_len;
  chunk_dof.push_back(tgt_col_h(0));
  chunk_start.push_back(my_start);
  chunk_len.push_back(1);
  int remap_min_dof = tgt_col_h(0);
  for (int ii=1; ii<my_chunk; ii++) {
    remap_min_dof = std::min(tgt_col_h(ii),remap_min_dof);
    if (tgt_col_h(ii) == chunk_dof.back()) {
      // Then we add one to the length for this chunk.
      chunk_len.back() ++;
    } else {
      // Start a new chunk for a new DOF
      chunk_dof.push_back(tgt_col_h(ii));
      chunk_start.push_back(my_start+ii);
      chunk_len.push_back(1);
    }
  }
  // Pass chunk information among all ranks so they can be consolidated.
  int  num_chunks          = chunk_dof.size();
  std::vector<int> num_chunks_per_rank(num_ranks), chunk_displacement(num_ranks);
  int  total_num_chunks;
  int  global_remap_min_dof;
  m_comm.all_gather(&num_chunks, num_chunks_per_rank.data(),1);
  m_comm.all_reduce(&remap_min_dof,&global_remap_min_dof,1,MPI_MIN);
  chunk_displacement[0] = 0;
  total_num_chunks = num_chunks_per_rank[0];
  for (int ii=1; ii<num_ranks; ii++) {
    chunk_displacement[ii] = total_num_chunks;
    total_num_chunks += num_chunks_per_rank[ii];
  }
  std::vector<int> buff_dof(total_num_chunks), buff_sta(total_num_chunks), buff_len(total_num_chunks);
  MPI_Allgatherv(chunk_dof.data(),  chunk_dof.size(),MPI_INT,buff_dof.data(),num_chunks_per_rank.data(),chunk_displacement.data(),MPI_INT,m_comm.mpi_comm());
  MPI_Allgatherv(chunk_start.data(),chunk_dof.size(),MPI_INT,buff_sta.data(),num_chunks_per_rank.data(),chunk_displacement.data(),MPI_INT,m_comm.mpi_comm());
  MPI_Allgatherv(chunk_len.data(),  chunk_dof.size(),MPI_INT,buff_len.data(),num_chunks_per_rank.data(),chunk_displacement.data(),MPI_INT,m_comm.mpi_comm());
  // Step 3: Now that all of the ranks are aware of all of the "sets" of source -> target mappings we
  //         construct and add segments for just the DOF's this rank cares about.
  std::vector<int> seg_dof, seg_start, seg_length;
  var_dof.clear();
  auto dofs_gids_h = Kokkos::create_mirror_view(m_dofs_gids);
  Kokkos::deep_copy(dofs_gids_h,m_dofs_gids);
  for (int ii=0; ii<total_num_chunks; ii++) {
    // Search dofs to see if this chunk matches dofs on this rank
    for (int jj=0; jj<dofs_gids_h.extent_int(0); jj++) {
      if (buff_dof[ii]-global_remap_min_dof == dofs_gids_h(jj)) {
        std::vector<int> var_tmp(buff_len[ii]);
        std::iota(var_tmp.begin(),var_tmp.end(),buff_sta[ii]);
        seg_dof.push_back(buff_dof[ii]);
        seg_start.push_back(var_dof.size()); 
        seg_length.push_back(buff_len[ii]);
        var_dof.insert(var_dof.end(),var_tmp.begin(),var_tmp.end());
      } 
    }
  }
  // Now that we know which parts of the remap file this rank cares about we can construct segments
  view_1d<int>  col("col",var_dof.size()); 
  view_1d<Real> S("S",var_dof.size()); 
  auto col_h = Kokkos::create_mirror_view(col);
  auto S_h = Kokkos::create_mirror_view(S);
  vec_of_dims = {"n_s"};
  i_decomp = std::string("int-col-n_s-") + std::to_string(var_dof.size());
  std::string r_decomp = std::string("Real-S-n_s-") + std::to_string(var_dof.size());
  scorpio::register_file(remap_filename,scorpio::Read);
  scorpio::get_variable(remap_filename, "col", "col", vec_of_dims, "int", i_decomp);
  scorpio::get_variable(remap_filename, "S", "S", vec_of_dims, "real", r_decomp);
  scorpio::set_dof(remap_filename,"col",var_dof.size(),var_dof.data());
  scorpio::set_dof(remap_filename,"S",var_dof.size(),var_dof.data());
  scorpio::set_decomp(remap_filename);
  scorpio::grid_read_data_array(remap_filename,"col",0,col_h.data(),col_h.size()); 
  scorpio::grid_read_data_array(remap_filename,"S",0,S_h.data(),S_h.size()); 
  scorpio::eam_pio_closefile(remap_filename);
  Kokkos::deep_copy(col,col_h);
  Kokkos::deep_copy(S,S_h);
  // Construct segments based on data just read from file
  for (size_t ii=0; ii<seg_dof.size(); ii++) {
    int seglength = seg_length[ii];
    int segstart  = seg_start[ii];
    view_1d<gid_type> source_dofs("",seglength);
    view_1d<Real>     weights("",seglength);
    Kokkos::parallel_for("", seglength, KOKKOS_LAMBDA (const int& jj) {
      int idx = segstart + jj;
      source_dofs(jj) = col(idx)-global_remap_min_dof;  // Offset to zero based dofs
      weights(jj)     = S(idx);
    });
    HorizontalMapSegment seg(seg_dof[ii]-global_remap_min_dof,seglength,source_dofs,weights);
    add_remap_segment(seg);
  }
  stop_timer("EAMxx::HorizontalMap::set_remap_segments_from_file");
}
/*-----------------------------------------------------------------------------------------------*/
// This function is used to set the internal set of degrees of freedom (dof) this map is responsible for.
// We use the global dofs, offset by the minimum global dof to make everything zero based.  Note, when
// gathering remap parameters from a file, depending on the algorithm that made the file the dof
// indices may be 1-based or 0-based.  By offsetting everything to 0-based we avoid potential bugs.
void HorizontalMap::set_dof_gids(const view_1d<const gid_type>& dofs_gids, const gid_type min_dof)
{
  start_timer("EAMxx::HorizontalMap::set_dof_gids");
  EKAT_REQUIRE(dofs_gids.size()>0);
  m_dofs_gids = view_1d<gid_type>("",dofs_gids.size());
  const auto l_dofs_gids = m_dofs_gids;
  m_num_dofs = m_dofs_gids.extent(0);
  Kokkos::parallel_for("", m_num_dofs, KOKKOS_LAMBDA (const int& ii) {
    l_dofs_gids(ii) = dofs_gids(ii)-min_dof;
  });
  m_dofs_set = true;
  stop_timer("EAMxx::HorizontalMap::set_dof_gids");
}
/*-----------------------------------------------------------------------------------------------*/
// This function adds a remap segment to a HorizontalMap, note, we want each segment to represent a full
// remapping.  This function also checks if a segment already exists for the degree of freedom
// in question.  If it does then instead of add the segment to the end, this function finds that
// segment and combines them into a new comprehensive segment.
void HorizontalMap::add_remap_segment(const HorizontalMapSegment& seg)
{
  // First determine if a segment already exists in this map for the seg_dof.
  gid_type seg_dof = seg.get_dof();
  int match_loc = -999;
  for (int iseg=0; iseg<m_num_segments; iseg++) {
    auto seg_ii = m_map_segments[iseg];
    if (seg_ii.get_dof()==seg_dof) {
      // Match has been found
      match_loc = iseg;
      break;
    }
  }
  // Now add the segment or create a new combined segment
  if (match_loc == -999) {
    // Match never found, just add this one to the end
    m_map_segments.push_back(seg);
  } else {
    auto seg_match = m_map_segments[match_loc];
    // Combine this segment with the one already found
    int new_length = seg.get_length() + seg_match.get_length();
    HorizontalMapSegment new_seg(seg_dof,new_length);
    auto new_source_dofs = new_seg.get_source_dofs();
    auto new_source_idx  = new_seg.get_source_idx();
    auto new_weights     = new_seg.get_weights();
    // First copy old seg values:
    int match_len = seg_match.get_length();
    auto match_source_dofs = seg_match.get_source_dofs();
    auto match_source_idx  = seg_match.get_source_idx();
    auto match_weights     = seg_match.get_weights();
    Kokkos::parallel_for("", match_len, KOKKOS_LAMBDA (const int& ii) {
      new_source_dofs(ii) = match_source_dofs(ii);
      new_source_idx(ii)  = match_source_idx(ii);
      new_weights(ii)     = match_weights(ii);
    });
    // Next add the new segment
    int seg_len = seg.get_length();
    auto seg_source_dofs = seg.get_source_dofs();
    auto seg_source_idx  = seg.get_source_idx();
    auto seg_weights     = seg.get_weights();
    Kokkos::parallel_for("", seg_len, KOKKOS_LAMBDA (const int& ii) {
      new_source_dofs(ii+match_len) = seg_source_dofs(ii);
      new_source_idx(ii+match_len)  = seg_source_idx(ii);
      new_weights(ii+match_len)     = seg_weights(ii);
    });
    // Replace old segment with this new one
    m_map_segments[match_loc] = new_seg;
  }

  m_num_segments = m_map_segments.size();
  if (m_unique_set) {
    // We need to reset the unique set of source columns taking
    // into account the new segment.
    set_unique_source_dofs();
  }
} 
/*-----------------------------------------------------------------------------------------------*/
// This function defines the unqiue set of dofs associated with this HorizontalMap.  This is important for
// knowing just exactly which dofs in the source data this HorizontalMap intends to use.  By just focusing on
// the dofs the HorizontalMap will use we can decrease the memory footprint.
void HorizontalMap::set_unique_source_dofs()
{
  start_timer("EAMxx::HorizontalMap::set_unique_dofs");
  // Create a temporary vector of dofs which will be used for lookup
  std::vector<gid_type> unique_dofs;
  // Check all segments and add unique dofs.  Done on HOST so we can use a vector, only done once
  // per map so it's alright to not be optimized for performance.
  int min_gid = INT_MAX;
  int max_gid = INT_MIN;
  for (int iseg=0; iseg<m_num_segments; iseg++) {
    auto seg = m_map_segments[iseg];
    const auto& src_dofs = seg.get_source_dofs();
    const auto& src_dofs_h = Kokkos::create_mirror_view(src_dofs);
    Kokkos::deep_copy(src_dofs_h,src_dofs);
    for (int ii=0; ii<seg.get_length(); ii++) {
      auto idx = std::find(unique_dofs.begin(), unique_dofs.end(), src_dofs_h(ii));
      if (idx == unique_dofs.end()) {
        unique_dofs.push_back(src_dofs_h(ii));
        min_gid = std::min(min_gid,src_dofs_h(ii));
        max_gid = std::max(max_gid,src_dofs_h(ii));
      }
    }
  }
  std::sort(unique_dofs.begin(), unique_dofs.end());
  // Assign unique source dofs to internal view
  m_num_unique_dofs = unique_dofs.size();
  m_unique_dofs = view_1d<gid_type>("",m_num_unique_dofs);
  auto temp_h = Kokkos::create_mirror_view(m_unique_dofs);
  for (int ii=0; ii<m_num_unique_dofs; ii++) {
    temp_h(ii) = unique_dofs[ii];
  }
  Kokkos::deep_copy(m_unique_dofs,temp_h);
  m_unique_set = true;

  // While we are here we can construct the lookup of source_dofs to unique_source_dofs
  // to be used with segment mapping.  Also done on HOST so we can take advantage of the
  // find function for vectors.
  for (int iseg=0; iseg<m_num_segments; iseg++)
  {
    auto& seg = m_map_segments[iseg];
    // Lookup for each segment w.r.t. m_dofs_gids on this map
    const auto seg_dof = seg.get_dof();
    bool seg_dof_found = false;
    auto dofs_gids_h = Kokkos::create_mirror_view(m_dofs_gids);
    Kokkos::deep_copy(dofs_gids_h,m_dofs_gids);
    for (int ii=0; ii<m_num_dofs; ii++)
    {
      if (seg_dof == dofs_gids_h(ii))
      {
        seg.set_dof_idx(ii);
        seg_dof_found = true;
        break;
      }
    }
    EKAT_REQUIRE(seg_dof_found);
    // Lookup for source data w.r.t. unique dofs
    auto src_dofs  = seg.get_source_dofs();
    auto src_dofs_h = Kokkos::create_mirror_view(src_dofs);
    Kokkos::deep_copy(src_dofs_h,src_dofs);
    auto src_idx   = seg.get_source_idx();
    auto src_idx_h = Kokkos::create_mirror_view(src_idx);
    Kokkos::deep_copy(src_idx_h,src_idx);
    for (int ii=0; ii<seg.get_length(); ii++) {
      auto dof = src_dofs_h(ii);
      // First we find where in the set of unique columns a particular
      // source dof is.
      auto idx = std::find(unique_dofs.begin(),unique_dofs.end(),dof);
      // Then assign the src_idx, which is a map from source dof to where in the set
      // of unique dofs a source dof is.
      int idx_ii = idx - unique_dofs.begin();
      src_idx_h(ii) = idx_ii;
    }
    Kokkos::deep_copy(src_idx,src_idx_h);
    // Sync to Host
    seg.sync_to_host();
  }
  stop_timer("EAMxx::HorizontalMap::set_unique_dofs");
}
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
void HorizontalMap::check() const
{
  EKAT_REQUIRE_MSG(m_dofs_set,"Error in HorizontalMap " + m_name + " on rank " + std::to_string(m_comm.rank()) +" - Global DOFS not yet set, need to call set_dof_gids."); 
  for (const auto& seg : m_map_segments) {
    EKAT_REQUIRE_MSG(seg.check(),"Error in HorizontalMap " + m_name + " on rank " + std::to_string(m_comm.rank()) +" - problem with a remapping segment for dof = " + std::to_string(seg.get_dof()) + ".");
  }
}
/*-----------------------------------------------------------------------------------------------*/
void HorizontalMap::print() const
{
    // TODO - Make this print everything on ROOT, probably need to pass map info
    //        via MPI first.  That will make the print out information easier to
    //        parse, for debugging.
    printf("\n=============================================\n");
    printf("Printing map information for map: %s\n",m_name.c_str());
    for (int ii=0; ii<m_num_segments; ii++) {
      auto& seg = m_map_segments[ii];
      seg.print();
    }
 
    printf(" Unique dofs info\n");
    if (m_unique_set) {
      auto unique_dofs_h = Kokkos::create_mirror_view(m_unique_dofs);
      Kokkos::deep_copy(unique_dofs_h,m_unique_dofs);
      for (size_t ii=0; ii<m_unique_dofs.extent(0); ii++) {
        printf("%10zu: %10d\n",ii, unique_dofs_h(ii));
      }
    } else {
      if (m_comm.am_i_root()) {
        printf("  WARNING - Unique DOFs have not been set yet\n");
      }
    }

    printf(" dofs_gids\n");
    auto dofs_gids_h = Kokkos::create_mirror_view(m_dofs_gids);
    Kokkos::deep_copy(dofs_gids_h,m_dofs_gids);
    for (size_t ii=0; ii<dofs_gids_h.extent(0); ii++) {
      printf("%10zu: %10d\n",ii,dofs_gids_h(ii));
    } 
    printf("\n=============================================\n");
}
/*-----------------------------------------------------------------------------------------------*/
// This overload of apply remap assumes a single horizontal slice of source data being mapped onto
// a horizontal slice of remapped data.  The assumption is that there are no levels in this data.
void HorizontalMap::apply_remap(const view_1d<const Real>& source_data, const view_1d<Real>& remapped_data) {
  start_timer("EAMxx::HorizontalMap::apply_remap_1d");
  if (m_num_dofs==0) { return; } // This HorizontalMap has nothing to do for this rank.
  auto remapped_data_h = Kokkos::create_mirror_view(remapped_data);
  auto source_data_h = Kokkos::create_mirror_view(source_data);
  Kokkos::deep_copy(source_data_h,source_data);
  Kokkos::deep_copy(remapped_data_h,0.0);
  for (int iseg=0; iseg<m_num_segments; iseg++) {
    auto seg = m_map_segments[iseg];
    auto seg_dof_idx  = seg.get_dof_idx();
    auto seg_length   = seg.get_length();
    auto source_idx_h = seg.get_source_idx_on_host();
    auto weights_h    = seg.get_weights_on_host();
    for (int ii=0; ii<seg_length; ii++) {
      int idx = source_idx_h(ii);
      remapped_data_h(seg_dof_idx) += source_data_h(idx)*weights_h(ii);
    }
  }
  Kokkos::deep_copy(remapped_data,remapped_data_h);
  stop_timer("EAMxx::HorizontalMap::apply_remap_1d");
}
/*-----------------------------------------------------------------------------------------------*/
// This overload of apply remap assumes a set of horizontal slices of source data being mapped onto
// a set horizontal slices of remapped data.  The assumption is that the second dimension is number
// of levels
void HorizontalMap::apply_remap(const view_2d<const Real>& source_data, const view_2d<Real>& remapped_data) {
  start_timer("EAMxx::HorizontalMap::apply_remap_2d");
  if (m_num_dofs==0) { return; } // This HorizontalMap has nothing to do for this rank.
  int num_levs = source_data.extent(1);
  auto remapped_data_h = Kokkos::create_mirror_view(remapped_data);
  auto source_data_h = Kokkos::create_mirror_view(source_data);
  Kokkos::deep_copy(source_data_h,source_data);
  Kokkos::deep_copy(remapped_data_h,0.0);
  for (int iseg=0; iseg<m_num_segments; iseg++) {
    auto seg = m_map_segments[iseg];
    auto seg_dof_idx  = seg.get_dof_idx();
    auto seg_length   = seg.get_length();
    auto source_idx_h = seg.get_source_idx_on_host();
    auto weights_h    = seg.get_weights_on_host();
    for (int ii=0; ii<seg_length; ii++) {
      for (int kk=0; kk<num_levs; kk++) {
        int idx = source_idx_h(ii);
        remapped_data_h(seg_dof_idx,kk) += source_data_h(idx,kk)*weights_h(ii);
      }
    }
  }
  Kokkos::deep_copy(remapped_data,remapped_data_h);
  stop_timer("EAMxx::HorizontalMap::apply_remap_2d");
}
/*-----------------------------------------------------------------------------------------------*/
// This overload of apply remap assumes a set of horizontal slices of source data being mapped onto
// a set of horizontal slices of remapped data.  The assumption is that there are levels and one other
// dimension for this data.
void HorizontalMap::apply_remap(const view_3d<const Real>& source_data, const view_3d<Real>& remapped_data) {
  start_timer("EAMxx::HorizontalMap::apply_remap_3d");
  if (m_num_dofs==0) { return; } // This HorizontalMap has nothing to do for this rank.
  int num_levs = source_data.extent(2);
  int num_bands = source_data.extent(1);
  auto remapped_data_h = Kokkos::create_mirror_view(remapped_data);
  auto source_data_h = Kokkos::create_mirror_view(source_data);
  Kokkos::deep_copy(source_data_h,source_data);
  Kokkos::deep_copy(remapped_data_h,0.0);
  for (int iseg=0; iseg<m_num_segments; iseg++) {
    auto seg = m_map_segments[iseg];
    auto seg_dof_idx  = seg.get_dof_idx();
    auto seg_length   = seg.get_length();
    auto source_idx_h = seg.get_source_idx_on_host();
    auto weights_h    = seg.get_weights_on_host();
    for (int ii=0; ii<seg_length; ii++) {
      for (int nn=0; nn<num_bands; nn++) {
        for (int kk=0; kk<num_levs; kk++) {
          int idx = source_idx_h(ii);
          remapped_data_h(seg_dof_idx,nn,kk) += source_data_h(idx,nn,kk)*weights_h(ii);
        }
      }
    }
  }
  Kokkos::deep_copy(remapped_data,remapped_data_h);
  stop_timer("EAMxx::HorizontalMap::apply_remap_3d");
}
/*-----------------------------------------------------------------------------------------------*/
HorizontalMapSegment::HorizontalMapSegment(const gid_type dof_gid, const int length)
  : m_dof    (dof_gid)
  , m_length (length)
{
  m_source_dofs  = view_1d<gid_type>("",m_length);
  m_source_idx   = view_1d<int>("",m_length);
  m_weights      = view_1d<Real>("",m_length);
  m_source_idx_h = Kokkos::create_mirror_view(m_source_idx);
  m_weights_h    = Kokkos::create_mirror_view(m_weights);
}
/*-----------------------------------------------------------------------------------------------*/
HorizontalMapSegment::HorizontalMapSegment(const gid_type dof_gid, const int length,  const view_1d<const gid_type>& source_dofs, const view_1d<const Real>& weights)
  : m_dof         (dof_gid)
  , m_length      (length)
{
  m_source_dofs = view_1d<gid_type>("",m_length);
  m_weights     = view_1d<Real>("",m_length);
  m_source_idx  = view_1d<int>("",m_length);
  Kokkos::deep_copy(m_source_dofs,source_dofs);
  Kokkos::deep_copy(m_weights,weights);
  m_source_idx_h = Kokkos::create_mirror_view(m_source_idx);
  m_weights_h    = Kokkos::create_mirror_view(m_weights);
}
/*-----------------------------------------------------------------------------------------------*/
void HorizontalMapSegment::sync_to_host()
{
  Kokkos::deep_copy(m_source_idx_h,m_source_idx);
  Kokkos::deep_copy(m_weights_h,m_weights);
}
/*-----------------------------------------------------------------------------------------------*/
bool HorizontalMapSegment::check() const
{
  // Basic check for bounds for arrays
  EKAT_REQUIRE_MSG(m_source_dofs.extent_int(0)==m_length,"Error remap segment for DOF: " + std::to_string(m_dof) + ", source_dofs view not the correct length");
  EKAT_REQUIRE_MSG(m_weights.extent_int(0)==m_length,"Error remap segment for DOF: " + std::to_string(m_dof) + ", weightss view not the correct length");
  EKAT_REQUIRE_MSG(m_source_idx.extent_int(0)==m_length,"Error remap segment for DOF: " + std::to_string(m_dof) + ", source_idx view not the correct length");
  // Check that the segment weight adds up to 1.0 
  Real wgt = 0.0;
  const auto weights = m_weights;
  Kokkos::parallel_reduce("", m_length, KOKKOS_LAMBDA (const int& ii, Real& lsum) {
    lsum += weights(ii);
  },wgt);
  Real tol = std::numeric_limits<Real>::epsilon() * 100.0;
  if (std::abs(wgt-1.0)>=tol) {
    printf("ERROR: HorizontalMap: checking remap segment for DOF = %d, total weight = %e.\n",m_dof,wgt);
    return false;
  }

  // If we made it this far, things all passed
  return true;
}
/*-----------------------------------------------------------------------------------------------*/
void HorizontalMapSegment::print() const
{
    printf("\n--------------------\n");
    printf("Printing information for segment with DOF = %d, DOF_idx for local decomp = %d\n",m_dof,m_dof_idx);
    printf("  length = %d\n",m_length);

    auto source_dofs_h = Kokkos::create_mirror_view(m_source_dofs);
    auto source_idx_h  = Kokkos::create_mirror_view(m_source_idx);
    auto weights_h     = Kokkos::create_mirror_view(m_weights);
    Kokkos::deep_copy(source_dofs_h,m_source_dofs);
    Kokkos::deep_copy(source_idx_h ,m_source_idx );
    Kokkos::deep_copy(weights_h    ,m_weights    );
    Real total_wgt = 0.0;
    printf("%10s: %10s, %10s, %s\n","ii","source dof","source idx","weight");
    for (int ii=0; ii<m_length; ii++) {
      printf("%10d: %10d, %10d, %e\n",ii, source_dofs_h(ii), source_idx_h(ii), weights_h(ii));
      total_wgt += weights_h(ii);
    }
    printf("%36s%e\n","",total_wgt);
    printf("\n--------------------\n");
}
/*-----------------------------------------------------------------------------------------------*/

} // namespace scream
