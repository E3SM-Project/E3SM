#include "share/grid/remap/horizontal_remap_utility.hpp"

namespace scream {

/*-----------------------------------------------------------------------------------------------*/
GSMap::GSMap(const ekat::Comm& comm)
  : m_comm (comm)
{
  m_dofs_set = false;
}
/*-----------------------------------------------------------------------------------------------*/
GSMap::GSMap(const ekat::Comm& comm, const std::string& map_name)
  : m_comm (comm)
  , m_name (map_name)
{
  m_dofs_set = false;
}
/*-----------------------------------------------------------------------------------------------*/
GSMap::GSMap(const ekat::Comm& comm, const std::string& map_name, const view_1d<gid_type>& dofs_gids, const gid_type min_dof)
  : m_comm (comm)
  , m_name (map_name)
{
  set_dof_gids(dofs_gids,min_dof);
}
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
// This function is used to set the internal set of degrees of freedom (dof) this map is responsible for.
// We use the global dofs, offset by the minimum global dof to make everything zero based.  Note, when
// gathering remap parameters from a file, depending on the algorithm that made the file the dof
// indices may be 1-based or 0-based.  By offsetting everything to 0-based we avoid potential bugs.
void GSMap::set_dof_gids(const view_1d<gid_type>& dofs_gids, const gid_type min_dof)
{
  EKAT_REQUIRE(dofs_gids.size()>0);
  m_dofs_gids = view_1d<gid_type>(dofs_gids);
  m_num_dofs = m_dofs_gids.extent(0);
  Kokkos::parallel_for("", m_num_dofs, KOKKOS_LAMBDA (const int& ii) {
    m_dofs_gids(ii) = dofs_gids(ii)-min_dof;
  });
  m_dofs_set = true;
}
/*-----------------------------------------------------------------------------------------------*/
// This function adds a remap segment to a GSMap, note, we want each segment to represent a full
// remapping.  This function also checks if a segment already exists for the degree of freedom
// in question.  If it does then instead of add the segment to the end, this function finds that
// segment and combines them into a new comprehensive segment.
void GSMap::add_remap_segment(const GSSegment& seg)
{
  // First determine if a segment already exists in this map for the seg_dof.
  gid_type seg_dof = seg.get_dof();
  Int match_loc = -999;
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
    Int new_length = seg.get_length() + seg_match.get_length();
    GSSegment new_seg(seg_dof,new_length);
    auto new_source_dofs = new_seg.get_source_dofs();
    auto new_source_idx  = new_seg.get_source_idx();
    auto new_weights     = new_seg.get_weights();
    // First copy old seg values:
    Int match_len = seg_match.get_length();
    auto match_source_dofs = seg_match.get_source_dofs();
    auto match_source_idx  = seg_match.get_source_idx();
    auto match_weights     = seg_match.get_weights();
    Kokkos::parallel_for("", match_len, KOKKOS_LAMBDA (const int& ii) {
      new_source_dofs(ii) = match_source_dofs(ii);
      new_source_idx(ii)  = match_source_idx(ii);
      new_weights(ii)     = match_weights(ii);
    });
    // Next add the new segment
    Int seg_len = seg.get_length();
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
// This function defines the unqiue set of dofs associated with this GSMap.  This is important for
// knowing just exactly which dofs in the source data this GSMap intends to use.  By just focusing on
// the dofs the GSMap will use we can decrease the memory footprint.
void GSMap::set_unique_source_dofs()
{
  // Create a temporary vector of dofs which will be used for lookup
  std::vector<gid_type> unique_dofs;
  // Check all segments and add unique dofs.  Done on HOST so we can use a vector, only done once
  // per map so it's alright to not be optimized for performance.
  Int min_gid = INT_MAX;
  Int max_gid = INT_MIN;
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
  m_unique_dofs = view_1d<gid_type>("",unique_dofs.size());
  auto temp_h = Kokkos::create_mirror_view(m_unique_dofs);
  for (int ii=0; ii<unique_dofs.size(); ii++) {
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
  }
}
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
void GSMap::check() const
{
  for (int iseg=0; iseg<m_map_segments.size(); iseg++) {
    auto seg = m_map_segments[iseg];
    bool seg_check = seg.check();
    EKAT_REQUIRE_MSG(seg_check,"Error in GSMap " + m_name + " - problem with a remapping segment for dof = " + std::to_string(seg.get_dof()) + ".");
  }
}
/*-----------------------------------------------------------------------------------------------*/
void GSMap::print_map() const
{
    printf("\n=============================================\n");
    printf("Printing map information for map: %s\n",m_name.c_str());
    for (int ii=0; ii<m_num_segments; ii++) {
      auto& seg = m_map_segments[ii];
      seg.print_seg();
    }
 
    printf(" Unique dofs info\n");
    if (m_unique_set) {
      auto unique_dofs_h = Kokkos::create_mirror_view(m_unique_dofs);
      Kokkos::deep_copy(unique_dofs_h,m_unique_dofs);
      for (int ii=0; ii<m_unique_dofs.extent(0); ii++) {
        printf("%10d: %10d\n",ii, unique_dofs_h(ii));
      }
    } else {
      printf("  WARNING - Unique DOFs have not been set yet\n");
    }

    printf(" dofs_gids\n");
    auto dofs_gids_h = Kokkos::create_mirror_view(m_dofs_gids);
    Kokkos::deep_copy(dofs_gids_h,m_dofs_gids);
    for (int ii=0; ii<dofs_gids_h.extent(0); ii++) {
      printf("%10d: %10d\n",ii,dofs_gids_h(ii));
    } 
    printf("\n=============================================\n");
}
/*-----------------------------------------------------------------------------------------------*/
void GSMap::apply_remap(const view_1d<Real>& source_data, view_1d<Real>& remapped_data) {
  if (m_num_dofs==0) { return; } // This GSMap has nothing to do for this rank.
  auto remapped_data_h = Kokkos::create_mirror_view(remapped_data);
  for (int iseg=0; iseg<m_num_segments; iseg++) {
    auto seg = m_map_segments[iseg];
    auto seg_dof_idx = seg.get_dof_idx();
    Real remap_val = seg.apply_segment(source_data);
    remapped_data_h(seg_dof_idx) = remap_val;
  }
  Kokkos::deep_copy(remapped_data,remapped_data_h);
}
//ASDtemplate <typename ScalarT>
//ASDvoid GSMap::apply_remap(const view_2d<ScalarT>& source_data, const view_2d<ScalarT>& remapped_data)
//ASD{
//ASD  // do nothing
//ASD}
//ASD
//ASDtemplate <typename ScalarT>
//ASDvoid GSMap::apply_remap(const view_3d<ScalarT>& source_data, const view_3d<ScalarT>& remapped_data)
//ASD{
//ASD  // do nothing
//ASD}
/*-----------------------------------------------------------------------------------------------*/
/*-----------------------------------------------------------------------------------------------*/
GSSegment::GSSegment(const gid_type dof_gid, const Int length)
  : m_dof    (dof_gid)
  , m_length (length)
{
  m_source_dofs = view_1d<gid_type>("",m_length);
  m_source_idx  = view_1d<Int>("",m_length);
  m_weights     = view_1d<Real>("",m_length);
}
/*-----------------------------------------------------------------------------------------------*/
Real GSSegment::apply_segment(const view_1d<Real>& source_data)
{
  Real ret;
  Kokkos::parallel_reduce("", m_length, KOKKOS_LAMBDA (const int& ii, Real& loc) {
    Int idx = m_source_idx(ii);
    loc += source_data(idx)*m_weights(ii);
  },ret);
  return ret;
}
/*-----------------------------------------------------------------------------------------------*/
GSSegment::GSSegment(const gid_type dof_gid, const Int length,  const view_1d<gid_type>& source_dofs, const view_1d<Real>& weights)
  : m_dof         (dof_gid)
  , m_length      (length)
  , m_source_dofs (source_dofs)
  , m_weights     (weights)
{
  m_source_idx  = view_1d<Int>("",m_length);
}
/*-----------------------------------------------------------------------------------------------*/
bool GSSegment::check() const
{
  // Basic check for bounds for arrays
  EKAT_REQUIRE_MSG(m_source_dofs.extent(0)==m_length,"Error remap segment for DOF: " + std::to_string(m_dof) + ", source_dofs view not the correct length");
  EKAT_REQUIRE_MSG(m_weights.extent(0)==m_length,"Error remap segment for DOF: " + std::to_string(m_dof) + ", weightss view not the correct length");
  EKAT_REQUIRE_MSG(m_source_idx.extent(0)==m_length,"Error remap segment for DOF: " + std::to_string(m_dof) + ", source_idx view not the correct length");
  // Check that the segment weight adds up to 1.0 
  Real wgt = 0.0;
  Kokkos::parallel_reduce("", m_length, KOKKOS_LAMBDA (const int& ii, Real& lsum) {
    lsum += m_weights(ii);
  },wgt);
  Real tol = std::numeric_limits<Real>::epsilon() * 100.0;
  if (std::abs(wgt-1.0)>=tol) {
    printf("ERROR: GSMAP: checking remap segment for DOF = %d, total weight = %e.\n",m_dof,wgt);
    return false;
  }

  // If we made it this far, things all passed
  return true;
}
/*-----------------------------------------------------------------------------------------------*/
void GSSegment::print_seg() const
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
