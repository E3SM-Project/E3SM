#ifndef EAMXX_HORIZONTAL_REMAP_UTILITY_HPP
#define EAMXX_HORIZONTAL_REMAP_UTILITY_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/io/scorpio_input.hpp"
#include "share/scream_types.hpp"
#include "ekat/kokkos/ekat_subview_utils.hpp"

#include <numeric>

namespace scream {

/*===============================================================================================*/
/*                                        HorizontalMapSegment
 * Lightweight structure to represent a single remapping of source data to a single target column.
 *     Y_target = sum_(n=1)^N ( w_n * Y_source_n )
 * See description of HorizontalMap structure for more details on the mapping.
 *
 * This structure is used to organize the overall horizontal remapping stored in the HorizontalMap structure.
 * Each segment represents a single target column in the HorizontalMap.
 * --------------------------------------
 *  A.S. Donahue (LLNL): 2022-09-07
 *===============================================================================================*/
struct HorizontalMapSegment {

  using gid_type = AbstractGrid::gid_type;
  using KT = KokkosTypes<DefaultDevice>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  
  template <typename S>
  using view_1d_host = typename KT::template view_1d<S>::HostMirror;
  
  // Constructors/Destructor
  HorizontalMapSegment() {};
  HorizontalMapSegment(const gid_type dof_gid, const int length);
  HorizontalMapSegment(const gid_type dof_gid, const int length, const view_1d<const gid_type>& source_dofs, const view_1d<const Real>& weights);

  // Helper Functions
  bool check() const;      // Check if this segment is valid
  void print() const;  // Useful for debugging, print the segment mapping info

  // Setter Functions
  void set_dof_idx(const int idx) { m_dof_idx = idx; }
  void sync_to_host();

  // Getter Functions
  gid_type get_dof()     const { return m_dof; }
  int      get_dof_idx() const { return m_dof_idx; }
  int      get_length()  const { return m_length; }
  view_1d<gid_type>  get_source_dofs()        const { return m_source_dofs; }
  view_1d<int>       get_source_idx()         const { return m_source_idx; }
  view_1d<Real>      get_weights()            const { return m_weights; }
  view_1d_host<int>  get_source_idx_on_host() const { return m_source_idx_h; }
  view_1d_host<Real> get_weights_on_host()    const { return m_weights_h; }
 
  // TODO: Not sure why, but this can't be set as private, otherwise `set_dof_idx` doesn't work. 
  int      m_dof_idx = -999; // The degree of freedom w.r.t. to the local index for this map
private:

  // Remap views
  view_1d<gid_type>  m_source_dofs;
  view_1d<int>       m_source_idx;
  view_1d<Real>      m_weights;
  view_1d_host<int>  m_source_idx_h;
  view_1d_host<Real> m_weights_h;

  // Segment ID
  gid_type m_dof;     // The global degree of freedom this segment maps to
  int      m_length;
}; // HorizontalMapSegment

/*===============================================================================================*/
/*                                        HorizontalMap
 * Structure which can be used to setup and control a horizontal remapping.  This structure
 * follows the basic premise that there are a set of source columns >=1 that will map to a
 * single target column with a specific set of weights.  Mapping follows the expression:
 *   Y_target = sum_(n=1)^N ( w_n * Y_source_n )
 * where,
 *   Y_target:   Is the remapped value on the target column.
 *   w_n:        Is the n'th weight
 *   Y_source_n: Is the n'th column in the source data
 *   N:          Is the total number of source columns mapping to the target column (N>=1)
 *
 * This structure follows the format used by the component coupler.
 * --------------------------------------
 *  A.S. Donahue (LLNL): 2022-09-07
 *===============================================================================================*/

class HorizontalMap {
  // Note: The name used for mapping in the component coupler is HorizontalMap.  We could adopt a different
  // name if desired.
  using gid_type = AbstractGrid::gid_type;
  using KT = KokkosTypes<DefaultDevice>;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;

  
public:
  // Constructors/Destructor
  ~HorizontalMap() = default;
  HorizontalMap() {};
  explicit HorizontalMap(const ekat::Comm& comm);
  HorizontalMap(const ekat::Comm& comm, const std::string& map_name);
  HorizontalMap(const ekat::Comm& comm, const std::string& map_name, const view_1d<const gid_type>& dofs_gids, const gid_type min_dof);
 
  // Main remap functions
  void apply_remap(const view_1d<const Real>& source_data, const view_1d<Real>& remapped_data);
  void apply_remap(const view_2d<const Real>& source_data, const view_2d<Real>& remapped_data);
  void apply_remap(const view_3d<const Real>& source_data, const view_3d<Real>& remapped_data);
 
  // Helper functions
  void check() const;      // A check to make sure the map is valid
  void print() const;  // Useful for debugging

  // Builder functions - used to build the HorizontalMap
  void set_dof_gids(const view_1d<const gid_type>& dofs_gids, const gid_type min_dof);
  void set_unique_source_dofs();
  void add_remap_segment(const HorizontalMapSegment& seg);
  void set_remap_segments_from_file(const std::string& remap_filename);

  // Getter functions
  view_1d<gid_type>      get_unique_source_dofs() const { return m_unique_dofs; }
  int                    get_num_unique_dofs() const { return m_num_unique_dofs; }
  int                    get_num_of_dofs() const { return m_num_dofs; }
  std::vector<HorizontalMapSegment> get_map_segments() const { return m_map_segments; }
  int                    get_num_of_segments() const { return m_num_segments; }
  
private:

  // Global degrees of freedom information on target grid
  view_1d<gid_type> m_dofs_gids;
  int               m_num_dofs = 0;
  // Global degrees of freedom information on source grid
  view_1d<gid_type> m_unique_dofs;
  int               m_num_unique_dofs;
  bool              m_unique_set = false;
  // HorizontalMap data
  std::string            m_name = "";
  ekat::Comm             m_comm;
  bool                   m_dofs_set = false;
  std::vector<HorizontalMapSegment> m_map_segments;
  int                    m_num_segments = 0;

}; // struct HorizontalMap

/*===============================================================================================*/

} //namespace scream

#endif // EAMXX_HORIZONTAL_REMAP_UTILITY_HPP
