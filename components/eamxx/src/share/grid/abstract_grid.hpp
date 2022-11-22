#ifndef SCREAM_ABSTRACT_GRID_HPP
#define SCREAM_ABSTRACT_GRID_HPP

#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "share/grid/grid_utils.hpp"
#include "share/field/field_layout.hpp"
#include "share/scream_types.hpp"

#include "ekat/mpi//ekat_comm.hpp"

#include <map>
#include <list>
#include <memory>

namespace scream
{

/*
 * An interface base class for Grid objects
 *
 * An abstract grid represents a 2d grid, vertically extruded. It needs to be able
 * to provide the following information:
 *   - number of local degrees of freedom (dofs) along the horizontal direction
 *   - number of vertical levels
 *   - gids for all local (2d) dofs
 *   - type and name of the grid (see grid_utils.hpp for types)
 *   - the layout of a (2d) dof on this grid (see field_layout.hpp)
 *   - the gid of a (2d) dof given indices (i_1,...,i_N), with N being the
 *     rank of the scalar field layout, and i_k less than the i-th dimension
 *     in the scalar field layout
 *
 * The methods get_Xd_Y_layout, with X=2,3, and Y=scalar,vector, will return the
 * FieldLayout of a 2d/3d scalar/vector field on this grid.
 * Having the grid exposing these layout allows downstream classes to avoid
 * assumptions or switches on the particular grid, and simply rely on common
 * grid interfaces to retrieve the correct layout, based on known field properties,
 * like spatial dimension or "physical" rank.
 */

class AbstractGrid : public ekat::enable_shared_from_this<AbstractGrid>
{
public:
  using gid_type            = int;           // TODO: use int64_t? int? template class on gid_type?
  using device_type         = DefaultDevice; // TODO: template class on device type
  using kokkos_types        = KokkosTypes<device_type>;
  using kokkos_types_host   = KokkosTypes<HostDevice>;

  template<typename T>
  using view_1d = kokkos_types::view_1d<T>;
  template<typename T>
  using hview_1d = kokkos_types_host::view_1d<T>;
  template<typename T>
  using view_2d = kokkos_types::view_2d<T>;

  using geo_view_type       = view_1d<Real>;
  using geo_view_h_type     = hview_1d<Real>;
  using geo_view_map_type   = std::map<std::string,geo_view_type>;
  using geo_view_h_map_type = std::map<std::string,geo_view_h_type>;

  // The list of all dofs' gids
  using dofs_list_type   = view_1d<gid_type>;
  using dofs_list_h_type = hview_1d<gid_type>;

  // Row i of this 2d view gives the indices of the ith local dof
  // in the native 2d layout
  using lid_to_idx_map_type = view_2d<int>;

  // Constructor(s) & Destructor
  AbstractGrid (const std::string& name,
                const GridType type,
                const int num_local_dofs,
                const int num_vertical_lev,
                const ekat::Comm& comm);

  virtual ~AbstractGrid () = default;

  // Grid description utilities
  GridType type () const { return m_type; }
  const std::string& name () const { return m_name; }
  const ekat::Comm& get_comm () const { return m_comm; }
  const std::vector<std::string>& aliases () const { return m_aliases; }
  void add_alias (const std::string& alias);

  // Native layout of a dof. This is the natural way to index a dof in the grid.
  // E.g., for a scalar 2d field on a SE grid, this will be (nelem,np,np),
  //       for a vector 3d field on a Point grid it will be (ncols,vector_dim,nlevs)
  virtual FieldLayout get_2d_scalar_layout () const = 0;
  virtual FieldLayout get_2d_vector_layout (const FieldTag vector_tag, const int vector_dim) const = 0;
  virtual FieldLayout get_3d_scalar_layout (const bool midpoints) const = 0;
  virtual FieldLayout get_3d_vector_layout (const bool midpoints, const FieldTag vector_tag, const int vector_dim) const = 0;

  int get_num_vertical_levels () const { return m_num_vert_levs; }

  // Whether this grid contains unique dof GIDs
  bool is_unique () const;

  // When running with multiple ranks, fields are partitioned across ranks along this FieldTag
  virtual FieldTag get_partitioned_dim_tag () const = 0;

  // The extent of the partitioned dimension, locally and globally
  virtual int get_partitioned_dim_local_size () const = 0;
  virtual int get_partitioned_dim_global_size () const = 0;

  // The number of dofs on this MPI rank
  int get_num_local_dofs  () const { return m_num_local_dofs;  }
  gid_type get_num_global_dofs () const { return m_num_global_dofs; }
  gid_type get_global_min_dof_gid () const { return m_global_min_dof_gid; }
  gid_type get_global_max_dof_gid () const { return m_global_max_dof_gid; }

  // Set the dofs list
  // NOTE: this method calls valid_dofs_list, which may contain collective
  //       operations over the stored communicator.
  void set_dofs (const dofs_list_type& dofs);
  void set_global_min_dof_gid ();
  void set_global_max_dof_gid ();

  // Get a 1d view containing the dof gids
  const dofs_list_type& get_dofs_gids () const;
  const dofs_list_h_type& get_dofs_gids_host () const;

  // Set the the map dof_lid->dof_indices, where the indices are the ones used
  // to access the dof in the layout returned by get_2d_scalar_layout().
  // NOTE: this method calls valid_lid_to_idx_map, which may contain collective
  //       operations over the stored communicator.
  void set_lid_to_idx_map (const lid_to_idx_map_type& lid_to_idx);

  // Get a 2d view, where (i,j) entry contains the j-th coordinate of
  // the i-th dof in the native dof layout.
  const lid_to_idx_map_type& get_lid_to_idx_map () const;

  // Set/get geometric views.
  void set_geometry_data (const std::string& name, const geo_view_type& data);
  const geo_view_type& get_geometry_data (const std::string& name) const;
  const geo_view_h_type& get_geometry_data_host (const std::string& name) const;

  bool has_geometry_data (const std::string& name) const {
    return m_geo_views.find(name)!=m_geo_views.end();
  }

  // Get list of currently stored geometry data views
  std::list<std::string> get_geometry_data_names () const;

  // Creates a copy of this grid. If shallow=true, the copy shares views with
  // *this, otherwise each stored array is deep copied
  virtual std::shared_ptr<AbstractGrid> clone (const std::string& clone_name,
                                               const bool shallow) const = 0;

  // Allows to change the number of vertical levels associated with this grid.
  void reset_num_vertical_lev (const int num_vertical_lev);

  // Get a list of GIDs that are unique across all ranks in the grid comm. That is,
  // if a dof is present on 2+ ranks, it will (globally) appear just once in the
  // view returned by this method.
  dofs_list_type get_unique_gids () const;

  // For each entry in the input list of GIDs, retrieve the process id that owns it
  hview_1d<int> get_owners (const hview_1d<const gid_type>& gids) const;

  // Handy version of the above method, to allow passing a std::vector
  hview_1d<int> get_owners (const std::vector<gid_type>& gids) const {
    hview_1d<const gid_type> gids_v(gids.data(),gids.size());
    return get_owners(gids_v);
  }
protected:

  // Derived classes can override these methods, which are called inside the
  // set_dofs, set_lid_to_idx_map, and set_geometry_data methods respectively, to verify that the
  // views have been set to something that satisfies any requirement of the grid type.
  // This class already checks the extents of the view, but derived classes can add
  // some extra consistency check.
  virtual bool valid_dofs_list (const dofs_list_type& /*dofs_gids*/)      const { return true; }
  virtual bool valid_lid_to_idx_map (const lid_to_idx_map_type& /*lid_to_idx*/) const { return true; }

  void copy_views (const AbstractGrid& src, const bool shallow = true);

private:

  // The grid name and type
  GridType     m_type;
  std::string  m_name;

  std::vector<std::string> m_aliases;

  // Counters
  int m_num_local_dofs;
  int m_num_global_dofs;
  int m_num_vert_levs;

  // Whether the dofs have been set
  bool m_dofs_set = false;
  // Whether the lid->idx map has been set
  bool m_lid_to_idx_set = false;

  // The global ID of each dof
  dofs_list_type        m_dofs_gids;
  dofs_list_h_type      m_dofs_gids_host;
  gid_type              m_global_min_dof_gid;
  gid_type              m_global_max_dof_gid;

  // The map lid->idx
  lid_to_idx_map_type   m_lid_to_idx;

  geo_view_map_type     m_geo_views;
  geo_view_h_map_type   m_geo_views_host;

  // The MPI comm containing the ranks across which the global mesh is partitioned
  ekat::Comm            m_comm;
};

} // namespace scream

#endif // SCREAM_ABSTRACT_GRID_HPP
