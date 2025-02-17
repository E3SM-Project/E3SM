#ifndef SCREAM_ABSTRACT_GRID_HPP
#define SCREAM_ABSTRACT_GRID_HPP

#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "share/grid/grid_utils.hpp"
#include "share/field/field_layout.hpp"
#include "share/field/field.hpp"

#include "ekat/mpi//ekat_comm.hpp"

#include <map>
#include <list>
#include <memory>
#include <mutex>

namespace scream
{

/*
 * An interface base class for Grid objects
 *
 * An abstract grid represents a 2d grid, vertically extruded. It needs to be able
 * to provide the following information:
 *   - number of local degrees of freedom (dofs) along the horizontal direction
 *   - number of vertical levels
 *   - gids for all local 2d dofs
 *   - type and name of the grid (see grid_utils.hpp for types)
 *   - the layout of a 2d/3d field on this grid (see field_layout.hpp)
 *   - the mapping of dof lid to a set of indices (i_1,...,i_N), with N being the
 *     rank of the scalar field layout, and i_k less than the i-th dimension
 *     in the scalar field layout. This is how a dof is 'naturally' indexed on this grid.
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
  // TODO: use int64_t? int? template class on gid_type?
  // So far, 32 bits seem enough for 2d dofs numbering
  using gid_type = int;
  using gid_view_h = Field::view_host_t<const gid_type*>;

  // Constructor(s) & Destructor
  AbstractGrid (const std::string& name,
                const GridType type,
                const int num_local_dofs,
                const int num_vertical_lev,
                const ekat::Comm& comm);

  AbstractGrid (const std::string& name,
                const GridType type,
                const int num_local_dofs,
                const int num_global_dofs,
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
  FieldLayout get_vertical_layout (const bool midpoints) const;
  FieldLayout get_vertical_layout (const bool midpoints,
                                   const int vector_dim,
                                   const std::string& vec_dim_name = e2str(FieldTag::Component)) const;
  virtual FieldLayout get_2d_scalar_layout () const = 0;
  virtual FieldLayout get_2d_vector_layout (const int vector_dim, const std::string& vec_dim_name) const = 0;
  virtual FieldLayout get_2d_tensor_layout (const std::vector<int>& cmp_dims,
                                            const std::vector<std::string>& cmp_dims_names) const = 0;
  virtual FieldLayout get_3d_scalar_layout (const bool midpoints) const = 0;
  virtual FieldLayout get_3d_vector_layout (const bool midpoints, const int vector_dim,
                                            const std::string& vec_dim_name) const = 0;
  virtual FieldLayout get_3d_tensor_layout (const bool midpoints,
                                            const std::vector<int>& cmp_dims,
                                            const std::vector<std::string>& cmp_dims_names) const = 0;

  // Some shortcut versions of the above ones, where the name of the vector/tensor
  // components are all equal to e2str(CMP)
  FieldLayout get_2d_vector_layout (const int vector_dim) const;
  FieldLayout get_2d_tensor_layout (const std::vector<int>& cmp_dims) const;

  FieldLayout get_3d_vector_layout (const bool midpoints, const int vector_dim) const;
  FieldLayout get_3d_tensor_layout (const bool midpoints, const std::vector<int>& cmp_dims) const;

  int get_num_vertical_levels () const { return m_num_vert_levs; }

  // Whether this grid contains unique dof GIDs
  bool is_unique () const;

  // Check if the input layout is compatible with this grid
  bool is_valid_layout (const FieldLayout& layout) const;

  // When running with multiple ranks, fields are partitioned across ranks along this FieldTag
  virtual FieldTag get_partitioned_dim_tag () const = 0;

  // The extent of the partitioned dimension, locally and globally
  virtual int get_partitioned_dim_local_size () const = 0;
  virtual int get_partitioned_dim_global_size () const = 0;

  // The number of dofs on this MPI rank
  int get_num_local_dofs  () const { return m_num_local_dofs;  }
  gid_type get_num_global_dofs () const { return m_num_global_dofs; }
  gid_type get_global_min_dof_gid () const;
  gid_type get_global_max_dof_gid () const;
  gid_type get_global_min_partitioned_dim_gid () const;
  gid_type get_global_max_partitioned_dim_gid () const;

  // Get a Field storing 1d data (the dof gids)
  Field get_dofs_gids () const;
  Field get_dofs_gids ();

  // Get Field storing the gids that this process owns along the partitioned dim
  // NOTE: for some grids, this is the same as get_dofs_gids. The SEGrid is a counterexample:
  //       the dofs are the GLL dofs, but the partitioned dim is the element dimension
  Field get_partitioned_dim_gids ();
  Field get_partitioned_dim_gids () const;

  // Get a Field storing 2d data, where (i,j) entry contains the j-th coordinate of
  // the i-th dof in the native dof layout. Const verison returns a read-only field
  Field get_lid_to_idx_map () const;
  Field get_lid_to_idx_map ();

  // Get geometry-related fields
  Field get_geometry_data (const std::string& name) const;

  // Create geometry data, throws if already existing. Returns writable field
  Field create_geometry_data (const FieldIdentifier& fid, const int pack_size = 1);
  Field create_geometry_data (const std::string& name, const FieldLayout& layout,
                              const ekat::units::Units& units = ekat::units::Units::invalid(),
                              const DataType data_type = DataType::RealType,
                              const int pack_size = 1) {
    return create_geometry_data(FieldIdentifier(name,layout,units,this->name(),data_type),pack_size);
  }

  // Sets pre-existing field as geometry data.
  // NOTE: setter is const, since we do allow adding new data even if grid is const
  //       E.g., this allows atm procs to define coordinate vars for dimensions
  //       peculiar to that process
  void set_geometry_data (const Field& f) const;
  void delete_geometry_data (const std::string& name);

  bool has_geometry_data (const std::string& name) const {
    return m_geo_fields.find(name)!=m_geo_fields.end();
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
  std::vector<gid_type> get_unique_gids () const;

  // For each entry in the input list of GIDs, retrieve the process id that owns it
  std::vector<int> get_owners (const gid_view_h& gids) const;
  std::vector<int> get_owners (const std::vector<gid_type>& gids) const {
    gid_view_h gids_v(gids.data(),gids.size());
    return get_owners(gids_v);
  }

  void get_remote_pids_and_lids (const gid_view_h& gids,
                                 std::vector<int>& pids,
                                 std::vector<int>& lids) const;
  void get_remote_pids_and_lids (const std::vector<gid_type>& gids,
                                 std::vector<int>& pids,
                                 std::vector<int>& lids) const {
    gid_view_h gids_v(gids.data(),gids.size());
    get_remote_pids_and_lids(gids_v,pids,lids);
  }

  // Derived classes can override these methods to verify that the
  // dofs have been set to something that satisfies any requirement of the grid type.
  virtual bool check_valid_dofs()        const { return true; }
  virtual bool check_valid_lid_to_idx () const { return true; }

  void reset_field_tag_name (const FieldTag t, const std::string& s) { m_special_tag_names[t] = s; }
  bool has_special_tag_name (const FieldTag t) const { return m_special_tag_names.count(t)==1; }
  std::string get_special_tag_name (const FieldTag t) const { return m_special_tag_names.at(t); }

  // This member is used mostly by IO: if a field exists on multiple grids
  // with the same name, IO can use this as a suffix to diambiguate the fields in
  // the IO file, by appending each grid's suffix to the fields names.
  // NOTE: we'd need setter/getter for this, so we might as well make it public
  std::string m_short_name = "";

  int get_unique_grid_id () const { return m_unique_grid_id; }

  const std::map<gid_type,int>& get_gid2lid_map () const;

protected:

  void copy_data (const AbstractGrid& src, const bool shallow = true);

  // Note: this method must be called from the derived classes,
  //       since it calls get_2d_scalar_layout.
  void create_dof_fields (const int scalar2d_layout_rank);

  // The grid name and type
  GridType     m_type;
  std::string  m_name;

  int m_unique_grid_id;

  std::vector<std::string> m_aliases;

  std::map<FieldTag, std::string> m_special_tag_names;

  // Counters
  int m_num_local_dofs;
  int m_num_global_dofs;
  int m_num_vert_levs;

  // The global ID of each dof
  Field     m_dofs_gids;

  // The global ID of the owned entries of the partitioned dimension (if any)
  Field     m_partitioned_dim_gids;

  // The max/min dof GID across all ranks. Mutable, to allow for lazy calculation
  mutable gid_type  m_global_min_dof_gid =  std::numeric_limits<gid_type>::max();
  mutable gid_type  m_global_max_dof_gid = -std::numeric_limits<gid_type>::max();
  // Same as above, but for partitioned dim gids
  mutable gid_type  m_global_min_partitioned_dim_gid =  std::numeric_limits<gid_type>::max();
  mutable gid_type  m_global_max_partitioned_dim_gid = -std::numeric_limits<gid_type>::max();

  // The fcn is_unique is expensive, so we lazy init this at the first call.
  mutable bool m_is_unique;
  mutable bool m_is_unique_computed = false;

  // The map lid->idx
  Field     m_lid_to_idx;

  mutable std::map<std::string,Field>  m_geo_fields;

  // Mutable, for lazy calculation
  mutable std::map<gid_type,int> m_gid2lid;

  // For thread safety in modifying mutable items (just in case someone ever runs this code in threaded regions)
  mutable std::mutex m_mutex;

  // The MPI comm containing the ranks across which the global mesh is partitioned
  ekat::Comm            m_comm;
};

} // namespace scream

#endif // SCREAM_ABSTRACT_GRID_HPP
