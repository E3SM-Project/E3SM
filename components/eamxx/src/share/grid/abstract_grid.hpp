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
  gid_type get_global_min_dof_gid () const;
  gid_type get_global_max_dof_gid () const;

  // Get a Field storing 1d data (the dof gids)
  Field get_dofs_gids () const;
  Field get_dofs_gids ();

  // Get a Field storing 2d data, where (i,j) entry contains the j-th coordinate of
  // the i-th dof in the native dof layout. Const verison returns a read-only field
  Field get_lid_to_idx_map () const;
  Field get_lid_to_idx_map ();

  // Get geometry-related fields
  Field get_geometry_data (const std::string& name) const;

  // Create geometry data, throws if already existing. Returns writable field
  Field create_geometry_data (const FieldIdentifier& fid);
  Field create_geometry_data (const std::string& name, const FieldLayout& layout,
                              const ekat::units::Units& units = ekat::units::Units::invalid(),
                              const DataType data_type = DataType::RealType) {
    return create_geometry_data(FieldIdentifier(name,layout,units,this->name(),data_type));
  }

  // Sets pre-existing field as geometry data.
  void set_geometry_data (const Field& f);

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

  // Derived classes can override these methods to verify that the
  // dofs have been set to something that satisfies any requirement of the grid type.
  virtual bool check_valid_dofs()        const { return true; }
  virtual bool check_valid_lid_to_idx () const { return true; }

  // This member is used mostly by IO: if a field exists on multiple grids
  // with the same name, IO can use this as a suffix to diambiguate the fields in
  // the IO file, by appending each grid's suffix to the fields names.
  // NOTE: we'd need setter/getter for this, so we might as well make it public
  std::string m_short_name = "";

protected:

  void copy_data (const AbstractGrid& src, const bool shallow = true);

  // Note: this method must be called from the derived classes,
  //       since it calls get_2d_scalar_layout.
  void create_dof_fields (const int scalar2d_layout_rank);

private:

  // The grid name and type
  GridType     m_type;
  std::string  m_name;

  std::vector<std::string> m_aliases;

  // Counters
  int m_num_local_dofs;
  int m_num_global_dofs;
  int m_num_vert_levs;

  // The global ID of each dof
  Field     m_dofs_gids;

  // The max/min dof GID across all ranks. Mutable, to allow for lazy calculation
  mutable gid_type  m_global_min_dof_gid =  std::numeric_limits<gid_type>::max();
  mutable gid_type  m_global_max_dof_gid = -std::numeric_limits<gid_type>::max();

  // The map lid->idx
  Field     m_lid_to_idx;

  std::map<std::string,Field>  m_geo_fields;

  // The MPI comm containing the ranks across which the global mesh is partitioned
  ekat::Comm            m_comm;
};

} // namespace scream

#endif // SCREAM_ABSTRACT_GRID_HPP
