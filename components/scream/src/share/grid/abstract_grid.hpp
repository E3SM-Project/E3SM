#ifndef SCREAM_ABSTRACT_GRID_HPP
#define SCREAM_ABSTRACT_GRID_HPP

#include "share/grid/grid_utils.hpp"
#include "share/scream_types.hpp"

namespace scream
{

class AbstractGrid
{
public:
  using gid_type       = long;          // TODO: template class on gid type?
  using device_type    = DefaultDevice; // TODO: template class on device type
  using kokkos_types   = KokkosTypes<device_type>;
  using dofs_list_type = kokkos_types::view_1d<gid_type>;

  virtual ~AbstractGrid () = default;

  virtual GridType type () const = 0;

  virtual const std::string& name () const = 0;

  virtual int get_num_local_dofs () const = 0;

  virtual const dofs_list_type& get_dofs_gids () const = 0;
};

} // namespace scream

#endif // SCREAM_ABSTRACT_GRID_HPP
