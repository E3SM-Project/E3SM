#ifndef SCREAM_ABSTRACT_GRID_HPP
#define SCREAM_ABSTRACT_GRID_HPP

#include "share/grid/grid_utils.hpp"
#include "share/util/factory.hpp"
#include "share/parameter_list.hpp"
#include "share/scream_types.hpp"
#include "share/util/string_utils.hpp"

namespace scream
{

class AbstractGrid
{
public:
  // TODO: template everything on the device type
  using device_type   = DefaultDevice;
  using kokkos_types  = KokkosTypes<device_type>;
  using dofs_map_type = kokkos_types::view<int*[4]>; // elem, igp, jgp, col_gid

  virtual ~AbstractGrid () = default;

  virtual GridType type () const = 0;

  virtual std::string name () const = 0;

  virtual int num_dofs () const = 0;

  virtual dofs_map_type get_dofs_map () const = 0;
};

} // namespace scream

#endif // SCREAM_ABSTRACT_GRID_HPP
