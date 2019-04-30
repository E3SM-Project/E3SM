#ifndef SCREAM_DEFAULT_GRID_HPP
#define SCREAM_DEFAULT_GRID_HPP

#include "share/grid/abstract_grid.hpp"

namespace scream
{

template<GridType gridType>
class DefaultGrid : public AbstractGrid
{
public:
  // TODO: template everything on the device type
  using device_type   = DefaultDevice;
  using kokkos_types  = KokkosTypes<device_type>;
  using dofs_map_type = kokkos_types::view<int*[4]>; // elem, igp, jgp, col_gid

  DefaultGrid (const std::string& name)
   : m_num_dofs (0)
   , m_name     (name)
  {
    // Nothing to do here
  }

  DefaultGrid (dofs_map_type col_to_elgp, const std::string& name)
   : m_col_to_elgp (col_to_elgp)
   , m_num_dofs    (m_col_to_elgp.extent_int(0))
   , m_name        (name)
  {
    // Nothing to do here
  }

  virtual ~DefaultGrid () = default;

  GridType type () const { return gridType; }

  std::string name () const { return m_name; }

  int num_dofs () const { return m_num_dofs; }

  dofs_map_type get_dofs_map () const { return m_col_to_elgp; }

protected:
  dofs_map_type   m_col_to_elgp;
  int             m_num_dofs;
  std::string     m_name;
};

} // namespace scream

#endif // SCREAM_DEFAULT_GRID_HPP
