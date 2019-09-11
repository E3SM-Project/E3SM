#ifndef SCREAM_SE_GRID_HPP
#define SCREAM_SE_GRID_HPP

#include "share/grid/abstract_grid.hpp"
#include "share/scream_types.hpp"
#include "share/scream_assert.hpp"

namespace scream
{

class SEGrid : public AbstractGrid
{
public:
  // TODO: template everything on the device type
  using device_type   = DefaultDevice;
  using kokkos_types  = KokkosTypes<device_type>;
  using dofs_map_type = kokkos_types::view<int*[4]>; // elem, igp, jgp, col_gid

  SEGrid (const std::string& name, const GridType type)
   : m_num_dofs (0)
   , m_name     (name)
   , m_type     (type)
  {
    scream_require_msg(type==GridType::SE_CellBased || type==GridType::SE_NodeBased,
                       "Error! Grid type not (yet) supported by DefaultGrid.\n");
  }

  SEGrid (dofs_map_type col_to_elgp, const std::string& name, const GridType type)
   : m_col_to_elgp (col_to_elgp)
   , m_num_dofs    (m_col_to_elgp.extent_int(0))
   , m_name        (name)
   , m_type        (type)
  {
    scream_require_msg(type==GridType::SE_CellBased || type==GridType::SE_NodeBased,
                       "Error! Grid type not (yet) supported by DefaultGrid.\n");
  }

  virtual ~SEGrid () = default;

  GridType type () const override { return m_type; }

  std::string name () const override { return m_name; }

  int get_num_dofs () const override { return m_num_dofs; }

  dofs_map_type get_dofs_map () const { return m_col_to_elgp; }

protected:
  dofs_map_type   m_col_to_elgp;
  int             m_num_dofs;
  std::string     m_name;

  GridType        m_type;
};

} // namespace scream

#endif // SCREAM_SE_GRID_HPP
