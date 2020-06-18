#ifndef SCREAM_SE_GRID_HPP
#define SCREAM_SE_GRID_HPP

#include "share/grid/abstract_grid.hpp"
#include "ekat/scream_types.hpp"
#include "ekat/scream_assert.hpp"

namespace scream
{

class SEGrid : public AbstractGrid
{
public:
  using base_type     = AbstractGrid;
  using dofs_map_type = kokkos_types::view<int*[3]>; // elem, igp, jgp

  SEGrid (const std::string& grid_name,
          const GridType type)
   : SEGrid(grid_name,type,0)
  {
    // Nothing to do here
  }

  SEGrid (const std::string& grid_name,
          const GridType type,
          const int num_local_dofs)
   : SEGrid(dofs_list_type("",num_local_dofs),grid_name,type)
  {
    // Put some invalid number in the list of dofs
    Kokkos::deep_copy(m_dofs_gids,-1);
  }

  SEGrid (const dofs_list_type& dofs_gids,
          const std::string& grid_name,
          const GridType type)
   : SEGrid (dofs_map_type("dof_to_elgp",dofs_gids.extent_int(0)),
             dofs_gids,grid_name,type)
  {
    // Put some invalid number in the dof_to_elgp view
    Kokkos::deep_copy(m_dof_to_elgp,-1);
  }

  SEGrid (const dofs_map_type&  dof_to_elgp,
          const dofs_list_type& dofs_gids,
          const std::string& grid_name,
          const GridType type)
   : m_grid_name      (grid_name)
   , m_type           (type)
   , m_num_local_dofs (dofs_gids.extent_int(0))
   , m_dofs_gids      (dofs_gids)
   , m_dof_to_elgp    (dof_to_elgp)
  {
    scream_require_msg(type==GridType::SE_CellBased || type==GridType::SE_NodeBased,
                       "Error! Grid type not (yet) supported by SEGrid.\n");
    scream_require_msg(dofs_gids.extent_int(0)==dof_to_elgp.extent_int(0),
                       "Error! Dofs gids and dofs map views have mismatching dimensions.\n");
  }

  virtual ~SEGrid () = default;

  // Getters
  GridType type () const { return m_type; }

  const std::string& name () const { return m_grid_name; }

  int get_num_local_dofs () const { return m_num_local_dofs; }

  const dofs_list_type& get_dofs_gids () const { return m_dofs_gids; }

  dofs_map_type get_dofs_map () const { return m_dof_to_elgp; }

protected:

  const std::string   m_grid_name;
  const GridType      m_type;

  int                 m_num_local_dofs;
  dofs_list_type      m_dofs_gids;
  dofs_map_type       m_dof_to_elgp;
};

} // namespace scream

#endif // SCREAM_SE_GRID_HPP
