#ifndef SCREAM_ABSTRACT_GRID_HPP
#define SCREAM_ABSTRACT_GRID_HPP

#include "share/grid/grid_utils.hpp"

namespace scream
{

class AbstractGrid
{
public:
  virtual ~AbstractGrid () = default;

  virtual GridType type () const = 0;

  virtual std::string name () const = 0;

  virtual int get_num_dofs () const = 0;

  // TODO: perhaps store a View with the GID of the dofs?
  //       Sounds legit, but it may duplicate storage if the
  //       derived classes already store GIDs, but with a peculiar
  //       format/layout.
};

} // namespace scream

#endif // SCREAM_ABSTRACT_GRID_HPP
