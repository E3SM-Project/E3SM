/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_MPI_CONTEXT_HPP
#define HOMMEXX_MPI_CONTEXT_HPP

#include <string>
#include <map>
#include <memory>

namespace Homme {

class BuffersManager;
class Comm;
class Connectivity;

/* A MpiContext manages resources previously treated as singletons. MpiContext is
 * meant to have two roles. First, a MpiContext singleton is the only singleton in
 * the program. Second, a context need not be a singleton, and each MpiContext
 * object can have different Elements, Derivative, etc., objects. (That
 * probably isn't needed, but MpiContext immediately supports it.)
 *
 * Finally, MpiContext has two singleton functions: singleton(), which returns
 * MpiContext&, and finalize_singleton(). The second is called in a unit test exe
 * main before Kokkos::finalize().
 */
class MpiContext {
public:
  using BMMap = std::map<int,std::shared_ptr<BuffersManager>>;

private:
  // Note: using uniqe_ptr disables copy construction
  std::unique_ptr<Comm>                   comm_;
  std::shared_ptr<Connectivity>           connectivity_;
  std::shared_ptr<BMMap>                  buffers_managers_;

  // Clear the objects MpiContext manages.
  void clear();

public:
  MpiContext();
  virtual ~MpiContext();

  // Getters for each managed object.
  void create_comm(const int f_comm);
  Comm& get_comm();
  std::shared_ptr<Connectivity> get_connectivity();
  BMMap& get_buffers_managers();
  std::shared_ptr<BuffersManager> get_buffers_manager(short int exchange_type);

  // Exactly one singleton.
  static MpiContext& singleton();

  static void finalize_singleton();
};

} // namespace Homme

#endif // HOMMEXX_MPI_CONTEXT_HPP
