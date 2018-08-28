/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_VERTICAL_REMAP_MANAGER_HPP
#define HOMMEXX_VERTICAL_REMAP_MANAGER_HPP

#include <memory>

namespace Homme {

struct VerticalRemapManager {
  VerticalRemapManager();

  void run_remap(int np1, int np1_qdp, double dt) const;

private:
  struct Impl;
  std::shared_ptr<Impl> p_;
};

}

#endif
