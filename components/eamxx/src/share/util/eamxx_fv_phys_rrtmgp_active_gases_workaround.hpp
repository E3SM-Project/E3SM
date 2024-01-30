#include "share/grid/remap/abstract_remapper.hpp"

#include "ekat/ekat_parameter_list.hpp"

namespace scream {

// TODO [rrtmgp active gases] This is to address issue #1782. It supports option
// 1 in that issue. These fv_phys_rrtmgp_active_gases_* routines can be removed
// once rrtmgp active_gases initialization is treated properly.

struct TraceGasesWorkaround
{
private:

  TraceGasesWorkaround() { restart = false; }

  bool restart;
  std::shared_ptr<AbstractRemapper> remapper;
  std::vector<std::string> active_gases; // other than h2o

public:

  static TraceGasesWorkaround& singleton() {
    static TraceGasesWorkaround self;
    return self;
  }

  void set_restart (const bool is_restart) {
    restart = is_restart;
  }
  void set_remapper (const std::shared_ptr<AbstractRemapper>& remap_ptr) {
    remapper = remap_ptr;
  }
  void erase_remapper () { remapper = nullptr; }
  void add_active_gas (const std::string gas_name) {
    active_gases.push_back(gas_name);
  }

  bool is_restart() const { return restart; }
  std::shared_ptr<AbstractRemapper> get_remapper() const { return remapper; }
  std::vector<std::string> get_active_gases() const { return active_gases; }
};

extern bool fvphyshack;
void fv_phys_rrtmgp_active_gases_set_restart(const bool restart);

} // namespace scream
