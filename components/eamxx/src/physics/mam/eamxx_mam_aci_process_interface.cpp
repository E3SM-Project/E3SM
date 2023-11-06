#include <physics/mam/eamxx_mam_aci_process_interface.hpp>

namespace scream
{

MAMAci::MAMAci(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params){
}

} // namespace scream
