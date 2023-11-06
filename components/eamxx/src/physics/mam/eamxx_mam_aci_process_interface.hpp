#ifndef EAMXX_MAM_ACI_HPP
#define EAMXX_MAM_ACI_HPP

#include <share/atm_process/atmosphere_process.hpp>

namespace scream
{

class MAMAci final : public scream::AtmosphereProcess {

public:
  // Constructor
  MAMAci(const ekat::Comm& comm, const ekat::ParameterList& params);
}; // MAMAci

} // namespace scream


#endif // EAMXX_MAM_ACI_HPP
