#ifndef SCREAM_CLD_FRAC_NET_HPP
#define SCREAM_CLD_FRAC_NET_HPP

#include "share/atm_process/atmosphere_process.hpp"

namespace scream
{

/*
 * An ML emulator for the CldFraction process
 *
 * This process is NOT to be used in real runs, and is exclusively meant
 * to be an example of how to wrap a torch model in an eamxx atm process
*/

class CldFracNet : public AtmosphereProcess
{
public:
  CldFracNet (const ekat::Comm& comm, const ekat::ParameterList& params);

  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }
  std::string name () const { return "cld_frac_net"; }

  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

protected:

  void initialize_impl (const RunType run_type);
  void run_impl        (const double dt);
  void finalize_impl   ();
};

} // namespace scream

#endif // SCREAM_CLD_FRAC_NET_HPP
