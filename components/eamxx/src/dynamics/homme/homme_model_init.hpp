#ifndef EAMXX_HOMME_MODEL_INIT_HPP
#define EAMXX_HOMME_MODEL_INIT_HPP

#include "share/data_managers/model_init.hpp"

namespace scream
{

// A ModelInit for runs where Homme is the dynamics package.
//
// On top of the generic init of the atm input fields (done by ModelInit), this class takes
// care of initializing Homme itself, so that, by the time ModelInit::run returns, both the
// atm fields and Homme's internal state are fully valid:
//  - complete Homme's init1 phase, and set the views of the dyn-grid fields into Homme's C++
//    data structures (see HommeInstance);
//  - for an initial run, init Homme's dyn state from the GLL initial conditions; for a
//    restart run, finalize the restart of the dyn state (read from the restart file by
//    ModelInit::run) and restore Homme's time step counter;
//  - complete Homme's model init.
// This implies that this class requires the physics grid to be the GLL one, and that the
// HommeInstance has been setup by the HommeDynamics atm process (during its create_requests
// and set_computed_group_impl calls), and that the atm processes buffers have been allocated.
class HommeModelInit : public ModelInit {
public:
  using ModelInit::ModelInit;

  void run (const std::shared_ptr<FieldManager>& fm,
            const util::TimeStamp& t0,
            const RunType run_type) override;
};

} // namespace scream

#endif // EAMXX_HOMME_MODEL_INIT_HPP
