#ifndef EAMXX_MODEL_INIT_HOMME_HPP
#define EAMXX_MODEL_INIT_HOMME_HPP

#include "share/data_managers/model_init.hpp"

namespace scream
{

// Base ModelInit for runs using Homme as the dynamics package, regardless of
// whether the physics grid is the dynamics-native GLL grid (ModelInitGLL) or
// Homme's PG2 finite-volume grid (ModelInitPG2). On top of the usual
// constant/copy/file field initialization (handled by the base ModelInit
// class), a startup run also needs Homme's dynamics-grid state (the
// "dynamics" grid's prognostic fields) populated from the "physics_gll"
// grid's (by then already-inited) initial condition: this is what
// init_homme_dyn_state_from_gll_ic (declared in
// dynamics/homme/eamxx_homme_process_interface.hpp) does, operating on the
// dyn_state_fields Homme has, by then, made available via
// EamxxHommeContext::singleton() -- the same Field objects HommeDynamics
// itself uses, so this class needs no knowledge of HommeDynamics (or of any
// atmosphere process) to do its job.
class ModelInitHomme : public ModelInit {
public:
  using ModelInit::ModelInit;

  void run (const std::shared_ptr<FieldManager>& fm,
            const util::TimeStamp& t0,
            const RunType run_type) override;
};

// A ModelInitHomme for runs where the physics grid is Homme's own GLL grid
// (i.e., no PG2 finite-volume remap involved): the base ModelInit's straight
// FieldManager/file interaction is enough to select which fields need
// initialization, so this class adds no override on top of ModelInitHomme.
class ModelInitGLL : public ModelInitHomme {
public:
  using ModelInitHomme::ModelInitHomme;
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_HOMME_HPP
