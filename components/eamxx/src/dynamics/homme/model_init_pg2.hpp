#ifndef EAMXX_MODEL_INIT_PG2_HPP
#define EAMXX_MODEL_INIT_PG2_HPP

#include "share/data_managers/model_init.hpp"

namespace scream
{

// A ModelInit for runs using Homme's PG2 (finite-volume, 2x2 cells per
// element) physics grid.
//
// Background: when PG2 is active, HommeDynamics additionally registers a
// (Required) "physics_gll" copy of each atm-input state field it owns
// (T_mid, horiz_winds, ps, phis, the tracers group), alongside the normal
// (Updated) "physics_pg2"-grid copy that the rest of the atm procs actually
// use. "physics_gll" exists *only* for this purpose -- no other atm proc
// registers fields on it -- and is a transient IC-staging grid: Homme uses
// it to derive the dynamics-grid state, and from there (via its own
// internal, lossy, conservative remap) the physics_pg2-grid state, during
// its own initialization. ModelInit must therefore:
//   - init the "physics_gll" copies normally (constant/copy/file, exactly
//     like any other field): this is what Homme's remap needs as input;
//   - NOT attempt to init any "physics_pg2" field from the main IC file:
//     state fields with a physics_gll companion are Homme's to remap, and
//     the IC file (GLL-native) has no data for any of the rest either
//     (even though it may, incidentally, carry same-named variables for
//     some other configuration) -- they start at their default value and
//     get filled in by their owning process on its first call;
//   - NOT read "physics_gll"'s own RESTART data (that grid is never part
//     of a restart file, and does not even survive past atm-proc init);
//   - skip "physics_pg2"'s topography fields that have a physics_gll
//     companion (currently, just phis): Homme computes those via remap
//     too. Fields with no companion (sgh, sgh30) are physics_pg2-only, and
//     are read from the topography file normally.
//
// Since "physics_gll" is exclusively Homme's own staging grid, "does a
// same-named field also exist on physics_gll" is an *exact* (not
// heuristic) test for "this is one of Homme's remapped state fields" --
// no hardcoded field list is needed here, and this class needs no
// knowledge of which fields Homme happens to remap.
class ModelInitPG2 : public ModelInit {
public:
  using ModelInit::ModelInit;

protected:
  std::vector<Field>
  get_leaf_fields (const std::shared_ptr<FieldManager>& fm,
                   const std::string& group_name,
                   const std::string& grid_name) override;

  std::vector<Field>
  get_fields (const std::shared_ptr<FieldManager>& fm,
              const std::string& group_name,
              const std::string& grid_name) override;
};

} // namespace scream

#endif // EAMXX_MODEL_INIT_PG2_HPP
