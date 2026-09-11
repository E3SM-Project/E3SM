#include "dynamics/homme/model_init_pg2.hpp"

namespace scream
{

std::vector<Field>
ModelInitPG2::
get_leaf_fields (const std::shared_ptr<FieldManager>& fm,
                 const std::string& group_name,
                 const std::string& grid_name)
{
  if (group_name=="STARTUP" and grid_name!="physics_gll") {
    // No field on physics_pg2 is read from the main IC file: state fields
    // with a physics_gll companion get it via Homme's GLL->PG2 remap
    // instead (see class doc), while physics_pg2-only fields (e.g. some
    // physics packages' own "previous step"/tendency-tracking fields) start
    // at their default (zero) value and get filled in by their owning
    // process on its first call. This matches the IC file possibly also
    // carrying (unused, for this configuration) data for such fields.
    return {};
  }
  return ModelInit::get_leaf_fields(fm,group_name,grid_name);
}

std::vector<Field>
ModelInitPG2::
get_fields (const std::shared_ptr<FieldManager>& fm,
            const std::string& group_name,
            const std::string& grid_name)
{
  if (group_name=="RESTART" and grid_name=="physics_gll") {
    // physics_gll is a transient IC-staging grid: it is never part of a
    // restart file (nor of the FieldManager, past atm-proc init).
    return {};
  }

  auto fields = ModelInit::get_fields(fm,group_name,grid_name);

  if (group_name=="TOPOGRAPHY" and grid_name!="physics_gll") {
    // Same rule as get_leaf_fields: a field also present on physics_gll
    // (e.g. phis) is Homme's to remap, not ours to read from file. Fields
    // with no physics_gll companion (e.g. sgh/sgh30) flow through as-is.
    std::vector<Field> filtered;
    for (auto& f : fields) {
      if (not fm->has_field(f.name(),"physics_gll")) {
        filtered.push_back(f);
      }
    }
    return filtered;
  }
  return fields;
}

} // namespace scream
