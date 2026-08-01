#include "share/data_managers/model_init.hpp"

namespace scream
{

DefaultModelInit::
DefaultModelInit (const ekat::ParameterList& params)
 : m_params(params)
{
  // nothing to do
}

void DefaultModelInit::
run (const ModelInitContext& ctx)
{
  if (ctx.run_type==RunType::Restart) {
    ctx.restart_model();
  } else {
    ctx.set_initial_conditions();
  }
}

std::shared_ptr<ModelInit>
create_model_init (const ekat::ParameterList& params)
{
  return std::make_shared<DefaultModelInit>(params);
}

void register_default_model_init ()
{
  static bool registered = false;
  if (not registered) {
    auto& f = ModelInitFactory::instance();
    f.register_product("default",&create_model_init);
    registered = true;
  }
}

} // namespace scream
