#include "eamxx_gw_emulator_process.hpp"

#include <ekat_assert.hpp>
#include <ekat_units.hpp>

namespace scream
{

GWEmulator::
GWEmulator(const ekat::Comm& comm, const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params)
{
  // Nothing to do here
}

void GWEmulator::
set_grids(const std::shared_ptr<const GridsManager> grids_manager)
{
  m_grid = grids_manager->get_grid("physics");

  auto 
}

void GWEmulator::initialize_impl (const RunType /* run_type */)
{
#ifdef EAMXX_ENABLE_PYTHON
  if (has_py_module()) {
    int nlevs = m_grid->get_num_vertical_levels();
    int nn_layers = m_params.get("num_hidden_layers",1);
    const auto& state_dict_file = m_params.get<std::string>("state_dict_file");
    py_module_call("init",nlevs,nn_layers,state_dict_file);
  }
#endif
}

void GWEmulator::run_impl (const double dt)
{

}

void GWEmulator::finalize_impl ()
{

}

} // namespace scream
