#include <physics/mam/eamxx_mam_aci_process_interface.hpp>

namespace scream
{
MAMAci::MAMAci(
    const ekat::Comm& comm,
    const ekat::ParameterList& params)
  : AtmosphereProcess(comm, params){
}


//Return type of the process
AtmosphereProcessType MAMAci::type() const {
  return AtmosphereProcessType::Physics;
}

//return name of the process
std::string MAMAci::name() const{
  return "mam4_aci";
  }

//set grid for all the inputs and outputs
void MAMAci::set_grids(const std::shared_ptr<const GridsManager> grids_manager) {
  m_atm_logger->log(ekat::logger::LogLevel::info,"ACI set grid = ");
}

void MAMAci::initialize_impl(const RunType run_type) {
   m_atm_logger->log(ekat::logger::LogLevel::info,"ACI init = ");
}

void MAMAci::run_impl(const double dt) {
  m_atm_logger->log(ekat::logger::LogLevel::info,"ACI run = ");
}

void MAMAci::finalize_impl(){
  m_atm_logger->log(ekat::logger::LogLevel::info,"ACI final = ");
}

} // namespace scream