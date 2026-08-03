#include "share/data_managers/model_init.hpp"

namespace scream
{

void ModelInit::
setup (const ekat::Comm& comm,
       const ekat::ParameterList& atm_params,
       const std::shared_ptr<FieldManager>& field_mgr,
       const std::shared_ptr<GridsManager>& grids_mgr,
       const std::shared_ptr<AtmosphereProcessGroup>& atm_procs,
       const std::shared_ptr<IOPDataManager>& iop_data_mgr,
       const std::shared_ptr<ekat::logger::LoggerBase>& logger,
       const util::TimeStamp& current_ts,
       RunType run_type)
{
  m_comm         = comm;
  m_atm_params   = atm_params;
  m_field_mgr    = field_mgr;
  m_grids_mgr    = grids_mgr;
  m_atm_procs    = atm_procs;
  m_iop_data_mgr = iop_data_mgr;
  m_logger       = logger;
  m_current_ts   = current_ts;
  m_run_type     = run_type;
}

} // namespace scream
