#ifndef SCREAM_MODEL_INIT_HPP
#define SCREAM_MODEL_INIT_HPP

#include "share/core/eamxx_types.hpp"
#include "share/util/eamxx_time_stamp.hpp"

#include <ekat_factory.hpp>
#include <ekat_parameter_list.hpp>
#include <ekat_comm.hpp>

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace ekat {
namespace logger { class LoggerBase; }
}

namespace scream
{

// Forward declarations
class FieldManager;
class GridsManager;
class AtmosphereProcessGroup;
class IOPDataManager;

/*
 * Base class for model initialization.
 *
 * Concrete derived classes implement restart_model() and
 * set_initial_conditions(), then run() dispatches to the appropriate
 * one based on the run type provided to setup().
 *
 * Usage:
 *   auto mi = ModelInitFactory::instance().create(type, params);
 *   mi->setup(comm, atm_params, field_mgr, grids_mgr, atm_procs,
 *             iop_mgr, logger, current_ts, run_type);
 *   mi->run();
 */
class ModelInit {
public:
  using strmap_t = std::map<std::string, std::vector<std::string>>;

  virtual ~ModelInit () = default;

  // Provide the data needed by run(). Must be called before run().
  void setup (const ekat::Comm& comm,
              const ekat::ParameterList& atm_params,
              const std::shared_ptr<FieldManager>& field_mgr,
              const std::shared_ptr<GridsManager>& grids_mgr,
              const std::shared_ptr<AtmosphereProcessGroup>& atm_procs,
              const std::shared_ptr<IOPDataManager>& iop_data_mgr,
              const std::shared_ptr<ekat::logger::LoggerBase>& logger,
              const util::TimeStamp& current_ts,
              RunType run_type);

  // Dispatch to restart_model() or set_initial_conditions() based on run type.
  virtual void run () = 0;

  // Fields that were initialized (for DAG processing).
  const strmap_t& get_fields_inited () const { return m_fields_inited; }

protected:
  virtual void restart_model () = 0;
  virtual void set_initial_conditions () = 0;

  ekat::Comm                                      m_comm;
  ekat::ParameterList                             m_atm_params;
  std::shared_ptr<FieldManager>                   m_field_mgr;
  std::shared_ptr<GridsManager>                   m_grids_mgr;
  std::shared_ptr<AtmosphereProcessGroup>         m_atm_procs;
  std::shared_ptr<IOPDataManager>                 m_iop_data_mgr;
  std::shared_ptr<ekat::logger::LoggerBase>       m_logger;
  util::TimeStamp                                 m_current_ts;
  RunType                                         m_run_type;
  strmap_t                                        m_fields_inited;
};

using ModelInitFactory
    = ekat::Factory<ModelInit,
                    ekat::CaseInsensitiveString,
                    std::shared_ptr<ModelInit>,
                    const ekat::ParameterList&>;

} // namespace scream

#endif // SCREAM_MODEL_INIT_HPP
