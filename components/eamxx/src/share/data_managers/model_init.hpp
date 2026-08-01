#ifndef SCREAM_MODEL_INIT_HPP
#define SCREAM_MODEL_INIT_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_factory.hpp>
#include <ekat_parameter_list.hpp>

#include <functional>
#include <memory>

namespace scream
{

struct ModelInitContext {
  RunType run_type;
  std::function<void()> restart_model;
  std::function<void()> set_initial_conditions;
};

class ModelInit {
public:
  virtual ~ModelInit () = default;

  virtual void run (const ModelInitContext& ctx) = 0;
};

class DefaultModelInit : public ModelInit {
public:
  explicit DefaultModelInit (const ekat::ParameterList& params);

  void run (const ModelInitContext& ctx) override;

private:
  ekat::ParameterList m_params;
};

std::shared_ptr<ModelInit>
create_model_init (const ekat::ParameterList& params);

void register_default_model_init ();

using ModelInitFactory
    = ekat::Factory<ModelInit,
                    ekat::CaseInsensitiveString,
                    std::shared_ptr<ModelInit>,
                    const ekat::ParameterList&>;

} // namespace scream

#endif // SCREAM_MODEL_INIT_HPP
