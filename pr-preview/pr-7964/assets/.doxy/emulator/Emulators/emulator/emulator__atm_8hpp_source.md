

# File emulator\_atm.hpp

[**File List**](files.md) **>** [**components**](dir_409f97388efe006bc3438b95e9edef48.md) **>** [**emulator\_comps**](dir_cd6ef227c082afa5b90fe3621cc9f093.md) **>** [**eatm**](dir_54689134e1a693092e83f56806593839.md) **>** [**src**](dir_1c3b735e18de9b9534f50214e18facf2.md) **>** [**emulator\_atm.hpp**](emulator__atm_8hpp.md)

[Go to the documentation of this file](emulator__atm_8hpp.md)


```C++


#ifndef EMULATOR_ATM_HPP
#define EMULATOR_ATM_HPP

#include "../../common/src/coupling_fields.hpp"
#include "../../common/src/emulator_comp.hpp"
#include "../../common/src/emulator_config.hpp"
#include "../../common/src/emulator_output_manager.hpp"
#include "../../common/src/inference/inference_backend.hpp"
#include "impl/atm_coupling.hpp"
#include "impl/atm_field_data_provider.hpp"
#include "impl/atm_field_manager.hpp"
#include <fstream>
#include <memory>
#include <string>
#include <vector>

namespace emulator {

class EmulatorAtm : public EmulatorComp {
public:
  EmulatorAtm();
  ~EmulatorAtm() override = default;

  void init_coupling_indices(const std::string &export_fields,
                             const std::string &import_fields);

  void set_inference_config(const inference::InferenceConfig &config);

  const inference::InferenceConfig &get_inference_config() const {
    return m_inference_config;
  }

  void set_log_file(const std::string &filename);

protected:
  void init_impl() override;

  void run_impl(int dt) override;

  void final_impl() override;

  void run_inference(const std::vector<double> &inputs,
                     std::vector<double> &outputs) override;

private:
  // =========================================================================
  // Coupling and Fields
  // =========================================================================
  CouplingFieldsBase m_coupling_fields;    
  impl::AtmCouplingIndices m_coupling_idx; 
  impl::AtmFieldManager m_fields;          

  // =========================================================================
  // Configuration and Inference
  // =========================================================================
  EmulatorConfig m_config;                       
  inference::InferenceConfig m_inference_config; 
  std::unique_ptr<inference::InferenceBackend>
      m_inference;                        
  EmulatorOutputManager m_output_manager; 
  std::unique_ptr<impl::AtmFieldDataProvider>
      m_field_provider; 

  // =========================================================================
  // Helper Methods
  // =========================================================================

  void import_coupling_fields();

  void export_coupling_fields();

  void prepare_inputs();

  void process_outputs();

  bool read_initial_conditions(const std::string &filename);
};

} // namespace emulator

#endif // EMULATOR_ATM_HPP
```


