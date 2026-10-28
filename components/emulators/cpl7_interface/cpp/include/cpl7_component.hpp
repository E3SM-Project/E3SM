#ifndef E3SM_EMULATOR_CPL7_COMPONENT_HPP
#define E3SM_EMULATOR_CPL7_COMPONENT_HPP

#include "field_registry.hpp"
#include "process.hpp"
#include <coupler.hpp>
#include <stdexcept>
#include <string_view>
namespace e3sm::emulator {
using namespace e3sm::coupler;

class CPL7Component {
public:
  static constexpr std::string_view name() { return "cpl7"; }
  void populate_registry(FieldRegistry&);

  void configure_coupling(const ActiveCouplingFields& plan) {
    coupling_.configure(plan);
  }

  void export_fields(std::span<ExportBuffer> buffers) const {
    coupling_.export_fields(buffers);
  }

  void import_fields(std::span<const ImportBuffer> buffers) {
    coupling_.import_fields(buffers);
  }

  void run(const ProcessRunOpts opts) {
    throw std::runtime_error("DON'T CALL RUN ON CPL7Component");
  }

private:
  CouplingAdapter coupling_;
  std::span<double> external_to_emulator_;
  std::span<double> emulator_to_external_;
};

} // namespace e3sm::emulator

#endif
