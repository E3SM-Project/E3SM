#ifndef E3SM_EMULATOR_CPL7_COMPONENT_HPP
#define E3SM_EMULATOR_CPL7_COMPONENT_HPP

#include "field_registry.hpp"
#include "process.hpp"
#include <coupler.hpp>
#include <stdexcept>
#include <string_view>
namespace e3sm::emulator {
using namespace e3sm::coupler;

struct ExternalField {
  FieldRole role;
  std::string label;
  std::size_t extent;
  std::size_t offset;
};

class CPL7Component {
public:
  //
  // Satisfy AnyComponent
  //
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

  void attach_fields(std::span<const ExternalField> requested_fields,
                     std::span<double> import_storage,
                     std::span<double> export_storage);

private:
  struct BoundField {
    RegisteredField field;
    std::span<double> storage;
  };
  std::vector<BoundField> fields_;
  CouplingAdapter coupling_;
  bool initialized_ = false;
};

} // namespace e3sm::emulator

#endif
