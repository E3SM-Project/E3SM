#include "diagnostics/diagnostics_interface.hpp"

namespace scream {

// =========================================================================================
DiagnosticInterface::DiagnosticInterface()
{
  // Nothing to do
}
// =========================================================================================
void DiagnosticInterface::add_diagnostic(const std::string name)
{
  printf("Add diagnostic %s\n",name.c_str());
}

} //namespace scream
