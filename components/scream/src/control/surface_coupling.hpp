#ifndef SCREAM_SURFACE_COUPLING_HPP
#define SCREAM_SURFACE_COUPLING_HPP

#include "share/field/field.hpp"
#include "share/scream_types.hpp"

#include <unordered_map>

namespace scream {
namespace control {

// This class is responsible to import/export fields from/to the rest of
// E3SM, via the mct coupler.
class SurfaceCoupling {
public:

  using host_device_type = HostDevice;
  using device_type      = DefaultDevice;

  using import_field_type  = Field<      Real>;
  using export_field_type  = Field<const Real>;

  // MCT always passes data as double
  using export_data_ptr_type =       double*;
  using import_data_ptr_type = const double*;

  // Almost a default c-tor. Only needs to set the repo state to Open
  SurfaceCoupling ();

  // Import host fields from the component coupler to device fields in the AD
  void do_import ();

  // Export device fields from the AD to host fields in the component coupler
  void do_export ();

  // Register import/export scream fields
  void register_import_field (const import_field_type& field);
  void register_export_field (const export_field_type& field);

  // Register import/export host pointers
  void register_import_data_ptr (const std::string& fname, import_data_ptr_type ptr);
  void register_export_data_ptr (const std::string& fname, export_data_ptr_type ptr);

  // Marks the end of the registration phase. Here, we check that there are no
  // field-data pairs where only one (field or data) has been set.
  void registration_ends ();

protected:

  // I could use std:;pair, but then you need to remember which one is p.first and p.second.
  // This way we can be more verbose.
  struct ImportPair {
    import_data_ptr_type  data;
    import_field_type     field;
  };
  struct ExportPair {
    export_data_ptr_type  data;
    export_field_type     field;
  };

  std::unordered_map<std::string, ImportPair>    m_imports;
  std::unordered_map<std::string, ExportPair>    m_exports;

  RepoState     m_state;
};

} // namespace control
} // namespace scream

#endif // SCREAM_SURFACE_COUPLING_HPP
