#include "control/surface_coupling.hpp"

namespace scream {
namespace control {

SurfaceCoupling::
SurfaceCoupling ()
 : m_state (RepoState::Open)
{
  // Nothing to do here
}

void SurfaceCoupling::do_import () {
}

void SurfaceCoupling::do_export () {
}

void SurfaceCoupling::register_import_field (const import_field_type& field) {
  EKAT_REQUIRE_MSG (m_state==RepoState::Open, "Error! Registration phase has already ended.\n");
  const auto& fname = field.get_header().get_identifier().name();
  m_imports[fname].field = field;
}

void SurfaceCoupling::register_export_field (const export_field_type& field) {
  EKAT_REQUIRE_MSG (m_state==RepoState::Open, "Error! Registration phase has already ended.\n");
  const auto& fname = field.get_header().get_identifier().name();
  m_exports[fname].field = field;
}

void SurfaceCoupling::register_import_data_ptr (const std::string& fname, import_data_ptr_type data) {
  EKAT_REQUIRE_MSG (m_state==RepoState::Open, "Error! Registration phase has already ended.\n");
  m_imports[fname].data = data;
}

void SurfaceCoupling::register_export_data_ptr (const std::string& fname, export_data_ptr_type data) {
  EKAT_REQUIRE_MSG (m_state==RepoState::Open, "Error! Registration phase has already ended.\n");
  m_exports[fname].data = data;
}

void SurfaceCoupling::registration_ends () {
  // TODO: should I simply return instead?
  EKAT_REQUIRE_MSG (m_state==RepoState::Open, "Error! Registration phase has already ended.\n");

  // Loop over import/exports, and make sure both data and field are set
  for (const auto& it : m_imports) {
    EKAT_REQUIRE_MSG (it.second.data!=nullptr,
                        "Error! No data pointer set for '" + it.first + "'.\n");
    EKAT_REQUIRE_MSG (it.second.field.get_view().data()!=nullptr,
                        "Error! No field set for '" + it.first + "'.\n");
  }

  // Finally, mark registration as completed.
  m_state = RepoState::Closed;
}

} // amespace control
}  // namespace scream
