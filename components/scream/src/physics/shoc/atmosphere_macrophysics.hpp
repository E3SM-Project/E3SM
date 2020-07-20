#ifndef SCREAM_SHOC_MACROPHYSICS_HPP
#define SCREAM_SHOC_MACROPHYSICS_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/scream_parameter_list.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the atmosphere microphysics
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 *
 *  Note: for now, scream is only going to accommodate SHOC as macrophysics
*/

class SHOCMacrophysics : public scream::AtmosphereProcess
{
public:
  using field_type       = Field<      Real,device_type>;
  using const_field_type = Field<const Real,device_type>;

  // Constructors
  SHOCMacrophysics (const Comm& comm, const ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Macrophysics"; }

  // The communicator used by subcomponent
  const Comm& get_comm () const { return m_shoc_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert("Physics");
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // The three main interfaces for the subcomponent
  void initialize (const util::TimeStamp& t0);
  void run        (const Real dt);
  void finalize   ();

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const;

  // Get the set of required/computed fields
  const std::set<FieldIdentifier>& get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>& get_computed_fields () const { return m_computed_fields; }

protected:

  // Setting the fields in the atmospheric process
  void set_required_field_impl (const Field<const Real, device_type>& f);
  void set_computed_field_impl (const Field<      Real, device_type>& f);

  std::set<FieldIdentifier> m_required_fields;
  std::set<FieldIdentifier> m_computed_fields;

  std::map<std::string,const_field_type>  m_shoc_fields_in;
  std::map<std::string,field_type>        m_shoc_fields_out;

  util::TimeStamp   m_current_ts;
  Comm              m_shoc_comm;

}; // class SHOCMacrophysics

} // namespace scream

#endif // SCREAM_SHOC_MACROPHYSICS_HPP
