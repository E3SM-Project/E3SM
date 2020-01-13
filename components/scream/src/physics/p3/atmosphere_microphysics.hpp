#ifndef SCREAM_P3_MICROPHYSICS_HPP
#define SCREAM_P3_MICROPHYSICS_HPP

#include "share/atmosphere_process.hpp"
#include "share/parameter_list.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the atmosphere microphysics
 * 
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
 * 
 *  Note: for now, scream is only going to accommodate P3 as microphysics 
*/

class P3Microphysics : public scream::AtmosphereProcess
{
public:
  using field_type       = Field<      Real,device_type>;
  using const_field_type = Field<const Real,device_type>;

  // Constructors
  P3Microphysics (const Comm& comm, const ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Microphysics"; }
  
  // The communicator used by subcomponent
  const Comm& get_comm () const { return m_p3_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_p3_params.get<std::string>("Grid"));
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

  std::map<std::string,const_field_type>  m_p3_fields_in;
  std::map<std::string,field_type>        m_p3_fields_out;

  std::map<std::string,const_field_type::view_type::HostMirror>  m_p3_host_views_in;
  std::map<std::string,field_type::view_type::HostMirror>        m_p3_host_views_out;

  std::map<std::string,const Real*>  m_raw_ptrs_in;
  std::map<std::string,Real*>        m_raw_ptrs_out;

  util::TimeStamp   m_current_ts;
  Comm              m_p3_comm;

  ParameterList     m_p3_params;

}; // class P3Microphysics

/*
 * This class is microphysics stub meant to be used for initialization of
 * stand-alone tests.  It is similar to P3Microphysics except that it won't
 * actually calculate anything, it will only initialize the input fields needed
 * by atmophseric microphysics.
 * TODO: Switch to a general stand-alone class for initialization of any stand-alone
 * test.
*/
class SAMicrophysics : public scream::AtmosphereProcess
{
public:
  using field_type       = Field<      Real,device_type>;
  using const_field_type = Field<const Real,device_type>;

  // Constructors
  SAMicrophysics (const Comm& comm, const ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "Stand-Alone"; }
  
  // The communicator used by subcomponent
  const Comm& get_comm () const { return m_sa_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_sa_params.get<std::string>("Grid"));
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

  std::map<std::string,const_field_type>  m_sa_fields_in;
  std::map<std::string,field_type>        m_sa_fields_out;

  std::map<std::string,const_field_type::view_type::HostMirror>  m_sa_host_views_in;
  std::map<std::string,field_type::view_type::HostMirror>        m_sa_host_views_out;

  std::map<std::string,const Real*>  m_raw_ptrs_in;
  std::map<std::string,Real*>        m_raw_ptrs_out;

  util::TimeStamp   m_current_ts;
  Comm              m_sa_comm;

  ParameterList     m_sa_params;

}; // class SAMicrophysics

} // namespace scream

#endif // SCREAM_P3_MICROPHYSICS_HPP
