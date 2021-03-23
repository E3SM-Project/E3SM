#ifndef SCREAM_CLDFRACTION_HPP
#define SCREAM_CLDFRACTION_HPP

#include "physics/cld_fraction/cldfraction_functions.hpp"
#include "share/atm_process/atmosphere_process.hpp"
#include "ekat/ekat_parameter_list.hpp"

#include <string>

namespace scream
{

/*
 * The class responsible to handle the calculation of the subgrid cloud fractions
 *
 * The AD should store exactly ONE instance of this class stored
 * in its list of subcomponents (the AD should make sure of this).
*/

  using namespace cldfrac;
  using cldfracF = Functions<Real, DefaultDevice>;
  using Spack    = cldfracF::Spack;
  using Smask    = cldfracF::Smask;
  using Pack     = ekat::Pack<Real,Spack::n>;

class CldFraction : public AtmosphereProcess
{
public:
  using field_type       = Field<      Real>;
  using const_field_type = Field<const Real>;

  // Constructors
  CldFraction (const ekat::Comm& comm, const ekat::ParameterList& params);

  // The type of subcomponent
  AtmosphereProcessType type () const { return AtmosphereProcessType::Physics; }

  // The name of the subcomponent
  std::string name () const { return "CldFraction"; }

  // The communicator used by subcomponent
  const ekat::Comm& get_comm () const { return m_cldfraction_comm; }

  // Get the required grid for subcomponent
  std::set<std::string> get_required_grids () const {
    static std::set<std::string> s;
    s.insert(m_cldfraction_params.get<std::string>("Grid"));
    return s;
  }

  // Set the grid
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real>& field_repo) const;

  // Get the set of required/computed fields
  const std::set<FieldIdentifier>& get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>& get_computed_fields () const { return m_computed_fields; }

protected:

  // The three main overrides for the subcomponent
  void initialize_impl (const util::TimeStamp& t0);
  void run_impl        (const Real dt);
  void finalize_impl   ();

  // Setting the fields in the atmospheric process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

  std::set<FieldIdentifier> m_required_fields;
  std::set<FieldIdentifier> m_computed_fields;

  std::map<std::string,const_field_type>  m_cldfraction_fields_in;
  std::map<std::string,field_type>        m_cldfraction_fields_out;

  template<typename T>
  using view_type = field_type::view_type<T*>;

  template<typename T>
  using host_view_type = field_type::get_view_type<view_type<T>,Host>;

  using host_view_in_type   = host_view_type<const_field_type::RT>;
  using host_view_out_type  = host_view_type<      field_type::RT>;

  std::map<std::string,host_view_in_type>   m_cldfraction_host_views_in;
  std::map<std::string,host_view_out_type>  m_cldfraction_host_views_out;

  std::map<std::string,const Real*>  m_raw_ptrs_in;
  std::map<std::string,Real*>        m_raw_ptrs_out;

  util::TimeStamp     m_current_ts;
  ekat::Comm          m_cldfraction_comm;
  ekat::ParameterList m_cldfraction_params;

  // Keep track of field dimensions and the iteration count
  Int m_num_cols; 
  Int m_num_levs;
  Int m_nk_pack;

}; // class CldFraction

} // namespace scream

#endif // SCREAM_CLDFRACTION_HPP
