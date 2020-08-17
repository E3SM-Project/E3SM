#ifndef SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP
#define SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP

#include "share/grid/remap/abstract_remapper.hpp"
#include "share/atm_process/atmosphere_process.hpp"

#include "ekat/ekat_parameter_list.hpp"

#include <string>
#include <list>

namespace scream
{

/*
 *  A class representing a group of atmosphere processes as a single process.
 *
 *  This class allows to create nested lists of processes, by representing a group
 *  of processes as a single atmosphere process.
 *  All the calls to setup/run methods are simply forwarded to the stored list of
 *  atm processes, and the stored list of required/computed fields is simply a
 *  concatenation of the correspong lists in the underlying atm processes.
 */

class AtmosphereProcessGroup : public AtmosphereProcess
{
public:
  using atm_proc_type     = AtmosphereProcess;
  using remapper_type     = AbstractRemapper<Real, device_type>;
  using remapper_ptr_type = std::shared_ptr<remapper_type>;

  // Constructor(s)
  explicit AtmosphereProcessGroup (const ekat::Comm& comm, const ekat::ParameterList& params);

  virtual ~AtmosphereProcessGroup () = default;

  // The type of the block (e.g., dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Group; }

  // The type of grids on which the process is defined
  std::set<std::string> get_required_grids () const { return m_required_grids; }

  // The name of the block
  std::string name () const { return m_group_name; }

  // The communicator associated with this atm process
  const ekat::Comm& get_comm () const { return m_comm; }

  // Grab the proper grid from the grids manager
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // The initialization, run, and finalization methods
  void initialize (const TimeStamp& t0);
  void run        (const Real dt);
  void finalize   (/* what inputs? */);

  void final_setup ();

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real, device_type>& field_repo) const;

  // The methods used to query the process for its inputs/outputs
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_computed_fields; }
  const std::set<FieldIdentifier>&  get_internal_fields () const { return m_internal_fields; }

  // --- Methods specific to AtmosphereProcessGroup --- //

  int get_num_processes () const { return m_atm_processes.size(); }

  std::shared_ptr<const atm_proc_type> get_process (const int i) const {
    return m_atm_processes.at(i);
  }

  void setup_remappers (const FieldRepository<Real, device_type>& field_repo);

  const std::vector<std::map<std::string,remapper_ptr_type>>&
  get_inputs_remappers () const { return m_inputs_remappers; }

  const std::vector<std::map<std::string,remapper_ptr_type>>&
  get_outputs_remappers () const { return m_outputs_remappers; }

  ScheduleType get_schedule_type () const { return m_group_schedule_type; }

#ifdef SCREAM_DEBUG
  void set_field_repos (const FieldRepository<Real, device_type>& repo,
                        const FieldRepository<Real, device_type>& bkp_repo);
#endif

  void set_internal_field (const Field<Real, device_type>& f);

protected:

  void run_sequential (const Real dt);
  void run_parallel   (const Real dt);

  // The methods to set the fields in the process
  void set_required_field_impl (const Field<const Real, device_type>& f);
  void set_computed_field_impl (const Field<      Real, device_type>& f);

  // Method to build the identifier of a field on the reference grid given
  // an identifier on a different grid
  FieldIdentifier create_ref_fid (const FieldIdentifier& fid,
                                  const remapper_ptr_type& remapper);

  // The communicator that each process in this group uses
  ekat::Comm        m_comm;

  // The name of the group. This is usually a concatenation of the names of the individual processes
  std::string       m_group_name;
  int               m_group_size;

  // The list of atm processes in this group
  std::vector<std::shared_ptr<atm_proc_type>>  m_atm_processes;

  // The grids required by this process
  std::set<std::string>  m_required_grids;

  // The reference grid name.
  std::string m_ref_grid_name;

  // The schedule type: Parallel vs Sequential
  ScheduleType   m_group_schedule_type;

  // The cumulative set of required/computed fields of the atm processes in the group
  std::set<FieldIdentifier>      m_required_fields;
  std::set<FieldIdentifier>      m_computed_fields;

  // This is useful only for sequential splitting. Internal fields are fields that are
  // required fields for some atm proc, but that are computed by some other PREVIOUS
  // atm proc in the group. These fields may not be used outside of the group.
  std::set<FieldIdentifier>      m_internal_fields;

  // The remappers are to map output/input fields to/from the reference grid
  // Note: the i-th entry of the vector is a map of remappers needed by the i-th process.
  //       the map's key is the name of the non-reference grid.
  //       Notice that, as of today (07/2019), only AtmosphereProcessGroup (APG)
  //       can have more than one grid associated with it (the union of the
  //       grids of all the stored processes), while 'normal' atm processes
  //       should act on *only one grid*. Furthermore, we don't need to perform
  //       remapping for a stored APG, since the APG will already do it itself,
  //       so we only need remappers for 'normal' processes.
  //       It would then be tempting to replace std::map<std::string,remapper_ptr_type>
  //       with simply a remapper_ptr_type. This would be probably fine.
  //       However, if this assumption (only one grid per atm proc) becomes
  //       wrong, it may be hard for someone to jump in and adapt the code.
  //       At the very least, it would be harder than it is for me to write
  //       the code in a more general fashion right from the beginning.
  //       There's no real performance issue in coding for a more general
  //       scenario, and may prevent headaches in the future.
  std::vector<std::map<std::string,remapper_ptr_type>> m_inputs_remappers;
  std::vector<std::map<std::string,remapper_ptr_type>> m_outputs_remappers;

#ifdef SCREAM_DEBUG
  using field_type = Field<Real, device_type>;
  bool views_are_equal (const field_type& v1, const field_type& v2);

  const FieldRepository<Real, device_type>*   m_field_repo;
  const FieldRepository<Real, device_type>*   m_bkp_field_repo;

#endif

  TimeStamp                                   m_current_ts;
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP
