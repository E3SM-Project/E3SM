#ifndef SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP
#define SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP

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

  void final_setup ();

  void set_required_group (const FieldGroup<const Real>& group);
  void set_updated_group (const FieldGroup<Real>& group);

  // Register all fields in the given repo
  void register_fields (FieldRepository<Real>& field_repo) const;

  // The methods used to query the process for its inputs/outputs
  const std::set<FieldIdentifier>&  get_required_fields () const { return m_required_fields; }
  const std::set<FieldIdentifier>&  get_computed_fields () const { return m_computed_fields; }

  std::set<GroupRequest> get_required_groups () const;
  std::set<GroupRequest> get_updated_groups () const;

  // --- Methods specific to AtmosphereProcessGroup --- //

  int get_num_processes () const { return m_atm_processes.size(); }

  std::shared_ptr<const atm_proc_type> get_process (const int i) const {
    return m_atm_processes.at(i);
  }

  ScheduleType get_schedule_type () const { return m_group_schedule_type; }

#ifdef SCREAM_DEBUG
  void set_field_repos (const FieldRepository<Real>& repo,
                        const FieldRepository<Real>& bkp_repo);
#endif

protected:

  // Adds fid to the list of required/computed fields of the group (as a whole).
  void process_required_field (const FieldIdentifier& fid);

  // The initialization, run, and finalization methods
  void initialize_impl (const TimeStamp& t0);
  void run_impl        (const Real dt);
  void finalize_impl   (/* what inputs? */);

  void run_sequential (const Real dt);
  void run_parallel   (const Real dt);

  // The methods to set the fields in the process
  void set_required_field_impl (const Field<const Real>& f);
  void set_computed_field_impl (const Field<      Real>& f);

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

#ifdef SCREAM_DEBUG
  using field_type = Field<Real>;

  const FieldRepository<Real>*   m_field_repo;
  const FieldRepository<Real>*   m_bkp_field_repo;
#endif
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP
