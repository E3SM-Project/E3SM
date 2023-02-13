#ifndef SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP
#define SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP

#include "share/atm_process/atmosphere_process.hpp"
#include "share/property_checks/mass_and_energy_column_conservation_check.hpp"
#include "control/surface_coupling_utils.hpp"

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
 *  The only caveat is required fields in sequential scheduling: if an atm proc
 *  requires a field that is computed by a previous atm proc in the group,
 *  that field is not exposed as a required field of the group.
 */

class AtmosphereProcessGroup : public AtmosphereProcess
{
public:
  using atm_proc_type     = AtmosphereProcess;

  // Constructor(s)
  AtmosphereProcessGroup (const ekat::Comm& comm, const ekat::ParameterList& params);

  virtual ~AtmosphereProcessGroup () = default;

  // The type of the block (e.g., dynamics or physics)
  AtmosphereProcessType type () const { return AtmosphereProcessType::Group; }

  // The name of the block
  std::string name () const { return m_group_name; }

  // Grab the proper grid from the grids manager
  void set_grids (const std::shared_ptr<const GridsManager> grids_manager);

  // --- Methods specific to AtmosphereProcessGroup --- //
  int get_num_processes () const { return m_atm_processes.size(); }

  std::shared_ptr<const atm_proc_type> get_process (const int i) const {
    return m_atm_processes.at(i);
  }

  std::shared_ptr<atm_proc_type> get_process_nonconst (const int i) const {
    return m_atm_processes.at(i);
  }

  ScheduleType get_schedule_type () const { return m_group_schedule_type; }

  // Computes total number of bytes needed for local variables
  size_t requested_buffer_size_in_bytes () const;

  // Set local variables using memory provided by
  // the ATMBufferManager
  void init_buffers(const ATMBufferManager& buffer_manager);

  // The APG class needs to perform special checks before establishing whether
  // a required group/field is indeed a required group for this APG
  void set_required_field (const Field& field);
  void set_required_group (const FieldGroup& group);

  // Gather internal fields from all processes in the group
  // NOTE: this method *must* be called before any attempt to query this atm proc group
  //       for its internal fields, otherwise it will appear as if this atm proc group
  //       stores ZERO internal fields. In other words, this method populates the list
  //       of internal fields of the group.
  void gather_internal_fields ();

  // Returns true if any internal processes enables
  // the mass and energy conservation checks.
  bool are_column_conservation_checks_enabled () const;

  // Adds the mass and energy conservation
  // checks to appropriate physics processes.
  void setup_column_conservation_checks (
      const std::shared_ptr<MassAndEnergyColumnConservationCheck>& conservation_check,
      const CheckFailHandling                                      fail_handling_type) const;

  // Add nan checks after each non-group process, for each computed field.
  // If checks fail, we print all input and output fields of that process
  // (that are on the same grid) at the location of the fail.
  void add_postcondition_nan_checks () const;

protected:

  // Adds fid to the list of required/computed fields of the group (as a whole).
  void process_required_field (const FieldRequest& req);
  void process_required_group (const GroupRequest& req);

  // The initialization, run, and finalization methods
  void initialize_impl(const RunType run_type);
  void initialize_impl ();
  void run_impl        (const double dt);
  void finalize_impl   (/* what inputs? */);

  void run_sequential (const double dt);
  void run_parallel   (const double dt);

  // The methods to set the fields/groups in the right processes of the group
  void set_required_field_impl (const Field& f);
  void set_computed_field_impl (const Field& f);
  void set_required_group_impl (const FieldGroup& group);
  void set_computed_group_impl (const FieldGroup& group);

  // The name of the group. This is usually a concatenation of the names of the individual processes
  std::string       m_group_name;
  int               m_group_size;

  // The list of atm processes in this group
  std::vector<std::shared_ptr<atm_proc_type>>  m_atm_processes;

  // The schedule type: Parallel vs Sequential
  ScheduleType   m_group_schedule_type;

  // This is only needed to be able to access grids objects later on
  std::shared_ptr<const GridsManager>   m_grids_mgr;
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_GROUP_HPP
