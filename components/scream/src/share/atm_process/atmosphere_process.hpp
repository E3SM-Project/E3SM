#ifndef SCREAM_ATMOSPHERE_PROCESS_HPP
#define SCREAM_ATMOSPHERE_PROCESS_HPP

#include "share/atm_process/atmosphere_process_utils.hpp"
#include "share/atm_process/ATMBufferManager.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_manager.hpp"
#include "share/field/field_request.hpp"
#include "share/field/field.hpp"
#include "share/field/field_group.hpp"
#include "share/grid/grids_manager.hpp"

#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_factory.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"

#include <string>
#include <set>
#include <list>

namespace scream
{

/*
 *  The abstract interface of a process of the atmosphere (AP)
 *
 *  The process will handle a particular part of the atmosphere component.
 *  This includes physics (i.e., parametrizations), dynamics, as well
 *  as the surface coupling.
 *  The atmosphere driver (AD) will take care of calling init/run/finalize
 *  methods of each process, in an order that the AD establishes
 *  based on the user input options.
 *  A concrete process must implement all the purely virtual
 *  methods of this interface; for instance, it must provide a list of fields
 *  that it needs as input, together with a list of fields that are computed.
 *
 *  Notes to developers:
 *   - If an AP 'updates' a field (i.e., f = f + delta), then it should make
 *     sure that f is listed both as required and computed field. This helps
 *     the AD to check that all the AP's dependencies are met.
 *   - An AP can claim to require or compute a group of fields. This can be useful
 *     if the AP performs the *same* action on a bunch of fields, with no
 *     knowledge of what fields are (e.g., advect them, or apply fix/limiter).
 *   - Fields and groups must be requested via FieldRequest and GroupRequest
 *     respectively (see field_request.hpp). To add a request, use the methods
 *     add_field<RT>(..) and add_group<RT>(..), with RT=Required, Updated, or
 *     Computed (Updated = Required + Computed).
 *   - If the same group is needed on multiple grids, the AP will issue a separate
 *     request for each grid.
 *   - Notice that it is unlikely that an AP computes a group, without requiring
 *     it as an input (it should probably know what's in the group that it computes).
 *     Nevertheless, to simplify the code (treat fields and groups similarly),
 *     we expose required and computed groups, just like fields.
 */

class AtmosphereProcess : public ekat::enable_shared_from_this<AtmosphereProcess>
{
public:
  using TimeStamp      = util::TimeStamp;
  using ci_string      = ekat::CaseInsensitiveString;
  using       field_type = Field<      Real>;
  using const_field_type = Field<const Real>;
  using       group_type = FieldGroup<      Real>;
  using const_group_type = FieldGroup<const Real>;

  template<typename T>
  using str_map = std::map<std::string,T>;

  // Base constructor to set MPI communicator and params
  AtmosphereProcess (const ekat::Comm& comm, const ekat::ParameterList& params);

  virtual ~AtmosphereProcess () = default;

  // The type of the process (e.g., dynamics or physics)
  virtual AtmosphereProcessType type () const = 0;

  // The type of grids needed by the process
  virtual std::set<std::string> get_required_grids () const = 0;

  // The name of the process
  virtual std::string name () const = 0;

  // Give the grids manager to the process, so it can grab its grid
  // Upon return, the atm proc should have a valid and complete list
  // of in/out/inout FieldRequest and GroupRequest.
  virtual void set_grids (const std::shared_ptr<const GridsManager> grids_manager) = 0;

  // These are the three main interfaces:
  //   - the initialize method sets up all the stuff the process needs in order to run,
  //     including arrays/views, parameters, and precomputed stuff.
  //   - the run method time-advances the process by one time step, provided in seconds.
  //   - the finalize method makes sure, if necessary, that all resources are freed.
  // TODO: should we check that initialize/finalize are called once per simulation?
  //       Whether it's a needed check depends on what resources we init/free, and how.
  // TODO: should we check that initialize has been called, when calling run/finalize?
  void initialize (const TimeStamp& t0);
  void run (const int dt);
  void finalize   (/* what inputs? */);

  // Return the MPI communicator
  const ekat::Comm& get_comm () const { return m_comm; }

  // Return the parameter list
  const ekat::ParameterList& get_params () const { return m_params; }

  // These methods set fields/groups in the atm process. The fields/groups are stored
  // in a list (with some helpers maps that can be used to quickly retrieve them).
  // If derived class need additional bookkeping/checks, they can override the
  // corresponding set_xyz_impl method(s).
  // Note: this method will be called *after* set_grids, but *before* initialize.
  //       You are *guaranteed* that the views in the field/group are allocated by now.
  // Note: you are *unlikely* to need to override these methods. In 99.99% of the cases,
  //       overriding the corresponding _impl method should be enough. The class
  //       AtmosphereProcessGroup is the big exception to this, since it needs
  //       to perform some extra action *before* setting the field/group.
  virtual void set_required_field (const Field<const Real>& f);
  virtual void set_computed_field (const Field<Real>& f);
  virtual void set_required_group (const FieldGroup<const Real>& group);
  virtual void set_computed_group (const FieldGroup<Real>& group);

  // These methods checks that all the in/out fields of this atm process are valid.
  // For each field, these routines run all the property checks stored in the field.
  // Derived classes can perform additional checks (or repairs) by overriding the
  // corresponding check_xyz_fields_impl method(s).
  // Note: We don't want inputs to be 'repaired', since they might break assumptions
  //       in other atm procs (e.g., some atm procs might have a "backup" copy).
  //       However, we allow computed fields to be repaired, so check_computed_fields
  //       is not a const method.  If a process wants to repair computed fields
  //       it can do so by overriding the "check_computed_fields_impl" routine.
  void check_required_fields () const;
  void check_computed_fields ();

  // These methods allow the AD to figure out what each process needs, with very fine
  // grain detail. See field_request.hpp for more info on what FieldRequest and GroupRequest
  // are, and field_group.hpp for what groups of fields are.
  const std::set<FieldRequest>& get_required_field_requests () const { return m_required_field_requests; }
  const std::set<FieldRequest>& get_computed_field_requests () const { return m_computed_field_requests; }
  const std::set<GroupRequest>& get_required_group_requests () const { return m_required_group_requests; }
  const std::set<GroupRequest>& get_computed_group_requests () const { return m_computed_group_requests; }

  // These sets allow to get all the actual in/out fields stored by the atm proc
  // Note: if an atm proc requires a group, then all the fields in the group, as well as
  //       the bundled field (if present) will be added as required fields for this atm proc.
  //       See field_group.hpp for more info about groups of fields.
  const std::list<const_field_type>& get_fields_in  () const { return m_fields_in;  }
  const std::list<      field_type>& get_fields_out () const { return m_fields_out; }
  const std::list<const_group_type>& get_groups_in  () const { return m_groups_in;  }
  const std::list<      group_type>& get_groups_out () const { return m_groups_out; }

  // Whether this atm proc requested the field/group as in/out/inout, via a FieldRequest/GroupRequest.
  bool requires_field (const FieldIdentifier& id) const;
  bool computes_field (const FieldIdentifier& id) const;
  bool requires_group (const std::string& name, const std::string& grid) const;
  bool computes_group (const std::string& name, const std::string& grid) const;

  // Computes total number of bytes needed for local variables
  virtual int requested_buffer_size_in_bytes () const { return 0; }

  // Set local variables using memory provided by
  // the ATMBufferManager
  virtual void init_buffers(const ATMBufferManager& /*buffer_manager*/) {}

protected:

  enum RequestType {
    Required,
    Computed,
    Updated   // For convenience, triggers Required+Computed
  };

  // Derived classes can used these method, so that if we change how fields/groups
  // requirement are stored (e.g., change the std container), they don't need to change
  // their implementation.
  // We provide plenty of overloads to make it simple to add field request.
  // E.g., provide a FID directly vs provide its ctor args; or provide groups list and
  // no pack size vs single group and pack size.

  // Field requests
  template<RequestType RT>
  void add_field (const std::string& name, const FieldLayout& layout,
                  const ekat::units::Units& u, const std::string& grid_name,
                  const int ps = 1)
  { add_field<RT>(FieldIdentifier(name,layout,u,grid_name),ps); }

  template<RequestType RT>
  void add_field (const std::string& name, const FieldLayout& layout,
                  const ekat::units::Units& u, const std::string& grid_name,
                  const std::string& group, const int ps = 1)
  { add_field<RT>(FieldIdentifier(name,layout,u,grid_name),group,ps); }

  template<RequestType RT>
  void add_field (const std::string& name, const FieldLayout& layout,
                  const ekat::units::Units& u, const std::string& grid_name,
                  const std::list<std::string>& groups, const int ps = 1)
  { add_field<RT>(FieldIdentifier(name,layout,u,grid_name),groups,ps); }

  template<RequestType RT>
  void add_field (const FieldIdentifier& fid, const std::string& group, const int ps = 1)
  { add_field<RT>(FieldRequest(fid,group,ps)); }

  template<RequestType RT>
  void add_field (const FieldIdentifier& fid, const std::list<std::string>& groups)
  { add_field<RT>(FieldRequest(fid,groups)); }

  template<RequestType RT>
  void add_field (const FieldIdentifier& fid, const int ps = 1)
  { add_field<RT>(FieldRequest(fid,ps)); }

  template<RequestType RT>
  void add_field (const FieldIdentifier& fid, const std::list<std::string>& groups, const int ps)
  { add_field<RT>(FieldRequest(fid,groups,ps)); }

  template<RequestType RT>
  void add_field (const FieldRequest& req)
  {
    // Since we use C-style enum, let's avoid invalid integers casts
    static_assert(RT==Required || RT==Computed || RT==Updated,
                  "Error! Invalid request type in call to add_field.\n");

    if (RT==Updated) {
      add_field<Required>(req);
      add_field<Computed>(req);
    } else {
      auto& fields = RT==Required ? m_required_field_requests : m_computed_field_requests;
      fields.emplace(req);
    }
  }

  // Group requests
  template<RequestType RT>
  void add_group (const std::string& name, const std::string& grid, const int ps, const Bundling b,
                  const DerivationType t, const std::string& src_name, const std::string& src_grid,
                  const std::list<std::string>& excl)
  { add_group<RT>(GroupRequest(name,grid,ps,b,t,src_name,src_grid,excl)); }

  template<RequestType RT>
  void add_group (const std::string& name, const std::string& grid_name,
                  const Bundling b = Bundling::NotNeeded)
  { add_group<RT> (GroupRequest(name,grid_name,b)); }

  template<RequestType RT>
  void add_group (const std::string& name, const std::string& grid_name,
                  const int pack_size, const Bundling b = Bundling::NotNeeded)
  { add_group<RT> (GroupRequest(name,grid_name,pack_size,b)); }

  template<RequestType RT>
  void add_group (const GroupRequest& req)
  {
    // Since we use C-style enum, let's avoid invalid integers casts
    static_assert(RT==Required || RT==Updated || RT==Computed,
        "Error! Invalid request type in call to add_group.\n");
    if (RT==Updated) {
      add_group<Required>(req);
      add_group<Computed>(req);
    } else {
      auto& groups = RT==Required ? m_required_group_requests : m_computed_group_requests;
      groups.emplace(req);
    }
  }

  // Override this method to initialize the derived
  virtual void initialize_impl() = 0;

  // Override this method to define how the derived runs forward one step
  // (of size dt). This method is called before the timestamp is updated.
  virtual void run_impl(const int dt) = 0;

  // Override this method to finalize the derived class
  virtual void finalize_impl(/* what inputs? */) = 0;

  // This provides access to this process's timestamp.
  const TimeStamp& timestamp() const { return m_time_stamp; }

  // These three methods modify the FieldTracking of the input field (see field_tracking.hpp)
  void update_time_stamps ();
  void add_me_as_provider (const Field<Real>& f);
  void add_me_as_customer (const Field<const Real>& f);

  // The base class already registers the required/computed/updated fields/groups in
  // the set_required/computed_field and set_required/computed_group routines.
  // These impl methods provide a way for derived classes to add more specialized
  // actions, such as extra fields bookkeeping, extra checks, or create copies.
  // Since most derived classes do not need to perform additional actions,
  // we provide empty implementations.
  virtual void set_required_field_impl (const Field<const Real>& /* f */) {}
  virtual void set_computed_field_impl (const Field<      Real>& /* f */) {}
  virtual void set_required_group_impl (const FieldGroup<const Real>& /* group */) {}
  virtual void set_computed_group_impl (const FieldGroup<      Real>& /* group */) {}

  // The Base class already runs all registered field checks for all fields.
  // Similar to the set_required/computed_field_impl comment above, it is
  // not necessary for a derived class to overwrite these two routines.
  // An example of a case where a derived class may want to override these is
  // if the derived class is able to 'repair' a field in case the check fails.
  // See field_property_check.hpp for more info on check/repair.
  virtual void check_required_fields_impl () const {}
  virtual void check_computed_fields_impl () {}

  // Convenience function to retrieve input/output fields from the field/group (and grid) name.
  // Note: the version without grid name only works if there is only one copy of the field/group.
  //       In that case, the single copy is returned, regardless of the associated grid name.
  Field<const Real>& get_field_in(const std::string& field_name, const std::string& grid_name);
  Field<const Real>& get_field_in(const std::string& field_name);
  Field<Real>& get_field_out(const std::string& field_name, const std::string& grid_name);
  Field<Real>& get_field_out(const std::string& field_name);

  FieldGroup<const Real>& get_group_in(const std::string& group_name, const std::string& grid_name);
  FieldGroup<const Real>& get_group_in(const std::string& group_name);
  FieldGroup<Real>& get_group_out(const std::string& group_name, const std::string& grid_name);
  FieldGroup<Real>& get_group_out(const std::string& group_name);

  // These methods set up an extra pointer in the m_[fields|groups]_[in|out]_pointers,
  // for convenience of use (e.g., use a short name for a field/group).
  // Note: these methods do *not* create a copy of the field/group. Also, notice that
  //       these methods need to be created *after* set_fields_and_groups_pointers().
  void alias_field_in (const std::string& field_name,
                       const std::string& grid_name,
                       const std::string& alias_name);

  void alias_field_out (const std::string& field_name,
                        const std::string& grid_name,
                        const std::string& alias_name);

  void alias_group_in (const std::string& group_name,
                       const std::string& grid_name,
                       const std::string& alias_name);

  void alias_group_out (const std::string& group_name,
                        const std::string& grid_name,
                        const std::string& alias_name);

  // MPI communicator
  ekat::Comm m_comm;

  // Parameter list
  ekat::ParameterList m_params;

private:
  // Called from initialize, this method creates the m_[fields|groups]_[in|out]_pointers
  // maps, which are used inside the get_[field|group]_[in|out] methods.
  void set_fields_and_groups_pointers ();

  // Store input/output fields and groups.
  std::list<const_group_type>  m_groups_in;
  std::list<      group_type>  m_groups_out;
  std::list<const_field_type>  m_fields_in;
  std::list<      field_type>  m_fields_out;

  // These maps help to retrieve a field/group stored in the lists above. E.g.,
  //   auto ptr = m_field_in_pointers[field_name][grid_name];
  // then *ptr is a field in m_fields_in, with name $field_name, on grid $grid_name.
  str_map<str_map< const_group_type* >> m_groups_in_pointers;
  str_map<str_map<       group_type* >> m_groups_out_pointers;
  str_map<str_map< const_field_type* >> m_fields_in_pointers;
  str_map<str_map<       field_type* >> m_fields_out_pointers;

  // The list of in/out/inout field/group requests.
  std::set<FieldRequest>   m_required_field_requests;
  std::set<FieldRequest>   m_computed_field_requests;
  std::set<GroupRequest>   m_required_group_requests;
  std::set<GroupRequest>   m_computed_group_requests;

  // This process's copy of the timestamp, which is set on initialization and
  // updated during stepping.
  TimeStamp m_time_stamp;
};

// A short name for the factory for atmosphere processes
// WARNING: you do not need to write your own creator function to register your atmosphere process in the factory.
//          You could, but there's no need. You can simply register the common one right below, using your
//          atmosphere process class name as templated argument. If you roll your own creator function, you
//          *MUST* ensure that it correctly sets up the self pointer after creating the shared_ptr.
//          This is *necessary* until we can safely switch to std::enable_shared_from_this.
//          For more details, see the comments at the top of ekat_std_enable_shared_from_this.hpp.
using AtmosphereProcessFactory =
    ekat::Factory<AtmosphereProcess,
                  ekat::CaseInsensitiveString,
                  std::shared_ptr<AtmosphereProcess>,
                  const ekat::Comm&,const ekat::ParameterList&>;

// Create an atmosphere process, and correctly set up the (weak) pointer to self.
template <typename AtmProcType>
inline std::shared_ptr<AtmosphereProcess>
create_atmosphere_process (const ekat::Comm& comm, const ekat::ParameterList& p) {
  auto ptr = std::make_shared<AtmProcType>(comm,p);
  ptr->setSelfPointer(ptr);
  return ptr;
}

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_HPP
