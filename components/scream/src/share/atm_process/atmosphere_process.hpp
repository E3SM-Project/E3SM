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

#include "ekat/ekat_assert.hpp"
#include "ekat/mpi/ekat_comm.hpp"
#include "ekat/ekat_parameter_list.hpp"
#include "ekat/util/ekat_factory.hpp"
#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/std_meta/ekat_std_enable_shared_from_this.hpp"
#include "ekat/std_meta/ekat_std_utils.hpp"

#include <string>
#include <set>
#include <type_traits>

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
 *   - Most AP's will require a single grid. However, the special concrete
 *     class AtmosphereProcessGroup can store AP's running on different grids.
 *   - When a derived class implements the set methods for required/computed
 *     fields, it should make sure to set itself as customer/provider of
 *     the field. The methods add_me_as_provider/customer can be used on the
 *     input field to achieve this result.
 *   - An AP can claim to require or update a group of fields. This can be useful
 *     if the AP performs the *same* action on a bunch of fields, with no
 *     knowledge of what fields are (e.g., advect them, or apply fix/limiter).
 *     For each group, the AP must specify the group name and the grid name.
 *     If the same group is needed on multiple grids, the AP will issue a separate
 *     request for each grid. An AP that needs this feature must override the
 *     get_X_groups and the set_X_groups (with X=required and/or updated).
 *     Notice that it makes no sense to have computed_groups: if an AP computes
 *     a group of fields that are not an input (i.e., not updated fields), then
 *     it must know the names and layouts of the fields, which means it can handle
 *     all these fields directly in [set,get]_computed_fields.
 */

class AtmosphereProcess : public ekat::enable_shared_from_this<AtmosphereProcess>
{
public:
  using TimeStamp      = util::TimeStamp;
  using ci_string      = ekat::CaseInsensitiveString;

  virtual ~AtmosphereProcess () = default;

  // The type of the process (e.g., dynamics or physics)
  virtual AtmosphereProcessType type () const = 0;

  // The type of grids needed by the process
  virtual std::set<std::string> get_required_grids () const = 0;

  // The name of the process
  virtual std::string name () const = 0;

  // The communicator associated with this atm process
  virtual const ekat::Comm& get_comm () const = 0;

  // Give the grids manager to the process, so it can grab its grid
  // IMPORTANT: the process is *expected* to have valid field ids for
  //            required/computed fields once this method returns.
  //            This means that all tags/dims must be set upon return.
  virtual void set_grids (const std::shared_ptr<const GridsManager> grids_manager) = 0;

  // These are the three main interfaces:
  //   - the initialize method sets up all the stuff the process needs to run,
  //     including arrays/views, parameters, and precomputed stuff.
  //   - the run method time-advances the process by one time step.
  //     We could decide whether we want to assume that other process may
  //     be running at the same time, or whether each process can assume
  //     that no other process is currently inside a call to 'run'.
  //   - the finalize method makes sure, if necessary, that all resources are freed.
  // NOTE: You don't override these methods. Override the protected methods
  // NOTE: initialize_impl, run_impl, and finalize_impl instead.
  // The initialize/finalize method should be called just once per simulation (should
  // we enforce that? It depends on what resources we init/free, and how), while the
  // run method can (and usually will) be called multiple times.
  // We should put asserts to verify that the process has been init-ed, when
  // run/finalize is called.
  void initialize (const TimeStamp& t0) {
    t_ = t0;
    initialize_impl(t_);
  }
  void run        (const Real dt) {
    // Call the subclass's run method and update it afterward.
    run_impl(dt);
    t_ += dt;
  }
  void finalize   (/* what inputs? */) {
    finalize_impl(/* what inputs? */);
  }

  // These methods set fields in the atm process. Fields live on the default
  // device and they are all 1d.
  // If the process *needs* to store the field as n-dimensional field, use the
  // template function 'get_reshaped_view' (see field.hpp for details).
  // Note: this method will be called *after* set_grids, but *before* initialize.
  //       You are *guaranteed* that the view in the Field f is allocated by now.
  // Note: it would be tempting to add this process as provider/customer right here.
  //       However, if this process is of type Group, we don't really want to add it
  //       as provider/customer. The group is just a 'design layer', and the stored
  //       processes are the actuall providers/customers.
  void set_required_field (const Field<const Real>& f) {
    EKAT_REQUIRE_MSG (requires_field(f.get_header().get_identifier()),
        "Error! This atmosphere process does not require\n  " +
        f.get_header().get_identifier().get_id_string() +
        "\nSomething is wrong up the call stack. Please, contact developers.\n");
    set_required_field_impl (f);
  }
  void set_computed_field (const Field<Real>& f) {
    EKAT_REQUIRE_MSG (computes_field(f.get_header().get_identifier()),
        "Error! This atmosphere process does not compute\n  " +
        f.get_header().get_identifier().get_id_string() +
        "\nSomething is wrong up the call stack. Please, contact developers.\n");
    set_computed_field_impl (f);
  }

  // Note: for the following (unlike set_required/computed_field, we do provide an
  //       implementation, since requiring a group is "rare".
  // Note: from the group, derived class can extract individual fields, and,
  //       if needed, they can check that the group was allocated as a bundle,
  //       and if so, and if desired, they can access the bundled field directly.
  //       See field_group.hpp for more details.
  virtual void set_required_group (const FieldGroup<const Real>& /* group */) {
    EKAT_ERROR_MSG (
      "Error! This atmosphere process does not require a group of fields, meaning\n"
      "       that 'get_required_groups' was not overridden in this class, or that\n"
      "       its override returns an empty set.\n"
      "       If you override 'get_required_groups' to return a non-empty set,\n"
      "       then you must also override 'set_required_group' in your derived class.\n"
    );
  }
  virtual void set_updated_group (const FieldGroup<Real>& /* group */) {
    EKAT_ERROR_MSG (
      "Error! This atmosphere process does not update a group of fields, meaning\n"
      "       that 'get_updated_groups' was not overridden in this class, or that\n"
      "       its override returns an empty set.\n"
      "       If you override 'get_updated_groups' to return a non-empty set,\n"
      "       then you must also override 'set_updated_group' in your derived class.\n"
    );
  }

  // These two methods allow the driver to figure out what process need
  // a given field and what process updates a given field.
  const std::set<FieldRequest>& get_required_fields () const { return m_required_fields; }
  const std::set<FieldRequest>& get_computed_fields () const { return m_computed_fields; }

  // If needed, an Atm Proc can claim to need/update a whole group of fields, without really knowing
  // a priori how many they are, or even what they are. Each entry of the returned set is a pair
  // of strings, where the 1st is the group name, and the 2nd the grid name. If the same group is
  // needed on multiple grids, two separate entries are needed.
  const std::set<GroupRequest>& get_required_groups () const { return m_required_groups; }
  const std::set<GroupRequest>& get_updated_groups  () const { return m_updated_groups; }

  // NOTE: C++20 will introduce the method 'contains' for std::set. Till then, use our util free function
  bool requires_field (const FieldIdentifier& id) const {
    for (const auto& it : m_required_fields) {
      if (it.fid==id) {
        return true;
      }
    }
    return false;
  }
  bool computes_field (const FieldIdentifier& id) const {
    for (const auto& it : m_computed_fields) {
      if (it.fid==id) {
        return true;
      }
    }
    return false;
  }

  bool requires_group (const std::string& name, const std::string& grid) const {
    for (const auto& it : m_required_groups) {
      if (it.name==name && it.grid==grid) {
        return true;
      }
    }
    return false;
  }
  bool updates_group (const std::string& name, const std::string& grid) const {
    for (const auto& it : m_updated_groups) {
      if (it.name==name && it.grid==grid) {
        return true;
      }
    }
    return false;
  }

  // Computes total number of bytes needed for local variables
  virtual int requested_buffer_size_in_bytes () const { return 0; }

  // Set local variables using memory provided by
  // the ATMBufferManager
  virtual void init_buffers(const ATMBufferManager& /*buffer_manager*/) {}

protected:

  enum RequestType {
    Required,
    Computed,
    Updated
  };

  // Derived classes can used these method, so that if we change how fields/groups
  // requirement are stored (e.g., std::set -> std::list), they don't need to change
  // their implementation.

  // Field requests. Provide plenty of overloads to make it simple to add field request.
  // E.g., provide a FID directly vs provide its ctor args; or provide groups list and
  // no pack size vs single group and pack size.
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
      auto& fields = RT==Required ? m_required_fields : m_computed_fields;
      fields.emplace(req);
    }
  }

  // Group requests
  template<RequestType RT>
  void add_group (const std::string& name, const std::string& grid, const int ps, const Bundling b,
                  const GroupRequest* p, const Relationship t, const std::list<std::string>& excl)
  { add_group<RT>(GroupRequest(name,grid,ps,b,p,t,excl)); }

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
    static_assert(RT==Required || RT==Updated,
        "Error! Invalid request type in call to add_group.\n");
    auto& groups = RT==Required ? m_required_groups : m_updated_groups;
    groups.emplace(req);
  }

  // Override this method to initialize your subclass.
  virtual void initialize_impl(const TimeStamp& t0) = 0;

  // Override this method to define how your subclass runs forward one step
  // (of size dt). This method is called before the timestamp is updated.
  virtual void run_impl(const Real dt) = 0;

  // Override this method to finalize your subclass.
  virtual void finalize_impl(/* what inputs? */) = 0;

  // This provides access to this process's timestamp.
  const TimeStamp& timestamp() const {
    return t_;
  }

  void add_me_as_provider (const Field<Real>& f) {
    f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
  }

  void add_me_as_customer (const Field<const Real>& f) {
    f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
  }

  virtual void set_required_field_impl (const Field<const Real>& f) = 0;
  virtual void set_computed_field_impl (const Field<      Real>& f) = 0;

private:

  std::set<FieldRequest>   m_required_fields;
  std::set<FieldRequest>   m_computed_fields;

  std::set<GroupRequest>   m_required_groups;
  std::set<GroupRequest>   m_updated_groups;

  // This process's copy of the timestamp, which is set on initialization and
  // updated during stepping.
  TimeStamp t_;
};

// A short name for the factory for atmosphere processes
// WARNING: you do not need to write your own creator function to register your atmosphere process in the factory.
//          You could, but there's no need. You can simply register the common one right below, using your
//          atmosphere process class name as templated argument. If you roll your own creator function, you
//          *MUST* ensure that it correctly sets up the self pointer after creating the shared_ptr.
//          This is *necessary* until we can safely switch to std::enable_shared_from_this.
//          For more details, see the comments at the top of util/scream_std_enable_shared_from_this.hpp.
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
