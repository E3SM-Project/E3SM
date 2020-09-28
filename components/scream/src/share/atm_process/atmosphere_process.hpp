#ifndef SCREAM_ATMOSPHERE_PROCESS_HPP
#define SCREAM_ATMOSPHERE_PROCESS_HPP

#include "share/atm_process/atmosphere_process_utils.hpp"
#include "share/field/field_identifier.hpp"
#include "share/field/field_repository.hpp"
#include "share/field/field.hpp"
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
 */

class AtmosphereProcess : public ekat::enable_shared_from_this<AtmosphereProcess>
{
public:
  using device_type = DefaultDevice; // may need to template class on this
  using TimeStamp = util::TimeStamp;

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

  // These methods set fields in the atm process. Fields live on device and they are all 1d.
  // If the process *needs* to store the field as n-dimensional field, use the
  // template function 'get_reshaped_view' (see field.hpp for details).
  // Note: this method will be called *after* set_grid, but *before* initialize.
  //       You are *guaranteed* that the view in the Field f is allocated by now.
  // Note: it would be tempting to add this process as provider/customer right here.
  //       However, if this process is of type Group, we don't really want to add it
  //       as provider/customer. The group is just a 'design layer', and the stored processes
  //       are the actuall providers/customers.
  void set_required_field (const Field<const Real, device_type>& f) {
    ekat::error::runtime_check(requires(f.get_header().get_identifier()),
                         "Error! This atmosphere process does not require this field. "
                         "Something is wrong up the call stack. Please, contact developers.\n");
    set_required_field_impl (f);
  }
  void set_computed_field (const Field<Real, device_type>& f) {
    ekat::error::runtime_check(computes(f.get_header().get_identifier()),
                         "Error! This atmosphere process does not compute this field. "
                         "Something is wrong up the call stack. Please, contact developers.\n");
    set_computed_field_impl (f);
  }

  // Register required/computed fields in the field repo
  virtual void register_fields (FieldRepository<Real, device_type>& field_repo) const = 0;

  // These two methods allow the driver to figure out what process need
  // a given field and what process updates a given field.
  virtual const std::set<FieldIdentifier>& get_required_fields () const = 0;
  virtual const std::set<FieldIdentifier>& get_computed_fields () const = 0;

  // NOTE: C++20 will introduce the method 'contains' for std::set. Till then, use our util free function
  bool requires (const FieldIdentifier& id) const { return ekat::contains(get_required_fields(),id); }
  bool computes (const FieldIdentifier& id) const { return ekat::contains(get_computed_fields(),id); }

protected:

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

  void add_me_as_provider (const Field<Real, device_type>& f) {
    f.get_header_ptr()->get_tracking().add_provider(weak_from_this());
  }

  void add_me_as_customer (const Field<const Real, device_type>& f) {
    f.get_header_ptr()->get_tracking().add_customer(weak_from_this());
  }

  virtual void set_required_field_impl (const Field<const Real, device_type>& f) = 0;
  virtual void set_computed_field_impl (const Field<      Real, device_type>& f) = 0;

private:

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
