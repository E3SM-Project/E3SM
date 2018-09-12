#ifndef SCREAM_ATMOSPHERE_PROCESS_HPP
#define SCREAM_ATMOSPHERE_PROCESS_HPP

#include <string>
#include <list>

namespace scream
{

enum class AtmosphereProcessType {
  Dynamics,
  Physics
};

/*
 *  The abstract interface of a process of the atmosphere 
 *
 *  The process will handle a particular part of the atmosphere component.
 *  This includes both physics (i.e., parametrizations) and dynamics.
 *  The atmosphere driver will take care of calling init/run/finalize
 *  methods of each process, in an order that the driver
 *  establishes. A process must provide a list of fields
 *  that it needs as input, together with a list of fields that
 *  are computed.
 */

class AtmosphereProcess
{
public:
  // The type of the block (dynamics or physics)
  virtual AtmosphereProcessType type () const = 0;

  // The name of the block 
  virtual std::string name () const = 0;

  // These are the three main interfaces:
  //   - the initialize method sets up all the stuff the process needs to run,
  //     including arrays/views, parameters, and precomputed stuff.
  //   - the run method time-advances the process by one time step.
  //     We could decide whether we want to assume that other process may
  //     be running at the same time, or whether each process can assume
  //     that no other process is currently inside a call to 'run'.
  //   - the finalize method makes sure, if necessary, that all resources are freed.
  // The initialize/finalize method should be called just once per simulation (should
  // we enforce that? It depends on what resources we init/free, and how), while the
  // run method can (and usually will) be called multiple times.
  // We should put asserts to verify that the process has been init-ed, when
  // run/finalize is called.
  virtual void initialize (/* what inputs? */) = 0;
  virtual void run        (/* what inputs? */) = 0;
  virtual void finalize   (/* what inputs? */) = 0;

  // These two methods allow the driver to figure out what process need
  // a given field and what process updates a given field.
  virtual const std::list<std::string>&  get_required_fields () const = 0;
  virtual const std::list<std::string>&  get_computed_fields () const = 0;
};

} // namespace scream

#endif // SCREAM_ATMOSPHERE_PROCESS_HPP
