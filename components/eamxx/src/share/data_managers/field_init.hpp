#ifndef EAMXX_FIELD_INIT_HPP
#define EAMXX_FIELD_INIT_HPP

#include "share/field/field.hpp"

#include <map>
#include <string>
#include <vector>

namespace scream
{

// A recipe to initialize a field from other fields, for the initial conditions.
//
// An atm process that "owns" a field (e.g., a tracer that it introduces, and that is not
// in the IC file) can provide a FieldInit for it (see AtmosphereProcess::add_field_init).
// ModelInit uses it for startup runs, after the IC file has been read, and only for fields
// that are not overridden, not in the IC file, and without a default constant. In particular,
// this happens before any atm process is initialized, so it can init fields that the dynamics
// needs (e.g., tracers) as if they were read from the IC file.
//
// A FieldInit is a self-contained recipe: it must not hold on to the atm process that creates
// it, nor to any of its fields. Anything it needs must be copied in at construction. The reason
// is that the field to init (and its inputs) do not necessarily live on the physics grid,
// since the ICs may be on a different grid than the one the process uses (e.g., for pg2 runs).
// Therefore, a FieldInit should be column-local: the value of a column of the output
// must only depend on the same column of the inputs.
class FieldInit {
public:
  FieldInit (const std::string& name,
             const std::vector<std::string>& inputs = {})
   : m_name (name), m_inputs (inputs)
  { /* Nothing to do here */ }

  virtual ~FieldInit () = default;

  // The name of the field this recipe inits
  const std::string& name () const { return m_name; }

  // The names of the fields needed to init it. They are inited before the call to init_field.
  const std::vector<std::string>& inputs () const { return m_inputs; }

  // Init the (allocated) field y, given the inputs xs (by name). The inputs are on the same grid as y.
  virtual void init_field (Field& y, const std::map<std::string,Field>& xs) const = 0;

private:
  std::string               m_name;
  std::vector<std::string>  m_inputs;
};

} // namespace scream

#endif // EAMXX_FIELD_INIT_HPP
