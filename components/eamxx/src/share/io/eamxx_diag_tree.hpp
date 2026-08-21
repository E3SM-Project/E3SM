#ifndef SCREAM_EAMXX_DIAG_TREE_HPP
#define SCREAM_EAMXX_DIAG_TREE_HPP

#include "share/io/eamxx_diag_spec.hpp"

#include "share/diagnostics/abstract_diagnostic.hpp"
#include "share/data_managers/field_manager.hpp"
#include "share/field/field.hpp"
#include "share/grid/abstract_grid.hpp"

#include <map>
#include <memory>
#include <string>
#include <vector>

namespace scream {

/*
 * Instantiation of a DiagSpec tree (see eamxx_diag_spec.hpp).
 *
 * Walks the spec tree in post-order: each node's operands are built first, then
 * handed to the node as input fields. Two things fall out of that for free:
 *  - a diagnostic always appears in the evaluation order *after* everything it
 *    depends on, with no need to discover the ordering by repeated scanning;
 *  - a node's inputs are actual Fields, so nothing has to be looked up by name
 *    in the FieldManager, nor re-parsed from a name.
 */

using diag_ptr_t = std::shared_ptr<AbstractDiagnostic>;

struct DiagTree {
  // The root's field, named after the root's 'registered' name.
  Field top;

  // Every diag that must be evaluated, in an order that respects dependencies.
  std::vector<diag_ptr_t> eval_order;
};

/*
 * Build the diags described by 'spec'.
 *
 *  - 'fm' provides the model fields the expression bottoms out in. It is only
 *    read from: registering the resulting field(s) is up to the caller, which
 *    is also the only one that knows about DiagNames::write.
 *  - 'repo' caches diags by canonical name, so that a sub-expression appearing
 *    more than once (in this tree or in an earlier one) is computed once. The
 *    caller owns it, so that diags can be shared across output streams.
 *  - 'known' are fields the caller has already built and given a name to, and
 *    is consulted before 'fm'. This is what lets an expression refer to another
 *    expression by the name it was registered under.
 */
DiagTree build_diag_tree (const DiagSpec& spec,
                          const std::shared_ptr<const AbstractGrid>& grid,
                          const FieldManager& fm,
                          std::map<std::string,diag_ptr_t>& repo,
                          const std::map<std::string,Field>& known = {});

} // namespace scream

#endif // SCREAM_EAMXX_DIAG_TREE_HPP
