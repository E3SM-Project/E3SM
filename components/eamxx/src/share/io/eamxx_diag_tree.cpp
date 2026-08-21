#include "share/io/eamxx_diag_tree.hpp"

#include "share/io/eamxx_legacy_diag_names.hpp"

#include <ekat_assert.hpp>

#include <algorithm>

namespace scream {

namespace {

Field build (const DiagSpec& s,
             const std::shared_ptr<const AbstractGrid>& grid,
             const FieldManager& fm,
             std::map<std::string,diag_ptr_t>& repo,
             const std::map<std::string,Field>& known,
             std::vector<diag_ptr_t>& order)
{
  const auto& cname = s.names.canonical;

  if (s.is_leaf()) {
    // A name the caller has already built and registered. Checked first, so
    // that an expression can refer to another expression by name.
    if (auto it=known.find(cname); it!=known.end()) {
      return it->second;
    }
    // An unresolved leaf that the model provides is simply a field. A model
    // field wins over a diagnostic of the same name, as it always has.
    if (fm.has_field(cname)) {
      return fm.get_field(cname);
    }
    // Not a field, and not a diag registered under this very name: it may be a
    // legacy mangled name. Translating it here, rather than up front, is what
    // preserves the rule above at every level of the expression -- a model
    // field called 'foo_prev' was already returned above, and never gets here.
    if (not DiagnosticFactory::instance().has_product(cname)) {
      const auto expr = legacy_to_expr(cname);
      if (expr!=cname) {
        const auto sub = lower_to_diag_spec(expr,cname);
        // Alias back to the legacy name: that is what the parent's params, and
        // hence the consuming diag's input list, refer to it as.
        return build(sub,grid,fm,repo,known,order).alias(cname);
      }
    }
  }

  // A diag can be shared with an earlier tree, or show up twice in this one.
  // Either way it must appear exactly once in *this* tree's evaluation order.
  auto add_to_order = [&](const diag_ptr_t& d) {
    if (std::find(order.begin(),order.end(),d)==order.end()) {
      order.push_back(d);
    }
  };

  // Already built for this very sub-expression? Reuse it.
  if (auto it=repo.find(cname); it!=repo.end()) {
    add_to_order(it->second);
    return it->second->get().alias(cname);
  }

  // Operands first: this is what puts 'order' in dependency order.
  std::map<std::string,Field> operands;
  for (const auto& c : s.children) {
    operands[c.names.canonical] = build(c,grid,fm,repo,known,order);
  }

  const std::string key = s.is_leaf() ? cname : s.factory_key;
  auto& factory = DiagnosticFactory::instance();
  EKAT_REQUIRE_MSG (factory.has_product(key),
      "Error! Unrecognized diagnostic.\n"
      " - name      : " + key + "\n"
      " - expression: " + cname + "\n"
      " - note: this is neither a field computed by the model, nor a diagnostic\n"
      "         registered in the diagnostic factory.\n");

  // Lowering is grid-independent, so the grid is injected here. Diags that do
  // not care about it simply never read the parameter.
  auto params = s.params;
  params.set<std::string>("grid_name",grid->name());
  auto diag = factory.create(key,grid->get_comm(),params,grid);

  // Wire the inputs. Most of them are the operands we just built, but a diag may
  // also need inputs that never appear in the expression. Those are often model
  // fields (a dp-weighted vertical average needs 'pseudo_density'), but they can
  // be diagnostics in their own right (FieldAtHeight needs 'height_mid'), so
  // resolve them exactly as a leaf of the expression would be resolved.
  // NOTE: this recursion happens before the diag is added to 'order', so an
  //       implicit input still ends up ahead of the diag that consumes it.
  for (const auto& in_name : diag->get_input_fields_names()) {
    auto it = operands.find(in_name);
    if (it!=operands.end()) {
      diag->set_input_field(it->second);
    } else {
      const auto dep = build(lower_to_diag_spec(in_name,in_name),grid,fm,repo,known,order);
      diag->set_input_field(dep.alias(in_name));
    }
  }

  diag->initialize();
  repo[cname] = diag;
  add_to_order(diag);

  // The diag names its own output field, following the legacy convention. The
  // tree refers to sub-expressions by canonical name instead, so alias it.
  // Aliasing shares the data, the tracking, and the extra data, so a mask or a
  // timestamp set on the diag field is still visible through the alias.
  return diag->get().alias(cname);
}

} // anonymous namespace

DiagTree build_diag_tree (const DiagSpec& spec,
                          const std::shared_ptr<const AbstractGrid>& grid,
                          const FieldManager& fm,
                          std::map<std::string,diag_ptr_t>& repo,
                          const std::map<std::string,Field>& known)
{
  DiagTree tree;
  auto f = build(spec,grid,fm,repo,known,tree.eval_order);

  // Unlike the canonical name, the registered one is what the FieldManager and
  // the output file will know this by.
  tree.top = spec.names.registered.empty()
           ? f : f.alias(spec.names.registered);
  return tree;
}

} // namespace scream
