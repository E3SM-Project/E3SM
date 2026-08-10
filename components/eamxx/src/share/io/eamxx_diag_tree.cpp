#include "share/io/eamxx_diag_tree.hpp"

#include "share/io/eamxx_io_utils.hpp"
#include "share/io/eamxx_legacy_diag_names.hpp"

#include <ekat_assert.hpp>

#include <algorithm>
#include <stdexcept>

namespace scream {

namespace {

// The factory key to build 'name' with, when 'name' is a bare leaf. Usually the
// name itself, but a diag computing several fields is registered under none of
// the names it produces, so ask the registry for it. Returns 'name' unchanged if
// nothing claims it, so that the caller reports the failure as it always has.
std::string leaf_factory_key (const std::string& name)
{
  if (DiagnosticFactory::instance().has_product(name)) {
    return name;
  }
  const auto provider = DiagOutputsRegistry::instance().provider_of(name);
  return provider.empty() ? name : provider;
}

// The output of 'd' that this sub-expression asked for. Only a diag computing
// several fields is asked by name; anything else names its single output
// however it likes, and 'cname' is merely what we call it here.
Field diag_output (const diag_ptr_t& d, const std::string& cname)
{
  return d->has_output(cname) ? d->get(cname) : d->get();
}

// A diag is built once and shared by everything that asks for it, including
// other output streams. If a later request configures it differently from how
// it was built, the options it gave would simply be ignored, so say so instead.
void check_diag_params_match (AbstractDiagnostic& diag,
                              const std::string& key,
                              const std::string& cname,
                              const ekat::ParameterList& diag_params)
{
  if (not diag_params.isSublist(key)) {
    return;
  }
  const auto& wanted = diag_params.sublist(key);

  // What the diag actually got built with, restricted to the options at hand.
  ekat::ParameterList got (wanted.name());
  for (const auto& name : wanted.param_names()) {
    if (diag.get_params().isParameter(name)) {
      overlay_params(got,diag.get_params(),{name});
    }
  }

  const auto wanted_sig = params_signature(wanted);
  const auto got_sig    = params_signature(got);
  EKAT_REQUIRE_MSG (wanted_sig==got_sig,
      "Error! A diagnostic was already built with different options.\n"
      " - diag      : " + key + "\n"
      " - expression: " + cname + "\n"
      " - requested here : " + wanted_sig + "\n"
      " - built earlier with: " + got_sig + "\n"
      " - note: one diagnostic is shared by every output stream that asks for it,\n"
      "         so all of them must give it the same 'diag_params' options.\n");
}

Field build (const DiagSpec& s,
             const std::shared_ptr<const AbstractGrid>& grid,
             const FieldManager& fm,
             std::map<std::string,diag_ptr_t>& repo,
             const std::map<std::string,Field>& known,
             const ekat::ParameterList& diag_params,
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
    // Not a field, and not a name any diag computes: it may be a legacy mangled
    // name. Translating it here, rather than up front, is what preserves the
    // rule above at every level of the expression -- a model field called
    // 'foo_prev' was already returned above, and never gets here.
    if (not DiagnosticFactory::instance().has_product(cname) and
        DiagOutputsRegistry::instance().provider_of(cname).empty()) {
      const auto expr = legacy_to_expr(cname);
      if (expr!=cname) {
        const auto sub = lower_to_diag_spec(expr,cname);
        // Alias back to the legacy name: that is what the parent's params, and
        // hence the consuming diag's input list, refer to it as.
        return build(sub,grid,fm,repo,known,diag_params,order).alias(cname);
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

  const std::string key = s.is_leaf() ? leaf_factory_key(cname) : s.factory_key;

  // Already built for this very sub-expression? Reuse it.
  if (auto it=repo.find(cname); it!=repo.end()) {
    check_diag_params_match(*it->second,key,cname,diag_params);
    add_to_order(it->second);
    return diag_output(it->second,cname).alias(cname);
  }

  // Operands first: this is what puts 'order' in dependency order.
  std::map<std::string,Field> operands;
  for (const auto& c : s.children) {
    operands[c.names.canonical] = build(c,grid,fm,repo,known,diag_params,order);
  }

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

  // Options given in the output yaml for this kind of diag, if any. They are
  // overlaid last, so a knob set there wins over what the expression lowered
  // to. This is how a diag gets settings that its expression has no syntax for
  // (e.g. how often COSP should actually run).
  if (diag_params.isSublist(key)) {
    overlay_params(params,diag_params.sublist(key));
  }

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
      // Report a failure here in terms of the diag that wanted the input, not
      // in terms of the recursion. Otherwise a diag needing a field the model
      // does not compute is reported as an unrecognized *diagnostic* named after
      // that field, which sends the reader looking in the wrong place.
      Field dep;
      try {
        dep = build(lower_to_diag_spec(in_name,in_name),grid,fm,repo,known,diag_params,order);
      } catch (const std::exception& e) {
        EKAT_ERROR_MSG (
            "Error! Could not resolve an input of a diagnostic.\n"
            " - diag      : " + diag->name() + "\n"
            " - input     : " + in_name + "\n"
            " - expression: " + cname + "\n"
            " - note: the diag needs this field, but it is neither computed by the\n"
            "         model nor obtainable as a diagnostic. Note that a field only\n"
            "         exists if some atm process asks for it, so a diagnostic can\n"
            "         need an input that no process in this run provides.\n"
            " - underlying error:\n" + std::string(e.what()) + "\n");
      }
      diag->set_input_field(dep.alias(in_name));
    }
  }

  diag->initialize();
  repo[cname] = diag;

  // A diag computing several fields must be found again when a *sibling* output
  // is requested, or a second instance would be built and would recompute the
  // whole suite. The canonical name above only covers the one output that
  // happened to be asked for first.
  const auto out_names = diag->get_output_names();
  if (out_names.size()>1) {
    for (const auto& o : out_names) {
      repo[o] = diag;
    }
  }

  add_to_order(diag);

  // The diag names its own output field, following the legacy convention. The
  // tree refers to sub-expressions by canonical name instead, so alias it.
  // Aliasing shares the data, the tracking, and the extra data, so a mask or a
  // timestamp set on the diag field is still visible through the alias.
  return diag_output(diag,cname).alias(cname);
}

} // anonymous namespace

DiagTree build_diag_tree (const DiagSpec& spec,
                          const std::shared_ptr<const AbstractGrid>& grid,
                          const FieldManager& fm,
                          std::map<std::string,diag_ptr_t>& repo,
                          const std::map<std::string,Field>& known,
                          const ekat::ParameterList& diag_params)
{
  DiagTree tree;
  auto f = build(spec,grid,fm,repo,known,diag_params,tree.eval_order);

  // Unlike the canonical name, the registered one is what the FieldManager and
  // the output file will know this by.
  tree.top = spec.names.registered.empty()
           ? f : f.alias(spec.names.registered);
  return tree;
}

} // namespace scream
