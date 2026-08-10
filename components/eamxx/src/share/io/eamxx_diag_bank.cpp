#include "share/io/eamxx_diag_bank.hpp"

#include "share/io/eamxx_legacy_diag_names.hpp"

#include <ekat_assert.hpp>
#include <ekat_string_utils.hpp>

#include <algorithm>

namespace scream {

namespace {

// ekat::split does not trim, so ' a := b ' would otherwise register a field
// whose name literally starts and ends with a space. Whitespace is meaningless
// on both sides of ':=', so drop it.
std::string trim (const std::string& s)
{
  const auto beg = s.find_first_not_of(" \t");
  if (beg==std::string::npos) return "";
  const auto end = s.find_last_not_of(" \t");
  return s.substr(beg,end-beg+1);
}

// The bare names an expression bottoms out in. Any of them may be another
// request, which then has to be built first.
void collect_leaves (const DiagSpec& s, std::vector<std::string>& out)
{
  if (s.is_leaf()) {
    out.push_back(s.names.canonical);

    // A legacy name keeps its operands hidden until it is translated, so look
    // through it: 'new_stuff_horiz_avg' depends on 'new_stuff', and we cannot
    // see that until the name is desugared. Without this, the dependency is
    // only discovered while walking the tree, by which time it is too late to
    // build it, and the request would have to be listed first to work.
    // NOTE: the translation strips one suffix at a time, so this terminates.
    const auto expr = legacy_to_expr(s.names.canonical);
    if (expr!=s.names.canonical) {
      try {
        collect_leaves(lower_to_diag_spec(expr),out);
      } catch (...) {
        // Not a valid expression after all. Say nothing here: the walk will
        // report it properly, with the full context.
      }
    }
  }
  for (const auto& c : s.children) {
    collect_leaves(c,out);
  }
}

} // anonymous namespace

DiagBank::
DiagBank (const std::shared_ptr<const AbstractGrid>& grid,
          const FieldManager& fm)
 : DiagBank(grid,fm,m_own_repo)
{
  // Nothing else to do
}

DiagBank::
DiagBank (const std::shared_ptr<const AbstractGrid>& grid,
          const FieldManager& fm,
          std::map<std::string,diag_ptr_t>& repo)
 : m_grid(grid)
 , m_fm(fm)
 , m_repo(repo)
{
  EKAT_REQUIRE_MSG (grid, "[DiagBank] Error! Invalid grid pointer.\n");
}

std::string DiagBank::add (const std::string& request, const bool write)
{
  const auto tokens = ekat::split(request,":=");
  EKAT_REQUIRE_MSG (tokens.size()==1 or tokens.size()==2,
      "Error! Invalid diagnostic request. Should be 'expr' or 'name:=expr'.\n"
      " - request: " + request + "\n");

  Request r;
  std::string registered;
  if (tokens.size()==2) {
    registered = trim(tokens[0]);
    r.expr     = trim(tokens[1]);
    EKAT_REQUIRE_MSG (not registered.empty() and not r.expr.empty(),
        "Error! Invalid diagnostic request. Should be 'name:=expr'.\n"
        " - request: " + request + "\n");
  } else {
    // No name given, so the request doubles as its own name. This is what makes
    // the legacy names work unchanged: they are already legal nc var names.
    registered = trim(request);
    r.expr     = registered;
  }
  r.write = write;

  EKAT_REQUIRE_MSG (m_requests.count(registered)==0,
      "Error! The same name is requested twice.\n"
      " - name       : " + registered + "\n"
      " - first  expr: " + m_requests[registered].expr + "\n"
      " - second expr: " + r.expr + "\n");

  m_requests[registered] = r;
  m_added.push_back(registered);

  return registered;
}

const std::string& DiagBank::expr_of (const std::string& registered) const
{
  auto it = m_requests.find(registered);
  EKAT_REQUIRE_MSG (it!=m_requests.end(),
      "Error! No diagnostic was requested under this name.\n"
      " - name: " + registered + "\n");
  return it->second.expr;
}

void DiagBank::build ()
{
  for (const auto& name : m_added) {
    build_one(name);
  }
}

const Field& DiagBank::build_one (const std::string& registered)
{
  if (auto it=m_entries.find(registered); it!=m_entries.end()) {
    return it->second.field;
  }

  EKAT_REQUIRE_MSG (m_building.count(registered)==0,
      "Error! Circular dependency between diagnostic requests.\n"
      " - name: " + registered + "\n"
      " - expr: " + m_requests.at(registered).expr + "\n");
  m_building.insert(registered);

  const auto& r = m_requests.at(registered);
  const auto spec = lower_to_diag_spec(r.expr,registered,r.write);

  // A leaf naming another request must be built first, so that the walk has
  // something to bottom out on. This is why requests can be given in any order.
  std::vector<std::string> leaves;
  collect_leaves(spec,leaves);
  for (const auto& l : leaves) {
    // Skip the request's own name: 'T_mid' asks for the field T_mid, and
    // 'foo:=foo' renames it. Neither is a request depending on itself.
    if (l!=registered and m_requests.count(l)==1 and m_entries.count(l)==0) {
      build_one(l);
    }
  }

  const auto tree = build_diag_tree(spec,m_grid,m_fm,m_repo,m_known,m_diag_params);

  // A diag shared with an earlier request is already in the order, ahead of
  // everything that consumes it, so it must not be added again.
  for (const auto& d : tree.eval_order) {
    if (std::find(m_eval_order.begin(),m_eval_order.end(),d)==m_eval_order.end()) {
      m_eval_order.push_back(d);
    }
  }

  m_building.erase(registered);
  m_known[registered]   = tree.top;

  // If the request bottomed out on a model field, no diag was built for it, and
  // this entry is nothing but a rename.
  const auto diag = tree.eval_order.empty() ? nullptr : tree.eval_order.back();
  m_entries[registered] = Entry{tree.top,r.write,diag};

  return m_entries[registered].field;
}

const Field& DiagBank::field (const std::string& registered) const
{
  auto it = m_entries.find(registered);
  EKAT_REQUIRE_MSG (it!=m_entries.end(),
      "Error! No diagnostic was built under this name.\n"
      " - name: " + registered + "\n");
  return it->second.field;
}

} // namespace scream
