#include "share/io/eamxx_diag_bank.hpp"

#include <ekat_assert.hpp>
#include <ekat_string_utils.hpp>

#include <algorithm>

namespace scream {

namespace {

// The bare names an expression bottoms out in. Any of them may be another
// request, which then has to be built first.
void collect_leaves (const DiagSpec& s, std::vector<std::string>& out)
{
  if (s.is_leaf()) {
    out.push_back(s.names.canonical);
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

void DiagBank::add (const std::string& request, const bool write)
{
  const auto tokens = ekat::split(request,":=");
  EKAT_REQUIRE_MSG (tokens.size()==1 or tokens.size()==2,
      "Error! Invalid diagnostic request. Should be 'expr' or 'name:=expr'.\n"
      " - request: " + request + "\n");

  Request r;
  std::string registered;
  if (tokens.size()==2) {
    EKAT_REQUIRE_MSG (not tokens[0].empty() and not tokens[1].empty(),
        "Error! Invalid diagnostic request. Should be 'name:=expr'.\n"
        " - request: " + request + "\n");
    registered = tokens[0];
    r.expr     = tokens[1];
  } else {
    // No name given, so the request doubles as its own name. This is what makes
    // the legacy names work unchanged: they are already legal nc var names.
    registered = request;
    r.expr     = request;
  }
  r.write = write;

  EKAT_REQUIRE_MSG (m_requests.count(registered)==0,
      "Error! The same name is requested twice.\n"
      " - name       : " + registered + "\n"
      " - first  expr: " + m_requests[registered].expr + "\n"
      " - second expr: " + r.expr + "\n");

  m_requests[registered] = r;
  m_added.push_back(registered);
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

  const auto tree = build_diag_tree(spec,m_grid,m_fm,m_repo,m_known);

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
