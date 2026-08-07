#ifndef SCREAM_EAMXX_DIAG_BANK_HPP
#define SCREAM_EAMXX_DIAG_BANK_HPP

#include "share/io/eamxx_diag_tree.hpp"

#include <map>
#include <set>
#include <string>
#include <vector>

namespace scream {

/*
 * The bank of names for one set of diagnostic requests.
 *
 * A diagnostic name has three jobs, and the whole point of this class is that
 * it keeps them apart (see DiagNames in eamxx_diag_spec.hpp):
 *
 *  - the *canonical* name identifies a sub-expression. It is what makes two
 *    requests share a diag when they ask for the same thing, and it is how a
 *    node refers to its operands. Internal, unambiguous, never re-parsed.
 *
 *  - the *registered* name is what a request is known by: what the FieldManager
 *    and the output file call it, and what other requests can refer to it as.
 *    It comes from the 'name:=expr' syntax, or defaults to the request itself.
 *
 *  - the *write* flag says whether it reaches the nc file. It is false for the
 *    intermediates declared in the yaml 'aliases' section, which must exist so
 *    that other requests can use them, but must not be written.
 *
 * Requests may refer to each other by registered name, in any order; the bank
 * works out what to build first, and rejects cycles.
 */
class DiagBank {
public:
  struct Entry {
    Field field;   // named after the request's registered name
    bool  write;
  };

  DiagBank (const std::shared_ptr<const AbstractGrid>& grid,
            const FieldManager& fm);

  // Add a request, spelled as it is in the output yaml: either 'expr' or
  // 'name:=expr'. 'write' is false for entries in the 'aliases' section.
  void add (const std::string& request, const bool write = true);

  // Build every request that was added.
  void build ();

  // All the diags to evaluate, with inputs always ahead of their consumers.
  const std::vector<diag_ptr_t>& eval_order () const { return m_eval_order; }

  // What was built, by registered name.
  const std::map<std::string,Entry>& entries () const { return m_entries; }

  const Field& field (const std::string& registered) const;

private:
  const Field& build_one (const std::string& registered);

  struct Request {
    std::string expr;
    bool        write;
  };

  std::shared_ptr<const AbstractGrid> m_grid;
  const FieldManager&                 m_fm;

  std::map<std::string,Request> m_requests;
  std::vector<std::string>      m_added;      // requests, in the order added

  std::map<std::string,Entry>      m_entries;    // registered -> field
  std::map<std::string,Field>      m_known;      // registered -> field, for the walk
  std::map<std::string,diag_ptr_t> m_repo;       // canonical  -> diag, for reuse
  std::vector<diag_ptr_t>          m_eval_order;
  std::set<std::string>            m_building;   // to catch cycles
};

} // namespace scream

#endif // SCREAM_EAMXX_DIAG_BANK_HPP
