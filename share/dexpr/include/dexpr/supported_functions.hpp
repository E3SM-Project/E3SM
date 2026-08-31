/**
 * @file supported_functions.hpp
 * @brief The set of callables a component makes available to expressions.
 *
 * dexpr owns the grammar; a component owns the vocabulary. The parser accepts
 * any call syntactically -- `foo(a, b=c)` parses whether or not `foo` exists --
 * and a component says what it can evaluate by filling a FunctionRegistry and
 * running validate_calls() over the AST. Keeping that out of the parser lets
 * one expression be checked against several components' vocabularies.
 */
#ifndef DEXPR_SUPPORTED_FUNCTIONS_HPP
#define DEXPR_SUPPORTED_FUNCTIONS_HPP

#include <dexpr/ast.hpp>

#include <map>
#include <ostream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <vector>

namespace dexpr {

// A keyword argument, e.g. the `dims` of `derivative(dx, dims=['col'])`.
struct ParamSpec {
  std::string name;
  bool required = false;
};

struct FunctionSpec {
  std::string name;
  std::string desc;

  // Positional arity, not counting keyword arguments. For a method call `x.f()`
  // the receiver is not positional, so `x.mean('lev')` has one.
  int min_positional = 0;
  int max_positional = 0;

  std::vector<ParamSpec> keywords;

  // How the call is written, e.g. "T_mid.mean('lev')". Documentation, and what
  // validate_registry() checks the spec against.
  std::string example;

  // "name(arg, kw=..)\n--- desc\n--- e.g. example", the `dexpr` tool listing.
  std::string to_string() const;
};

inline std::ostream& operator<<(std::ostream& os, const FunctionSpec& function) {
  return os << function.to_string() << '\n';
}

// Ordered, so listings and error messages are stable rather than hash-order.
class FunctionRegistry {
public:
  // Throws std::invalid_argument if the spec is malformed (no name, impossible
  // arity, a nameless or repeated keyword) or the name is already registered.
  void add(FunctionSpec spec);

  const FunctionSpec* find(std::string_view name) const;
  bool contains(std::string_view name) const { return find(name) != nullptr; }

  std::vector<std::string> names() const;

  auto begin() const { return fns_.begin(); }
  auto end() const { return fns_.end(); }

private:
  // std::less<> so lookups take a string_view without allocating.
  std::map<std::string, FunctionSpec, std::less<>> fns_;
};

// Generic callables, offered as a seed. Nothing consults it implicitly.
const FunctionRegistry& builtin_functions();

class ValidationError : public std::runtime_error {
public:
  explicit ValidationError(const std::vector<std::string>& errors);
};

// Checks every call in the tree: callee is a plain name, name is registered,
// positional arity fits, keywords are declared and present as required.
// Collects every problem and throws once.
void validate_calls(const ast::Expression& root, const FunctionRegistry& reg);

// Checks a vocabulary against itself: every spec's example must parse, pass
// validate_calls() against `reg`, and call the function that declared it. Run
// it after registering a function, from a test or from the registry itself.
void validate_registry(const FunctionRegistry& reg);

} // namespace dexpr

#endif
