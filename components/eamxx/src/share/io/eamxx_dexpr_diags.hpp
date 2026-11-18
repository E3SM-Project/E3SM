#ifndef SCREAM_EAMXX_DEXPR_DIAGS_HPP
#define SCREAM_EAMXX_DEXPR_DIAGS_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"
#include "share/grid/abstract_grid.hpp"

#include <memory>
#include <string>
#include <vector>

namespace scream {

/*
 * The expression front end for diagnostics.
 *
 * Parses with share/dexpr and drives the diagnostic factory from the AST.
 * create_diagnostic() sends every request here first; only a plain name
 * (a field, or a diag class like LiqWaterPath) comes back unresolved.
 *
 * Composition is by field name, as always. A diag declares its operands by
 * their canonical expression string, and the customer resolving dependencies
 * feeds those back through create_diagnostic(): '(qc+qv)*p_mid' gives
 * BinaryOp{arg1="(qc+qv)", arg2="p_mid"}, and '(qc+qv)' comes back here.
 *
 * Diags built here carry two extra params:
 *   - 'output_field_name': the expression, so the output field is named after
 *     the request rather than the diag's own param concatenation.
 *   - 'from_expression': marks a name as resolved by this front end. An
 *     expression is not a usable NetCDF variable name, so a customer writing
 *     one to file must require 'name := expr'.
 *
 * NOTE: exposes no dexpr types, keeping dexpr private to eamxx_io.
 */

// Returns nullptr if 'expr' is a plain identifier: that is a diag class name
// (or a typo) for the caller to resolve as before. Throws if 'expr' is an
// expression we cannot parse, validate, or translate.
std::shared_ptr<AbstractDiagnostic>
dexpr_create_diagnostic (const std::string& expr,
                         const std::shared_ptr<const AbstractGrid>& grid);

} // namespace scream

#endif // SCREAM_EAMXX_DEXPR_DIAGS_HPP
