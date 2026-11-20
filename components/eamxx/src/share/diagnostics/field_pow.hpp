#ifndef EAMXX_FIELD_POW_DIAG_HPP
#define EAMXX_FIELD_POW_DIAG_HPP

#include "share/diagnostics/abstract_diagnostic.hpp"

namespace scream {

/*
 * Raises a field X to a constant power p, elementwise.
 *
 * Parameters:
 *  - field_name: the input field X (Real data type, any layout)
 *  - exponent:   p, as a string: either a fraction 'n/d' (two integers, taken
 *                exactly: '1/3'), or a decimal number ('0.5', '-2'). It cannot be 0.
 *  - invalid:    what an entry outside the domain of X^p becomes:
 *                  nan   (default) whatever the math library returns (NaN or
 *                        inf). With floating point exceptions enabled, this
 *                        traps: use it when X is not expected to be invalid.
 *                  mask  the entry is masked out (and set to the fill value)
 *                The domain is X>0; X==0 if p>0; and X<0 if p is an integer.
 *  - units:      check   (default) the output has units u^p, with u the units of
 *                        X. That needs p to be a fraction n/d (see below), unless
 *                        X is dimensionless, which stays so for any p.
 *                ignore  X's units are not checked, and the output is unitless
 *
 * Units and decimal exponents: a unit is raised to a rational power, but a
 * decimal like 0.5 is not one. So we look for a fraction n/d (with d up to
 * 12) that matches it, to 1e-12. 0.5 and -1.5 do; 0.2857 does not, and needs
 * either units='ignore' or a dimensionless X. Writing the exponent as 'n/d'
 * skips this search, and has no limit on d.
 *
 * Squares, cubes and inverses (p=2, 3, -1) use multiplications and a division
 * rather than the (much more expensive) general power.
 *
 * If X has a valid mask, the output has one too, and masked entries of X are
 * not computed (they may hold the fill value). With invalid=nan, the output
 * shares X's mask.
 */

class FieldPow : public AbstractDiagnostic {
public:
  enum class Kind { General, Square, Cube, Inverse };

  FieldPow (const ekat::Comm &comm, const ekat::ParameterList &params,
            const std::shared_ptr<const AbstractGrid>& grid);

  // The name of the diagnostic CLASS (not the computed field)
  std::string name() const override { return "FieldPow"; }

protected:
#ifdef KOKKOS_ENABLE_CUDA
public:
#endif
  void compute_impl() override;

  void initialize_impl() override;

  std::string m_field_name;
  std::string m_exponent_str;
  Real        m_exponent = 0;
  Kind        m_kind = Kind::General;

  // The exponent as the fraction num/den, if it is known to be one
  bool        m_has_fraction = false;
  long long   m_num = 0;
  long long   m_den = 1;

  bool        m_mask_invalid = false;
  bool        m_out_masked   = false;
};

} // namespace scream

#endif // EAMXX_FIELD_POW_DIAG_HPP
