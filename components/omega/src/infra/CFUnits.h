#ifndef OMEGA_CFUNITS_H
#define OMEGA_CFUNITS_H
//===-- infra/CFUnits.h - units algebra for CF units strings ----*- C++ -*-===//
//
/// \file
/// \brief Defines the CFUnits class for deriving units strings
///
/// Omega fields carry their units as a string attribute in the plain form
/// used by the CF standard name table: unit symbols separated by single
/// spaces, each followed by an optional signed integer exponent (m s-1,
/// kg m-3, W m-2) and the dimensionless quantity written as 1. When a field
/// is derived from others, for example by an area-weighted sum or a product,
/// its units must be derived too. CFUnits parses a units string into its
/// factors, supports the product, quotient and integer power of two such
/// expressions, and formats the result back in the plain form.
///
/// A units string may start with a positive number, the scale factor of a
/// scaled unit in udunits (1e+06 m3 s-1 is a million cubic meters per
/// second). Rescaling a field's values, for example to report a transport in
/// units of 10^6 m3 s-1, rescales its units the other way, so CFUnits carries
/// the scale factor through every operation.
///
/// A field with no units attribute has unknown units, represented by an empty
/// string. Unknown units propagate: any expression involving them is unknown,
/// so a derived field never claims units its inputs did not have. Parsing
/// also accepts the udunits spellings with a caret before the exponent (m^2)
/// and a slash for division (m/s), but always emits the plain form. Anything
/// else is a parse error, since a malformed units string is a programming
/// error in the field definition and should stop the run at initialization.
///
//===----------------------------------------------------------------------===//

#include "Error.h"

#include <string>
#include <utility>
#include <vector>

namespace OMEGA {

/// A units expression as a product of unit symbols raised to integer powers,
/// or the unknown units of a field that has no units attribute
class CFUnits {

 public:
   //---------------------------------------------------------------------------
   /// Default constructor creates unknown units
   CFUnits();

   //---------------------------------------------------------------------------
   /// Parses a units string. An empty string gives unknown units and "1" gives
   /// dimensionless units. Returns a Fail error if the string is not a product
   /// of unit symbols with optional integer exponents, in the CF plain form or
   /// the udunits caret and slash spellings.
   static Error parse(const std::string &UnitsStr, ///< [in] units string
                      CFUnits &Result ///< [out] parsed units expression
   );

   //---------------------------------------------------------------------------
   /// Returns the product of two units strings in the CF plain form, or an
   /// empty string if either is unknown. Aborts if either string is invalid.
   static std::string multiply(const std::string &UnitsA, ///< [in] first units
                               const std::string &UnitsB  ///< [in] second units
   );

   /// Returns the quotient of two units strings in the CF plain form, or an
   /// empty string if either is unknown. Aborts if either string is invalid.
   static std::string divide(const std::string &UnitsA, ///< [in] numerator
                             const std::string &UnitsB  ///< [in] denominator
   );

   /// Returns a units string raised to an integer power in the CF plain form,
   /// or an empty string if it is unknown. Aborts if the string is invalid.
   static std::string power(const std::string &UnitsA, ///< [in] units
                            int Exponent               ///< [in] power
   );

   /// Returns a units string multiplied by a positive scale factor in the CF
   /// plain form (scale("m3 s-1", 1e6) is "1e+06 m3 s-1"), or an empty string
   /// if it is unknown. Aborts if the string is invalid or the factor is not
   /// positive and finite.
   static std::string scale(const std::string &UnitsA, ///< [in] units
                            double Factor ///< [in] positive scale factor
   );

   //---------------------------------------------------------------------------
   /// True if the units are unknown (the field has no units attribute)
   bool isUnknown() const;

   /// True if the units are known and dimensionless (a pure number, which
   /// may be scaled)
   bool isDimensionless() const;

   /// Formats the units in the CF plain form: the scale factor if it is not
   /// 1, then symbols with positive exponents, then those with negative
   /// exponents, each in order of first appearance; "1" if dimensionless and
   /// unscaled; an empty string if unknown
   std::string str() const;

   //---------------------------------------------------------------------------
   /// Product of two units expressions; unknown if either is unknown
   CFUnits operator*(const CFUnits &Other) const;

   /// Quotient of two units expressions; unknown if either is unknown
   CFUnits operator/(const CFUnits &Other) const;

   /// Units expression raised to an integer power; unknown if unknown
   CFUnits pow(int Exponent) const;

   /// Units expression multiplied by a positive scale factor; unknown if
   /// unknown. Aborts if the factor is not positive and finite.
   CFUnits scaled(double Factor) const;

   /// True if both are unknown, or both are known with the same factors in
   /// any order and the same scale factor to within rounding
   bool operator==(const CFUnits &Other) const;
   bool operator!=(const CFUnits &Other) const;

 private:
   /// False for the unknown units of a field with no units attribute
   bool Known;

   /// Unit symbols and their non-zero exponents, in order of first appearance
   std::vector<std::pair<std::string, int>> Factors;

   /// Positive scale factor multiplying the unit symbols, 1 if unscaled
   double Scale;

   /// Parses a units string, aborting on failure, for the string interfaces
   static CFUnits parseOrAbort(const std::string &UnitsStr);

   /// Multiplies in one symbol raised to a power, merging with an existing
   /// factor of the same symbol and dropping the factor if its exponent
   /// becomes zero
   void multiplyFactor(const std::string &Symbol, int Exponent);

   /// Parses one whitespace-delimited token (a symbol with an optional
   /// exponent, or the digit 1) and multiplies it in with the given sign on
   /// the exponent. Returns false if the token is malformed.
   bool parseFactor(const std::string &Token, int Sign);

   /// Parses a token that starts with a digit or a decimal point as a scale
   /// factor. Returns false unless the whole token is a positive, finite
   /// number.
   static bool parseScale(const std::string &Token, double &Value);

}; // end class CFUnits

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
#endif // OMEGA_CFUNITS_H
