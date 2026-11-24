//===-- infra/CFUnits.cpp - units algebra for CF units strings --*- C++ -*-===//
//
// Implementation of the CFUnits class. A units string is split on slashes
// into a numerator and denominators, each of those on whitespace into
// factors, and each factor into a unit symbol and an optional integer
// exponent. Factors are kept in order of first appearance so that the
// formatted result reads like the conventional spelling (kg m-3, not
// m-3 kg). A leading number is the scale factor of a scaled unit.
//
//===----------------------------------------------------------------------===//

#include "CFUnits.h"
#include "Logging.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <sstream>

namespace OMEGA {

namespace {

// Relative tolerance for comparing scale factors, which pick up rounding
// error in products and powers (1/1e-6 is not exactly 1e6)
constexpr double ScaleTol = 1.0e-12;

bool sameScale(double A, double B) {
   return std::abs(A - B) <= ScaleTol * std::max(std::abs(A), std::abs(B));
}

// Formats a scale factor with the fewest significant digits that reproduce
// it to within rounding, so 1/1e-6 is written 1e+06, not 999999.9999999999,
// and without an exponent when that is no longer (1000, not 1e+03)
std::string formatScale(double Value) {
   char Buffer[32];
   int Digits = 1;
   for (; Digits < 17; ++Digits) {
      std::snprintf(Buffer, sizeof(Buffer), "%.*g", Digits, Value);
      if (sameScale(std::strtod(Buffer, nullptr), Value))
         break;
   }
   std::string Shortest = Buffer;

   // %g uses an exponent once the integer part has more digits than the
   // precision; asking for all of them gives the fixed spelling instead
   int Exponent =
       static_cast<int>(std::floor(std::log10(std::strtod(Buffer, nullptr))));
   if (Exponent >= Digits) {
      std::snprintf(Buffer, sizeof(Buffer), "%.*g", Exponent + 1, Value);
      if (std::string(Buffer).size() <= Shortest.size())
         return Buffer;
   }
   return Shortest;
}

} // namespace

//------------------------------------------------------------------------------
// Default constructor creates unknown units
CFUnits::CFUnits() : Known(false), Scale(1.0) {}

//------------------------------------------------------------------------------
// Parses a units string into its factors
Error CFUnits::parse(const std::string &UnitsStr, // [in] units string
                     CFUnits &Result // [out] parsed units expression
) {

   Error Err;

   Result = CFUnits();

   // Trim surrounding whitespace; an empty string is the unknown units of a
   // field with no units attribute
   auto First = UnitsStr.find_first_not_of(" \t");
   if (First == std::string::npos)
      return Err;
   auto Last           = UnitsStr.find_last_not_of(" \t");
   std::string Trimmed = UnitsStr.substr(First, Last - First + 1);

   Result.Known = true;

   // A trailing slash has no denominator after it; getline would skip it
   if (Trimmed.back() == '/') {
      Result = CFUnits();
      RETURN_ERROR(Err, ErrorCode::Fail,
                   "CFUnits: cannot parse units string '{}': "
                   "empty denominator",
                   UnitsStr);
   }

   // Each slash starts a denominator: every factor after it divides
   int Sign = 1;
   std::string Segment;
   std::istringstream Segments(Trimmed);
   bool FirstSegment = true;
   bool FirstToken   = true;
   while (std::getline(Segments, Segment, '/')) {

      // A slash with nothing before or after it is malformed
      std::istringstream Tokens(Segment);
      std::string Token;
      bool AnyToken = false;
      while (Tokens >> Token) {
         AnyToken = true;

         // A number other than 1 is a scale factor, which udunits allows
         // only before the unit symbols
         bool IsNumber = std::isdigit(static_cast<unsigned char>(Token[0])) ||
                         Token[0] == '.';
         if (IsNumber && Token != "1") {
            double Value = 1.0;
            if (!FirstToken || !parseScale(Token, Value)) {
               Result = CFUnits();
               RETURN_ERROR(Err, ErrorCode::Fail,
                            "CFUnits: cannot parse units string '{}': "
                            "'{}' is not a positive scale factor before the "
                            "unit symbols",
                            UnitsStr, Token);
            }
            Result.Scale = Value;
            FirstToken   = false;
            continue;
         }
         FirstToken = false;

         if (!Result.parseFactor(Token, Sign)) {
            Result = CFUnits();
            RETURN_ERROR(Err, ErrorCode::Fail,
                         "CFUnits: cannot parse units string '{}': "
                         "unexpected factor '{}'",
                         UnitsStr, Token);
         }
      }
      if (!AnyToken) {
         Result = CFUnits();
         RETURN_ERROR(Err, ErrorCode::Fail,
                      "CFUnits: cannot parse units string '{}': "
                      "empty {}",
                      UnitsStr, FirstSegment ? "numerator" : "denominator");
      }

      FirstSegment = false;
      Sign         = -1;
   }

   return Err;

} // end parse

//------------------------------------------------------------------------------
// Parses one factor: a unit symbol of letters and underscores followed by an
// optional integer exponent, with or without a caret (m2, s-1, m^2, s^-1), or
// the digit 1 which contributes nothing
bool CFUnits::parseFactor(const std::string &Token, // [in] one factor
                          int Sign                  // [in] +1 or -1
) {

   if (Token == "1")
      return true;

   // Symbol: one or more letters or underscores
   std::size_t Pos = 0;
   while (Pos < Token.size() &&
          (std::isalpha(static_cast<unsigned char>(Token[Pos])) ||
           Token[Pos] == '_'))
      ++Pos;
   if (Pos == 0)
      return false;
   std::string Symbol = Token.substr(0, Pos);

   // Exponent: optional caret, optional minus sign, then digits to the end
   int Exponent = 1;
   if (Pos < Token.size()) {
      if (Token[Pos] == '^')
         ++Pos;
      int ExpSign = 1;
      if (Pos < Token.size() && Token[Pos] == '-') {
         ExpSign = -1;
         ++Pos;
      }
      if (Pos >= Token.size())
         return false;
      int Magnitude = 0;
      for (; Pos < Token.size(); ++Pos) {
         if (!std::isdigit(static_cast<unsigned char>(Token[Pos])))
            return false;
         Magnitude = 10 * Magnitude + (Token[Pos] - '0');
      }
      Exponent = ExpSign * Magnitude;
   }

   multiplyFactor(Symbol, Sign * Exponent);
   return true;

} // end parseFactor

//------------------------------------------------------------------------------
// Parses a scale factor: the whole token must be a positive, finite number
bool CFUnits::parseScale(const std::string &Token, // [in] numeric token
                         double &Value             // [out] scale factor
) {

   char *End = nullptr;
   Value     = std::strtod(Token.c_str(), &End);
   return End == Token.c_str() + Token.size() && std::isfinite(Value) &&
          Value > 0.0;

} // end parseScale

//------------------------------------------------------------------------------
// Multiplies in one symbol raised to a power
void CFUnits::multiplyFactor(const std::string &Symbol, // [in] unit symbol
                             int Exponent               // [in] its power
) {

   if (Exponent == 0)
      return;

   for (auto Iter = Factors.begin(); Iter != Factors.end(); ++Iter) {
      if (Iter->first == Symbol) {
         Iter->second += Exponent;
         if (Iter->second == 0)
            Factors.erase(Iter);
         return;
      }
   }
   Factors.emplace_back(Symbol, Exponent);

} // end multiplyFactor

//------------------------------------------------------------------------------
// Parses a units string, aborting on failure
CFUnits CFUnits::parseOrAbort(const std::string &UnitsStr // [in] units string
) {
   CFUnits Result;
   Error Err = parse(UnitsStr, Result);
   CHECK_ERROR_ABORT(Err, "CFUnits: invalid units string");
   return Result;
}

//------------------------------------------------------------------------------
// String interfaces for the product, quotient and power
std::string CFUnits::multiply(const std::string &UnitsA, // [in] first units
                              const std::string &UnitsB  // [in] second units
) {
   return (parseOrAbort(UnitsA) * parseOrAbort(UnitsB)).str();
}

std::string CFUnits::divide(const std::string &UnitsA, // [in] numerator
                            const std::string &UnitsB  // [in] denominator
) {
   return (parseOrAbort(UnitsA) / parseOrAbort(UnitsB)).str();
}

std::string CFUnits::power(const std::string &UnitsA, // [in] units
                           int Exponent               // [in] power
) {
   return parseOrAbort(UnitsA).pow(Exponent).str();
}

std::string CFUnits::scale(const std::string &UnitsA, // [in] units
                           double Factor // [in] positive scale factor
) {
   return parseOrAbort(UnitsA).scaled(Factor).str();
}

//------------------------------------------------------------------------------
// Queries
bool CFUnits::isUnknown() const { return !Known; }

bool CFUnits::isDimensionless() const { return Known && Factors.empty(); }

//------------------------------------------------------------------------------
// Formats the units in the CF plain form
std::string CFUnits::str() const {

   if (!Known)
      return "";

   bool Unscaled = sameScale(Scale, 1.0);
   if (Factors.empty())
      return Unscaled ? "1" : formatScale(Scale);

   std::string Result = Unscaled ? "" : formatScale(Scale);
   for (int Pass = 0; Pass < 2; ++Pass) {
      for (const auto &Factor : Factors) {
         bool Positive = Factor.second > 0;
         if (Positive != (Pass == 0))
            continue;
         if (!Result.empty())
            Result += " ";
         Result += Factor.first;
         if (Factor.second != 1)
            Result += std::to_string(Factor.second);
      }
   }
   return Result;

} // end str

//------------------------------------------------------------------------------
// Algebra
CFUnits CFUnits::operator*(const CFUnits &Other) const {

   if (!Known || !Other.Known)
      return CFUnits();

   CFUnits Result = *this;
   Result.Scale *= Other.Scale;
   for (const auto &Factor : Other.Factors)
      Result.multiplyFactor(Factor.first, Factor.second);
   return Result;
}

CFUnits CFUnits::operator/(const CFUnits &Other) const {
   return *this * Other.pow(-1);
}

CFUnits CFUnits::pow(int Exponent) const {

   if (!Known)
      return CFUnits();

   CFUnits Result;
   Result.Known = true;
   Result.Scale = std::pow(Scale, Exponent);
   for (const auto &Factor : Factors)
      Result.multiplyFactor(Factor.first, Factor.second * Exponent);
   return Result;
}

CFUnits CFUnits::scaled(double Factor) const {

   OMEGA_REQUIRE(std::isfinite(Factor) && Factor > 0.0,
                 "CFUnits: scale factor {} is not positive and finite", Factor);

   if (!Known)
      return CFUnits();

   CFUnits Result = *this;
   Result.Scale *= Factor;
   return Result;
}

//------------------------------------------------------------------------------
// Comparison, ignoring the order of factors
bool CFUnits::operator==(const CFUnits &Other) const {

   if (Known != Other.Known)
      return false;
   std::map<std::string, int> Mine(Factors.begin(), Factors.end());
   std::map<std::string, int> Theirs(Other.Factors.begin(),
                                     Other.Factors.end());
   return Mine == Theirs && sameScale(Scale, Other.Scale);
}

bool CFUnits::operator!=(const CFUnits &Other) const {
   return !(*this == Other);
}

} // end namespace OMEGA

//===----------------------------------------------------------------------===//
