#ifndef EAMXX_LEVEL_INDEX_UTILS_HPP
#define EAMXX_LEVEL_INDEX_UTILS_HPP

#include "share/core/eamxx_types.hpp"

#include <ekat_assert.hpp>

#include <cstdio>
#include <limits>
#include <string>

// Helpers shared by VerticalLevelIndex and FieldVerticalInterp, so that they agree on
// the name of the index diag, and on how a vertical target is interpreted.

namespace scream
{
namespace level_index_utils
{

// Canonical spelling of a vertical target value, used to name the range masks, so that
// equivalent requests (e.g., 900hPa and 900mb) share the same mask. digits10 significant
// digits absorb the rounding error of unit conversions (e.g., 850.3hPa to Pa and back), and
// trailing zeros are dropped. No exponent, since it would end up in a variable name.
inline std::string value_string (const Real v)
{
  char buf[64];
  std::snprintf(buf,sizeof(buf),"%.*g",std::numeric_limits<Real>::digits10,static_cast<double>(v));
  std::string s = buf;
  if (s.find('e')!=std::string::npos) {
    std::snprintf(buf,sizeof(buf),"%.12f",static_cast<double>(v));
    s = buf;
    while (s.back()=='0') {
      s.pop_back();
    }
    if (s.back()=='.') {
      s.pop_back();
    }
  }
  return s;
}

// Convert a pressure value, expressed in the given units ('mb', 'hPa', 'Pa'), to Pa
inline Real pressure_in_Pa (const std::string& value, const std::string& units)
{
  EKAT_REQUIRE_MSG (units=="mb" or units=="hPa" or units=="Pa",
      "Error! Invalid pressure units.\n"
      " - input units: " + units + "\n"
      " - valid units: 'mb', 'hPa', 'Pa'\n");
  const Real v = std::stod(value);
  return (units=="Pa") ? v : v*100;
}

// A vertical target, as requested (e.g., 500 hPa, or 2 m above the surface), along with
// what we need to know about the vertical coordinate it refers to.
// The coordinate is one of
//  - "p": the pressure (p_mid/p_int), increasing with the level index,
//  - "z": the height above sea level (z_mid/z_int), decreasing with the level index,
//  - "height": the height above the surface (height_mid/height_int), decreasing with
//    the level index.
struct VerticalTarget
{
  VerticalTarget (const std::string& coordinate,
                  const std::string& value,
                  const std::string& units_)
   : coord (coordinate)
   , value_str (value)
   , units (units_)
  {
    EKAT_REQUIRE_MSG (coord=="p" or coord=="z" or coord=="height",
        "Error! Invalid vertical coordinate.\n"
        " - input value : " + coord + "\n"
        " - valid values: 'p', 'z', 'height'\n");
    if (coord=="p") {
      target = pressure_in_Pa(value_str,units);
    } else {
      EKAT_REQUIRE_MSG (units=="m",
          "Error! Invalid height units.\n"
          " - input units: " + units + "\n"
          " - valid units: m\n");
      target = std::stod(value_str);
    }
    // Pressure grows with the level index, heights decrease with it
    descending = coord!="p";
    // An "above surface" height is always well defined near the surface (every column
    // has one), so we extrapolate a target below the bottom entry. An "above sealevel"
    // height is not (e.g., a target elevation below a mountain's surface), and neither
    // is a pressure higher than the bottom one, so those are flagged as invalid.
    extrapolate_bottom = coord=="height";
  }

  // Name of the coordinate field, for a given layer ('mid' or 'int')
  std::string coord_field (const std::string& layer) const
  {
    return coord + "_" + layer;
  }

  // Name of the VerticalLevelIndex diag for this target, for a given layer. It is the
  // expression that builds it (see eamxx_dexpr_diags.cpp), so it must be spelled like
  // the expression front end spells it.
  std::string index_name (const std::string& layer) const
  {
    return coord_field(layer) + ".searchsorted(" + value_str + ",units='" + units + "')";
  }

  // Name of the range mask for this target, for a given layer. It only depends on the
  // level, not on how the request spells it (e.g., 900hPa and 900mb), so that all equivalent
  // requests share the same mask name, and IO the same valid count. The pressure is spelled in hPa,
  // which is what the requests normally use, so that their mask names (and the IO count names)
  // are the ones of the request.
  std::string mask_name (const std::string& layer) const
  {
    const std::string lvl = coord=="p"
                          ? value_string(target/100) + "hPa"
                          : value_string(target) + "m_above_" + (coord=="z" ? "sealevel" : "surface");
    return lvl + (layer=="int" ? "_int" : "") + "_mask";
  }

  // Suffix of the name of the field interpolated at this target
  std::string field_suffix () const
  {
    return "_at_" + value_str + units + (coord=="p" ? "" : (coord=="z" ? "_above_sealevel" : "_above_surface"));
  }

  std::string coord;
  std::string value_str;
  std::string units;

  Real target;              // In Pa for the pressure, in m for the heights
  bool descending;          // Whether the coordinate decreases with the level index
  bool extrapolate_bottom;  // Whether a target below the bottom entry is valid
};

} // namespace level_index_utils
} // namespace scream

#endif // EAMXX_LEVEL_INDEX_UTILS_HPP
