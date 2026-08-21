#include "share/io/eamxx_legacy_diag_names.hpp"

#include <array>
#include <string>
#include <vector>

namespace scream {

namespace {

bool ends_with (const std::string& s, const std::string& suffix)
{
  return s.size()>suffix.size() and s.compare(s.size()-suffix.size(),suffix.size(),suffix)==0;
}

// A non-negative integer or fixed-point number, i.e. what '\d+(\.\d+)?' matched.
bool is_number (const std::string& s)
{
  if (s.empty()) return false;
  bool seen_dot = false, seen_digit_after_dot = true;
  for (size_t i=0; i<s.size(); ++i) {
    if (s[i]=='.') {
      if (seen_dot or i==0) return false;
      seen_dot = true;
      seen_digit_after_dot = false;
    } else if (std::isdigit(static_cast<unsigned char>(s[i]))) {
      seen_digit_after_dot = true;
    } else {
      return false;
    }
  }
  return seen_digit_after_dot;
}

// Position of the last occurrence of 'sep' in 's', or npos.
size_t rfind_sep (const std::string& s, const std::string& sep)
{
  return s.rfind(sep);
}

// ---------------------------------------------------------------------------
//  '_at_' family: level, pressure, height
// ---------------------------------------------------------------------------

// Try to read 's' as the tail of a legacy '_at_' name. Returns the argument
// list to pass to .at(), or "" if the tail does not match this flavor.
std::string at_level_args (const std::string& tail)
{
  if (tail=="model_top") return "lev='top'";
  if (tail=="model_bot") return "lev='bot'";
  if (tail.rfind("lev_",0)==0) {
    const auto lev = tail.substr(4);
    if (not lev.empty() and lev.find_first_not_of("0123456789")==std::string::npos) {
      return "lev=" + lev;
    }
  }
  return "";
}

std::string at_pressure_args (const std::string& tail)
{
  for (const auto& units : {std::string("hPa"),std::string("mb"),std::string("Pa")}) {
    if (ends_with(tail,units)) {
      const auto value = tail.substr(0,tail.size()-units.size());
      if (is_number(value)) {
        return "p=" + value + ", units='" + units + "'";
      }
    }
  }
  return "";
}

std::string at_height_args (const std::string& tail)
{
  for (const auto& ref : {std::string("sealevel"),std::string("surface")}) {
    const auto suffix = "m_above_" + ref;
    if (ends_with(tail,suffix)) {
      const auto value = tail.substr(0,tail.size()-suffix.size());
      if (is_number(value)) {
        return "z=" + value + ", units='m', ref='" + ref + "'";
      }
    }
  }
  return "";
}

// The legacy regexes anchor at the end and use a greedy prefix, so the operand
// runs as far right as it can: try the rightmost '_at_' first. Each flavor is
// tried over all positions before moving on to the next, which is what testing
// three separate regexes in sequence did.
std::string translate_at (const std::string& name)
{
  using args_fn = std::string (*)(const std::string&);
  for (const args_fn flavor : {&at_level_args,&at_pressure_args,&at_height_args}) {
    size_t pos = name.size();
    while ((pos=name.rfind("_at_",pos==0 ? 0 : pos-1))!=std::string::npos) {
      const auto args = flavor(name.substr(pos+4));
      if (not args.empty()) {
        return name.substr(0,pos) + ".at(" + args + ")";
      }
      if (pos==0) break;
    }
  }
  return "";
}

// ---------------------------------------------------------------------------
//  Conditional sampling and binary ops
// ---------------------------------------------------------------------------

std::string translate_where (const std::string& name)
{
  const auto wpos = rfind_sep(name,"_where_");
  if (wpos==std::string::npos or wpos==0) return "";

  const auto field = name.substr(0,wpos);
  const auto cond  = name.substr(wpos+7);

  // The lhs is greedy too, so pick the rightmost comparison token.
  constexpr std::array<std::pair<const char*,const char*>,6> cmps {{
    {"_gt_",">"}, {"_ge_",">="}, {"_eq_","=="},
    {"_ne_","!="}, {"_lt_","<"}, {"_le_","<="},
  }};
  size_t best_pos = std::string::npos;
  const char* best_op = nullptr;
  size_t best_len = 0;
  for (const auto& [tok,op] : cmps) {
    const auto p = cond.rfind(tok);
    if (p!=std::string::npos and (best_pos==std::string::npos or p>best_pos)) {
      best_pos = p;
      best_op  = op;
      best_len = std::string(tok).size();
    }
  }
  if (best_pos==std::string::npos or best_pos==0) return "";

  const auto lhs = cond.substr(0,best_pos);
  const auto rhs = cond.substr(best_pos+best_len);
  if (lhs.empty() or rhs.empty()) return "";

  return field + ".where(" + lhs + best_op + rhs + ")";
}

std::string translate_binary_op (const std::string& name)
{
  constexpr std::array<std::pair<const char*,const char*>,4> ops {{
    {"_plus_","+"}, {"_minus_","-"}, {"_times_","*"}, {"_over_","/"},
  }};
  size_t best_pos = std::string::npos;
  const char* best_op = nullptr;
  size_t best_len = 0;
  for (const auto& [tok,op] : ops) {
    const auto p = name.rfind(tok);
    if (p!=std::string::npos and (best_pos==std::string::npos or p>best_pos)) {
      best_pos = p;
      best_op  = op;
      best_len = std::string(tok).size();
    }
  }
  if (best_pos==std::string::npos or best_pos==0) return "";

  const auto lhs = name.substr(0,best_pos);
  const auto rhs = name.substr(best_pos+best_len);
  if (lhs.empty() or rhs.empty()) return "";

  return lhs + best_op + rhs;
}

// ---------------------------------------------------------------------------
//  Suffix-only patterns
// ---------------------------------------------------------------------------

// If 'name' ends with 'suffix', return what precedes it, else "".
std::string strip (const std::string& name, const std::string& suffix)
{
  return ends_with(name,suffix) ? name.substr(0,name.size()-suffix.size()) : std::string();
}

std::string translate_zonal_avg (const std::string& name)
{
  const auto head = strip(name,"_bins");
  if (head.empty()) return "";

  const auto pos = head.rfind("_zonal_avg_");
  if (pos==std::string::npos or pos==0) return "";

  const auto bins = head.substr(pos+11);
  if (bins.empty() or bins.find_first_not_of("0123456789")!=std::string::npos) return "";

  return head.substr(0,pos) + ".zonal_avg(bins=" + bins + ")";
}

std::string translate_histogram (const std::string& name)
{
  const auto pos = name.rfind("_histogram_");
  if (pos==std::string::npos or pos==0) return "";

  const auto config = name.substr(pos+11);
  std::vector<std::string> bins;
  size_t beg = 0;
  while (beg<=config.size()) {
    const auto end = config.find('_',beg);
    const auto tok = config.substr(beg,end==std::string::npos ? std::string::npos : end-beg);
    if (not is_number(tok)) return "";
    bins.push_back(tok);
    if (end==std::string::npos) break;
    beg = end+1;
  }
  if (bins.empty()) return "";

  std::string list;
  for (const auto& b : bins) {
    if (not list.empty()) list += ", ";
    list += b;
  }
  return name.substr(0,pos) + ".histogram([" + list + "])";
}

} // anonymous namespace

std::string legacy_to_expr (const std::string& name)
{
  std::string head, out;

  // The order below mirrors, one for one, the order of the regexes in the old
  // create_diagnostic. Do not reshuffle it; see the header for why.

  // Built-in alias, expanded before anything else.
  if (not (head=strip(name,"_atm_backtend")).empty()) return head + ".tend()";

  // _at_ family (level, then pressure, then height).
  if (not (out=translate_at(name)).empty()) return out;

  // NOTE: the bare-name diags (LiqWaterPath, precip_liq_surf_mass_flux, z_mid,
  //       dz, ...) need no translation: they are already valid expressions, and
  //       the lowering recognizes them. They sit here in the original order.

  // _over_dt, before the binary ops so that it is not read as X over dt.
  if (not (head=strip(name,"_over_dt")).empty()) return head + "/dt";

  if (not (head=strip(name,"_horiz_avg")).empty()) return head + ".horiz_avg()";

  // Vertical contraction, with optional weighting.
  for (const auto& method : {std::string("avg"),std::string("sum")}) {
    for (const auto& w : {std::string("dp"),std::string("dz")}) {
      if (not (head=strip(name,"_vert_" + method + "_" + w + "_weighted")).empty()) {
        return head + ".vert_" + method + "(weight='" + w + "')";
      }
    }
    if (not (head=strip(name,"_vert_" + method)).empty()) {
      return head + ".vert_" + method + "()";
    }
  }

  for (const auto& dx : {std::string("p"),std::string("z")}) {
    if (not (head=strip(name,"_" + dx + "vert_derivative")).empty()) {
      return head + ".derivative(dx='" + dx + "')";
    }
  }

  if (not (out=translate_zonal_avg(name)).empty()) return out;
  if (not (out=translate_where(name)).empty())     return out;
  if (not (out=translate_binary_op(name)).empty()) return out;

  // _prev, after the binary ops so that X_minus_X_prev is a subtraction.
  if (not (head=strip(name,"_prev")).empty()) return head + ".prev()";

  if (not (out=translate_histogram(name)).empty()) return out;

  return name;
}

} // namespace scream
