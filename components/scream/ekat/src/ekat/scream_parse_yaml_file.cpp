#include "ekat/scream_parse_yaml_file.hpp"
#include "ekat/scream_assert.hpp"
#include "ekat/scream_assert.hpp"

#include <yaml-cpp/yaml.h>

#include <iostream>

namespace scream {

template<YAML::NodeType::value Type>
void parse_node (const YAML::Node& node,
                 const std::string& key,
                 ParameterList& list);

// Helpers
bool is_bool (const std::string& key) {
  return key=="true" || key=="false" ||
         key=="TRUE" || key=="FALSE";
}

bool str2bool (const std::string& key) {
  return (key=="true" || key=="TRUE");
}

bool is_int (const std::string& key) {
  const auto& f = std::use_facet<std::ctype<char>>(std::locale());
  constexpr auto digit = std::ctype_base::digit;
  bool first_char = true;
  for (auto ch : key) {
    if (ch=='+' || ch=='-') {
      // A sign is only allowed at the beginning
      if (!first_char) {
        return false;
      }
    } else if (!f.is(digit,ch)) {
      return false;
    }
    first_char = false;
  }
  return true;
}

bool is_double (const std::string& key) {
  const auto& f = std::use_facet<std::ctype<char>>(std::locale());
  constexpr auto digit = std::ctype_base::digit;
  bool point_found = false;
  bool exp_found = false;
  bool exp_symbol_found = false;
  bool exp_symbol_just_found = false;
  bool first_char = true;
  for (auto ch : key) {
    if (ch=='+' || ch=='-') {
      // Ok to add +/- in front of a number
      // Ok to add +/- after 'e'/'E', for the exponent
      // Any other case is wrong
      if (!first_char || !exp_symbol_just_found) {
        return false;
      }
      exp_symbol_just_found = false;
    } else if (ch=='.') {
      if(point_found || exp_symbol_found) {
        // Only one decimal point allowed, and it cannot be in the exponent
        return false;
      }
      // Mark that we found the decimal point
      point_found = true;
      exp_symbol_just_found = false;
    } else if (ch=='e' || ch=='E') {
      // Mark that we found the exp letter (e or E)
      exp_symbol_found = true;
      exp_symbol_just_found = true;
    } else if (f.is(digit,ch)) {
      // If e/E was found, we found the exponent
      exp_found = exp_symbol_found;
      exp_symbol_just_found = false;
    } else {
      return false;
    }
    first_char = false;
  }

  // If we found 'e' or 'E', we need to have found an exponent
  return exp_found==exp_symbol_found;
}

// ---------- IMPLEMENTATION -------------- // 

template<>
void parse_node<YAML::NodeType::Scalar> (
    const YAML::Node& node,
    const std::string& key,
    ParameterList& list)
{
  scream_require_msg (node.Type()==YAML::NodeType::Scalar,
                      "Error! Actual node type incompatible with template parameter.\n");

  // Extract scalar as string, then try some casts
  // First int, then double, and finally fall back to string
  std::string str = node.as<std::string>();
  if (is_int(str)) {
    list.set(key,std::stoi(str));
  } else if (is_bool(str)) {
    list.set(key,str2bool(str));
  } else if (is_double(str)) {
    list.set(key,std::stod(str));
  } else {
    list.set(key,str);
  }
}

template<>
void parse_node<YAML::NodeType::Sequence> (
    const YAML::Node& node,
    const std::string& key,
    ParameterList& list)
{
  scream_require_msg (node.Type()==YAML::NodeType::Sequence,
                      "Error! Actual node type incompatible with template parameter.\n");

  constexpr int str_type = 0;
  constexpr int dbl_type = 1;
  constexpr int int_type = 2;
  constexpr int bool_type = 3;
  int seq_type = -1;
  int n = node.size();
  for (int i=0; i<n; ++i) {
    std::string str = node[i].as<std::string>();
    int ith_type = is_int(str)
                     ? int_type
                     : (is_double(str)
                          ? dbl_type
                          : (is_bool(str) ? bool_type : str_type));

    scream_require_msg(seq_type==-1 || seq_type==ith_type,
                       "Error! Found a squence with entries of mixed types.\n"
                       "       Common case is a sequence with int and doubles, such as\n"
                       "          [ 2.3 1 ]\n");

    seq_type = ith_type;
  }

  if (seq_type==int_type) {
    std::vector<int> vec(n);
    for (int i=0; i<n; ++i) {
      std::string str = node[i].as<std::string>();
      vec[i] = std::stoi(str);
    }
    list.set(key,vec);
  } else if (seq_type==bool_type) {
    // NOTE: std::vector<bool> is a BAD thing to use.
    std::vector<int> vec(n);
    for (int i=0; i<n; ++i) {
      std::string str = node[i].as<std::string>();
      vec[i] = str2bool(str) ? 1 : 0;
    }
    list.set(key,vec);
  } else if (seq_type==dbl_type) {
    std::vector<double> vec(n);
    for (int i=0; i<n; ++i) {
      std::string str = node[i].as<std::string>();
      vec[i] = std::stod(str);
    }
    list.set(key,vec);
  } else {
    std::vector<std::string> vec(n);
    for (int i=0; i<n; ++i) {
      vec[i] = node[i].as<std::string>();
    }
    list.set(key,vec);
  }
}


template<>
void parse_node<YAML::NodeType::Map> (
    const YAML::Node& node,
    const std::string& key,
    ParameterList& list)
{
  using YNT = YAML::NodeType;
  scream_require_msg (node.Type()==YNT::Map,
                      "Error! Actual node type incompatible with template parameter.\n");

  ParameterList& sublist = list.sublist(key);
  for (auto it : node) {
    std::string item_key = it.first.as<std::string>();
    const auto& item = it.second;

    switch (it.second.Type()) {
      case YNT::Null:
        printf("Null node\n");
        break;
      case YNT::Scalar:
        parse_node<YNT::Scalar>(item, item_key, sublist);
        break;
      case YNT::Undefined:
        printf("Undefined node\n");
        break;
      case YNT::Map:
        parse_node<YNT::Map>(item, item_key, sublist);
        break ;
      case YAML::NodeType::Sequence:
        parse_node<YNT::Sequence>(item, item_key, sublist);
        break;
      default:
        printf("Unexpected node type\n");
    }
  }
}

ParameterList parse_yaml_file (const std::string& fname) {
  ParameterList params;
  parse_yaml_file(fname,params);
  return params;
}

void parse_yaml_file (const std::string& fname, ParameterList& params) {
  YAML::Node root;
  try {
    root = YAML::LoadFile(fname);
  } catch (YAML::BadFile&) {
    scream_error_msg ("Error! Something went wrong while opening file " + fname + "'.\n");
  }
  ParameterList temp(params.name());
  parse_node<YAML::NodeType::Map> (root, temp.name(), temp);
  params = temp.sublist(params.name());
}

} // namespace scream

