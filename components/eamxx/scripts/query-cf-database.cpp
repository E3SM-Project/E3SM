#include "ekat/util/ekat_string_utils.hpp"
#include "ekat/ekat_assert.hpp"

#include <yaml-cpp/yaml.h>
#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <map>

#define STRINGIFY2(X) #X
#define STRINGIFY(X) STRINGIFY2(X)

#ifndef CF_STANDARD_NAME_FILE
#error "Error! This file needs the macro CF_STANDARD_NAME_FILE to be defined\n"
#endif
#ifndef CF_SCREAM_NAME_FILE
#error "Error! This file needs the macro CF_SCREAM_NAME_FILE to be defined\n"
#endif

// Default string similarity threshold.
static const double SIMILARY_THRESH = 0.6;

// Default max number of results.
static const size_t MAX_RESULTS = 5;

enum Algo {
  Jaccard = 0,
  Jaro    = 1,
  All     = 2
};

void usage(const char* prog) {
  fprintf(stderr, "%s: queries a CF database for the name of a field.\n", prog);
  fprintf(stderr, "Usage:\n");
  fprintf(stderr, "%s [options] <field_name>\n\n", prog);
  fprintf(stderr, "Options:\n");
  fprintf(stderr, "-a [name]   Use only the specified algorithm.\n");
  fprintf(stderr, "            Valid options: jaccard, jaro (default=jaccard)\n");
  fprintf(stderr, "-d          Retrieve the description of the given field\n");
  fprintf(stderr, "-w          Matches input string as a whole. E.g, passing\n");
  fprintf(stderr, "            -w my_field won't match \"field_my\", \"my\", or \"field\".\n");
  fprintf(stderr, "            NOTE: not usable with '-a jaro'.\n");
  fprintf(stderr, "-s [thresh] Report names similar to the given field,\n");
  fprintf(stderr, "            with a given string similarity index threshold\n");
  fprintf(stderr, "            that determines the results (default: %g)\n", SIMILARY_THRESH);
  fprintf(stderr, "-b          Report all similar results, regardless of any\n");
  fprintf(stderr, "            similarity threshold\n");
  fprintf(stderr, "-n N        Report at most N similar results (default: %zu)\n", MAX_RESULTS);
  fprintf(stderr, "-v          Verbose mode: show similarity indices for matches.\n");
  exit(1);
}

// This stores options for searching. See the usage information above for
// a description of its members.
struct SearchOptions {
  bool show_description;
  double similarity_threshold;
  bool ignore_threshold;
  bool whole_string_atomic;
  size_t max_results;
  bool verbose;
  Algo algorithm;


  // Default constructor
  SearchOptions():
    show_description(false), similarity_threshold(SIMILARY_THRESH),
    ignore_threshold(false), whole_string_atomic(false),
    max_results(MAX_RESULTS), verbose(false), algorithm(All) {}
};

// Parses options, populating the arguments with values.
void parse_options(int argc, char* argv[], std::string& field_name,
                   SearchOptions& options) {
  for (int i = 1; i < argc; ++i) {
    if (argv[i][0] == '-') { // argument is a flag
      char flag = argv[i][1];
      if (flag == 'd') {
        options.show_description = true;
      } else if (flag == 's') {
        if (i < argc-1) {
          options.similarity_threshold = atof(argv[i+1]);
          i++;
        } else {
          usage(argv[0]);
        }
      } else if (flag == 'a') {
        if (i < argc-1) {
          std::string algo_str = argv[i+1];
          if (algo_str=="jaccard") {
            options.algorithm = Jaccard;
          } else if (algo_str=="jaro") {
            options.algorithm = Jaro;
          } else {
            printf ("unrecognized algorithm: %s\n",argv[i+1]);
            usage(argv[0]);
          }
          i++;
        } else {
          usage(argv[0]);
        }
      } else if (flag == 'w') {
        options.whole_string_atomic = true;
      } else if (flag == 'b') {
        options.ignore_threshold = true;
      } else if (flag == 'n') {
        if (i < argc-1) {
          options.max_results = atoi(argv[i+1]);
          i++;
        } else {
          usage(argv[0]);
        }
      } else if (flag == 'v') {
        options.verbose = true;
      } else {
        fprintf(stderr, "Unrecognized flag: %s\n", argv[i]);
        usage(argv[0]);
      }
    } else {
      field_name = std::string(argv[i]);
    }
  }

  // --- Fix/check options --- //

  if (options.ignore_threshold) {
    options.similarity_threshold = 0.0;
  }

  EKAT_REQUIRE_MSG (not options.whole_string_atomic or options.algorithm!=Jaro,
      "Error! The -w flag can only be used with the Jaccard algorithm.\n");

  // Did we get everything we needed?
  if (field_name.empty()) {
    usage(argv[0]);
  }
}

// Reads a YAML file containing a CF field name database, returning a
// map of field names to their descriptions.
std::map<std::string, std::string> read_cf_database(const char* filename) {
  std::map<std::string, std::string> db;

  // Try to load the file.
  try {
    auto root = YAML::LoadFile(filename);

    // Cool--it worked. Let's read the database into memory.
    for (auto iter = root.begin(); iter != root.end(); ++iter) {
      auto name = iter->first.as<std::string>();
      auto field = iter->second;
      // Fetch the field's description if it has one.
      std::string description;
      if (field["description"]) {
        description = field["description"].as<std::string>();
      }
      db[name] = description;

      // Handle aliases.
      if (field["aliases"]) {
        if (field["aliases"].IsSequence()) { // > 1 alias
          auto aliases = field["aliases"].as<std::vector<std::string> >();
          for (const auto& alias: aliases) {
            db[alias] = description;
          }
        } else if (field["aliases"].IsScalar()) { // only 1 alias
          // Aliases in the cf table can be listed as a "scalar" entry, as
          //  aliases: alias1 alias2 alias3
          // So split the entry using ' ' as a delimiter
          auto aliases = ekat::gather_tokens(field["aliases"].as<std::string>(),{' '});
          for (const auto& alias: aliases) {
            db[alias] = description;
          }
        }
      }
    }
  }
  catch (YAML::BadFile& e) {
    fprintf(stderr, "Error loading CF database %s:\n %s\n", filename, e.what());
    exit(1);
  }
  catch (YAML::ParserException& e) {
    fprintf(stderr, "Error parsing CF database %s:\n %s\n", filename, e.what());
    exit(1);
  }
  return db;
}

// This stores a matching result for the field name search.
struct SearchResult {
  std::string field_name;
  std::string description;
  double index;
  std::string index_type;
};

// Searches the given database for the given field name, subject to
// criteria defined in the search options. Matching results are appended
// to the results vector.
void search_jaccard(const std::map<std::string, std::string>& database,
                    const std::string& field_name,
                    const SearchOptions& options,
                    std::vector<SearchResult>& results) {
  if (options.verbose) {
    fprintf(stdout, "Searching for results similar to \"%s\" (similarity > %g) using jaccard algorithm.\n",
            field_name.c_str(), options.similarity_threshold);
  }

  for (auto iter = database.begin(); iter != database.end(); ++iter) {
    auto db_field_name = iter->first;
    auto db_field_desc = iter->second;

    // First compute the Jaccard similarity index, which is good at matching
    // phrases with the same words in different orders.
    double jaccard_index = ekat::jaccard_similarity(field_name, db_field_name,
                                                    {' ', '_'},
                                                    not options.whole_string_atomic);
    if ((options.verbose) and (jaccard_index > 0.25*options.similarity_threshold)) {
      fprintf(stdout, "  Jaccard similarity (\"%s\", \"%s\") = %g\n",
              field_name.c_str(), db_field_name.c_str(), jaccard_index);
    }
    if (jaccard_index >= options.similarity_threshold) {
      results.push_back({db_field_name, db_field_desc, jaccard_index, "jaccard"});
    }
  }
}

void search_jaro(const std::map<std::string, std::string>& database,
                 const std::string& field_name,
                 const SearchOptions& options,
                  std::vector<SearchResult>& results) {
  if (options.verbose) {
    fprintf(stdout, "Searching for results similar to \"%s\" (similarity > %g) using jaro algorithm.\n",
            field_name.c_str(), options.similarity_threshold);
  }

  for (auto iter = database.begin(); iter != database.end(); ++iter) {
    auto db_field_name = iter->first;
    auto db_field_desc = iter->second;

    double jaro_index = ekat::jaro_similarity(field_name, db_field_name);
    if ((options.verbose) and (jaro_index > options.similarity_threshold)) {
      fprintf(stdout, "  Jaro similarity (\"%s\", \"%s\") = %g\n",
              field_name.c_str(), db_field_name.c_str(), jaro_index);
    }
    if (jaro_index >= options.similarity_threshold) {
      results.push_back({db_field_name, db_field_desc, jaro_index, "jaro"});
    }
  }
}

void search_database(const std::map<std::string, std::string>& database,
                     const std::string& field_name,
                     const SearchOptions& options,
                     std::vector<SearchResult>& results) {
  if (options.verbose) {
    fprintf(stdout, "Searching for results similar to \"%s\" (similarity > %g)\n",
            field_name.c_str(), options.similarity_threshold);
  }

  switch (options.algorithm) {
    case Jaccard:
      search_jaccard(database,field_name,options,results);
      break;
    case Jaro:
      search_jaro(database,field_name,options,results);
      break;
    case All:
      search_jaccard(database,field_name,options,results);
      if (results.size()<options.max_results and not options.whole_string_atomic) {
        fprintf(stdout, "  ** Note: only %zu results found with Jaccard similarity test (max_results=%zu).\n"
                        "           Falling back to Jaro similarity test to find more.\n\n",
                        results.size(),options.max_results);
        search_jaro(database,field_name,options,results);
      }
      break;
  }

  // Now sort the results from best to worst matches.
  std::sort(results.begin(), results.end(),
            [](const SearchResult& a, const SearchResult& b) {
              if (a.index_type==b.index_type) {
                return a.index > b.index;
              } else {
                return a.index_type=="jaccard";
              }
            });

  if (options.verbose) {
    fprintf(stdout, "Total results found: %ld\n\n", results.size());
  }
}

// Prints matching results for the given field name, possibly filtered
// by options.
void report_results(const std::string& field_name,
                    const std::vector<SearchResult>& results,
                    const SearchOptions& options,
                    FILE* output_stream) {
  fprintf(output_stream, "   Field name searched: %s\n", field_name.c_str());

  // If tolerance is too strict,
  if (results.empty()) {
    fprintf(output_stream,
        " ** No matches found for the provided criteria. Possible reasons:\n"
        "\n"
        "  - the tolerance was too high. Try reducing it, using the -s flag.\n"
        "  - you used the -w flag. The script will NOT fall back to the Jaro algirithm,\n"
        "    in the case where no matches are found with the Jaccard algorithm.\n"
        "    To allow fall back to Jaro algorithm, try droping the -w flag (if present).\n"
        "\n");
    return;
  }

  fprintf(output_stream, " +--------------------------------------------------+\n");
  fprintf(output_stream, " | Best matches [name, similarity index, algorithm] |\n");
  fprintf(output_stream, " +--------------------------------------------------+\n");
  unsigned num_reported = 0;
  for (const auto& result: results) {
    fprintf(output_stream, "   %-80s %f, %15s\n", (result.field_name+",").c_str(),result.index, result.index_type.c_str());

    if (options.show_description) {
      fprintf(output_stream, "      %s\n", result.description.c_str());
    }
    num_reported++;
    if (num_reported >= options.max_results) {
      break;
    }
  }
}

int main(int argc, char* argv[]) {

  if (argc < 2) {
    usage(argv[0]);
  }

  // Handle command-line options.
  std::string field_name;
  SearchOptions options;
  parse_options(argc, argv, field_name, options);

  // Read the CF field names into a pair of dictionaries that map them
  // to their descriptions.

  // Location of CF database YAML file.
  const char* cf_standard_name_file = STRINGIFY(CF_STANDARD_NAME_FILE);
  const char* cf_scream_name_file   = STRINGIFY(CF_SCREAM_NAME_FILE);
  auto cf_std_fields = read_cf_database(cf_standard_name_file);
  auto cf_scr_fields = read_cf_database(cf_scream_name_file);

  // Search each dictionary for the given field name, building a set of
  // results.
  std::vector<SearchResult> results;
  search_database(cf_std_fields, field_name, options, results);
  if (results.size()<options.max_results) {
    fprintf(stdout, "  ** Note: only %zu results found in the cf-standard database (max_results=%zu).\n"
                    "           Falling back to scream database to find more.\n\n",
                    results.size(),options.max_results);
    search_database(cf_scr_fields, field_name, options, results);
  }
  // Report back.
  report_results(field_name, results, options, stdout);

  return 0;
}
