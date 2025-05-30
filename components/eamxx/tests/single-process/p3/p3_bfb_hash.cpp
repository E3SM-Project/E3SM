#include <catch2/catch.hpp>

#include <ekat_test_utils.hpp>
#include <ekat_assert.hpp>

#include <fstream>

namespace scream {

TEST_CASE("check_bfb_sha") {

  auto& ts = ekat::TestSession::get();

  auto log_file = ts.params.at("log-file");
  std::ifstream ifile(log_file);

  EKAT_REQUIRE_MSG (ifile.is_open(),
      "Error! Could not open log file '" << log_file << ".\n");

  std::string line;
  std::vector<std::string> collected_lines;
  while (std::getline(ifile, line)) {
    // Check if the line starts with "eamxx hash>" and ends with "naccum=N"
    if (line.find("eamxx hash>") == std::string::npos)
      continue;

    // Find the position of "naccum="
    size_t pos = line.find("naccum=");
    EKAT_REQUIRE_MSG (pos!=std::string::npos,
        "Error! Could not locate 'naccum' in the 'eamxx hash>' line.\n"
        "  - line: " + line + "\n");

    // Extract the value of N
    std::string naccum = line.substr(pos + 7); // Skip "naccum="
    int N = std::stoi(naccum); // Convert to integer

    // Store the current line
    collected_lines.push_back(line);

    // Read the next N lines
    for (int i = 0; i < N; ++i) {
      if (std::getline(ifile, line)) {
        collected_lines.push_back(line);
      } else {
        EKAT_ERROR_MSG (
            "Error! Could not find the declared number of accum lines.\n"
            " - expected accum lines: " << N << "\n"
            " - actual accum lines  : " << i << "\n");
      }
    }
  }

  auto output_fname = ts.params.at("output-file");
  std::ofstream ofile (output_fname);

  for (const auto& l : collected_lines) {
    ofile << l << "\n";
  }
}

} // e namespace
