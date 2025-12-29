/**
 * @file atm_field_data_provider.hpp
 * @brief FieldDataProvider adapter for AtmFieldManager.
 *
 * Provides the FieldDataProvider interface implementation that wraps
 * AtmFieldManager for use with the EmulatorOutputManager.
 */

#ifndef ATM_FIELD_DATA_PROVIDER_HPP
#define ATM_FIELD_DATA_PROVIDER_HPP

#include "../../../common/src/emulator_output_stream.hpp"
#include "atm_field_manager.hpp"
#include <regex>
#include <set>
#include <string>
#include <vector>

namespace emulator {
namespace impl {

/**
 * @brief Adapter implementing FieldDataProvider for AtmFieldManager.
 *
 * Wraps an AtmFieldManager instance and provides the FieldDataProvider
 * interface required by EmulatorOutputStream for diagnostic output.
 *
 * ## Features
 *
 * - Field access by name via AtmFieldManager::get_field_ptr()
 * - Field stacking for sliced variables (e.g., wind_0, wind_1 → wind_3d)
 * - Tracking of available field names
 *
 * ## Field Stacking
 * Detects sliced field patterns like "field_N" (where N is an integer)
 * and can combine them into stacked 3D fields with shape [ncols, nlevs].
 */
class AtmFieldDataProvider : public FieldDataProvider {
public:
  /**
   * @brief Construct adapter with reference to field manager.
   * @param fields Reference to AtmFieldManager
   * @param ncols_local Number of local columns
   */
  AtmFieldDataProvider(AtmFieldManager &fields, int ncols_local);

  ~AtmFieldDataProvider() override = default;

  /**
   * @brief Get pointer to field data by name.
   *
   * First checks for exact match, then checks for stacked fields.
   *
   * @param name Field name
   * @return Pointer to data vector, or nullptr if not found
   */
  const std::vector<double> *get_field(const std::string &name) const override;

  /**
   * @brief Get list of all available field names.
   *
   * Includes both direct fields and auto-detected stacked fields.
   */
  std::vector<std::string> get_field_names() const override;

  /**
   * @brief Get number of local columns.
   */
  int get_ncols() const override { return m_ncols; }

  /**
   * @brief Get number of vertical levels for a field.
   * @param name Field name
   * @return Number of levels (1 for 2D fields, >1 for stacked fields)
   */
  int get_field_nlevs(const std::string &name) const override;

  /**
   * @brief Scan fields and detect stackable slice patterns.
   *
   * Finds fields matching "basename_N" pattern and registers them
   * as stackable groups.
   */
  void detect_stacked_fields();

  /**
   * @brief Check if a field is a stacked field.
   */
  bool is_stacked_field(const std::string &name) const;

  /**
   * @brief Get the stacked data for a 3D field.
   *
   * Combines sliced fields [field_0, field_1, ...] into a single
   * contiguous buffer with shape [nlevs, ncols].
   *
   * @param basename Base field name (without _N suffix)
   * @return Stacked data buffer
   */
  const std::vector<double> &
  get_stacked_field(const std::string &basename) const;

private:
  AtmFieldManager &m_fields;
  int m_ncols;

  // Stacked field detection and caching
  // Maps basename → list of level indices found
  mutable std::map<std::string, std::vector<int>> m_stacked_field_levels;

  // Cache for stacked field data
  mutable std::map<std::string, std::vector<double>> m_stacked_cache;

  // Set of all known field names (for get_field_names)
  mutable std::set<std::string> m_all_field_names;
  mutable bool m_field_names_cached = false;

  /**
   * @brief Parse a field name for slice pattern.
   * @param name Full field name (e.g., "wind_0")
   * @param[out] basename Base name without suffix (e.g., "wind")
   * @param[out] level_idx Level index
   * @return true if pattern matched
   */
  bool parse_slice_pattern(const std::string &name, std::string &basename,
                           int &level_idx) const;

  /**
   * @brief Build stacked field cache for a given basename.
   */
  void build_stacked_cache(const std::string &basename) const;
};

} // namespace impl
} // namespace emulator

#endif // ATM_FIELD_DATA_PROVIDER_HPP
