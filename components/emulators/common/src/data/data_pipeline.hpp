/**
 * @file data_pipeline.hpp
 * @brief Bidirectional data pipeline orchestrator for emulators.
 *
 * DataPipeline composes DataSource → DataView → TensorAdapter / IOAdapter
 * into a convenient, optional orchestration layer. Emulators can use it
 * for a turnkey data flow, or skip it and compose the pieces manually.
 */

#ifndef E3SM_EMULATOR_DATA_PIPELINE_HPP
#define E3SM_EMULATOR_DATA_PIPELINE_HPP

#include "data_source.hpp"
#include "data_view.hpp"
#include "field_spec.hpp"
#include "io_adapter.hpp"
#include "tensor_adapter.hpp"

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace emulator {

/**
 * @brief Bidirectional data pipeline for emulator data flow.
 *
 * ## Lifecycle
 * 1. Construct with a name and npoints.
 * 2. Optionally set sources, adapters via set_*() methods.
 * 3. Register fields via add_import_field() / add_export_field(),
 *    or call auto_register_fields() to pull from the DataSource.
 * 4. Call initialize() to allocate DataViews.
 * 5. During run: call ingest() → process → egest().
 *
 * ## Minimal use (ML-only, no IO, no MCT)
 * ```cpp
 * DataPipeline pipe("ml_only", npoints);
 * pipe.add_import_field(FieldSpec("x"));
 * pipe.add_export_field(FieldSpec("y"));
 * pipe.set_tensor_adapter(std::make_unique<RawTensorAdapter>());
 * pipe.initialize();
 * // Fill import view manually, read tensor, etc.
 * ```
 *
 * ## Full use (MCT → ML → MCT)
 * ```cpp
 * DataPipeline pipe("atm", npoints);
 * pipe.set_import_source(std::make_unique<MctDataSource>(...));
 * pipe.set_export_source(std::make_unique<MctDataSource>(...));
 * pipe.set_tensor_adapter(std::make_unique<RawTensorAdapter>());
 * pipe.auto_register_fields();
 * pipe.initialize();
 *
 * // Each timestep:
 * pipe.ingest();                    // MCT → import_view
 * const double* in = pipe.tensor_data();
 * // ... run ML inference ...
 * pipe.write_tensor(ml_output);     // ML → export_view
 * pipe.egest();                     // export_view → MCT
 * ```
 */
class DataPipeline {
public:
  /**
   * @brief Construct a DataPipeline.
   *
   * @param name    Pipeline name (used as prefix for DataView names)
   * @param npoints Number of local grid points
   * @param layout  Internal DataView layout (default: POINT_MAJOR)
   */
  DataPipeline(std::string name, int npoints,
               DataLayout layout = DataLayout::POINT_MAJOR);

  ~DataPipeline();

  // Non-copyable, movable
  DataPipeline(const DataPipeline &) = delete;
  DataPipeline &operator=(const DataPipeline &) = delete;
  DataPipeline(DataPipeline &&) noexcept;
  DataPipeline &operator=(DataPipeline &&) noexcept;

  // ── Assembly (call during init) ─────────────────────────────────

  /** @brief Set the import data source (external → import_view). */
  void set_import_source(std::unique_ptr<DataSource> source);

  /** @brief Set the export data source (export_view → external). */
  void set_export_source(std::unique_ptr<DataSource> source);

  /** @brief Set the tensor adapter for ML data access. */
  void set_tensor_adapter(std::unique_ptr<TensorAdapter> adapter);

  /** @brief Set the IO adapter for file read/write. */
  void set_io_adapter(std::unique_ptr<IOAdapter> io);

  /** @brief Register a field on the import DataView. */
  void add_import_field(const FieldSpec &spec);

  /** @brief Register a field on the export DataView. */
  void add_export_field(const FieldSpec &spec);

  /**
   * @brief Auto-register fields from the attached DataSources.
   *
   * If import_source is set, its available_fields() are added to
   * the import view. Same for export_source → export view.
   */
  void auto_register_fields();

  /**
   * @brief Allocate both DataViews.
   *
   * Must be called after field registration.
   * @throws std::runtime_error if no fields on either view, or
   *         if already initialized
   */
  void initialize();

  /** @brief Whether initialize() has been called. */
  bool is_initialized() const { return m_initialized; }

  // ── Runtime (call during run) ───────────────────────────────────

  /**
   * @brief Pull data from the import source into the import view.
   * @throws std::runtime_error if no import source is set
   */
  void ingest();

  /**
   * @brief Push data from the export view to the export source.
   * @throws std::runtime_error if no export source is set
   */
  void egest();

  // ── Tensor access ──────────────────────────────────────────────

  /**
   * @brief Zero-copy pointer to the import view's data for ML input.
   * @throws std::runtime_error if no tensor adapter is set
   */
  const double *tensor_data() const;

  /**
   * @brief Tensor shape of the import view.
   * @throws std::runtime_error if no tensor adapter is set
   */
  std::vector<int64_t> tensor_shape() const;

  /**
   * @brief Write ML output into the export view.
   *
   * @param ml_output Source buffer (from ML inference)
   * @throws std::runtime_error if no tensor adapter is set
   */
  void write_tensor(const void *ml_output);

  // ── IO access ──────────────────────────────────────────────────

  /**
   * @brief Write the import view to a file.
   * @throws std::runtime_error if no IO adapter is set
   */
  void write_io(const std::string &filename);

  /**
   * @brief Read data from a file into the import view.
   * @throws std::runtime_error if no IO adapter is set
   */
  void read_io(const std::string &filename);

  // ── Direct access ──────────────────────────────────────────────

  DataView &import_view();
  const DataView &import_view() const;

  DataView &export_view();
  const DataView &export_view() const;

  const std::string &name() const { return m_name; }

private:
  std::string m_name;
  int m_npoints;
  DataLayout m_layout;
  bool m_initialized = false;

  std::unique_ptr<DataView> m_import_view;
  std::unique_ptr<DataView> m_export_view;

  std::unique_ptr<DataSource> m_import_source;
  std::unique_ptr<DataSource> m_export_source;
  std::unique_ptr<TensorAdapter> m_tensor_adapter;
  std::unique_ptr<IOAdapter> m_io_adapter;

  // Staged fields (before initialize)
  std::vector<FieldSpec> m_import_fields;
  std::vector<FieldSpec> m_export_fields;
};

} // namespace emulator

#endif // E3SM_EMULATOR_DATA_PIPELINE_HPP
