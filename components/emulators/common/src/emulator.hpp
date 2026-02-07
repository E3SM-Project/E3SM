/**
 * @file emulator.hpp
 * @brief Abstract base class for all E3SM emulators.
 */

#ifndef EMULATOR_HPP
#define EMULATOR_HPP

#include <mpi.h>
#include <string>
#include <vector>

namespace emulator {

/**
 * @brief Enumeration of emulator types in E3SM.
 */
enum class EmulatorType {
  ATM_COMP = 0, ///< Atmosphere component emulator
  OCN_COMP = 1, ///< Ocean component emulator
  ICE_COMP = 2, ///< Sea ice component emulator
  LND_COMP = 3  ///< Land component emulator
};

/**
 * @brief Abstract base class for all E3SM emulators.
 *
 * Provides the common infrastructure for emulators.
 * Derived classes implement the pure virtual methods for
 * emulator-specific behavior. Grid-agnostic: grids provided
 * externally via MCT or test harness.
 */
class Emulator {
public:
  /**
   * @brief Construct a new Emulator.
   *
   * @param type Emulator type
   * @param fcomm Fortran MPI communicator handle
   * @param id Emulator ID (-1 if unassigned)
   * @param name Emulator name (empty if unassigned)
   */
  explicit Emulator(EmulatorType type, int fcomm = MPI_Comm_c2f(MPI_COMM_NULL),
                    int id = -1, const std::string &name = "");
  virtual ~Emulator() = default;

  // Lifecycle methods
  void initialize();
  void run(int dt);
  void finalize();

  // Accessors
  EmulatorType type() const { return m_type; }
  int id() const { return m_id; }
  const std::string &name() const { return m_name; }
  bool is_initialized() const { return m_initialized; }
  int step_count() const { return m_step_count; }

  // MPI accessors
  MPI_Comm comm() const { return m_comm; }
  int rank() const { return m_rank; }
  int nprocs() const { return m_nprocs; }

  // Local decomposition size (set by derived class)
  int lsize() const { return m_lsize; }
  const std::vector<int> &gindex() const { return m_gindex; }

protected:
  // Virtual methods for derived classes
  virtual void init_impl() = 0;
  virtual void run_impl(int dt) = 0;
  virtual void final_impl() = 0;

  EmulatorType m_type;
  MPI_Comm m_comm;
  int m_rank;
  int m_nprocs;
  int m_id;
  std::string m_name;
  bool m_initialized = false;
  int m_step_count = 0;

  // Grid-agnostic decomposition (set by derived class or MCT layer)
  int m_lsize = 0;
  std::vector<int> m_gindex;
};

} // namespace emulator

#endif // EMULATOR_HPP
