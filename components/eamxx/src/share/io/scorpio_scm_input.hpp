#ifndef SCREAM_SCORPIO_SCM_INPUT_HPP
#define SCREAM_SCORPIO_SCM_INPUT_HPP

#include "share/grid/abstract_grid.hpp"

#include <ekat/logging/ekat_logger.hpp>

namespace scream
{

// Similar to AtmosphereInput, but reads in a single column from
// a file with N columns. A few assumptions:
//  - lat and lon variables are present in the file
//  - fields have layout <COL [, ...]>
class SCMInput
{
public:
  // --- Constructor(s) & Destructor --- //
  SCMInput (const std::string& filename,
            const double lat, const double lon,
            const std::vector<Field>& fields,
            const ekat::Comm& comm);

  ~SCMInput ();

  // Due to resource acquisition (in scorpio), avoid copies
  SCMInput (const SCMInput&) = delete;
  SCMInput& operator= (const SCMInput&) = delete;

  // Read fields that were required via parameter list.
  void read_variables (const int time_index = -1);

  // Option to add a logger
  void set_logger(const std::shared_ptr<ekat::logger::LoggerBase>& atm_logger) {
      m_atm_logger = atm_logger;
  }

#ifndef KOKKOS_ENABLE_CUDA
  // Cuda requires methods enclosing __device__ lambda's to be public
protected:
#endif
  void create_closest_col_info (double target_lat, double target_lon);
protected:

  struct ClosestColInfo {
    // MPI rank which owns the columns whose lat/lon pair is the closest to target lat/lon
    int mpi_rank;
    // Local column index of on rank=mpi_rank (-1 on all other ranks)
    int col_lid;
  };

  void create_io_grid ();
  void init_scorpio_structures ();
  void set_decompositions();

  // Internal variables
  ekat::Comm    m_comm;

  std::shared_ptr<AbstractGrid>   m_io_grid;

  std::string               m_filename;
  std::vector<Field>        m_fields;
  std::vector<Field>        m_io_fields;

  ClosestColInfo            m_closest_col_info;

  // The logger to be used throughout the ATM to log message
  std::shared_ptr<ekat::logger::LoggerBase> m_atm_logger;
};

} //namespace scream

#endif // SCREAM_SCORPIO_SCM_INPUT_HPP
