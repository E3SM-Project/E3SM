/**
 * @file eatm_atm_emulator.cpp
 * @brief Implementation of the AtmEmulator class.
 *
 * compute_decomposition() mirrors the logic in
 * components/xcpl_comps/xshare/dead_mod.F90::dead_setNewGrid
 * (1-D contiguous decomposition by latitude rows).
 */

#include "eatm_atm_emulator.hpp"

#include <cmath>
#include <iostream>
#include <mpi.h>

namespace eatm {

// ---------- constants matching dead_mod.F90 ----------
static constexpr double pi      = 3.14159265358979323846;
static constexpr double deg2rad = pi / 180.0;
static constexpr double re      = 6.37122e6;   // Earth radius (m)

// ---------- construction ----------
AtmEmulator::AtmEmulator(int fcomm, int compid, int nxg, int nyg)
    : emulator::Emulator(emulator::EmulatorType::ATM_COMP,
                         compid, "eatm"),
      m_fcomm(fcomm),
      m_comp_id(compid),
      m_nxg(nxg),
      m_nyg(nyg),
      m_lsize(0) {}

// ---------- grid decomposition ----------
void AtmEmulator::compute_decomposition() {
  // Convert Fortran communicator → C communicator
  MPI_Comm comm = MPI_Comm_f2c(m_fcomm);
  int rank, nprocs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &nprocs);

  const int gsize = m_nxg * m_nyg;

  // Simple contiguous 1-D decomposition
  int base  = gsize / nprocs;
  int extra = gsize % nprocs;
  m_lsize   = base + (rank < extra ? 1 : 0);
  int offset = rank * base + std::min(rank, extra);

  // Resize vectors
  m_gindex.resize(m_lsize);
  m_lons.resize(m_lsize);
  m_lats.resize(m_lsize);
  m_areas.resize(m_lsize);
  m_masks.resize(m_lsize);
  m_fracs.resize(m_lsize);

  // Fill synthetic grid data (reproduces dead_setNewGrid logic)
  const double dx = 360.0 / m_nxg * deg2rad;

  for (int n = 0; n < m_lsize; ++n) {
    int g  = offset + n;            // 0-based global index
    int ig = g % m_nxg;             // 0-based column
    int jg = g / m_nxg;             // 0-based row

    double ys = -90.0 + static_cast<double>(jg)     * 180.0 / m_nyg;
    double yc = -90.0 + (static_cast<double>(jg) + 0.5) * 180.0 / m_nyg;
    double yn = -90.0 + (static_cast<double>(jg) + 1.0) * 180.0 / m_nyg;
    double dy = std::sin(yn * deg2rad) - std::sin(ys * deg2rad);

    m_gindex[n] = g + 1;                       // 1-based for MCT
    m_lons[n]   = static_cast<double>(ig) * 360.0 / m_nxg;
    m_lats[n]   = yc;
    m_areas[n]  = dx * dy * re * re;
    m_masks[n]  = 1.0;
    m_fracs[n]  = 1.0;
  }
}

// ---------- Emulator lifecycle stubs ----------
void AtmEmulator::init_impl() {
  MPI_Comm comm = MPI_Comm_f2c(m_fcomm);
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    std::cout << "(eatm) AtmEmulator::init_impl  "
              << "grid " << m_nxg << "x" << m_nyg
              << "  lsize=" << m_lsize << std::endl;
  }
}

void AtmEmulator::run_impl(int dt) {
  // Skeleton: no-op.  Future work fills a2x fields here.
  (void)dt;
}

void AtmEmulator::final_impl() {
  MPI_Comm comm = MPI_Comm_f2c(m_fcomm);
  int rank;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    std::cout << "(eatm) AtmEmulator::final_impl  "
              << "steps=" << step_count() << std::endl;
  }
}

} // namespace eatm
