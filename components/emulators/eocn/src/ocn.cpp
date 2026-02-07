/**
 * @file ocn.cpp
 * @brief Implementation of the Ocn emulator class.
 *
 * Grid-agnostic: decomposition is provided externally.
 */

#include "ocn.hpp"

#include <algorithm>
#include <iostream>

namespace emulator {
namespace ocn {

// ---------- construction ----------
Ocn::Ocn(int fcomm, int compid)
    : Emulator(EmulatorType::OCN_COMP, fcomm, compid, "eocn") {}

// ---------- set decomposition from MCT ----------
void Ocn::set_decomposition(int lsize, const int *gindex) {
  m_lsize = lsize;
  m_gindex.assign(gindex, gindex + lsize);
}

// ---------- dummy decomposition for testing ----------
void Ocn::set_dummy_decomposition(int global_size) {
  // Simple contiguous 1-D decomposition
  int base = global_size / m_nprocs;
  int extra = global_size % m_nprocs;
  m_lsize = base + (m_rank < extra ? 1 : 0);
  int offset = m_rank * base + std::min(m_rank, extra);

  m_gindex.resize(m_lsize);
  for (int n = 0; n < m_lsize; ++n) {
    m_gindex[n] = offset + n + 1; // 1-based for MCT
  }
}

// ---------- Emulator lifecycle ----------
void Ocn::init_impl() {
  if (m_rank == 0) {
    std::cout << "(eocn) Ocn::init_impl  lsize=" << m_lsize << std::endl;
  }
}

void Ocn::run_impl(int dt) {
  // Skeleton: no-op. Future work fills o2x fields here.
  (void)dt;
}

void Ocn::final_impl() {
  if (m_rank == 0) {
    std::cout << "(eocn) Ocn::final_impl  steps=" << step_count() << std::endl;
  }
}

} // namespace ocn
} // namespace emulator
