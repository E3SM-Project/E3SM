/********************************************************************************
 * HOMMEXX 1.0: Copyright of Sandia Corporation
 * This software is released under the BSD license
 * See the file 'COPYRIGHT' in the HOMMEXX/src/share/cxx directory
 *******************************************************************************/

#ifndef HOMMEXX_PPM_REMAP_HPP
#define HOMMEXX_PPM_REMAP_HPP

#include "ErrorDefs.hpp"

#include "KernelVariables.hpp"
#include "Types.hpp"
#include "ExecSpaceDefs.hpp"
#include "utilities/LoopsUtils.hpp"
#include "utilities/MathUtils.hpp"
#include "utilities/SubviewUtils.hpp"
#include "utilities/SyncUtils.hpp"

#include "profiling.hpp"

namespace Homme {
namespace Remap {
namespace Ppm {

namespace _ppm_consts {
// TODO: Hide these values from users
using Kokkos::Impl::MEMORY_ALIGNMENT;
// If sizeof(Real) doesn't divide the memory alignment value for the
// architecture, we're in trouble regardless
static constexpr int Real_Alignment =
    max(int(Kokkos::Impl::MEMORY_ALIGNMENT / sizeof(Real)), 1);
static constexpr int Vector_Alignment = max(Real_Alignment / VECTOR_SIZE, 1);

static constexpr int gs = 2;

// Padding to improve memory access alignment
static constexpr int INITIAL_PADDING =
    lcm(gs, int(VECTOR_SIZE), Real_Alignment);
static constexpr int VECTOR_PADDING = INITIAL_PADDING / VECTOR_SIZE;

// ghost cells, length 2, on both boundaries
static constexpr int DPO_PHYSICAL_LEV = NUM_PHYSICAL_LEV + INITIAL_PADDING + gs;
static constexpr int DPO_LEV = DPO_PHYSICAL_LEV / VECTOR_SIZE;

// cumulative integral of source, 0 start, with extra level as absolute maximum
static constexpr int PIO_PHYSICAL_LEV = NUM_PHYSICAL_LEV + 2;
static constexpr int PIO_LEV = PIO_PHYSICAL_LEV / VECTOR_SIZE;

// cumulative integral of target, 0 start
static constexpr int PIN_PHYSICAL_LEV = NUM_PHYSICAL_LEV + 1;
static constexpr int PIN_LEV = PIN_PHYSICAL_LEV / VECTOR_SIZE;

static constexpr int PPMDX_PHYSICAL_LEV = NUM_PHYSICAL_LEV + 2;
static constexpr int PPMDX_LEV = PPMDX_PHYSICAL_LEV / VECTOR_SIZE;

// ghost cells, length 2, on both boundaries
static constexpr int AO_PHYSICAL_LEV = NUM_PHYSICAL_LEV + INITIAL_PADDING + gs;
static constexpr int AO_LEV = AO_PHYSICAL_LEV / VECTOR_SIZE;

static constexpr int MASS_O_PHYSICAL_LEV = NUM_PHYSICAL_LEV + 2;
static constexpr int MASS_O_LEV = MASS_O_PHYSICAL_LEV / VECTOR_SIZE;

static constexpr int DMA_PHYSICAL_LEV = NUM_PHYSICAL_LEV + 2;
static constexpr int DMA_LEV = DMA_PHYSICAL_LEV / VECTOR_SIZE;

static constexpr int AI_PHYSICAL_LEV = NUM_PHYSICAL_LEV + 1;
static constexpr int AI_LEV = AI_PHYSICAL_LEV / VECTOR_SIZE;

} // namespace _ppm_consts

struct PpmBoundaryConditions {};

// Corresponds to remap alg = 1
struct PpmMirrored : public PpmBoundaryConditions {
  static constexpr int fortran_remap_alg = 1;

  KOKKOS_INLINE_FUNCTION
  static void apply_ppm_boundary(
      ExecViewUnmanaged<const Real[_ppm_consts::AO_PHYSICAL_LEV]> /* cell_means */,
      ExecViewUnmanaged<Real[3][NUM_PHYSICAL_LEV]> /* parabola_coeffs */)
  {
    // Nothing to do here
  }

  KOKKOS_INLINE_FUNCTION
  static void fill_cell_means_gs(
      KernelVariables &kv,
      const ExecViewUnmanaged<Real[_ppm_consts::DPO_PHYSICAL_LEV]>&,
      ExecViewUnmanaged<Real[_ppm_consts::AO_PHYSICAL_LEV]> cell_means) {
    const int gs = _ppm_consts::gs;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, gs),
                         [&](const int &k_0) {
      cell_means(_ppm_consts::INITIAL_PADDING - 1 - k_0 - 1 + 1) =
          cell_means(k_0 + _ppm_consts::INITIAL_PADDING);

      cell_means(NUM_PHYSICAL_LEV + _ppm_consts::INITIAL_PADDING - gs + k_0 +
                 1 + 1) =
          cell_means(NUM_PHYSICAL_LEV + _ppm_consts::INITIAL_PADDING - gs + 1 -
                     k_0 - 1 + 1);
    }); // end ghost cell loop
  }

  static constexpr const char *name() { return "Mirrored PPM"; }
};

// Corresponds to remap alg = 10
struct PpmLimitedExtrap : public PpmBoundaryConditions {
  static constexpr int fortran_remap_alg = 10;

  KOKKOS_INLINE_FUNCTION static void apply_ppm_boundary (
    const ExecViewUnmanaged<const Real[_ppm_consts::AO_PHYSICAL_LEV]>&,
    const ExecViewUnmanaged<Real[3][NUM_PHYSICAL_LEV]>&)
  {
    // Nothing to do here
  }

  KOKKOS_INLINE_FUNCTION static void linextrap (
    const Real& dx1, const Real& dx2, const Real& dx3, const Real& dx4,
    const Real& y1, const Real& y2, Real& y3, Real& y4,
    const Real& lo, const Real& hi)
  {
    const auto den = (dx1 + dx2)/2;
    auto num = den + (dx2 + dx3)/2;
    auto a = num/den;
    y3 = (1-a)*y1 + a*y2;

    num = num + (dx3 + dx4)/2;
    a = num/den;
    y4 = (1-a)*y1 + a*y2;

    y3 = max(lo, min(hi, y3));
    y4 = max(lo, min(hi, y4));
  }

  KOKKOS_INLINE_FUNCTION static void fill_cell_means_gs (
    KernelVariables& kv, const ExecViewUnmanaged<Real[_ppm_consts::DPO_PHYSICAL_LEV]>& dpo,
    const ExecViewUnmanaged<Real[_ppm_consts::AO_PHYSICAL_LEV]>& ao)
  {
    using Kokkos::parallel_reduce;
    using Kokkos::Min;
    using Kokkos::Max;

    const int ip = _ppm_consts::INITIAL_PADDING, plev = NUM_PHYSICAL_LEV;

    Real lo, hi; {
      const auto tvr = Kokkos::ThreadVectorRange(kv.team, plev);
      parallel_reduce(tvr, [&] (const int k, Real& lo) { lo = min(lo, ao(ip+k)); },
                      Min<Real>(lo));
      parallel_reduce(tvr, [&] (const int k, Real& hi) { hi = max(hi, ao(ip+k)); },
                      Max<Real>(hi));
    }

    linextrap(dpo(ip+1), dpo(ip), dpo(ip-1), dpo(ip-2),
              ao(ip+1), ao(ip), ao(ip-1), ao(ip-2),
              lo, hi);
    linextrap(dpo(ip+plev-2), dpo(ip+plev-1), dpo(ip+plev), dpo(ip+plev+1),
              ao(ip+plev-2), ao(ip+plev-1), ao(ip+plev), ao(ip+plev+1),
              lo, hi);
  }

  static constexpr const char* name () { return "PPM with limited extrapolation"; }
};

// Piecewise Parabolic Method stencil
template <typename boundaries>
struct PpmVertRemap : public VertRemapAlg {
  static_assert(std::is_base_of<PpmBoundaryConditions, boundaries>::value,
                "PpmVertRemap requires a valid PPM "
                "boundary condition");
  const int gs = _ppm_consts::gs;

  explicit PpmVertRemap(const int num_elems, const int num_remap)
      : m_dpo("dpo", num_elems)
      , m_pio("pio", num_elems)
      , m_pin("pin", num_elems)
      , m_ppmdx("ppmdx", num_elems)
      , m_z2("z2", num_elems)
      , m_kid("kid", num_elems)
      , m_ppm_tu(get_default_team_policy<ExecSpace>(num_elems * num_remap))
      , m_ao("a0", m_ppm_tu.get_num_ws_slots())
      , m_mass_o("mass_o",m_ppm_tu.get_num_ws_slots())
      , m_dma("dma", m_ppm_tu.get_num_ws_slots())
      , m_ai("ai", m_ppm_tu.get_num_ws_slots())
      , m_parabola_coeffs("Coefficients for the interpolating parabola", m_ppm_tu.get_num_ws_slots())
  {
    // Nothing to do here
  }

  KOKKOS_INLINE_FUNCTION
  void compute_grids_phase(
      KernelVariables &kv,
      ExecViewUnmanaged<const Scalar[NP][NP][NUM_LEV]> src_layer_thickness,
      ExecViewUnmanaged<const Scalar[NP][NP][NUM_LEV]> tgt_layer_thickness)
      const {
    compute_partitions(kv, src_layer_thickness, tgt_layer_thickness);
    compute_integral_bounds(kv);
  }

  KOKKOS_INLINE_FUNCTION
  void compute_remap_phase(KernelVariables &kv,
                           ExecViewUnmanaged<Scalar[NP][NP][NUM_LEV]> remap_var)
      const {
    // From here, we loop over tracers for only those portions which depend on
    // tracer data, which includes PPM limiting and mass accumulation
    // More parallelism than we need here, maybe break it up?
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                           [&](const int k) {
        const int ilevel = k / VECTOR_SIZE;
        const int ivector = k % VECTOR_SIZE;
        m_ao(kv.team_idx, igp, jgp, k + _ppm_consts::INITIAL_PADDING) =
            remap_var(igp, jgp, ilevel)[ivector] /
            m_dpo(kv.ie, igp, jgp, k + _ppm_consts::INITIAL_PADDING);
      });

      boundaries::fill_cell_means_gs(kv, Homme::subview(m_dpo, kv.ie, igp, jgp),
                                     Homme::subview(m_ao, kv.team_idx, igp, jgp));

      Dispatch<ExecSpace>::parallel_scan(
          kv.team, NUM_PHYSICAL_LEV,
          [=](const int &k, Real &accumulator, const bool last) {
            // Accumulate the old mass up to old grid cell interface locations
            // to simplify integration during remapping. Also, divide out the
            // grid spacing so we're working with actual tracer values and can
            // conserve mass.
            const int ilevel = k / VECTOR_SIZE;
            const int ivector = k % VECTOR_SIZE;
            accumulator += remap_var(igp, jgp, ilevel)[ivector];
            if (last) {
              m_mass_o(kv.team_idx, igp, jgp, k + 1) = accumulator;
            }
      });

      // Computes a monotonic and conservative PPM reconstruction
      compute_ppm(kv,
                  Homme::subview(m_ao, kv.team_idx, igp, jgp),
                  Homme::subview(m_ppmdx, kv.ie, igp, jgp),
                  Homme::subview(m_dma, kv.team_idx, igp, jgp),
                  Homme::subview(m_ai, kv.team_idx, igp, jgp),
                  Homme::subview(m_parabola_coeffs, kv.team_idx, igp, jgp));

      compute_remap(kv,
                    Homme::subview(m_kid, kv.ie, igp, jgp),
                    Homme::subview(m_z2, kv.ie, igp, jgp),
                    Homme::subview(m_parabola_coeffs, kv.team_idx, igp, jgp),
                    Homme::subview(m_mass_o, kv.team_idx, igp, jgp),
                    Homme::subview(m_dpo, kv.ie, igp, jgp),
                    Homme::subview(remap_var, igp, jgp));
    }); // End team thread range
    kv.team_barrier();
  }

  KOKKOS_FORCEINLINE_FUNCTION
  Real compute_mass(const Real sq_coeff, const Real lin_coeff,
                    const Real const_coeff, const Real prev_mass,
                    const Real prev_dp, const Real x2) const {
    // This remapping assumes we're starting from the left interface of an
    // old grid cell
    // In fact, we're usually integrating very little or almost all of the
    // cell in question
    const Real x1 = -0.5;
    const Real integral =
        integrate_parabola(sq_coeff, lin_coeff, const_coeff, x1, x2);
    const Real mass = prev_mass + integral * prev_dp;
    return mass;
  }

  template <typename ExecSpaceType = ExecSpace>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<!Homme::OnGpu<ExecSpaceType>::value, void>::type
  compute_remap(KernelVariables &/* kv */,
      ExecViewUnmanaged<const int[NUM_PHYSICAL_LEV]> k_id,
      ExecViewUnmanaged<const Real[NUM_PHYSICAL_LEV]> integral_bounds,
      ExecViewUnmanaged<const Real[3][NUM_PHYSICAL_LEV]> parabola_coeffs,
      ExecViewUnmanaged<Real[_ppm_consts::MASS_O_PHYSICAL_LEV]> mass,
      ExecViewUnmanaged<const Real[_ppm_consts::DPO_PHYSICAL_LEV]> prev_dp,
      ExecViewUnmanaged<Scalar[NUM_LEV]> remap_var) const {
    // Compute tracer values on the new grid by integrating from the old cell
    // bottom to the new cell interface to form a new grid mass accumulation.
    // Store the mass in the integral bounds for that level
    // Then take the difference between accumulation at successive interfaces
    // gives the mass inside each cell. Since Qdp is supposed to hold the full
    // mass this needs no normalization.
    Real mass1 = 0;
    Real mass2;
    ExecViewUnmanaged<Real[NUM_PHYSICAL_LEV]> rvar(reinterpret_cast<Real*>(remap_var.data()));
    for (int k=0; k<NUM_PHYSICAL_LEV; ++k) {
      const int kk_cur_lev = k_id(k);
      assert(kk_cur_lev < parabola_coeffs.extent_int(1));

      const Real x2_cur_lev = integral_bounds(k);
      // Repurpose the mass buffer to store the new mass.
      // WARNING: This may not be thread safe in future architectures which
      //          use this level of parallelism!!!
      mass2 = compute_mass(
          parabola_coeffs(2, kk_cur_lev), parabola_coeffs(1, kk_cur_lev),
          parabola_coeffs(0, kk_cur_lev), mass(kk_cur_lev),
          prev_dp(kk_cur_lev + _ppm_consts::INITIAL_PADDING), x2_cur_lev);
      rvar(k) = mass2 - mass1;
      mass1 = mass2;
    }
  }

  template <typename ExecSpaceType = ExecSpace>
  KOKKOS_INLINE_FUNCTION
  typename std::enable_if<Homme::OnGpu<ExecSpaceType>::value, void>::type
  compute_remap(KernelVariables &kv,
      ExecViewUnmanaged<const int[NUM_PHYSICAL_LEV]> k_id,
      ExecViewUnmanaged<const Real[NUM_PHYSICAL_LEV]> integral_bounds,
      ExecViewUnmanaged<const Real[3][NUM_PHYSICAL_LEV]> parabola_coeffs,
      ExecViewUnmanaged<Real[_ppm_consts::MASS_O_PHYSICAL_LEV]> prev_mass,
      ExecViewUnmanaged<const Real[_ppm_consts::DPO_PHYSICAL_LEV]> prev_dp,
      ExecViewUnmanaged<Scalar[NUM_LEV]> remap_var) const {
    // This duplicates work, but the parallel gain on CUDA is >> 2
    assert(VECTOR_SIZE==1);
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                         [&](const int k) {
      const Real mass_1 =
          (k > 0)
              ? compute_mass(
                    parabola_coeffs(2, k_id(k - 1)),
                    parabola_coeffs(1, k_id(k - 1)),
                    parabola_coeffs(0, k_id(k - 1)), prev_mass(k_id(k - 1)),
                    prev_dp(k_id(k - 1) + _ppm_consts::INITIAL_PADDING),
                    integral_bounds(k - 1))
              : 0.0;

      const Real x2_cur_lev = integral_bounds(k);

      const int kk_cur_lev = k_id(k);
      assert(kk_cur_lev < parabola_coeffs.extent_int(1));

      const Real mass_2 = compute_mass(
          parabola_coeffs(2, kk_cur_lev), parabola_coeffs(1, kk_cur_lev),
          parabola_coeffs(0, kk_cur_lev), prev_mass(kk_cur_lev),
          prev_dp(kk_cur_lev + _ppm_consts::INITIAL_PADDING), x2_cur_lev);

      remap_var(k)[0] = mass_2 - mass_1;
    }); // k loop
  }

  KOKKOS_INLINE_FUNCTION
  void compute_grids(KernelVariables &kv,
      const ExecViewUnmanaged<const Real[_ppm_consts::DPO_PHYSICAL_LEV]> dx,
      const ExecViewUnmanaged<Real[10][_ppm_consts::PPMDX_PHYSICAL_LEV]> grids) const
  {
    constexpr int dpo_offset = _ppm_consts::INITIAL_PADDING - _ppm_consts::gs;
    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,
                                                   NUM_PHYSICAL_LEV + 2),
                         [&](const int j) {
      grids(0, j) = dx(j + 1 + dpo_offset) /
                    (dx(j + dpo_offset) + dx(j + 1 + dpo_offset) +
                     dx(j + 2 + dpo_offset));

      grids(1, j) = (2.0 * dx(j + dpo_offset) + dx(j + 1 + dpo_offset)) /
                    (dx(j + 1 + dpo_offset) + dx(j + 2 + dpo_offset));

      grids(2, j) = (dx(j + 1 + dpo_offset) + 2.0 * dx(j + 2 + dpo_offset)) /
                    (dx(j + dpo_offset) + dx(j + 1 + dpo_offset));
    });

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,
                                                   NUM_PHYSICAL_LEV + 1),
                         [&](const int j) {
      grids(3, j) = dx(j + 1 + dpo_offset) /
                    (dx(j + 1 + dpo_offset) + dx(j + 2 + dpo_offset));

      grids(4, j) = 1.0 / (dx(j + dpo_offset) + dx(j + 1 + dpo_offset) +
                           dx(j + 2 + dpo_offset) + dx(j + 3 + dpo_offset));

      grids(5, j) = (2.0 * dx(j + 1 + dpo_offset) * dx(j + 2 + dpo_offset)) /
                    (dx(j + 1 + dpo_offset) + dx(j + 2 + dpo_offset));

      grids(6, j) = (dx(j + dpo_offset) + dx(j + 1 + dpo_offset)) /
                    (2.0 * dx(j + 1 + dpo_offset) + dx(j + 2 + dpo_offset));

      grids(7, j) = (dx(j + 3 + dpo_offset) + dx(j + 2 + dpo_offset)) /
                    (2.0 * dx(j + 2 + dpo_offset) + dx(j + 1 + dpo_offset));

      grids(8, j) = dx(j + 1 + dpo_offset) *
                    (dx(j + dpo_offset) + dx(j + 1 + dpo_offset)) /
                    (2.0 * dx(j + 1 + dpo_offset) + dx(j + 2 + dpo_offset));

      grids(9, j) = dx(j + 2 + dpo_offset) *
                    (dx(j + 2 + dpo_offset) + dx(j + 3 + dpo_offset)) /
                    (dx(j + 1 + dpo_offset) + 2.0 * dx(j + 2 + dpo_offset));
    });
  }

  KOKKOS_INLINE_FUNCTION
  void compute_ppm(KernelVariables &kv,
      // input  views
      ExecViewUnmanaged<const Real[_ppm_consts::AO_PHYSICAL_LEV]> cell_means,
      ExecViewUnmanaged<const Real[10][_ppm_consts::PPMDX_PHYSICAL_LEV]> dx,
      // buffer views
      ExecViewUnmanaged<Real[_ppm_consts::DMA_PHYSICAL_LEV]> dma,
      ExecViewUnmanaged<Real[_ppm_consts::AI_PHYSICAL_LEV]> ai,
      // result view
      ExecViewUnmanaged<Real[3][NUM_PHYSICAL_LEV]> parabola_coeffs) const
  {
    const auto INITIAL_PADDING = _ppm_consts::INITIAL_PADDING;

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team,
                                                   NUM_PHYSICAL_LEV + 2),
                         [&](const int j) {
      if ((cell_means(j + INITIAL_PADDING) -
           cell_means(j + INITIAL_PADDING - 1)) *
              (cell_means(j + INITIAL_PADDING - 1) -
               cell_means(j + INITIAL_PADDING - gs)) >
          0.0) {
        Real da =
            dx(0, j) * (dx(1, j) * (cell_means(j + INITIAL_PADDING) -
                                    cell_means(j + INITIAL_PADDING - 1)) +
                        dx(2, j) * (cell_means(j + INITIAL_PADDING - 1) -
                                    cell_means(j + INITIAL_PADDING - gs)));

        dma(j) = min(fabs(da), 2.0 * fabs(cell_means(j + INITIAL_PADDING - 1) -
                                          cell_means(j + INITIAL_PADDING - gs)),
                     2.0 * fabs(cell_means(j + INITIAL_PADDING) -
                                cell_means(j + INITIAL_PADDING - 1))) *
                 copysign(1.0, da);
      } else {
        dma(j) = 0.0;
      }
    });

    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV + 1),
        [&](const int j) {
          ai(j) = cell_means(j + INITIAL_PADDING - 1) +
                  dx(3, j) * (cell_means(j + INITIAL_PADDING) -
                              cell_means(j + INITIAL_PADDING - 1)) +
                  dx(4, j) * (dx(5, j) * (dx(6, j) - dx(7, j)) *
                                  (cell_means(j + INITIAL_PADDING) -
                                   cell_means(j + INITIAL_PADDING - 1)) -
                              dx(8, j) * dma(j + 1) + dx(9, j) * dma(j));
        });

    Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                         [&](const int j_prev) {
      const int j = j_prev + 1;
      Real al = ai(j - 1);
      Real ar = ai(j);
      if ((ar - cell_means(j + INITIAL_PADDING - 1)) *
              (cell_means(j + INITIAL_PADDING - 1) - al) <=
          0.) {
        al = cell_means(j + INITIAL_PADDING - 1);
        ar = cell_means(j + INITIAL_PADDING - 1);
      }
      if ((ar - al) * (cell_means(j + INITIAL_PADDING - 1) - (al + ar) / 2.0) >
          (ar - al) * (ar - al) / 6.0) {
        al = 3.0 * cell_means(j + INITIAL_PADDING - 1) - 2.0 * ar;
      }
      if ((ar - al) * (cell_means(j + INITIAL_PADDING - 1) - (al + ar) / 2.0) <
          -(ar - al) * (ar - al) / 6.0) {
        ar = 3.0 * cell_means(j + INITIAL_PADDING - 1) - 2.0 * al;
      }

      // Computed these coefficients from the edge values
      // and cell mean in Maple. Assumes normalized
      // coordinates: xi=(x-x0)/dx

      assert(parabola_coeffs.data() != nullptr);
      assert(j - 1 < parabola_coeffs.extent_int(1));
      assert(2 < parabola_coeffs.extent_int(0));

      parabola_coeffs(0, j - 1) =
          1.5 * cell_means(j + INITIAL_PADDING - 1) - (al + ar) / 4.0;
      parabola_coeffs(1, j - 1) = ar - al;
      parabola_coeffs(2, j - 1) =
          3.0 * (-2.0 * cell_means(j + INITIAL_PADDING - 1) + (al + ar));
    });

    Kokkos::single(Kokkos::PerThread(kv.team), [&]() {
      boundaries::apply_ppm_boundary(cell_means, parabola_coeffs);
    });
  }

  KOKKOS_INLINE_FUNCTION
  void compute_partitions(
      KernelVariables &kv,
      ExecViewUnmanaged<const Scalar[NP][NP][NUM_LEV]> src_layer_thickness,
      ExecViewUnmanaged<const Scalar[NP][NP][NUM_LEV]> tgt_layer_thickness)
      const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      ExecViewUnmanaged<Real[_ppm_consts::PIO_PHYSICAL_LEV]> pt_pio =
          Homme::subview(m_pio, kv.ie, igp, jgp);
      ExecViewUnmanaged<Real[_ppm_consts::PIN_PHYSICAL_LEV]> pt_pin =
          Homme::subview(m_pin, kv.ie, igp, jgp);
      ExecViewUnmanaged<const Scalar[NUM_LEV]> pt_src_thickness =
          Homme::subview(src_layer_thickness, igp, jgp);
      ExecViewUnmanaged<const Scalar[NUM_LEV]> pt_tgt_thickness =
          Homme::subview(tgt_layer_thickness, igp, jgp);

      Dispatch<ExecSpace>::parallel_scan(
          kv.team, NUM_PHYSICAL_LEV,
          [&](const int &level, Real &accumulator, const bool last) {
            if (last) {
              pt_pio(level) = accumulator;
            }
            const int ilev = level / VECTOR_SIZE;
            const int vlev = level % VECTOR_SIZE;
            accumulator += pt_src_thickness(ilev)[vlev];
          });
      Dispatch<ExecSpace>::parallel_scan(
          kv.team, NUM_PHYSICAL_LEV,
          [&](const int &level, Real &accumulator, const bool last) {
            if (last) {
              pt_pin(level) = accumulator;
            }
            const int ilev = level / VECTOR_SIZE;
            const int vlev = level % VECTOR_SIZE;
            accumulator += pt_tgt_thickness(ilev)[vlev];
          });

      Kokkos::single(Kokkos::PerThread(kv.team), [&]() {
        const int ilast = (NUM_PHYSICAL_LEV - 1) / VECTOR_SIZE;
        const int vlast = (NUM_PHYSICAL_LEV - 1) % VECTOR_SIZE;
        pt_pio(NUM_PHYSICAL_LEV) =
            pt_pio(NUM_PHYSICAL_LEV - 1) + pt_src_thickness(ilast)[vlast];
        pt_pin(NUM_PHYSICAL_LEV) =
            pt_pin(NUM_PHYSICAL_LEV - 1) + pt_tgt_thickness(ilast)[vlast];
        // This is here to allow an entire block of k
        // threads to run in the remapping phase. It makes
        // sure there's an old interface value below the
        // domain that is larger.
        assert(fabs(m_pio(kv.ie, igp, jgp, NUM_PHYSICAL_LEV) -
                    m_pin(kv.ie, igp, jgp, NUM_PHYSICAL_LEV)) < 1.0);

        pt_pio(_ppm_consts::PIO_PHYSICAL_LEV - 1) =
            pt_pio(_ppm_consts::PIO_PHYSICAL_LEV - 2) + 1.0;

        // The total mass in a column does not change.
        // Therefore, the pressure of that mass cannot
        // either.
        pt_pin(NUM_PHYSICAL_LEV) = pt_pio(NUM_PHYSICAL_LEV);
      });

      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                           [&](const int &k) {
        int ilevel = k / VECTOR_SIZE;
        int ivector = k % VECTOR_SIZE;
        m_dpo(kv.ie, igp, jgp, k + _ppm_consts::INITIAL_PADDING) =
            src_layer_thickness(igp, jgp, ilevel)[ivector];
      });
    });
    kv.team_barrier();
    // Fill in the ghost regions with mirrored values.
    // if vert_remap_q_alg is defined, this is of no
    // consequence.
    // Note that the range of k makes this completely parallel,
    // without any data dependencies
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, gs),
                           [&](const int &k) {
        m_dpo(kv.ie, igp, jgp, _ppm_consts::INITIAL_PADDING - 1 - k) =
            m_dpo(kv.ie, igp, jgp, k + _ppm_consts::INITIAL_PADDING);
        m_dpo(kv.ie, igp, jgp,
            NUM_PHYSICAL_LEV + _ppm_consts::INITIAL_PADDING + k) =
            m_dpo(kv.ie, igp, jgp,
                NUM_PHYSICAL_LEV + _ppm_consts::INITIAL_PADDING - 1 - k);
      });
    });
    kv.team_barrier();
  }

  KOKKOS_INLINE_FUNCTION
  void compute_integral_bounds(KernelVariables &kv) const {
    Kokkos::parallel_for(Kokkos::TeamThreadRange(kv.team, NP * NP),
                         [&](const int &loop_idx) {
      const int igp = loop_idx / NP;
      const int jgp = loop_idx % NP;
      Kokkos::parallel_for(Kokkos::ThreadVectorRange(kv.team, NUM_PHYSICAL_LEV),
                           [&](const int k) {
        // Compute remapping intervals once for all
        // tracers. Find the old grid cell index in which
        // the k-th new cell interface resides. Then
        // integrate from the bottom of that old cell to
        // the new interface location. In practice, the
        // grid never deforms past one cell, so the search
        // can be simplified by this. Also, the interval
        // of integration is usually of magnitude close to
        // zero or close to dpo because of minimial
        // deformation. Numerous tests confirmed that the
        // bottom and top of the grids match to machine
        // precision, so set them equal to each other.
        int kk = k;// + 1;
        // This reduces the work required to find the index where this
        // fails at, and is typically less than NUM_PHYSICAL_LEV^2 Since
        // the top bounds match anyway, the value of the coefficients
        // don't matter, so enforcing kk <= NUM_PHYSICAL_LEV doesn't
        // affect anything important
        //
        // Note that because we set
        // pio(:, :, :, NUM_PHYSICAL_LEV + 1) = pio(:, :, :,
        // NUM_PHYSICAL_LEV) + 1.0 and pin(:, :, :, NUM_PHYSICAL_LEV) =
        // pio(:, :, :, NUM_PHYSICAL_LEV) this loop ensures kk <=
        // NUM_PHYSICAL_LEV + 2 Furthermore, since we set pio(:, :, :,
        // 0) = 0.0 and pin(:, :, :, 0) = 0.0 kk must be incremented at
        // least once
        assert(m_pio(kv.ie, igp, jgp, _ppm_consts::PIO_PHYSICAL_LEV - 1) >
               m_pin(kv.ie, igp, jgp, k + 1));
        if (m_pio(kv.ie,igp,jgp,kk) <= m_pin(kv.ie,igp,jgp,k+1)) {
          while (m_pio(kv.ie, igp, jgp, kk) <= m_pin(kv.ie, igp, jgp, k + 1)) {
            kk++;
            assert(kk < m_pio.extent_int(3));
          }
          --kk;
        } else {
          binary_search(Homme::subview(m_pio,kv.ie,igp,jgp),m_pin(kv.ie,igp,jgp,k+1),kk);
        }

        // This is to keep the indices in bounds.
        if (kk == NUM_PHYSICAL_LEV) {
          --kk;
        }
        // kk is now the cell index we're integrating over.

        // Save kk for reuse
        m_kid(kv.ie, igp, jgp, k) = kk;
        // PPM interpolants are normalized to an independent coordinate
        // domain
        // [-0.5, 0.5].
        m_z2(kv.ie, igp, jgp, k) =
            (m_pin(kv.ie, igp, jgp, k + 1) -
             (m_pio(kv.ie, igp, jgp, kk) + m_pio(kv.ie, igp, jgp, kk+1)) * 0.5) /
            m_dpo(kv.ie, igp, jgp, kk + _ppm_consts::INITIAL_PADDING);
      });

      auto point_dpo   = Homme::subview(m_dpo, kv.ie, igp, jgp);
      auto point_ppmdx = Homme::subview(m_ppmdx, kv.ie, igp, jgp);
      compute_grids(kv, point_dpo, point_ppmdx);
    });
  }

  KOKKOS_FORCEINLINE_FUNCTION Real
  integrate_parabola(const Real sq_coeff, const Real lin_coeff,
                     const Real const_coeff, Real x1, Real x2) const {
    return (const_coeff * (x2 - x1) + lin_coeff * (x2 * x2 - x1 * x1) / 2.0) +
           sq_coeff * (x2 * x2 * x2 - x1 * x1 * x1) / 3.0;
  }

  ExecViewManaged<Real * [NP][NP][_ppm_consts::DPO_PHYSICAL_LEV]> m_dpo;
  // pio corresponds to the points in each layer of the source layer thickness
  ExecViewManaged<Real * [NP][NP][_ppm_consts::PIO_PHYSICAL_LEV]> m_pio;
  // pin corresponds to the points in each layer of the target layer thickness
  ExecViewManaged<Real * [NP][NP][_ppm_consts::PIN_PHYSICAL_LEV]> m_pin;
  ExecViewManaged<Real * [NP][NP][10][_ppm_consts::PPMDX_PHYSICAL_LEV]> m_ppmdx;
  ExecViewManaged<Real * [NP][NP][NUM_PHYSICAL_LEV]>  m_z2;
  ExecViewManaged<int * [NP][NP][NUM_PHYSICAL_LEV]>   m_kid;

  TeamUtils<ExecSpace> m_ppm_tu;
  ExecViewManaged<Real * [NP][NP][_ppm_consts::AO_PHYSICAL_LEV]> m_ao;
  ExecViewManaged<Real * [NP][NP][_ppm_consts::MASS_O_PHYSICAL_LEV]> m_mass_o;
  ExecViewManaged<Real * [NP][NP][_ppm_consts::DMA_PHYSICAL_LEV]> m_dma;
  ExecViewManaged<Real * [NP][NP][_ppm_consts::AI_PHYSICAL_LEV]> m_ai;
  ExecViewManaged<Real * [NP][NP][3][NUM_PHYSICAL_LEV]> m_parabola_coeffs;
};

} // namespace Ppm
} // namespace Remap
} // namespace Homme

#endif // HOMMEXX_PPM_REMAP_HPP
