// COMPOSE version 1.0: Copyright 2018 NTESS. This software is released under
// the BSD license; see LICENSE in the top-level directory.

#ifndef INCLUDE_CEDR_LOCAL_HPP
#define INCLUDE_CEDR_LOCAL_HPP

#include "cedr.hpp"
#include "cedr_kokkos.hpp"

namespace cedr {
namespace local {

// The following routines solve
//     min_x norm(x - y; w)
//      st   a'x = b
//           xlo <= x <= xhi,
// a > 0, w > 0.

// Minimize the weighted 2-norm. Return 0 on success and x == y, 1 on success
// and x != y, -1 if infeasible, -2 if max_its hit with no solution. See section
// 3 of Bochev, Ridzal, Shashkov, Fast optimization-based conservative remap of
// scalar fields through aggregate mass transfer.
KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp(const Int n, const Real* w, const Real* a, const Real b,
                    const Real* xlo, const Real* xhi,
                    const Real* y, Real* x, const Int max_its = 100);

KOKKOS_INLINE_FUNCTION
Int solve_1eq_bc_qp_2d(const Real* w, const Real* a, const Real b,
                       const Real* xlo, const Real* xhi,
                       const Real* y, Real* x,
                       const bool clip = true,
                       // The 2D algorithm doesn't need to terminate based on a
                       // tolerance, but by default it will exit early if the ND
                       // algorithm would. Set this to false to run the full 2D
                       // algorithm w/o checking the tolerance at start. If this
                       // is false, the feasibility check is also disabled.
                       const bool early_exit_on_tol = true);

// ClipAndAssuredSum. Minimize the 1-norm with w = 1s. Does not check for
// feasibility.
KOKKOS_INLINE_FUNCTION
void caas(const Int n, const Real* a, const Real b,
          const Real* xlo, const Real* xhi,
          const Real* y, Real* x,
          const bool clip = true);

struct Method { enum Enum { least_squares, caas }; };

// Solve
//     min_x norm(x - y; w)
//      st   a'x = b
//           x >= 0,
// a, w > 0. Return 0 on success and x == y, 1 on success and x != y, -1 if
// infeasible. w is used only if lcl_method = least_squares.
KOKKOS_INLINE_FUNCTION
Int solve_1eq_nonneg(const Int n, const Real* a, const Real b, const Real* y, Real* x,
                     const Real* w, const Method::Enum lcl_method);

Int unittest();

} // namespace local
} // namespace cedr

#include "cedr_local_inl.hpp"

#endif
