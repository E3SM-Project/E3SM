#ifndef OMEGA_HORZUTIL_H
#define OMEGA_HORZUTIL_H

#include "DataTypes.h"
#include "Logging.h"
namespace OMEGA {

KOKKOS_INLINE_FUNCTION Real distance(const Real x, const Real y, const Real z) {
   const Real dist = Kokkos::sqrt(x * x + y * y + z * z);
   return dist;
}

KOKKOS_INLINE_FUNCTION Real sphere_angle(const Real ax, const Real ay,
                                         const Real az, const Real bx,
                                         const Real by, const Real bz,
                                         const Real cx, const Real cy,
                                         const Real cz) {
   // Computes the angle between arcs AB and AC, given points A, B, and C
   // Equation numbers w.r.t.
   // http://mathworld.wolfram.com/SphericalTrigonometry.html
   const Real one     = 1._Real;
   const Real neg_one = -1._Real;
   const auto a =
       Kokkos::acos(Kokkos::max(Kokkos::min(bx * cx + by * cy + bz * cz, one),
                                neg_one)); // Eqn. (3)
   const auto b =
       Kokkos::acos(Kokkos::max(Kokkos::min(ax * cx + ay * cy + az * cz, one),
                                neg_one)); // Eqn. (2)
   const auto c =
       Kokkos::acos(Kokkos::max(Kokkos::min(ax * bx + ay * by + az * bz, one),
                                neg_one)); // Eqn. (1)
   const auto ABx = bx - ax;
   const auto ABy = by - ay;
   const auto ABz = bz - az;

   const auto ACx = cx - ax;
   const auto ACy = cy - ay;
   const auto ACz = cz - az;

   const auto Dx = (ABy * ACz) - (ABz * ACy);
   const auto Dy = -((ABx * ACz) - (ABz * ACx));
   const auto Dz = (ABx * ACy) - (ABy * ACx);

   const auto s         = 0.5_Real * (a + b + c);
   const auto sin_angle = Kokkos::sqrt(Kokkos::min(
       one, Kokkos::max(0._Real,
                        (Kokkos::sin(s - b) * Kokkos::sin(s - c)) /
                            (Kokkos::sin(b) * Kokkos::sin(c))))); // Eqn. (28)
   Real sa              = 2._Real *
             Kokkos::asin(Kokkos::max(Kokkos::min(sin_angle, one), neg_one));
   if ((Dx * ax + Dy * ay + Dz * az) < 0.0)
      sa *= neg_one;
   return sa;
}
KOKKOS_INLINE_FUNCTION Real arc_length(const Real ax, const Real ay,
                                       const Real az, const Real bx,
                                       const Real by, const Real bz) {
   // Returns the length of the great circle arc from A=(ax, ay, az) to
   // B=(bx, by, bz). It is assumed that both A and B lie on the surface of the
   // same sphere centered at the origin.
   const auto cx = bx - ax;
   const auto cy = by - ay;
   const auto cz = bz - az;
   const auto r  = Kokkos::sqrt(ax * ax + ay * ay + az * az);
   const auto c  = Kokkos::sqrt(cx * cx + cy * cy + cz * cz);
   const auto al = r * 2.0 * Kokkos::asin(c / (2.0 * r));
   return al;
}
template <int NA>
static KOKKOS_FUNCTION void mat_t_mul(Real AtA[NA][NA], const Real A[][NA],
                                      const I4 Rows) {
   // Compute AtA = A^t * A for the first N rows of A
   const I4 Cols = NA;
   for (int I = 0; I < Cols; ++I)
      for (int J = 0; J < Cols; ++J)
         for (int K = 0; K < Rows; ++K)
            AtA[I][J] += A[K][I] * A[K][J]; // = A^t(I,K) * A(K,j)
}
template <int NA>
KOKKOS_INLINE_FUNCTION void elgs(Real A[NA][NA], I4 Indx[NA]) {
   // subroutine to perform the partial-pivoting gaussian elimination.
   // a(n,n) is the original matrix in the input and transformed matrix
   // plus the pivoting element ratios below the diagonal in the output.
   // indx(n) records the pivoting order.  tao pang 2001.
   //
   // A is a square matrix NA by NA but only the first N rows,columns are
   // pivoted.
   const I4 N = NA;
   for (int I = 0; I < N; ++I)
      Indx[I] = I;
   // find the rescaling factors, one from each row
   Real C[NA] = {};
   for (int I = 0; I < N; ++I)
      for (int J = 0; J < N; ++J)
         C[I] = Kokkos::max(C[I], Kokkos::abs(A[I][J]));

   // search the pivoting (largest) element from each column
   for (int J = 0; J < N - 1; ++J) {
      Real pi1 = 0._Real;
      I4 K     = 0;
      for (int I = J; I < N; ++I) {
         const Real pi = Kokkos::abs(A[Indx[I]][J]) / C[Indx[I]];
         if (pi1 < pi) {
            pi1 = pi;
            K   = I;
         }
      }
      // interchange the rows via indx(n) to record pivoting order
      const I4 ITmp = Indx[J];
      Indx[J]       = Indx[K];
      Indx[K]       = ITmp;
      for (int I = J + 1; I < N; ++I) {
         const Real pj = A[Indx[I]][J] / A[Indx[J]][J];
         // record pivoting ratios below the diagonal
         A[Indx[I]][J] = pj;
         // modify other elements accordingly
         for (int K = J + 1; K < N; ++K)
            A[Indx[I]][K] -= pj * A[Indx[J]][K];
      }
   }
}
template <int NA>
KOKKOS_INLINE_FUNCTION void migs(Real A[NA][NA], Real X[NA][NA]) {
   // subroutine to invert matrix a(n,n) with the inverse stored
   // in x(n,n) in the output. tao pang 2001.
   //
   // Matrix sizes are NA by NA but only the first N rowx & columns are
   // inverted.
   I4 Indx[NA]    = {};
   Real B[NA][NA] = {};
   for (int I = 0; I < NA; ++I)
      B[I][I] = 1._Real;
   elgs<NA>(A, Indx);
   for (int I = 0; I < NA - 1; ++I)
      for (int J = I + 1; J < NA; ++J)
         for (int K = 0; K < NA; ++K)
            B[Indx[J]][K] -= A[Indx[J]][I] * B[Indx[I]][K];
   for (int I = 0; I < NA; ++I) {
      X[NA - 1][I] = B[Indx[NA - 1]][I] / A[Indx[NA - 1]][NA - 1];
      for (int J = NA - 2; 0 <= J; --J) {
         X[J][I] = B[Indx[J]][I];
         for (int K = J + 1; K < NA; ++K)
            X[J][I] -= A[Indx[J]][K] * X[K][I];
         X[J][I] /= A[Indx[J]][J];
      }
   }
}
template <int NA>
KOKKOS_INLINE_FUNCTION void mat_mul_t(Array2DReal C, const Real B[NA][NA],
                                      const Real A[][NA], const I4 Rows) {
   // Compute  C = B * A^t for first N columns of C
   // C assumed to be initialized to 0.
   // The number of rows in C is inferred from the number of column in A and B.
   const I4 Cols = NA;
   for (int I = 0; I < Cols; ++I)
      for (int J = 0; J < Rows; ++J)
         C(I, J) = 0;
   for (int I = 0; I < Cols; ++I)
      for (int J = 0; J < Rows; ++J)
         for (int K = 0; K < Cols; ++K)
            C(I, J) += B[I][K] * A[J][K]; // Transpose A
}
template <int NA>
KOKKOS_INLINE_FUNCTION void poly_fit_2(const Real A[][NA], Array2DReal B,
                                       const I4 Rows) {
   Real AtA[NA][NA]    = {};
   Real AtAInv[NA][NA] = {};
   mat_t_mul<NA>(AtA, A, Rows);
   migs<NA>(AtA, AtAInv);
   mat_mul_t<NA>(B, AtAInv, A, Rows);
}

KOKKOS_INLINE_FUNCTION void arc_bisect(const Real ax, const Real ay,
                                       const Real az, const Real bx,
                                       const Real by, const Real bz, Real &cx,
                                       Real &cy, Real &cz) {

   // Returns the point C=(cx, cy, cz) that bisects the great circle arc from
   //   A=(ax, ay, az) to B=(bx, by, bz). It is assumed that A and B lie on the
   //   surface of a sphere centered at the origin.

   cx = 0.5_Real * (ax + bx);
   cy = 0.5_Real * (ay + by);
   cz = 0.5_Real * (az + bz);
   if (cx == 0. && cy == 0. && cz == 0.) {
      printf("arc_bisect: A and B are diametrically opposite");
   } else {
      const Real d = Kokkos::sqrt(cx * cx + cy * cy + cz * cz);
      const Real r = Kokkos::sqrt(ax * ax + ay * ay + az * az);
      cx           = r * cx / d;
      cy           = r * cy / d;
      cz           = r * cz / d;
   }
}

} // namespace OMEGA
#endif
